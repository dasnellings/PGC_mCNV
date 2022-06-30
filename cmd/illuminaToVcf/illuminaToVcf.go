package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/PGC_mCNV/illumina"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"log"
	"path"
	"strings"
)

const debug int = 1

const headerInfo string = "##fileformat=VCFv4.2\n" +
	"##INFO=<ID=ALLELE_A,Number=1,Type=Integer,Description=\"A allele\">\n" +
	"##INFO=<ID=ALLELE_B,Number=1,Type=Integer,Description=\"B allele\">\n" +
	"##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC ratio content around the variant\">\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">\n" +
	"##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

func usage() {
	fmt.Print(
		"illuminaToVcf - Convert SNP array data from GenomeStudio report format to VCF format.\n" +
			"Usage:\n" +
			"./illuminaToVcf [options] -gsReport sample1,sample2 -manifest arrayManifest.csv -ref reference.fasta\n\n")
	flag.PrintDefaults()
}

func main() {
	gsReportFilename := flag.String("gsReport", "", "GenomeStudio summary report file. The file "+
		"should be named by sample and have tab seperated fields. May be a comma-seperated list of files.")
	manifestFilename := flag.String("manifest", "", "Manifest file for the array used (.csv)")
	fastaFilename := flag.String("ref", "", "Reference fasta file for the assembly used for the GenomeStudio report.")
	output := flag.String("o", "stdout", "Output VCF file")
	mapmode := flag.Bool("hash", false, "Hash map lookup for snp IDs. Use for out of order data.")
	silent := flag.Bool("suppress", false, "Prevent warning messages.")
	flag.Parse()

	if *gsReportFilename == "" || *manifestFilename == "" || *fastaFilename == "" {
		usage()
		log.Fatal("ERROR: GenomeStudio report, manifest, and reference fasta files are required (-gsReport, -manifest, -ref)")
	}

	if *mapmode {
		illuminaToVcfMap(strings.Split(*gsReportFilename, ","), *manifestFilename, *fastaFilename, *output, *silent)
	} else {
		illuminaToVcf(strings.Split(*gsReportFilename, ","), *manifestFilename, *fastaFilename, *output, *silent)
	}
}

func illuminaToVcf(gsReportFiles []string, manifestFile, fastaFile, output string, silent bool) {
	out := fileio.EasyCreate(output)
	ref := fasta.NewSeeker(fastaFile, fastaFile+".fai")
	var header vcf.Header
	header.Text = strings.Split(headerInfo, "\n")

	trimSamples := slices.Clone(gsReportFiles)
	for i := range trimSamples {
		trimSamples[i] = strings.TrimRight(path.Base(trimSamples[i]), ".gz")
	}
	header.Text[len(header.Text)-1] += "\t" + strings.Join(trimSamples, "\t")
	vcf.NewWriteHeader(out, header)

	gsReportChans := make([]<-chan illumina.GsReport, len(gsReportFiles))
	for i := range gsReportFiles {
		gsReportChans[i] = illumina.GoReadGsReportToChan(gsReportFiles[i])
	}
	manifestData := illumina.GoReadManifestToChan(manifestFile)

	var err error
	var curr vcf.Vcf
	var gs illumina.GsReport
	curr.Filter = "."
	curr.Format = []string{"GT", "BAF", "LRR"}
	sb := new(strings.Builder)
	var alleleAint, alleleBint int16
	var alleleA, alleleB, gsAllele1, gsAllele2 string
	var seqBefore, seqAfter []dna.Base
	var stringBefore, stringAfter string
	var refBase []dna.Base
	var altNeedsRevComp bool
	var samplesWritten int

	for m := range manifestData {
		if m.Chr == "XY" || m.Chr == "chrXY" { // SERIOUSLY ILLUMINA... SERIOUSLY
			m.Chr = "X"
		}
		curr.Chr = "chr" + strings.TrimLeft(m.Chr, "chr")
		curr.Pos = m.Pos
		curr.Id = m.Name
		refBase, err = fasta.SeekByName(ref, "chr"+strings.TrimLeft(m.Chr, "chr"), m.Pos-1, m.Pos)
		exception.PanicOnErr(err)
		curr.Ref = strings.ToUpper(dna.BaseToString(refBase[0]))

		seqBefore, err = fasta.SeekByName(ref, "chr"+strings.TrimLeft(m.Chr, "chr"), (m.Pos-1)-len(m.SeqBefore), m.Pos-1)
		exception.PanicOnErr(err)
		stringBefore = strings.ToUpper(dna.BasesToString(seqBefore))
		seqAfter, err = fasta.SeekByName(ref, "chr"+strings.TrimLeft(m.Chr, "chr"), m.Pos, m.Pos+len(m.SeqAfter))
		if err != nil && !silent {
			fmt.Println("WARNING", err)
		}
		stringAfter = strings.ToUpper(dna.BasesToString(seqAfter))

		// check one of the alleles matches ref
		altNeedsRevComp = false
		switch {
		case levenshtein(stringBefore, m.SeqBefore) <= 5 ||
			levenshtein(stringAfter, m.SeqAfter) <= 5: // this is a really weak match, but you would not believe the things I have seen...
			if !m.TopStrand {
				altNeedsRevComp = true
			}

			// only do partial check on rev comps since if snp is not directly in middle of probe then before/after lengths differ
		case levenshtein(revComp(stringBefore)[:5], m.SeqAfter[:5]) <= 1 ||
			levenshtein(revComp(stringAfter)[len(stringAfter)-5:], m.SeqBefore[len(m.SeqBefore)-5:]) <= 1:
			if m.TopStrand {
				altNeedsRevComp = true
			}

		default:
			if !silent {
				log.Printf("WARNING: Context sequences did not match reference:\n%s+%s\n%s+%s\n", stringBefore, stringAfter, m.SeqBefore, m.SeqAfter)
				log.Println(m.Name, m.Chr, m.Pos)
			}
		}

		if altNeedsRevComp {
			alleleA = revComp(m.AlleleA)
			alleleB = revComp(m.AlleleB)
		} else {
			alleleA = m.AlleleA
			alleleB = m.AlleleB
		}

		switch curr.Ref {
		case alleleA:
			alleleAint = 0
			alleleBint = 1
			curr.Alt = []string{alleleB}
		case alleleB:
			alleleAint = 1
			alleleBint = 0
			curr.Alt = []string{alleleA}
		default:
			if alleleA == alleleB {
				alleleAint = 1
				alleleBint = 1
				curr.Alt = []string{alleleA}
			} else {
				alleleAint = 1
				alleleBint = 2
				curr.Alt = []string{alleleA, alleleB}
			}
		}

		curr.Info = fmt.Sprintf("ALLELE_A=%d;ALLELE_B=%d;GC=%.4g", alleleAint, alleleBint, m.GC)
		curr.Samples = make([]vcf.Sample, len(gsReportChans))
		sb.Reset()
		samplesWritten = 0
		for i := range curr.Samples {
			for gs.Chrom == "" || gs.Chrom == "0" {
				gs = <-gsReportChans[i]
				switch gs.Chrom {
				case "xy":
					gs.Chrom = "x"
				case "XY":
					gs.Chrom = "X"
				case "MT":
					gs.Chrom = "M"
				}
			}

			if !matchesManifest(gs, m) {
				if i != 0 {
					log.Panicf("something went horibly wrong with sample %s\n%v", gsReportFiles[i], gs)
				}
				if !silent {
					log.Printf("WARNING: Manifest mismatch. See report and manifest data below\n%v\n%v\n", gs, m)
					log.Printf("waiting for %v", gs)
					log.Println("moving to next manifest record")
				}
				break
			}
			samplesWritten++
			gsAllele1 = gs.Allele1
			gsAllele2 = gs.Allele2
			if (gs.ReportedAsFwd && m.TopStrand != m.SrcTopStrand) || (!gs.ReportedAsFwd && !m.TopStrand) {
				gsAllele1 = revComp(gsAllele1)
				gsAllele2 = revComp(gsAllele2)
			}

			curr.Samples[i].FormatData = []string{"", fmt.Sprintf("%.4g", gs.BAlleleFreq), fmt.Sprintf("%.4g", gs.LogRRatio)}
			switch gsAllele1 {
			case m.AlleleA:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleAint)
			case m.AlleleB:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleBint)
			}
			switch gsAllele2 {
			case m.AlleleA:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleAint)
			case m.AlleleB:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleBint)
			}
			curr.Samples[i].Phase = make([]bool, len(curr.Samples[i].Alleles)) // leave as false for unphased
			gs.Chrom = ""
		}
		if samplesWritten > 0 && curr.Chr != "chrM" { // exclude chrM
			vcf.WriteVcf(out, curr)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
	err = ref.Close()
	exception.PanicOnErr(err)
}

func illuminaToVcfMap(gsReportFiles []string, manifestFile, fastaFile, output string, silent bool) {
	out := fileio.EasyCreate(output)
	ref := fasta.NewSeeker(fastaFile, fastaFile+".fai")
	var header vcf.Header
	header.Text = strings.Split(headerInfo, "\n")

	trimSamples := slices.Clone(gsReportFiles)
	for i := range trimSamples {
		trimSamples[i] = strings.TrimRight(path.Base(trimSamples[i]), ".gz")
	}
	header.Text[len(header.Text)-1] += "\t" + strings.Join(trimSamples, "\t")
	vcf.NewWriteHeader(out, header)

	gsReportChans := make([]<-chan illumina.GsReport, len(gsReportFiles))
	for i := range gsReportFiles {
		gsReportChans[i] = illumina.GoReadGsReportToChan(gsReportFiles[i])
	}
	mm := makeManifestMap(manifestFile)

	var err error
	var curr vcf.Vcf
	var gs illumina.GsReport
	curr.Filter = "."
	curr.Format = []string{"GT", "BAF", "LRR"}
	var alleleAint, alleleBint int16
	var alleleA, alleleB, gsAllele1, gsAllele2 string
	var seqBefore, seqAfter []dna.Base
	var stringBefore, stringAfter string
	var refBase []dna.Base
	var altNeedsRevComp, found bool
	var samplesWritten int
	var m illumina.Manifest

	for gs = range gsReportChans[0] {
		if debug > 0 {
			fmt.Println("debug: started - ", gs)
		}
		for gs.Chrom == "" || gs.Chrom == "0" {
			if debug > 0 {
				fmt.Println("debug: no chrom for - ", gs)
			}
			for i := 1; i < len(gsReportChans); i++ {
				if debug > 0 {
					fmt.Println("debug: burning - ", <-gsReportChans[i])
				} else {
					<-gsReportChans[i] // burn
				}
			}
			gs = <-gsReportChans[0]
		}
		switch gs.Chrom {
		case "xy":
			gs.Chrom = "x"
		case "XY":
			gs.Chrom = "X"
		case "MT":
			gs.Chrom = "M"
		}

		samplesWritten = 0
		m, found = mm[strings.ToLower(gs.Marker)]
		if !found {
			m.Chr = "NOT_FOUND"
			for i := 1; i < len(gsReportChans); i++ {
				<-gsReportChans[i] // burn
			}
			continue
		}
		if m.Chr == "XY" || m.Chr == "chrXY" { // SERIOUSLY ILLUMINA... SERIOUSLY
			m.Chr = "X"
		}
		curr.Chr = "chr" + strings.TrimLeft(m.Chr, "chr")
		curr.Pos = m.Pos
		curr.Id = m.Name
		refBase, err = fasta.SeekByName(ref, "chr"+strings.TrimLeft(m.Chr, "chr"), m.Pos-1, m.Pos)
		exception.PanicOnErr(err)
		curr.Ref = strings.ToUpper(dna.BaseToString(refBase[0]))

		seqBefore, err = fasta.SeekByName(ref, "chr"+strings.TrimLeft(m.Chr, "chr"), (m.Pos-1)-len(m.SeqBefore), m.Pos-1)
		exception.PanicOnErr(err)
		stringBefore = strings.ToUpper(dna.BasesToString(seqBefore))
		seqAfter, err = fasta.SeekByName(ref, "chr"+strings.TrimLeft(m.Chr, "chr"), m.Pos, m.Pos+len(m.SeqAfter))
		if err != nil && !silent {
			fmt.Println("WARNING", err)
		}
		stringAfter = strings.ToUpper(dna.BasesToString(seqAfter))

		// check one of the alleles matches ref
		altNeedsRevComp = false
		switch {
		case levenshtein(stringBefore, m.SeqBefore) <= 5 ||
			levenshtein(stringAfter, m.SeqAfter) <= 5: // this is a really weak match, but you would not believe the things I have seen...
			if !m.TopStrand {
				altNeedsRevComp = true
			}

			// only do partial check on rev comps since if snp is not directly in middle of probe then before/after lengths differ
		case levenshtein(revComp(stringBefore)[:5], m.SeqAfter[:5]) <= 1 ||
			levenshtein(revComp(stringAfter)[len(stringAfter)-5:], m.SeqBefore[len(m.SeqBefore)-5:]) <= 1:
			if m.TopStrand {
				altNeedsRevComp = true
			}

		default:
			log.Printf("WARNING: Context sequences did not match reference:\n%s+%s\n%s+%s\n", stringBefore, stringAfter, m.SeqBefore, m.SeqAfter)
			log.Println(m.Name, m.Chr, m.Pos)
		}

		if altNeedsRevComp {
			alleleA = revComp(m.AlleleA)
			alleleB = revComp(m.AlleleB)
		} else {
			alleleA = m.AlleleA
			alleleB = m.AlleleB
		}

		switch curr.Ref {
		case alleleA:
			alleleAint = 0
			alleleBint = 1
			curr.Alt = []string{alleleB}
		case alleleB:
			alleleAint = 1
			alleleBint = 0
			curr.Alt = []string{alleleA}
		default:
			if alleleA == alleleB {
				alleleAint = 1
				alleleBint = 1
				curr.Alt = []string{alleleA}
			} else {
				alleleAint = 1
				alleleBint = 2
				curr.Alt = []string{alleleA, alleleB}
			}
		}

		curr.Info = fmt.Sprintf("ALLELE_A=%d;ALLELE_B=%d;GC=%.4g", alleleAint, alleleBint, m.GC)
		curr.Samples = make([]vcf.Sample, len(gsReportChans))

		for i := 0; i < len(curr.Samples); i++ {
			if i > 0 {
				gs = <-gsReportChans[i]
				for gs.Chrom == "" || gs.Chrom == "0" {
					gs = <-gsReportChans[0]
					log.Println("skipped", gs)
				}
				switch gs.Chrom {
				case "xy":
					gs.Chrom = "x"
				case "XY":
					gs.Chrom = "X"
				case "MT":
					gs.Chrom = "M"
				}
			}
			if m.Chr == "NOT_FOUND" {
				continue
			}
			if strings.ToLower(m.Name) != strings.ToLower(gs.Marker) {
				log.Print(m)
				log.Print(gsReportFiles[i], gs)
				log.Panic("PANIC!!! DATA OUT OF ORDER")
			}

			if !matchesManifest(gs, m) && !silent {
				log.Printf("WARNING: Manifest mismatch. See report and manifest data below\n%v\n%v\n", gs, m)
			}
			samplesWritten++
			gsAllele1 = gs.Allele1
			gsAllele2 = gs.Allele2
			if (gs.ReportedAsFwd && m.TopStrand != m.SrcTopStrand) || (!gs.ReportedAsFwd && !m.TopStrand) {
				gsAllele1 = revComp(gsAllele1)
				gsAllele2 = revComp(gsAllele2)
			}

			curr.Samples[i].FormatData = []string{"", fmt.Sprintf("%.4g", gs.BAlleleFreq), fmt.Sprintf("%.4g", gs.LogRRatio)}
			switch gsAllele1 {
			case m.AlleleA:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleAint)
			case m.AlleleB:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleBint)
			}
			switch gsAllele2 {
			case m.AlleleA:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleAint)
			case m.AlleleB:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleBint)
			}
			curr.Samples[i].Phase = make([]bool, len(curr.Samples[i].Alleles)) // leave as false for unphased
			gs.Chrom = ""
		}
		if samplesWritten > 0 && curr.Chr != "chrM" { // exclude chrM
			vcf.WriteVcf(out, curr)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
	err = ref.Close()
	exception.PanicOnErr(err)
}

func makeManifestMap(manifest string) map[string]illumina.Manifest {
	var found bool
	m := make(map[string]illumina.Manifest)
	c := illumina.GoReadManifestToChan(manifest)
	for v := range c {
		if _, found = m[strings.ToLower(v.Name)]; !found {
			m[strings.ToLower(v.Name)] = v
		}
	}
	return m
}

func matchesManifest(gs illumina.GsReport, m illumina.Manifest) bool {
	if strings.ToUpper(gs.Marker) != strings.ToUpper(m.Name) {
		return false
	}
	if "chr"+strings.ToUpper(gs.Chrom) != "chr"+strings.TrimLeft(m.Chr, "chr") {
		return false
	}
	if gs.Pos != m.Pos {
		return false
	}
	return true
}

func revComp(base string) string {
	ans := make([]byte, len(base))
	var j int
	for i := len(base) - 1; i >= 0; i-- {
		switch base[i] {
		case 'A', 'a':
			ans[j] = 'T'
		case 'C', 'c':
			ans[j] = 'G'
		case 'G', 'g':
			ans[j] = 'C'
		case 'T', 't':
			ans[j] = 'A'
		default:
			//log.Panicf("ERROR: unrecognized base '%s'", base) // they throw Y's and shit in there
			ans[j] = base[i]
		}
		j++
	}
	return string(ans)
}

func levenshtein(s1, s2 string) int {
	if s1 == "" || s2 == "" {
		return numbers.Max(len(s1), len(s2))
	}
	s1len := len(s1)
	s2len := len(s2)
	column := make([]int, len(s1)+1)

	for y := 1; y <= s1len; y++ {
		column[y] = y
	}
	for x := 1; x <= s2len; x++ {
		column[0] = x
		lastkey := x - 1
		for y := 1; y <= s1len; y++ {
			oldkey := column[y]
			var incr int
			if s1[y-1] != s2[x-1] {
				incr = 1
			}

			column[y] = minimum(column[y]+1, column[y-1]+1, lastkey+incr)
			lastkey = oldkey
		}
	}
	return column[s1len]
}

func minimum(a, b, c int) int {
	if a < b {
		if a < c {
			return a
		}
	} else {
		if b < c {
			return b
		}
	}
	return c
}
