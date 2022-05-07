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
	"log"
	"strings"
)

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
	flag.Parse()

	if *gsReportFilename == "" || *manifestFilename == "" || *fastaFilename == "" {
		usage()
		log.Fatal("ERROR: GenomeStudio report, manifest, and reference fasta files are required (-gsReport, -manifest, -ref)")
	}

	illuminaToVcf(strings.Split(*gsReportFilename, ","), *manifestFilename, *fastaFilename, *output)
}

func illuminaToVcf(gsReportFiles []string, manifestFile, fastaFile, output string) {
	out := fileio.EasyCreate(output)
	ref := fasta.NewSeeker(fastaFile, fastaFile+".fai")
	var header vcf.Header
	header.Text = strings.Split(headerInfo, "\n")
	header.Text[len(header.Text)-1] += "\t" + strings.Join(gsReportFiles, "\t")
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
	var alleleA, alleleB string //, gsAllele1, gsAllele2 string
	var seqBefore, seqAfter []dna.Base
	var stringBefore, stringAfter string
	var refBase []dna.Base
	var altNeedsRevComp bool

	for m := range manifestData {
		if m.Chr == "XY" { // SERIOUSLY ILLUMINA... SERIOUSLY
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
		exception.PanicOnErr(err)
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

		curr.Info = fmt.Sprintf("ALLELE_A=%d;ALLELE_B=%d;GC=%g", alleleAint, alleleBint, m.GC)
		curr.Samples = make([]vcf.Sample, len(gsReportChans))

		sb.Reset()
		for i := range curr.Samples {
			gs = <-gsReportChans[i]
			if gs.Chrom == "xy" {
				gs.Chrom = "x"
			}
			if !matchesManifest(gs, m) {
				log.Panicf("ERROR: Manifest mismatch. See report and manifest data below\n%v\n%v\n", gs, m)
			}
			fmt.Println()
			fmt.Println(m.Name)
			fmt.Println("manifest:", m.AlleleA, m.AlleleB, m.TopStrand)
			fmt.Println("report  :", gs.Allele1, gs.Allele2)
			curr.Samples[i].FormatData = []string{"", fmt.Sprintf("%g", gs.BAlleleFreq), fmt.Sprintf("%g", gs.LogRRatio)}
			switch gs.Allele1 {
			case m.AlleleA:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleAint)
			case m.AlleleB:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleBint)
			}
			switch gs.Allele2 {
			case m.AlleleA:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleAint)
			case m.AlleleB:
				curr.Samples[i].Alleles = append(curr.Samples[i].Alleles, alleleBint)
			}
			curr.Samples[i].Phase = make([]bool, len(curr.Samples[i].Alleles)) // leave as false for unphased
		}
		vcf.WriteVcf(out, curr)
	}

	err = out.Close()
	exception.PanicOnErr(err)
	err = ref.Close()
	exception.PanicOnErr(err)
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
