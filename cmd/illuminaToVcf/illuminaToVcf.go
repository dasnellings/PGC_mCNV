package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/PGC_mCNV/illumina"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
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
	var seqBefore, seqAfter []dna.Base
	var stringBefore, stringAfter string
	var refBase []dna.Base
	var altNeedsRevComp bool

	for m := range manifestData {
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
		case percentMatch(stringBefore, m.SeqBefore, 0.9) &&
			percentMatch(stringAfter, m.SeqAfter, 0.9):
			if !m.TopStrand {
				altNeedsRevComp = true
			}

			// only do partial check on rev comps since if snp is not directl in middle of probe then before/after lengths differ
		case percentMatch(revComp(stringBefore)[:20], m.SeqAfter[:20], 0.9) &&
			percentMatch(revComp(stringAfter)[len(stringAfter)-20:], m.SeqBefore[len(m.SeqBefore)-20:], 0.9):
			if m.TopStrand {
				altNeedsRevComp = true
			}

		default:
			fmt.Println(stringBefore)
			fmt.Println(m.SeqBefore)
			fmt.Println()
			fmt.Println(revComp(stringBefore)[:20])
			fmt.Println(m.SeqAfter[:20])
			fmt.Println()
			fmt.Println(stringAfter)
			fmt.Println(m.SeqAfter)
			fmt.Println()
			fmt.Println(revComp(stringAfter)[len(stringAfter)-20:])
			fmt.Println(m.SeqBefore[len(m.SeqBefore)-20:])
			fmt.Println()
			log.Panicf("ERROR: Context sequences did not match reference:\n%s+%s\n%s+%s\n", stringBefore, stringAfter, m.SeqBefore, m.SeqAfter)
		}

		if altNeedsRevComp {
			m.AlleleA = revComp(m.AlleleA)
			m.AlleleB = revComp(m.AlleleB)
		}

		switch curr.Ref {
		case m.AlleleA:
			alleleAint = 0
			alleleBint = 1
			curr.Alt = []string{m.AlleleB}
		case m.AlleleB:
			alleleAint = 1
			alleleBint = 0
			curr.Alt = []string{m.AlleleA}
		default:
			log.Panicf("ERROR: alternate alleles did not match reference\n%v\n", m)
		}

		curr.Info = fmt.Sprintf("ALLELE_A=%d;ALLELE_B=%d;GC=%g", alleleAint, alleleBint, m.GC)
		curr.Samples = make([]vcf.Sample, len(gsReportChans))

		sb.Reset()
		for i := range curr.Samples {
			gs = <-gsReportChans[i]
			if !matchesManifest(gs, m) {
				log.Panicf("ERROR: Manifest mismatch. See report and manifest data below\n%v\n%v\n", gs, m)
			}
			curr.Samples[i].Phase = make([]bool, 1) // leave as false for unphased
			curr.Samples[i].FormatData = []string{fmt.Sprintf("%g", gs.BAlleleFreq), fmt.Sprintf("%g", gs.LogRRatio)}
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
			log.Panicf("ERROR: unrecognized base '%s'", base)
			return ""
		}
		j++
	}
	return string(ans)
}

func percentMatch(a, b string, percent float64) bool {
	if len(a) != len(b) {
		return false
	}
	var mismatch int
	for i := range a {
		if a[i] != b[i] {
			mismatch++
		}
	}
	if float64(mismatch)/float64(len(a)) >= percent {
		return true
	}
	return false
}
