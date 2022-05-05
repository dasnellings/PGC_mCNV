package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/PGC_mCNV/illumina"
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
	ref := fasta.ReadToString(fastaFile)
	var header vcf.Header
	header.Text = strings.Split(headerInfo, "\n")
	header.Text[len(header.Text)-1] += "\t" + strings.Join(gsReportFiles, "\t")
	vcf.NewWriteHeader(out, header)

	gsReportChans := make([]<-chan illumina.GsReport, len(gsReportFiles))
	for i := range gsReportFiles {
		gsReportChans[i] = illumina.GoReadGsReportToChan(gsReportFiles[i])
	}

	manifestData := illumina.GoReadManifestToChan(manifestFile)

	var curr vcf.Vcf
	var gs illumina.GsReport
	curr.Filter = "."
	curr.Format = []string{"GT", "BAF", "LRR"}
	sb := new(strings.Builder)
	var alleleAint, alleleBint int16
	for m := range manifestData {
		curr.Chr = m.Chr
		curr.Pos = m.Pos
		curr.Id = m.Name
		curr.Ref = strings.ToUpper(string(ref[m.Chr][m.Pos]))

		// check one of the alleles matches ref
		switch curr.Ref {
		case m.TopStrandAlleleA:
			curr.Alt = []string{m.TopStrandAlleleB}
			alleleAint = 0
			alleleBint = 1
		case m.TopStrandAlleleB:
			curr.Alt = []string{m.TopStrandAlleleA}
			alleleAint = 1
			alleleBint = 0
		default:
			log.Panicf("ERROR: Reference mismatch. Check '%s'\n", m.Name)
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

	err := out.Close()
	exception.PanicOnErr(err)
}

func matchesManifest(gs illumina.GsReport, m illumina.Manifest) bool {
	if strings.ToUpper(gs.Marker) != strings.ToUpper(m.Name) {
		return false
	}
	if "chr"+strings.ToUpper(gs.Chrom) != m.Chr {
		return false
	}
	if gs.Pos != m.Pos {
		return false
	}
	return true
}

func revComp(base string) string {
	switch base {
	case "A", "a":
		return "T"
	case "C", "c":
		return "G"
	case "G", "g":
		return "C"
	case "T", "t":
		return "A"
	default:
		log.Panicf("ERROR: unrecognized base '%s'", base)
		return ""
	}
}