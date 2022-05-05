package illumina

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
	"strings"
)

const expectedGsReportHeader string = "SNP Name\tChromosome\tPosition\tAl1Fwd\tAl2Fwd\tX\tY\tB Allele Freq\tLog R Ratio\n"

type GsReport struct {
	Marker      string
	Chrom       string
	Pos         int
	Allele1     string
	Allele2     string
	X           float64
	Y           float64
	BAlleleFreq float64
	LogRRatio   float64
}

func GoReadGsReportToChan(filename string) <-chan GsReport {
	ans := make(chan GsReport, 100)
	go readReportToChan(filename, ans)
	return ans
}

func readReportToChan(filename string, ans chan<- GsReport) {
	file := fileio.EasyOpen(filename)
	for line, done := fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, "SNP Name") {
			if line != expectedGsReportHeader {
				log.Fatalf("ERROR: unexpected report header. check file.\n%s\n%s", line, expectedGsReportHeader)
			}
			continue
		}
		ans <- processGsReportLine(line)
	}
	err := file.Close()
	exception.PanicOnErr(err)
	close(ans)
}

func processGsReportLine(s string) GsReport {
	var ans GsReport
	var err error
	fields := strings.Split(s, "\t")
	if len(fields) != 9 {
		log.Panicf("ERROR: following lines has unexpected number of columns:\n%s", s)
	}
	ans.Marker = fields[0]
	ans.Chrom = fields[1]
	ans.Pos, err = strconv.Atoi(fields[2])
	exception.PanicOnErr(err)
	ans.Allele1 = strings.ToUpper(fields[3])
	ans.Allele2 = strings.ToUpper(fields[4])
	ans.X, err = strconv.ParseFloat(fields[5], 64)
	exception.PanicOnErr(err)
	ans.Y, err = strconv.ParseFloat(fields[6], 64)
	exception.PanicOnErr(err)
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[7], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[8], 64)
	exception.PanicOnErr(err)
	return ans
}
