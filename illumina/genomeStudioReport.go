package illumina

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
	"strings"
)

const gsHeader1 string = "SNP Name\tChromosome\tPosition\tAl1Fwd\tAl2Fwd\tX\tY\tB Allele Freq\tLog R Ratio"
const gsHeader2 string = "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\tX Raw\tY Raw\tX\tY\tTheta\tB Allele Freq\tLog R Ratio\tR\tAllele1 - Top\tAllele2 - Top\tGC Score\tGT Score\tCluster Sep"
const gsHeader3 string = "SNP Name\tChr\tPosition\tAllele1 - Top\tAllele2 - Top\tX\tY\tLog R Ratio\tB Allele Freq"

type GsReport struct {
	Marker  string
	Chrom   string
	Pos     int
	Allele1 string
	Allele2 string
	//X           float64
	//Y           float64
	BAlleleFreq   float64
	LogRRatio     float64
	ReportedAsFwd bool
}

func GoReadGsReportToChan(filename string) <-chan GsReport {
	ans := make(chan GsReport, 100)
	go readReportToChan(filename, ans)
	return ans
}

func readReportToChan(filename string, ans chan<- GsReport) {
	file := fileio.EasyOpen(filename)
	var processFunc func(string) GsReport
	for line, done := fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		line = strings.TrimRight(line, "\t") // remove trailing tab
		if strings.HasPrefix(line, "SNP Name") {
			switch line {
			case gsHeader1:
				processFunc = processGsHeader1
			case gsHeader2:
				processFunc = processGsHeader2
			case gsHeader3:
				processFunc = processGsHeader3
			default:
				log.Fatalf("ERROR: unexpected report header. check file.\n%v\n%v", line)
			}
			continue
		}
		ans <- processFunc(line)
	}
	err := file.Close()
	exception.PanicOnErr(err)
	close(ans)
}

func processGsHeader1(s string) GsReport {
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
	//ans.X, err = strconv.ParseFloat(fields[5], 64)
	//exception.PanicOnErr(err)
	//ans.Y, err = strconv.ParseFloat(fields[6], 64)
	//exception.PanicOnErr(err)
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[7], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[8], 64)
	exception.PanicOnErr(err)
	ans.ReportedAsFwd = true
	return ans
}

func processGsHeader2(s string) GsReport {
	var ans GsReport
	var err error
	fields := strings.Split(s, "\t")
	if len(fields) != 18 {
		log.Panicf("ERROR: following lines has unexpected number of columns:\n%s", s)
	}
	ans.Marker = fields[0]
	ans.Chrom = fields[1]
	ans.Pos, err = strconv.Atoi(fields[2])
	exception.PanicOnErr(err)
	ans.Allele1 = strings.ToUpper(fields[3])
	ans.Allele2 = strings.ToUpper(fields[4])
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[10], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[11], 64)
	exception.PanicOnErr(err)
	ans.ReportedAsFwd = true
	return ans
}

func processGsHeader3(s string) GsReport {
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
	//ans.X, err = strconv.ParseFloat(fields[5], 64)
	//exception.PanicOnErr(err)
	//ans.Y, err = strconv.ParseFloat(fields[6], 64)
	//exception.PanicOnErr(err)
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[7], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[8], 64)
	exception.PanicOnErr(err)
	ans.ReportedAsFwd = false
	return ans
}
