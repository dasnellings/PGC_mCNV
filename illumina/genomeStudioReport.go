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
const gsHeader4 string = "SNP Name\tChr\tPosition\tAl1Fwd\tAl2Fwd\tX\tY\tB Allele Freq\tLog R Ratio"
const gsHeader5 string = "SNP Name\tChromosome\tPosition\tGC Score\tAllele1 - Top\tAllele2 - Top\tAllele1 - AB\tAllele2 - AB\tX\tY\tRaw X\tRaw Y\tR Illumina\tTheta Illumina\tB Allele Freq\tLog R Ratio"
const gsHeader6 string = "sample.id\tSNP\tchr\tpos\tA1.forward\tA2.forward\tX\tY\tB.Allele.Freq\tLogRRatio"

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
	var gs GsReport
	var processFunc func(string) GsReport
	for line, done := fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		line = strings.TrimRight(line, "\t") // remove trailing tab
		if strings.HasPrefix(line, "SNP Name") || strings.HasPrefix(line, "sample.id") {
			switch line {
			case gsHeader1:
				processFunc = processGsHeader1
			case gsHeader2:
				processFunc = processGsHeader2
			case gsHeader3:
				processFunc = processGsHeader3
			case gsHeader4:
				processFunc = processGsHeader1
			case gsHeader5:
				processFunc = processGsHeader5
			case gsHeader6:
				processFunc = processGsHeader6
			default:
				log.Fatalf("ERROR: unexpected report header. check file.\n%v\n%v", line)
			}
			continue
		}
		gs = processFunc(line)
		if strings.ToLower(gs.Chrom) == "mt" {
			gs.Chrom = "M"
		}
		ans <- gs
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
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[8], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[7], 64)
	exception.PanicOnErr(err)
	ans.ReportedAsFwd = false
	return ans
}

func processGsHeader5(s string) GsReport {
	var ans GsReport
	var err error
	fields := strings.Split(s, "\t")
	if len(fields) != 16 {
		log.Panicf("ERROR: following lines has unexpected number of columns:\n%s", s)
	}
	ans.Marker = fields[0]
	ans.Chrom = fields[1]
	ans.Pos, err = strconv.Atoi(fields[2])
	exception.PanicOnErr(err)
	ans.Allele1 = strings.ToUpper(fields[4])
	ans.Allele2 = strings.ToUpper(fields[5])
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[14], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[15], 64)
	exception.PanicOnErr(err)
	ans.ReportedAsFwd = false
	return ans
}

func processGsHeader6(s string) GsReport {
	var ans GsReport
	var err error
	fields := strings.Split(s, "\t")
	if len(fields) != 10 {
		log.Panicf("ERROR: following lines has unexpected number of columns:\n%s", s)
	}
	ans.Marker = fields[1]
	ans.Chrom = fields[2]
	ans.Pos, err = strconv.Atoi(fields[3])
	exception.PanicOnErr(err)
	ans.Allele1 = strings.ToUpper(fields[4])
	ans.Allele2 = strings.ToUpper(fields[5])
	if ans.Allele1 == "-" || ans.Allele2 == "-" || strings.Contains(s, "\tNA\t") {
		return ans
	}
	ans.BAlleleFreq, err = strconv.ParseFloat(fields[8], 64)
	exception.PanicOnErr(err)
	ans.LogRRatio, err = strconv.ParseFloat(fields[9], 64)
	exception.PanicOnErr(err)
	ans.ReportedAsFwd = true
	return ans
}
