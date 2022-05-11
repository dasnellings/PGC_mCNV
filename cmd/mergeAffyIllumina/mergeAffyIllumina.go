package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"strings"
)

func main() {
	affyFile := os.Args[0]
	outFile := os.Args[1]
	if affyFile == "" {
		log.Fatal("ERROR: arg 1 must be affy file")
	}
	if outFile == "" {
		log.Fatal("ERROR: arg 2 must be outfile")
	}
	out := fileio.EasyCreate(outFile)

	data, header := vcf.GoReadToChan(affyFile)
	vcf.NewWriteHeader(out, header)

	var v vcf.Vcf
	var queryResult [][]string
	var idx int
	var err error
	var found bool
	var s vcf.Sample
	for v = range data {
		queryResult, found = vcf.QueryString(v, header.Info["DBSNP_RS_ID"].Key)
		if !found {
			log.Println("WARNING:", err)
		} else {
			v.Id = queryResult[0][0]
		}
		v.Info = strings.Join(strings.Split(v.Info, ";")[:3], ";")
		if v.Format[0] != "GT" && v.Format[6] != "BAF" && v.Format[7] != "LRR" {
			log.Panicln("bad format:", v.String())
		}
		v.Format = []string{"GT", "BAF", "LRR"}
		for _, s := range v.Samples {
			s.FormatData = s.FormatData[5:]
			s.FormatData[0] = ""
		}
	}
}
