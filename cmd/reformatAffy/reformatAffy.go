package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"strings"
)

func main() {
	affyFile := os.Args[1]
	outFile := os.Args[2]
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
	var i int
	var err error
	var found bool
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
		for i = range v.Samples {
			v.Samples[i].FormatData = v.Samples[i].FormatData[5:]
			v.Samples[i].FormatData[0] = ""
		}
	}
}
