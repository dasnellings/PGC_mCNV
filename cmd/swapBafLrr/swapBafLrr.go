package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
)

func main() {
	file := os.Args[1]
	out := os.Args[2]
	o := fileio.EasyCreate(out)
	data, header := vcf.GoReadToChan(file)
	vcf.NewWriteHeader(o, header)
	var i int
	var v vcf.Vcf
	for v = range data {
		if v.Format[1] != "BAF" || v.Format[2] != "LRR" {
			log.Panicln("baf and lrr out of order or not found")
		}
		for i = range v.Samples {
			v.Samples[i].FormatData[1], v.Samples[i].FormatData[2] = v.Samples[i].FormatData[2], v.Samples[i].FormatData[1]
		}
		vcf.WriteVcf(o, v)
	}
	err := o.Close()
	exception.PanicOnErr(err)
}
