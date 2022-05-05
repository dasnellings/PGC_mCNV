package illumina

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
	"strings"
)

const expectedManifestHeader string = "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID," +
	"AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq," +
	"TopGenomicSeq,BeadSetId\n"

type Manifest struct {
	IlmnId           string
	Name             string
	Strand           bool
	AlleleA          string
	AlleleB          string
	TopStrandAlleleA string
	TopStrandAlleleB string
	GenomeBuild      string
	Chr              string
	Pos              int
	GC               float64
}

func GoReadManifestToChan(filename string) <-chan Manifest {
	ans := make(chan Manifest, 1000)
	go readManifestToChan(filename, ans)
	return ans
}

func readManifestToChan(filename string, ans chan<- Manifest) {
	file := fileio.EasyOpen(filename)
	var throughHeader bool
	for line, done := fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		if throughHeader {
			ans <- processManifestLine(line)
		}
		if !strings.HasPrefix(line, "IlmnId") {
			continue
		}
		if line != expectedManifestHeader {
			log.Fatalf("ERROR: unexpected manifest header. check file.\n%s\n%s", line, expectedManifestHeader)
		}
		throughHeader = true
	}
	err := file.Close()
	exception.PanicOnErr(err)
	close(ans)
}

func processManifestLine(s string) Manifest {
	var ans Manifest
	var err error
	fields := strings.Split(s, ",")
	if len(fields) != 19 {
		log.Panicf("ERROR: following lines has unexpected number of columns:\n%s", s)
	}
	ans.IlmnId = fields[0]
	ans.Name = fields[1]
	switch fields[2] {
	case "Top":
		ans.Strand = true
	case "Bot":
		ans.Strand = false
	default:
		log.Panicf("Unrecognized strand '%s' in below line\n%s\n", fields[2], s)
	}
	var alleles []string
	alleles = strings.Split(strings.TrimLeft(strings.TrimRight(fields[3], "]"), "["), "/")
	ans.AlleleA = alleles[0]
	ans.AlleleB = alleles[1]
	if ans.Strand {
		ans.TopStrandAlleleA = ans.AlleleA
		ans.TopStrandAlleleB = ans.AlleleB
	} else {
		ans.TopStrandAlleleA = revComp(ans.AlleleA)
		ans.TopStrandAlleleB = revComp(ans.AlleleB)
	}
	ans.GenomeBuild = fields[8]
	ans.Chr = fields[9]
	ans.Pos, err = strconv.Atoi(fields[10])
	exception.PanicOnErr(err)
	var seqContext string
	seqContext = fields[18][0:strings.Index(fields[18], "[")] + fields[18][strings.Index(fields[18], "]"):]
	seqContext = strings.ToUpper(seqContext)
	var gcCount int
	var totalCount int
	for _, base := range seqContext {
		switch base {
		case 'C', 'G':
			gcCount++
			totalCount++
		case 'T', 'A':
			totalCount++
		default:
			log.Panicf("ERROR: Unknown base in '%s'. Check following line\n%s\n", fields[18], s)
		}
	}
	ans.GC = float64(gcCount) / float64(totalCount)
	return ans
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