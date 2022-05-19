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
	"TopGenomicSeq,BeadSetId"
const expectedManifestHeader2 string = "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID," +
	"AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq," +
	"TopGenomicSeq,BeadSetID,Exp_Clusters,Intensity_Only,RefStrand"

type Manifest struct {
	IlmnId       string
	Name         string
	TopStrand    bool
	SrcTopStrand bool
	AlleleA      string
	AlleleB      string
	SeqBefore    string
	SeqAfter     string
	GenomeBuild  string
	Chr          string
	Pos          int
	GC           float64
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
		if strings.HasPrefix(line, "[Controls]") {
			break
		}
		if throughHeader {
			ans <- processManifestLine(line)
			continue
		}
		if !strings.HasPrefix(line, "IlmnID") {
			continue
		}
		if line != expectedManifestHeader && line != expectedManifestHeader2 {
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
	if len(fields) != len(strings.Split(expectedManifestHeader, ",")) && len(fields) != len(strings.Split(expectedManifestHeader2, ",")) {
		log.Panicf("ERROR: following lines has unexpected number of columns:\n%s", s)
	}
	ans.IlmnId = fields[0]
	ans.Name = fields[1]
	switch fields[2] {
	case "Top", "TOP":
		ans.TopStrand = true
	case "Bot", "BOT":
		ans.TopStrand = false
	default:
		log.Panicf("Unrecognized strand '%s' in below line\n%s\n", fields[2], s)
	}
	switch fields[15] {
	case "Top", "TOP":
		ans.SrcTopStrand = true
	case "Bot", "BOT":
		ans.SrcTopStrand = false
	default:
		log.Panicf("Unrecognized strand '%s' in below line\n%s\n", fields[2], s)
	}
	var alleles []string
	alleles = strings.Split(strings.TrimLeft(strings.TrimRight(fields[3], "]"), "["), "/")
	ans.AlleleA = alleles[0]
	ans.AlleleB = alleles[1]
	ans.GenomeBuild = fields[8]
	ans.Chr = fields[9]
	ans.Pos, err = strconv.Atoi(fields[10])
	exception.PanicOnErr(err)
	var seqContext string
	ans.SeqBefore = strings.ToUpper(strings.Split(fields[17], "[")[0])
	ans.SeqAfter = strings.ToUpper(strings.Split(fields[17], "]")[1])
	seqContext = ans.SeqBefore + ans.SeqAfter
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
			totalCount++
			//log.Panicf("ERROR: Unknown base in '%s'. Check following line\n%s\n", seqContext, s) // they throw Y's and shit in there
		}
	}
	ans.GC = float64(gcCount) / float64(totalCount)
	return ans
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
