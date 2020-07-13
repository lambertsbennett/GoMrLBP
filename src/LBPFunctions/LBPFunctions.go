package LBPFunctions
/*
Bennett Lambert: Multi-resolution local binary patterns go implementation

The general workflow and a lot code was converted from github.com/skouchaki/MrGBP.
See their associated paper: doi.org/10.1038/s41598-018-38197-9 (Kouchaki et al. Sci. Reps. 2019)

In general the workflow proceeds as follows.
1) Read in a FASTA file of sequences. The readFasta function parses files only in FASTA format. If you have FASTQ convert first with seqtk or sed if you're too cool for school.
Files can be gzipped or uncompressed.

2) Iterate through sequences and calculate LBP histograms.
This is done through a few methods:
IntRep() finds the integer representation of a sequence (following Kouchaki et al.)
FindLBP(wind) calculates LBP codes over a sequence segment size of 'wind'
NewHistogram(wind) calculates the LBP histogram corresponding to the LBP codes found above. This returns a histogram and does not modify the Sequence struct in place.
 This allows you to loop over segment sizes for multi-resolution LBP appending the different histograms together.
 */

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"github.com/xitongsys/parquet-go-source/local"
	"github.com/xitongsys/parquet-go/parquet"
	"github.com/xitongsys/parquet-go/writer"
	"gonum.org/v1/gonum/stat/combin"
	"log"
	. "math"
	"os"
	"strconv"
	"strings"
)

// Read in a fasta file saving each sequence and header to a unique Sequence type. Store sequences in a slice
// https://www.socketloop.com/tutorials/golang-bufio-newreader-readline-to-read-file-line-by-line
func ReadFasta(fname string) []Sequence {
	file, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}

	bReader := bufio.NewReader(file)
	testBytes, err := bReader.Peek(2)
	if err != nil {
		log.Fatal(err)
	}

	file.Close()

	if testBytes[0] == 31 && testBytes[1] == 139 {
		file, err := os.Open(fname)
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()
		gzipReader, err := gzip.NewReader(file)
		if err != nil{
			log.Fatal(err)
		}

		defer gzipReader.Close()
		scanner := bufio.NewScanner(gzipReader)
		sequenceList := make([]Sequence, 0)
		seq := new(Sequence)
		for scanner.Scan() {
			lstring := string(scanner.Text())
			if strings.Contains(lstring, ">") {
				if lstring == "" {
					fmt.Print("Empty line encountered!")
					break
				}
				seq.Header = lstring
			} else {
				if lstring == "" {
					fmt.Print("Empty line encountered!")
					break
				}
				seq.Seq = lstring
				sequenceList = append(sequenceList, *seq)
			}

		}

		if err := scanner.Err(); err != nil{
			log.Fatal(err)
		}
		return sequenceList


	}else {
		file, err := os.Open(fname)
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()
		reader := bufio.NewReader(file)
		scanner := bufio.NewScanner(reader)

		sequenceList := make([]Sequence, 0)
		seq := new(Sequence)
		for scanner.Scan() {

			lstring := string(scanner.Text())
			if strings.Contains(lstring, ">") {
				if lstring == "" {
					fmt.Print("Empty line encountered!")
					break
				}
				seq.Header = lstring
			} else {
				if lstring == "" {
					fmt.Print("Empty line encountered!")
					break
				}
				seq.Seq = lstring
				sequenceList = append(sequenceList, *seq)
			}

		}
		if err := scanner.Err(); err != nil{
			log.Fatal(err)
		}
		return sequenceList
	}
	return nil
}

func (sc *SequenceCollection) ToCSV(fname string){
	file,err := os.Create(fname)
	if err != nil {log.Fatal(err)}
	w := csv.NewWriter(file)
	defer file.Close()
	for _,s := range sc.Items {
		tmp := make([]string,0)
		tmp = append(tmp, s.Header)
		tmp = append(tmp," ")
		for _,f := range s.Hist{
			tmp = append(tmp,strconv.FormatFloat(f,'f',-1,64))
		}
		tmp = append(tmp," ")
		for _,f := range s.Svd{
			tmp = append(tmp,strconv.FormatFloat(f,'f',-1,64))
		}

		w.Write(tmp)

	}
	w.Flush()
}


//Write results to parquet file. This is the recommended output format.
func (sc *SequenceCollection) ToParquet (fname string){
	type tmpseq struct {
		Header    string  `parquet:"name=name, type=UTF8, encoding=PLAIN_DICTIONARY"`
		Hist     []float64  `parquet:"name=HIST, type=DOUBLE, repetitiontype=REPEATED"`
		SVD     []float64   `parquet:"name=SVD, type=DOUBLE, repetitiontype=REPEATED"`
	}

	fw, err := local.NewLocalFileWriter(fname)
	if err != nil {
		log.Println("Can't open file", err)
		return
	}
	pw, err := writer.NewParquetWriter(fw, new(tmpseq),4)
	if err != nil {
		log.Println("Can't create parquet writer", err)
		return
	}
	pw.RowGroupSize = 5 * 1024 * 1024 //5M
	pw.CompressionType = parquet.CompressionCodec_GZIP

	for _,s := range sc.Items {
		seq := tmpseq{
			Header: s.Header,
			Hist:   s.Hist,
			SVD:    s.Svd,
		}
		if err = pw.Write(seq); err != nil {
			log.Println("Write error", err)
		}
	}
	if err = pw.WriteStop(); err != nil {
		log.Println("WriteStop error", err)
		return
	}
	log.Println("Write Finished")
	fw.Close()

}


// Create integer representation for a sequence
func (seq *Sequence) IntRep() {

	for _, rne := range seq.Seq {

		switch rne {
		case 'A':
			seq.NRep = append(seq.NRep, 2.0)
		case 'T':
			seq.NRep = append(seq.NRep, -2.0)
		case 'G':
			seq.NRep = append(seq.NRep, 1.0)
		case 'C':
			seq.NRep = append(seq.NRep, -1.0)

		}
	}
}

//Threshold using signbit function and sliding window of size wind. Wind should be odd and > 0.
func (seq *Sequence) FindLBP (wind int) {
	for j := 0; j < (len(seq.NRep) - wind)
	{
		tmp := seq.NRep[j:j+wind]
		seq.LbpCodes = append(seq.LbpCodes,calcLBP(wind,tmp))
		j++
	}
}

// Adapted from github.com/MrGBP/blob/master/oned_lbp.h
// Need to run through a hand calculation to double-check behaviour
func calcLBP(wind int, tmp []float64) float64{
	s:=0.0
	cen := int(Ceil(float64(wind/2)))
	for i := 0; i<cen
	{
		s=s+(1.0-checkBitSetVar(Signbit(tmp[i]-tmp[cen])))*Pow(2,float64(i))+(1.0-checkBitSetVar(Signbit(tmp[i+cen+1]-tmp[cen])))*Pow(2,float64(i+cen))
		i++
	}
	return s
}

// Convenient way to convert bool to float64 needed for math ops above.
// From https://stackoverflow.com/questions/38627078/how-to-convert-bool-to-int8-in-golang
func checkBitSetVar(mybool bool) float64{
	if mybool{
		return 1.0
	}
	return 0.0
}


// Calculate the LBP histogram
func (seq *Sequence) NewHistogram(wind int) []float64 {
	bins := genBins(wind)
	hist := make([]float64, len(bins))
	for i,val := range bins{
		for _,lc := range seq.LbpCodes{
			if lc == val{
				hist[i] += 1
			}
		}
	}
	return hist
}

func genBins(wind int) []float64{
	nums := genNums(wind)
	bins := make([]float64,0)
	for i := 1; i <wind
	{
		c := combin.Combinations(wind-1,i)
		for _, arr := range c {
			bins = append(bins,sumCombos(arr,nums))
		}
		i++
	}
	return bins
}

func genNums(wind int) []float64{
	els := make([]float64,0)
	for i :=0; i<wind-1
	{
		els = append(els,Pow(2,float64(i)))
		i++
	}
	return els
}

func sumCombos(sub []int, nums []float64) float64{
	b:=0.0
	for _, s := range sub {
		b += nums[s]
	}
	return b
}


