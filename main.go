package main

import (
	"MLBP/SeqBase"
	"fmt"
	"os"
	"runtime"
	"sync"
	"time"
)


func main() {

	//filePtr := flag.String("file","","File containing contigs. Must be in Fasta format.")
	//threadsPtr := flag.Int("nthreads",2,"Set the number of threads to run the Mr-LBP algorithms. Default = 2")
	//windPtr := flag.Int("maxWind",9,"Set the maximum window size used to calculate LBP codes")
	//outPtr := flag.String("out", ".","Output file to write to.")
	//
	//flag.Parse()

// Processing functions. Input a list of sequences and get back the svd representation, LBP codes resulting from sequence contigs

runtime.GOMAXPROCS(8)


// Memory problem with >1,000,000 contigs if window size 9 included

ls := readFasta("/home/ben/Desktop/MrLBP/TestData/test.fa")
lsp := SequenceCollection{}

fmt.Println("Processing sequences")
start := time.Now()
var wg sync.WaitGroup
in := make(chan Sequence, len(ls))
//out := make(chan Sequence,len(ls))

if 7 % 2 == 0{
	fmt.Println("Error: window sizes must be an odd number.")
	os.Exit(1)
}
windows := make([]int,0)
for i:=3; i<=9; i=i+2{
	windows = append(windows,i)
}

	for i:=0; i< 100; i++ {
		wg.Add(1)
		go calcLBPHist(&wg,in,windows,&lsp)

	}

	fmt.Println("Spawned workers.")

	for _, s := range ls {
		in <- s
	}

	close(in)

	fmt.Println("Closed input channel.")

	PrintMemUsage()

	wg.Wait()

	fmt.Println("Workers done.")


t := time.Since(start)
fmt.Printf("%v sequences analysed in %s \n",len(ls),t)

d := SeqColToDense(lsp)
start = time.Now()
trunc := truncatedSVD(d.T(),40)
t = time.Since(start)
fmt.Printf("SVD completed in %s \n", t)

lsp.AddSVD(trunc)

lsp.toCSV("/home/ben/Desktop/MrLBP/TestData/test.csv")


}

func calcLBPHist(wg *sync.WaitGroup, in chan Sequence, windows []int,lsp *SequenceCollection) {
	for s:= range in {

		s.IntRep()
		for _, w := range windows {

			s.FindLBP(w)
			hist := s.NewHistogram(w)
			s.Hist = append(s.Hist, hist...)

		}
		rs := SeqBase.NewReducedSequence()
		rs.SpeciesID = s.SpeciesID
		rs.Header = s.Header
		rs.Hist = s.Hist
		lsp.Append(*rs)
	}
	wg.Done()
}

