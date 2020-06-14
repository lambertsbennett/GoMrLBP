package main

import (
	"github.com/lambertsbennett/GoMrLBP/src/LBPFunctions"
	"fmt"
	"os"
	"runtime"
	"sync"
	"time"
	"flag"
)


func main() {

var proc int
flag.IntVar(&proc,"n",2,"Number of processors or threads to leverage.")

var contigfile string
flag.StringVar(&contigfile,"file","","Contig file in fasta format.")

var out string
flag.StringVar(&out,"o",".","Output file.")

var max_win int
flag.IntVar(&max_win,"max-win",9,"Largest window to use computing LBP codes. Must be odd.")

var single bool
flag.BoolVar(&single,"single-win", false, "Use a single window to calculate LBP codes.")

runtime.GOMAXPROCS(proc)


// Memory problem with >1,000,000 contigs if window size 9 included

ls := LBPFunctions.ReadFasta(contigfile)
lsp := LBPFunctions.SequenceCollection{}

fmt.Println("Processing sequences")
start := time.Now()
var wg sync.WaitGroup
in := make(chan LBPFunctions.Sequence, len(ls))
//out := make(chan Sequence,len(ls))

if max_win % 2 == 0{
	fmt.Println("Error: window sizes must be an odd number.")
	os.Exit(1)
}

if single == true{

	windows := make([]int,0)
	windows = append(windows, max_win)

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

	LBPFunctions.PrintMemUsage()

	wg.Wait()

	fmt.Println("Workers done.")

}else {

	windows := make([]int, 0)
	for i := 3; i <= max_win; i = i + 2 {
		windows = append(windows, i)
	}

	for i := 0; i < 100; i++ {
		wg.Add(1)
		go calcLBPHist(&wg, in, windows, &lsp)
	}

	fmt.Println("Spawned workers.")

	for _, s := range ls {
		in <- s
	}

	close(in)

	fmt.Println("Closed input channel.")

	LBPFunctions.PrintMemUsage()

	wg.Wait()

	fmt.Println("Workers done.")

}

t := time.Since(start)
fmt.Printf("%v sequences analysed in %s \n",len(ls),t)

d := LBPFunctions.SeqColToDense(lsp)
start = time.Now()
trunc := LBPFunctions.TruncatedSVD(d.T(),40)
t = time.Since(start)
fmt.Printf("SVD completed in %s \n", t)

lsp.AddSVD(trunc)

lsp.ToCSV(out)


}

func calcLBPHist(wg *sync.WaitGroup, in chan LBPFunctions.Sequence, windows []int,lsp *LBPFunctions.SequenceCollection) {
	for s:= range in {

		s.IntRep()
		for _, w := range windows {

			s.FindLBP(w)
			hist := s.NewHistogram(w)
			s.Hist = append(s.Hist, hist...)

		}
		rs := LBPFunctions.NewReducedSequence()
		rs.SpeciesID = s.SpeciesID
		rs.Header = s.Header
		rs.Hist = s.Hist
		lsp.Append(*rs)
	}
	wg.Done()
}

