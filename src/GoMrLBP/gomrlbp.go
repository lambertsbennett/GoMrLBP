package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/lambertsbennett/GoMrLBP/src/LBPFunctions"
	"log"
	"runtime"
	"sync"
	"time"
)


func main() {

var proc int
flag.IntVar(&proc,"n",2,"Number of processors or threads to leverage.")

var contigfile string
flag.StringVar(&contigfile,"file","","Contig file in fasta format.")

var out string
flag.StringVar(&out,"o","./gomrlbpout.parquet","Output file.")

var max_win int
flag.IntVar(&max_win,"max-win",9,"Largest window to use computing LBP codes. Must be odd.")

var single bool
flag.BoolVar(&single,"single-win", false, "Use a single window to calculate LBP codes.")

flag.Parse()

runtime.GOMAXPROCS(proc)

<<<<<<< HEAD

=======
>>>>>>> 363db2f0e4c4aec473ed5e9b030679911b9ac817
ls := LBPFunctions.ReadFasta(contigfile)
lsp := LBPFunctions.SequenceCollection{}

fmt.Println("Processing sequences")
start := time.Now()
var wg sync.WaitGroup
in := make(chan LBPFunctions.Sequence, len(ls))

_, err := check_windows(max_win)
if err != nil{log.Fatal(err)}

if single == true{
	fmt.Println("Executing in single mode ...")
	windows := make([]int,0)
	windows = append(windows, max_win)

	for i:=0; i< 100; i++ {
		wg.Add(1)
		go calcLBPHist(&wg,in,windows,&lsp)
	}

	fmt.Println("Spawned workers ...")

	for _, s := range ls {
		in <- s
	}

	close(in)

	fmt.Println("Closed input channel ...")

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

	fmt.Println("Spawned workers ...")

	for _, s := range ls {
		in <- s
	}

	close(in)

	fmt.Println("Closed input channel ...")

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

lsp.ToParquet(out)
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
		rs.Header = s.Header
		rs.Hist = s.Hist
		lsp.Append(*rs)
	}
	wg.Done()
}

func check_windows(max_win int) (int, error){
	if max_win % 2 == 0{
		return 0, errors.New("Maximum window size should be an odd integer.")
	} else{
		return 0, nil
	}
}
