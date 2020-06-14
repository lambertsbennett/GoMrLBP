package LBPFunctions

import "sync"

type Sequence struct {
	Header string //contig id from fasta file
	Seq string // contig sequence
	NRep []float64 // numerical representation of sequence (Integer representation)
	LbpCodes []float64 //  lbp intermediate slice
	Hist []float64 // Un-normalized LBP histogram
	SpeciesID int // species ID. Useful for visualization purposes allows you to track contig fate
	Kmers string
}

type ReducedSequence struct {
	Header string //contig id from fasta file
	Hist []float64 // Un-normalized LBP histogram
	Svd []float64 // SVD reduced dimension vector
	SpeciesID int // species ID. Useful for visualization purposes allows you to track contig fate
	Kmers string
}

//http://dnaeon.github.io/concurrent-maps-and-slices-in-go/
type SequenceCollection struct {
	sync.RWMutex
	Items [] ReducedSequence
}

func (cs *SequenceCollection) Append(item ReducedSequence){
	cs.Lock()
	defer cs.Unlock()
	cs.Items = append(cs.Items, item)
}

// Create a new sequence.
// All slices are instantiated as nil slices and speciesID is default zero.
func NewSequence() *Sequence {
	s := Sequence{"", "", nil, nil, nil,0,""}
	return &s
}

func NewReducedSequence() *ReducedSequence {
	s := ReducedSequence{"", nil, nil,0,""}
	return &s
}

