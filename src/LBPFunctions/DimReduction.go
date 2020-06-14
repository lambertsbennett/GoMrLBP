package LBPFunctions

import (
	"fmt"
	"github.com/james-bowman/nlp"
	. "gonum.org/v1/gonum/mat"
	"os"
)

//Convert a SequenceCollection to a Dense Matrix for SVD

func SeqColToDense(sc SequenceCollection) *Dense {
	cols := len(sc.Items[0].Hist)
	rows := len(sc.Items)
	if rows == 0{
		fmt.Println("Cannot create matrix with no histogram bins present for sequences")
		os.Exit(1)
	}

	if cols == 0{
		fmt.Println("Cannot create a matrix from empty Sequence Collection")
		os.Exit(1)
	}

	seqMat := NewDense(rows,cols,nil)

	for i,seq := range sc.Items{
		seqMat.SetRow(i,seq.Hist)
	}

	return seqMat
}

func TruncatedSVD(data Matrix, k int) *Dense {
	svd := nlp.NewTruncatedSVD(k)
	mat, err := svd.FitTransform(data)
	if err != nil {
		fmt.Println("SVD factorization was not successful")
		os.Exit(69)
	}
	d := DenseCopyOf(mat.T())
	return d

}

func (sc *SequenceCollection) AddSVD(svd *Dense){
	for i,_ := range sc.Items{
		sc.Items[i].Svd = svd.RawRowView(i)
	}
}
