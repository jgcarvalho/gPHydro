package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"bitbucket.org/binet/go-gnuplot/pkg/gnuplot"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/gonum/stat"
)

type ScaleAA map[string]float64

func readFasta(fn string) (name string, seq string, err error) {
	fFasta, err := os.Open(fn)
	defer fFasta.Close()
	if err != nil {
		return "", "", err
	}
	t := linear.NewSeq("", nil, alphabet.Protein)
	reader := fasta.NewReader(fFasta, t)
	s, err := reader.Read()
	if err != nil {
		return "", "", err
	}
	sl := s.(*linear.Seq)
	return sl.Name(), sl.String(), nil
}

func readScaleAA(sc string, norm bool) (ScaleAA, error) {
	aa := []string{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
		"Q", "R", "S", "T", "V", "W", "Y"}
	aaScale := make(map[string]float64)
	var fn string
	inv := false
	switch sc {
	case "Argos":
		fn, _ = filepath.Abs("./scales/ARGP820101.dat")
	case "Black":
		fn, _ = filepath.Abs("./scales/BLAS910101.dat")
	case "Casari":
		fn, _ = filepath.Abs("./scales/CASG920101.dat")
	case "Cid":
		fn, _ = filepath.Abs("./scales/CIDH920105.dat")
	case "Eisenberg":
		fn, _ = filepath.Abs("./scales/EISD840101.dat")
	case "Engelman":
		fn, _ = filepath.Abs("./scales/ENGD860101.dat")
		inv = true
	case "Fasman":
		fn, _ = filepath.Abs("./scales/FASG890101.dat")
		inv = true
	case "Fauchere":
		fn, _ = filepath.Abs("./scales/FAUJ830101.dat")
	case "Goldsack":
		fn, _ = filepath.Abs("./scales/GOLD730101.dat")
	case "Hopp":
		fn, _ = filepath.Abs("./scales/HOPT810101.dat")
		inv = true
	case "Jones":
		fn, _ = filepath.Abs("./scales/JOND750101.dat")
	case "Kyte":
		fn, _ = filepath.Abs("./scales/KYTJ820101.dat")
	case "Levitt":
		fn, _ = filepath.Abs("./scales/LEVM760101.dat")
		inv = true
	case "Ponnuswamy":
		fn, _ = filepath.Abs("./scales/PONP930101.dat")
	case "Prabhakaran":
		fn, _ = filepath.Abs("./scales/PRAM900101.dat")
		inv = true
	case "Radzicka":
		fn, _ = filepath.Abs("./scales/RADA880108.dat")
	case "Rose":
		fn, _ = filepath.Abs("./scales/ROSG850102.dat")
	case "Wolfenden":
		fn, _ = filepath.Abs("./scales/WOLR790101.dat")
	case "Zimmerman":
		fn, _ = filepath.Abs("./scales/ZIMJ680101.dat")
	default:
		fn = ""
	}

	fScale, err := ioutil.ReadFile(fn)
	if err != nil {
		panic(err)
	}
	sScale := strings.Split(string(fScale), "\t")
	scale := make([]float64, len(sScale))
	for i := 0; i < len(scale); i++ {
		scale[i], err = strconv.ParseFloat(strings.TrimSpace(sScale[i]), 64)
		if err != nil {
			fmt.Println(err)
			return nil, err
		}
	}
	if norm {
		scale = rescale(scale, inv)
	}
	for i := 0; i < len(scale); i++ {
		aaScale[aa[i]] = scale[i]
	}
	return aaScale, err
}

func rescale(sc []float64, inv bool) []float64 {
	max := sc[0]
	min := sc[0]
	reScale := make([]float64, len(sc))
	for i := 1; i < len(sc); i++ {
		if sc[i] > max {
			max = sc[i]
		}
		if sc[i] < min {
			min = sc[i]
		}
	}
	if inv {
		for i := 0; i < len(sc); i++ {
			reScale[i] = 1.0 - ((sc[i] - min) / (max - min))
		}
	} else {
		for i := 0; i < len(sc); i++ {
			reScale[i] = (sc[i] - min) / (max - min)
		}
	}
	return reScale
}

func (sc ScaleAA) Apply(seq string, sw int) ([]float64, []float64) {
	// sw must be odd (test before)
	hydro := make([]float64, len(seq))
	hydrosw := make([]float64, len(seq))
	for i := 0; i < len(seq); i++ {
		hydro[i] = sc[string(seq[i])]
	}
	for i := 0; i < len(hydro); i++ {
		if (i >= (sw / 2)) && (i < (len(hydro) - (sw / 2))) {
			b, e := i-(sw/2), i+(sw/2)+1
			hydrosw[i] = stat.Mean(hydro[b:e], nil)
		} else if i < (sw / 2) {
			b, e := 0, i+(sw/2)+1
			hydrosw[i] = stat.Mean(hydro[b:e], nil)
		} else if i >= (len(hydro) - (sw / 2)) {
			b, e := i-(sw/2), len(hydro)
			hydrosw[i] = stat.Mean(hydro[b:e], nil)
		}
	}
	return hydro, hydrosw
}

func Plot(h []float64, hsw []float64) {
	fname := ""
	persist := true
	debug := true

	p, err := gnuplot.NewPlotter(fname, persist, debug)
	if err != nil {
		err_string := fmt.Sprintf("** err: %v\n", err)
		panic(err_string)
	}
	defer p.Close()
	p.SetXLabel("Residue number")
	p.SetYLabel("Hydrophobicity")
	p.CheckedCmd("set yrange [0:1]")
	p.CheckedCmd(fmt.Sprintf("set xrange [1:%d]", len(hsw)))
	// p.PlotX(h, "raw")
	p.SetStyle("lines")
	p.PlotX(hsw, "sw")

	p.CheckedCmd("set terminal pdf")
	p.CheckedCmd("set output 'plot002.pdf'")
	p.CheckedCmd("replot")

	p.CheckedCmd("q")
	return
}

func main() {
	fn := "/home/jgcarvalho/gocode/src/github.com/jgcarvalho/gPHydro/examples/2a1h.fasta"
	name, seq, err := readFasta(fn)
	if err != nil {
		panic(err)
	}
	fmt.Println(name)
	fmt.Println(seq)

	scale, err := readScaleAA("Rose", true)
	if err != nil {
		panic(err)
	}
	fmt.Println(scale)
	fmt.Println(len(scale))

	hydro, hydrosw := scale.Apply(seq, 7)
	fmt.Println(hydro)
	fmt.Println(hydrosw)
	Plot(hydro, hydrosw)
	fmt.Println(len(hydro))
}
