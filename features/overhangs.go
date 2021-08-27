package main

import (
	"fmt"
	"math/rand"
	"time"

	"github.com/Open-Science-Global/poly/finder"
	"github.com/Open-Science-Global/poly/io/fasta"
)

type Position struct {
	Start int
	End   int
}

func main() {
	enzymes := fasta.Read("../data/output/output_manually.fasta")
	rand.Seed(time.Now().UnixNano())
	var fastas []fasta.Fasta
	for _, enzyme := range enzymes {
		updated := fasta.Fasta{enzyme.Name, createCdsRemoveProblems(enzyme.Sequence)}
		fastas = append(fastas, updated)
	}

	fasta.Write(fastas, "../data/output/outputWithOverhangs.fasta")
}

// 15 random bp -> bbsi cut site forward GAAGAC -> 2bp -> bbsi overhang GGAG -> random 8bp ->
// bsai site forward -> bsai overhang -> main sequence -> bsai overhang 2 -> bsai site reverse -> random 8bp ->
// bbsi overhang CGCT -> 2bp -> bbsi cute site reverse GTCTTC -> 15 random bp
func addBbsiStructureFoward(internalOverhang string) string {
	bbsiFoward := "GAAGAC"
	bbsiOverhangFoward := "GGAG"

	randomFoward := createRandomDnaSequenceRemoveForbidden(15)
	twoRandomFoward := createRandomDnaSequenceRemoveForbidden(2)
	eightRandomFoward := createRandomDnaSequenceRemoveForbidden(8)

	return randomFoward + bbsiFoward + twoRandomFoward + bbsiOverhangFoward + eightRandomFoward + internalOverhang

}

func createCdsRemoveProblems(sequence string) string {

	var functions []func(string) []finder.Match
	functions = append(functions, finder.ForbiddenSequence(restrictionBindingSitesListOverhangs()))

	originalProblems := finder.Find(sequence, functions)
	fmt.Println("Original problems:", originalProblems)
	check := true
	var cds string
	for check {
		cds = createCdsPart(sequence)
		var functions []func(string) []finder.Match
		functions = append(functions, finder.ForbiddenSequence(restrictionBindingSitesListOverhangs()))

		problems := finder.Find(cds, functions)
		if len(problems) < len(originalProblems)+8 {
			break
		}
	}

	problems := finder.Find(cds, functions)
	fmt.Println("New problems:", problems)
	return cds
}

func createRandomDnaSequenceRemoveForbidden(size int) string {

	check := true

	var randomSequence string
	for check {
		randomSequence = randomDnaSequence(size)
		var functions []func(string) []finder.Match
		functions = append(functions, finder.ForbiddenSequence(restrictionBindingSitesListOverhangs()))

		problems := finder.Find(randomSequence, functions)

		if len(problems) == 0 {
			break
		}
	}
	return randomSequence
}

func addBbsiStructureReverse(internalOverhang string) string {
	bbsiReverse := "GTCTTC"
	bbsiOverhangReverse := "CGCT"
	randomReverse := createRandomDnaSequenceRemoveForbidden(15)
	twoRandomReverse := createRandomDnaSequenceRemoveForbidden(2)
	eightRandomReverse := createRandomDnaSequenceRemoveForbidden(8)

	return internalOverhang + eightRandomReverse + bbsiOverhangReverse + twoRandomReverse + bbsiReverse + randomReverse

}

func createCdsPart(sequence string) string {
	spacer := "T"
	bsaiFoward := "GGTCTC" + spacer
	bsaiReverse := spacer + "GAGACC"
	fiveOverhang := bsaiFoward + "A"
	threeOverhang := "GCTT" + bsaiReverse

	return addBbsiStructureFoward(fiveOverhang) + sequence + addBbsiStructureReverse(threeOverhang)
}

func randomDnaSequence(length int) string {
	var dnaAlphabet = []rune("ATCG")

	randomSequence := make([]rune, length)

	for basepair := range randomSequence {
		randomIndex := rand.Intn(len(dnaAlphabet))
		randomSequence[basepair] = dnaAlphabet[randomIndex]
	}
	return string(randomSequence)
}

func restrictionBindingSitesListOverhangs() []string {
	BsaI_bind_5prime := "GGTCTC" //5" GGTCTC N|      3"
	//3" CCAGAG N NNNN| 5"

	BbsI_bind_5prime := "GAAGAC" //5" GAAGAC NN|
	//3" CTTCTG NN NNNN|

	BtgzI_bind_5prime := "GCGATG" //5" GCGATG NNN NNN NNN N|      3"
	//3" CGCTAC NNN NNN NNN N NNNN| 5"

	SapI_bind_5prime := "GCTCTTC" //5" GCTCTTC N|     3"
	//3" CGAGAAG N NNN| 5"

	BsmbI_bind_5prime := "CGTCTC" //5" CGTCTC N|      3"
	//3" GCAGAG N NNNN| 5"

	AarI_bind_5prime := "CACCTGC" //5" CACCTGC N NNN|       3"
	//3" GTGGACG N NNN NNNN|  5"

	PmeI_bind := "GTTTAAAC"

	HindIII := "AAGCTT"
	PstI := "CTGCAG"
	XbaI := "TCTAGA"
	BamHI := "GGATCC"
	SmaI := "CCCGGG"
	KpnI := "GGTACC"
	SacI := "GAGCTC"
	SalI := "GTCGAC"
	EcoRI := "GAATTC"
	SphI := "GCATGC"
	AvrII := "CCTAGG"
	SwaI := "ATTTAAAT"
	AscI := "GGCGCGCC"
	FseI := "GGCCGGCC"
	PacI := "TTAATTAA"
	SpeI := "ACTAGT"
	NotI := "GCGGCCGC"
	SanDI_A := "GGGACCC"
	SanDI_T := "GGGTCCC"
	BglII := "AGATCT"
	XhoI := "CTCGAG"
	ClaI := "ATCGAT"

	forbiddenSeqList := []string{
		BsaI_bind_5prime,
		BbsI_bind_5prime,
		SapI_bind_5prime,
		BsmbI_bind_5prime,
		BtgzI_bind_5prime,
		AarI_bind_5prime,
		PmeI_bind,
		HindIII,
		PstI,
		XbaI,
		BamHI,
		SmaI,
		KpnI,
		SacI,
		SalI,
		EcoRI,
		SphI,
		AvrII,
		SacI,
		SalI,
		SwaI,
		AscI,
		FseI,
		PacI,
		SpeI,
		NotI,
		SanDI_A,
		SanDI_T,
		BglII,
		XhoI,
		ClaI,
		"AAAAAA",
		"CCCCCC",
		"TTTTTT",
		"GGGGGG",
	}

	return forbiddenSeqList
}
