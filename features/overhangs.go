package main

import (
	"math/rand"

	"github.com/Open-Science-Global/poly/io/fasta"
)

func main() {
	enzymes := fasta.Read("../data/output/output.fasta")

	var fastas []fasta.Fasta
	for _, enzyme := range enzymes {
		updated := fasta.Fasta{enzyme.Name, createCdsPart(enzyme.Sequence)}
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
	randomFoward := randomDnaSequence(15, int64(rand.Intn(10000)))
	twoRandomFoward := randomDnaSequence(2, int64(rand.Intn(10000)))
	eightRandomFoward := randomDnaSequence(8, int64(rand.Intn(10000)))

	return randomFoward + bbsiFoward + twoRandomFoward + bbsiOverhangFoward + eightRandomFoward + internalOverhang

}

func addBbsiStructureReverse(internalOverhang string) string {
	bbsiReverse := "GTCTTC"
	bbsiOverhangReverse := "CGCT"
	randomReverse := randomDnaSequence(15, int64(rand.Intn(10000)))
	twoRandomReverse := randomDnaSequence(2, int64(rand.Intn(10000)))
	eightRandomReverse := randomDnaSequence(8, int64(rand.Intn(10000)))

	return internalOverhang + eightRandomReverse + bbsiOverhangReverse + twoRandomReverse + bbsiReverse + randomReverse

}

func createCdsPart(sequence string) string {
	spacer := "T"
	bsaiFoward := "GGTCTC" + spacer
	bsaiReverse := spacer + "GAGACC"
	fiveOverhang := bsaiFoward + "AATG"
	threeOverhang := "GCTT" + bsaiReverse
	return addBbsiStructureFoward(fiveOverhang) + sequence + addBbsiStructureReverse(threeOverhang)
}

func randomDnaSequence(length int, seed int64) string {
	var dnaAlphabet = []rune("ATCG")
	rand.Seed(seed)

	randomSequence := make([]rune, length)

	for basepair := range randomSequence {
		randomIndex := rand.Intn(len(dnaAlphabet))
		randomSequence[basepair] = dnaAlphabet[randomIndex]
	}
	return string(randomSequence)
}
