package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/Open-Science-Global/poly/io/fasta"
	"github.com/Open-Science-Global/poly/io/genbank"
	"github.com/Open-Science-Global/poly/synthesis"
	"github.com/Open-Science-Global/poly/transform/codon"
)

func main() {
	// The goal is to take the list of enzymes and make CDS optimization using three different strategies:
	// 1. Create a codon table for Bacillus Subtilis strain. KO7 and generate a optimized CDS
	// 2. Create a codon table with starvation highly expressed genes in Bacillus Subtilis strain PY79 and optimize
	// 3. Codon optimize for two species; e.g Bacillus Subtilis KO7 and E. coli K12

	// To create a codon table we need a list of CDSs from the target organism and poly will take care of the rest for us
	// Taking a list of fasta files to generate the respective codon tables
	files := []string{"data/bsub-ko7-cdss.fasta", "data/bsub-py79-cdss-starvation.fasta", "data/bsub-py79-cdss.fasta", "data/ecoli-k12-cdss.fasta"}

	codonTables := make(map[string]codon.Table)

	for _, file := range files {
		fileName := filepath.Base(file)[0 : len(file)-11]

		fmt.Printf("Reading file %s...\n", fileName)
		cdsSequences := fasta.Read(file)

		// Create a single big string with all the CDSs
		var allCdssFromFile strings.Builder
		for _, cds := range cdsSequences {
			allCdssFromFile.WriteString(cds.Sequence)
		}
		codingRegions := allCdssFromFile.String()

		fmt.Printf("Creating table for %s...\n", fileName)
		codonTable := codon.GetCodonTable(11)

		fmt.Printf("Optimizing table for %s...\n", fileName)
		optimizationTable := codonTable.OptimizeTable(codingRegions)
		codonTables[fileName] = optimizationTable

		codon.WriteCodonJSON(optimizationTable, "data/codon-table/"+fileName+".json")
		fmt.Printf("\n")
	}

	//But before starting let's make a intersection between two codon tables, one from PY79 and another from K12
	fmt.Println("Creating a codon table for both species: Bacillus Subtilis PY79 and E. coli K12...")
	fmt.Printf("\n")
	bsubCodonTable := codonTables["bsub-py79-cdss"]
	ecoliCodonTable := codonTables["ecoli-k12-cdss"]

	bothSpeciesTable, _ := codon.CompromiseCodonTable(bsubCodonTable, ecoliCodonTable, 0.1)

	codon.WriteCodonJSON(bothSpeciesTable, "data/codon-table/bsub-ecoli.json")
	codonTables["bsub-ecoli"] = bothSpeciesTable

	starvationCodonTable := codonTables["bsub-py79-cdss-starvation"]
	specialBothSpeciesTable, _ := codon.CompromiseCodonTable(starvationCodonTable, ecoliCodonTable, 0.1)

	codon.WriteCodonJSON(specialBothSpeciesTable, "data/codon-table/starvation-ecoli.json")
	codonTables["starvation-ecoli"] = specialBothSpeciesTable

	fmt.Println("Tables created and optimized! You could find each one as json files inside data/codon-table folder.")
	fmt.Printf("\n")

	fmt.Println("Creating a Kmer Table from Host Genome...")
	fmt.Printf("\n")
	// Another important think that we need is a 20-mer Table of the Host Genome
	bsubGenomeFile := genbank.Read("data/bsub-py79-genome.gb")
	bsubGenome := bsubGenomeFile.Sequence
	hostGenomeKmerTable := getKmerTable(20, bsubGenome)

	fmt.Println("Host Kmer Table created!")
	fmt.Printf("\n")

	// Taking the list of enzymes to codon optimize for each STRATEGY and eliminate some problems
	enzymes := fasta.Read("data/enzymes.fasta")

	var output []fasta.Fasta
	for _, enzyme := range enzymes {
		fmt.Println("Codon Optimizing using Strategies 1, 2 and 3...")
		// Strategy #1:  Create a codon table for Bacillus Subtilis strain. KO7
		bsubKO7CodonTable := codonTables["bsub-ko7-cdss"]
		enzymeOptimizedBsub := CodonOptimization(enzyme.Sequence, bsubKO7CodonTable)

		// Strategy #2:  Create a codon table with starvation highly expressed genes
		bsubStarvationCodonTable := codonTables["bsub-py79-cdss-starvation"]
		enzymeOptimizedStarvation := CodonOptimization(enzyme.Sequence, bsubStarvationCodonTable)

		// Strategy #3:  Codon optimize for two species; e.g Bacillus Subtilis KO7 and E. coli K12
		bsubEcoliCodonTable := codonTables["bsub-ecoli"]
		enzymeOptimizedBsubEcoli := CodonOptimization(enzyme.Sequence, bsubEcoliCodonTable)

		// Strategy #3:  Codon optimize for two species; e.g Bacillus Subtilis KO7 and E. coli K12
		starvationEcoliCodonTable := codonTables["starvation-ecoli"]
		enzymeOptimizedStarvationEcoli := CodonOptimization(enzyme.Sequence, starvationEcoliCodonTable)

		fmt.Println("Fixing optimized sequences if it has any problems...")

		strategyOne := fixSequence(enzymeOptimizedBsub, bsubKO7CodonTable, hostGenomeKmerTable)
		//copyStrategyOne := TwoCdsWithoutRepetition(strategyOne, bsubKO7CodonTable)
		strategyTwo := fixSequence(enzymeOptimizedStarvation, bsubStarvationCodonTable, hostGenomeKmerTable)
		//copyStrategyTwo := TwoCdsWithoutRepetition(strategyTwo, bsubStarvationCodonTable)
		strategyThree := fixSequence(enzymeOptimizedBsubEcoli, bsubEcoliCodonTable, hostGenomeKmerTable)
		//copyStrategyThree := TwoCdsWithoutRepetition(strategyThree, bsubEcoliCodonTable)
		strategyFour := fixSequence(enzymeOptimizedStarvationEcoli, starvationEcoliCodonTable, hostGenomeKmerTable)

		output = append(output,
			fasta.Fasta{enzyme.Name + " | Codon Optimized By Strategy #1 Bacillus Subtilis KO7", strategyOne},
			//fasta.Fasta{enzyme.Name + " | Copy without repetititon of Codon Optimized By Strategy #1 Bacillus Subtilis KO7", copyStrategyOne},
			fasta.Fasta{enzyme.Name + " | Codon Optimized By Strategy #2 Bacillus Subtilis Starvation Genes", strategyTwo},
			//fasta.Fasta{enzyme.Name + " | Copy without repetititon of Codon Optimized By Strategy #2 Bacillus Subtilis Starvation Genes", copyStrategyTwo},
			fasta.Fasta{enzyme.Name + " | Codon Optimized By Strategy #3 Both species Bacillus Subtilis KO7 and E. coli K12", strategyThree},
			//fasta.Fasta{enzyme.Name + " | Copy without repetititon of Codon Optimized By Strategy #3 Both species Bacillus Subtilis KO7 and E. coli K12", copyStrategyThree},
			fasta.Fasta{enzyme.Name + " | Codon Optimized By Strategy #4 Bacillus Subtilis Starvation genes and E. coli K12", strategyFour},
		)

	}
	fmt.Println("Writing outputs...")
	fasta.Write(output, "data/output/output.fasta")

	fmt.Println("Finished! Check your fasta file with the results in data/output directory.")

}

func CodonOptimization(enzymeSequence string, codonTable codon.Table) string {
	// Poly generally makes Codon Optimization by receiving a list of protein sequences, but we actually have now CDSs
	// So we should first translate CDSs. We will be using the Eubacterial genetic code table 11, you could take a look
	// at this table in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	enzymeTranslated, _ := codon.Translate(enzymeSequence, codon.GetCodonTable(11))

	// Optimize sequence using the protein sequence and codon table
	optimizedSequence, _ := codon.Optimize(enzymeTranslated, codonTable)

	// Lets check if the codon optimization actually works by making some checks:
	// First one is if both codon sequences are different
	if optimizedSequence == enzymeSequence {
		fmt.Println("Both sequences are equal, some problem occur. They should be different because one is optimized. Checks what happened and run again.")
		os.Exit(0)
	}

	// Check if both translated sequences are equal
	protein, _ := codon.Translate(optimizedSequence, codon.GetCodonTable(11))
	if protein != enzymeTranslated {
		fmt.Println("These protein sequences aren't equal, some problem occur. They should be equal because codon optimization don't change any aminoacid.")
		os.Exit(0)
	}
	return optimizedSequence
}

// List of sequences that we should avoid in our software
func forbiddenSequencesList() []string {
	BsaI_bind_5prime := "GGTCTC" //5" GGTCTC N|      3"
	//3" CCAGAG N NNNN| 5"
	BsaI_bind_3prime := "GAGACC" //5" |NNNN N GAGACC 3"

	BbsI_bind_5prime := "GAAGAC" //5" GAAGAC NN|
	//3" CTTCTG NN NNNN|
	BbsI_bind_3prime := "GTCTTC" //5" |NNNN NN GTCTTC 3"

	BtgzI_bind_5prime := "GCGATG" //5" GCGATG NNN NNN NNN N|      3"
	//3" CGCTAC NNN NNN NNN N NNNN| 5"
	BtgzI_bind_3prime := "CATCGC" //5" |NNNN N NNN NNN NNN CATCGC 3"
	//3"      |N NNN NNN NNN GTAGCG 5"

	SapI_bind_5prime := "GCTCTTC" //5" GCTCTTC N|     3"
	//3" CGAGAAG N NNN| 5"
	SapI_bind_3prime := "GAAGAGC" //5" |NNN N GAAGAGC 3"
	//3"     |N CTTCTCG 5"

	BsmbI_bind_5prime := "CGTCTC" //5" CGTCTC N|      3"
	//3" GCAGAG N NNNN| 5"
	BsmbI_bind_3prime := "GAGACG" //5" |NNNN N GAGACG 3"
	//3"      |N CTCTGC 5"

	AarI_bind_5prime := "CACCTGC" //5" CACCTGC N NNN|       3"
	//3" GTGGACG N NNN NNNN|  5"
	AarI_bind_3prime := "GCAGGTG" //5" |NNNN NNN N GCAGGTG 3"
	//3"      |NNN N CGTCCAC 5"

	PmeI_bind := "GTTTAAAC" //5" GTTT|AAAC 3" blunt cutter

	forbiddenSeqList := []string{BsaI_bind_5prime,
		BsaI_bind_3prime,
		BbsI_bind_5prime,
		BbsI_bind_3prime,
		SapI_bind_5prime,
		SapI_bind_3prime,
		BsmbI_bind_5prime,
		BsmbI_bind_3prime,
		BtgzI_bind_5prime,
		BtgzI_bind_3prime,
		AarI_bind_5prime,
		AarI_bind_3prime,
		PmeI_bind,
		"AAAAA",
		"CCCCC",
		"GGGGG",
		"TTTTT"}
	return forbiddenSeqList
}

// getKmerTable receive a sequence string and a k int and generates a set of unique k-mers
func getKmerTable(k int, sequence string) map[string]bool {
	kmers := make(map[string]bool)
	for i := 0; i <= len(sequence)-k; i++ {
		kmers[strings.ToUpper(sequence[i:i+k])] = true
	}

	return kmers
}

func fixSequence(sequence string, codonTable codon.Table, genomeKmerTable map[string]bool) string {
	// Construct function that will remove unwanted properties in our sequences
	// Function#1: Remove unwanted sequences as restriction binding sites and homopolymers of length 5
	forbiddenSequences := forbiddenSequencesList()
	removeSequenceFunc := synthesis.RemoveSequence(forbiddenSequences)

	// Function#2: Remove secondary structures from begining to 1/3 of the sequence
	removeSecondaryFunc := synthesis.RemoveHairpin(20, 200)

	// Function#3: Remove repetition greater than 10 inside the sequence
	removeRepeatFunc := synthesis.RemoveRepeat(10)

	// Function#4: Remove repetitions between sequence and host genome
	genomeRemoveRepeatFunc := synthesis.GlobalRemoveRepeat(20, genomeKmerTable)

	// Added all those functions to a list and pass through FixCds function that will take care of our problems
	var functions []func(string, chan synthesis.DnaSuggestion, *sync.WaitGroup)
	functions = append(functions, removeSequenceFunc, removeRepeatFunc, genomeRemoveRepeatFunc, removeSecondaryFunc)

	fixedSeq, _, _ := synthesis.FixCds(":memory:", sequence, codonTable, functions)
	// Because FixCds actually remove stop codon we will concatenate it
	return fixedSeq + "TAA"
}

func TwoCdsWithoutRepetition(sequence string, codonTable codon.Table) string {

	forbiddenSequences := forbiddenSequencesList()
	removeSequenceFunc := synthesis.RemoveSequence(forbiddenSequences)

	removeRepeatFunc := synthesis.RemoveRepeat(10)

	globalRemoveRepeatFunc := synthesis.GlobalRemoveRepeat(20, getKmerTable(20, sequence))

	var functions []func(string, chan synthesis.DnaSuggestion, *sync.WaitGroup)
	functions = append(functions, removeSequenceFunc, removeRepeatFunc, globalRemoveRepeatFunc)

	fixedSeq, _, _ := synthesis.FixCds(":memory:", sequence, codonTable, functions)
	// Because FixCds actually remove stop codon we will concatenate it
	return fixedSeq + "TAA"
}
