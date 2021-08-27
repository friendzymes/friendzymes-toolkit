package main

import (
	"strconv"
	"strings"

	"github.com/Open-Science-Global/poly"
	"github.com/Open-Science-Global/poly/finder"
	"github.com/Open-Science-Global/poly/io/genbank"
	"github.com/Open-Science-Global/poly/transform"
)

func main() {

	base := "../data/"
	partFiles := []string{
		"pfu-sso7d.gb",
	}
	var parts []poly.Sequence
	for _, partFile := range partFiles {
		parts = append(parts, genbank.Read(base+partFile))
	}

	//Read Bsub PY79 genome and use this as input
	file := genbank.Read("../data/bsub-py79-genome.gb")
	hostGenome := strings.ToUpper(file.Sequence)
	// BsaI restriction binding sites
	for i, part := range parts {
		var functions []func(string) []finder.Match

		functions = append(functions, finder.ForbiddenSequence(restrictionBindingSitesListFindProblems()))
		functions = append(functions, finder.ForbiddenSequence(homologySequencesFindProblems()))
		functions = append(functions, finder.RemoveRepeat(10))
		functions = append(functions, finder.GlobalRemoveRepeat(20, hostGenome))
		functions = append(functions, AvoidHairpin(20, 200))

		problems := finder.Find(strings.ToUpper(part.Sequence), functions)

		part = finder.AddMatchesToSequence(problems, part)
		genbank.Write(part, "../data/output/"+partFiles[i])
	}
}

func AvoidHairpin(stemSize int, hairpinWindow int) func(string) []finder.Match {
	return func(sequence string) []finder.Match {
		var matches []finder.Match
		reverse := transform.ReverseComplement(sequence)
		for i := 0; i < len(sequence)-stemSize && len(sequence)-(i+hairpinWindow) >= 0; i++ {
			word := sequence[i : i+stemSize]
			rest := reverse[len(sequence)-(i+hairpinWindow) : len(sequence)-(i+stemSize)]
			if strings.Contains(rest, word) {
				location := strings.Index(rest, word)
				matches = append(matches, finder.Match{i, i + hairpinWindow - location - 1, "Harpin found in next " + strconv.Itoa(hairpinWindow) + "bp in reverse complementary sequence: " + word})
			}
		}
		return matches
	}
}

func homologySequencesFindProblems() []string {
	// I don't have to worry about TTTTTT and GGGGGG because I already try to find also by reverse complementary of each sequence in finder
	return []string{"AAAAAA", "CCCCCC"}
}

func restrictionBindingSitesListFindProblems() []string {
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
	}

	return forbiddenSeqList
}
