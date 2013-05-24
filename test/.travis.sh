#!/bin/bash

green="\033[1;32m"
red="\033[1;31m"
reset="\033[0m"

die () {
        echo -e $red"############"$reset
	echo -e $red$1$reset
	echo -e $red"Test failed!"$reset
	echo -e $red"############"$reset
	exit 1
}

echo -e $green"Single organism"$reset

./dape init || die "./dape init"
./dape add Rm1021 -c red || die "./dape add"
./dgenome add test/input/Rm1021.faa Rm1021 || die "./dgenome add"
./dgenome add-ko test/input/ko_Rm1021.tsv || die "./dgenome add-ko"
./dgenome start || die "./dgenome start"
./dgenome stats || die "./dgenome stats"
./dgenome export || die "./dgenome export"

./dphenome add test/input/Rm1021.csv Rm1021 || die "./dphenome add"
./dphenome zero || die "./dphenome zero"
./dphenome start -f || die "./dphenome start"
./dphenome purge keep-max || die "./dphenome purge"
./dphenome restore || die "./dphenome restore"
./dphenome plot || die "./dphenome plot"
./dphenome rings || die "./dphenome rings"
./dphenome stats || die "./dphenome stats"
./dphenome export || die "./dphenome export"

./dape start -p || die "./dape start"

./dape export || die "./dape export"

./dape clear --keep-org || die "./dape clear"
./dape import kegg.tsv || die "./dape import"

./dape clear --keep-kegg || die "./dape clear"

mv ductape.db single.db

echo -e $green"Deletion AND insertion mutants"$reset

./dape init || die "./dape init"
./dape import test/input/kegg.tsv || die "./dape import"
./dape add Rm1021 -c red || die "./dape add"
./dape add-mut -k deletion -c blue -m Rm1021 del || die "./dape add-mut"
./dape add-mut -k insertion -c green -m Rm1021 add || die "./dape add-mut"
./dgenome add test/input/Rm1021.faa Rm1021 || die "./dgenome add"
./dgenome add test/input/del.faa del || die "./dgenome add"
./dgenome add test/input/add.faa add || die "./dgenome add"
./dgenome add-ko test/input/ko_Rm1021.tsv || die "./dgenome add-ko"
./dgenome add-ko test/input/del.tsv || die "./dgenome add-ko"
./dgenome add-ko test/input/add.tsv || die "./dgenome add-ko"
./dgenome start || die "./dgenome start"
./dgenome stats || die "./dgenome stats"
./dgenome export || die "./dgenome export"

./dphenome add test/input/Rm1021.csv Rm1021 || die "./dphenome add"
./dphenome add test/input/del.csv del || die "./dphenome add"
./dphenome add test/input/add.csv add || die "./dphenome add"
./dphenome zero || die "./dphenome zero"
./dphenome start -f || die "./dphenome start"
./dphenome purge keep-max || die "./dphenome purge"
./dphenome restore || die "./dphenome restore"
./dphenome plot || die "./dphenome plot"
./dphenome rings || die "./dphenome rings"
./dphenome stats || die "./dphenome stats"
./dphenome export || die "./dphenome export"

./dape start -p || die "./dape start"

./dape export || die "./dape export"

./dape clear --keep-org || die "./dape clear"
./dape import kegg.tsv || die "./dape import"

./dape clear --keep-kegg || die "./dape clear"

mv ductape.db mut.db

echo -e $green"Pangenome"$reset

./dape init || die "./dape init"
./dape add-multi Rm1021 AK83 AK58 BL225C || die "./dape add-multi"
./dape import test/input/kegg.tsv || die "./dape import"

./dgenome add-dir test/input/pangenome || die "./dgenome add-dir"
./dgenome add-ko test/input/pangenome/ko.tab || die "./dgenome add-ko"
./dgenome add-orth test/input/pangenome/pangenome.tsv || die "./dgenome add-orth"
./dgenome start || die "./dgenome start"
./dgenome annotate || die "./dgenome annotate"
./dgenome deannotate || die "./dgenome deannotate"
./dgenome annotate || die "./dgenome annotate"
./dgenome stats || die "./dgenome stats"
./dgenome export || die "./dgenome export"

./dphenome add-dir test/input/pangenome || die "./dphenome add-dir"
./dphenome zero || die "./dphenome zero"
./dphenome start -f || die "./dphenome start"
./dphenome purge keep-max || die "./dphenome purge"
./dphenome restore PM03B || die "./dphenome restore"
./dphenome purge keep-min PM03B || die "./dphenome purge"
./dphenome plot || die "./dphenome plot"
./dphenome rings || die "./dphenome rings"
./dphenome stats || die "./dphenome stats"
./dphenome export || die "./dphenome export"

./dape start -p || die "./dape start"

./dape export || die "./dape export"
./dape clear --keep-org || die "./dape clear"
./dape import kegg.tsv || die "./dape import"

./dape clear --keep-kegg || die "./dape clear"

mv ductape.db pangenome.db

echo -e $green"All tests passed"$reset
