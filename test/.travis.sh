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

./dape -v init || die "./dape -v init"
./dape -v add Rm1021 -c red || die "./dape -v add"
./dgenome -v add test/input/Rm1021.faa Rm1021 || die "./dgenome -v add"
./dgenome -v add-ko test/input/ko_Rm1021.tsv || die "./dgenome -v add-ko"
./dgenome start > /dev/null || die "./dgenome start"
./dgenome -v stats || die "./dgenome -v stats"
./dgenome -v export || die "./dgenome -v export"

./dphenome -v add test/input/Rm1021.csv Rm1021 || die "./dphenome -v add"
./dphenome -v zero || die "./dphenome -v zero"
./dphenome start -f > /dev/null || die "./dphenome start"
./dphenome -v purge keep-max || die "./dphenome -v purge"
./dphenome -v restore || die "./dphenome -v restore"
./dphenome plot > /dev/null || die "./dphenome plot"
./dphenome -v rings || die "./dphenome -v rings"
./dphenome -v stats || die "./dphenome -v stats"
./dphenome -v export || die "./dphenome -v export"

./dape start -p > /dev/null || die "./dape start"

./dape -v export || die "./dape -v export"

./dape -v clear --keep-org || die "./dape -v clear"
sleep 5
./dape -v import kegg.tsv || die "./dape -v import"

rm ductape.db
rm -rf tmp

echo -e $green"Deletion AND insertion mutants"$reset

./dape -v init || die "./dape -v init"
./dape -v import kegg.tsv || die "./dape -v import"
./dape -v add Rm1021 -c red || die "./dape -v add"
./dape -v add-mut -k deletion -c blue -m Rm1021 del || die "./dape -v add-mut"
./dape -v add-mut -k insertion -c green -m Rm1021 add || die "./dape -v add-mut"
./dgenome -v add test/input/Rm1021.faa Rm1021 || die "./dgenome -v add"
./dgenome -v add test/input/del.faa del || die "./dgenome -v add"
./dgenome -v add test/input/add.faa add || die "./dgenome -v add"
./dgenome -v add-ko test/input/ko_Rm1021.tsv || die "./dgenome -v add-ko"
./dgenome -v add-ko test/input/del.tsv || die "./dgenome -v add-ko"
./dgenome -v add-ko test/input/add.tsv || die "./dgenome -v add-ko"
./dgenome start > /dev/null || die "./dgenome start"
./dgenome -v stats || die "./dgenome -v stats"
./dgenome -v export || die "./dgenome -v export"

./dphenome -v add test/input/Rm1021.csv Rm1021 || die "./dphenome -v add"
./dphenome -v add test/input/del.csv del || die "./dphenome -v add"
./dphenome -v add test/input/add.csv add || die "./dphenome -v add"
./dphenome -v zero || die "./dphenome -v zero"
./dphenome start -f > /dev/null || die "./dphenome start"
./dphenome -v purge keep-max || die "./dphenome -v purge"
./dphenome -v restore || die "./dphenome -v restore"
./dphenome plot > /dev/null || die "./dphenome plot"
./dphenome -v rings || die "./dphenome -v rings"
./dphenome -v stats || die "./dphenome -v stats"
./dphenome -v export || die "./dphenome -v export"

./dape start -p > /dev/null || die "./dape start"

./dape -v export || die "./dape -v export"

./dape -v clear --keep-org || die "./dape -v clear"
sleep 5
./dape -v import kegg.tsv || die "./dape -v import"

rm ductape.db
rm -rf tmp

echo -e $green"Pangenome"$reset

./dape -v init || die "./dape -v init"
./dape -v add-multi Rm1021 AK83 AK58 BL225C || die "./dape -v add-multi"
./dape -v import kegg.tsv || die "./dape -v import"

./dgenome -v add-dir test/input/pangenome || die "./dgenome -v add-dir"
./dgenome -v add-ko test/input/pangenome/ko.tab || die "./dgenome -v add-ko"
./dgenome -v add-orth test/input/pangenome/pangenome.tsv || die "./dgenome -v add-orth"
./dgenome start > /dev/null || die "./dgenome start"
./dgenome -v annotate || die "./dgenome -v annotate"
./dgenome -v deannotate || die "./dgenome -v deannotate"
./dgenome -v annotate || die "./dgenome -v annotate"
./dgenome -v stats || die "./dgenome -v stats"
./dgenome -v export || die "./dgenome -v export"

./dphenome -v add-dir test/input/pangenome || die "./dphenome -v add-dir"
./dphenome -v zero || die "./dphenome -v zero"
./dphenome start -f > /dev/null || die "./dphenome start"
./dphenome -v purge keep-max || die "./dphenome -v purge"
./dphenome -v restore PM03B || die "./dphenome -v restore"
./dphenome -v purge keep-min PM03B || die "./dphenome -v purge"
./dphenome plot > /dev/null || die "./dphenome plot"
./dphenome -v rings || die "./dphenome -v rings"
./dphenome -v stats || die "./dphenome -v stats"
./dphenome -v export || die "./dphenome -v export"

./dape start -p > /dev/null || die "./dape start"

./dape -v export || die "./dape -v export"
./dape -v clear --keep-org || die "./dape -v clear"
sleep 5
./dape -v import kegg.tsv || die "./dape -v import"

rm ductape.db

echo -e $green"Custom plates"$reset

./dape -v init || die "dape init"
./dphenome -v import-plates test/input/newplate.tsv || die "dphenome import-plates (good)"
sleep 5
./dphenome -v export || die "dphenome export"
./dphenome -v import-plates test/input/newplate_wrong.tsv && die "dphenome import-plates (wrong)"
./dape -v add Rm1021 -c red || die "dape add"
./dphenome -v add test/input/Rm1021newplate.csv Rm1021 || die "dphenome add"
./dphenome -v zero || die "dphenome zero"
./dphenome start -f > /dev/null || die "dphenome start"
sleep 5
./dphenome -v purge keep-min || die "dphenome purge"
./dphenome -v restore || die "dphenome restore"
./dphenome plot > /dev/null|| die "dphenome plot"
./dphenome -v rings || die "dphenome rings"
./dphenome -v stats || die "dphenome stats"
./dphenome -v export || die "dphenome export"

rm ductape.db

echo -e $green"All tests passed"$reset
