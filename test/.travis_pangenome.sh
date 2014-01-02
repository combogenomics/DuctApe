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

echo -e $green"Pangenome"$reset

./dape -v init || die "./dape -v init"
./dape -v add-multi Rm1021 AK83 AK58 BL225C || die "./dape -v add-multi"
./dape -v import kegg.tsv || die "./dape -v import"
sleep 5

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
./dphenome trim || die "./dphenome trim"
./dphenome start -f > /dev/null || die "./dphenome start"
./dphenome -v purge replica -r 2 || die "./dphenome purge replica"
./dphenome -v restore -r 2 || die "./dphenome restore replica"
./dphenome -v purge keep-max || die "./dphenome purge"
./dphenome -v restore PM03B || die "./dphenome restore"
./dphenome -v purge keep-min PM03B || die "./dphenome purge plate"
./dphenome plot > /dev/null || die "./dphenome plot"
./dphenome plot PM01 > /dev/null || die "./dphenome plot PM01"
./dphenome plot PM01 H12 > /dev/null || die "./dphenome plot PM01 H12"
./dphenome -v rings || die "./dphenome -v rings"
./dphenome -v stats || die "./dphenome -v stats"
./dphenome -v export || die "./dphenome -v export"

./dape start -p > /dev/null || die "./dape start"

./dape -v export || die "./dape -v export"
./dape -v clear --keep-org || die "./dape -v clear"
sleep 5
./dape -v import kegg.tsv || die "./dape -v import"

rm ductape.db
