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

python dape -v init || die "python dape -v init"
python dape -v add Rm1021 -c red || die "python dape -v add"
python dgenome -v add test/input/Rm1021.faa Rm1021 || die "python dgenome -v add"
python dgenome -v add-ko test/input/ko_Rm1021.tsv || die "python dgenome -v add-ko"
python dgenome start || die "python dgenome start"
python dgenome -v stats || die "python dgenome -v stats"
python dgenome -v export || die "python dgenome -v export"

python dphenome -v add test/input/Rm1021.csv Rm1021 || die "python dphenome -v add"
python dphenome -v zero || die "python dphenome -v zero"
python dphenome start -f > /dev/null || die "python dphenome start"
python dphenome -v purge keep-max || die "python dphenome -v purge"
python dphenome -v restore || die "python dphenome -v restore"
python dphenome plot > /dev/null || die "python dphenome plot"
python dphenome -v rings || die "python dphenome -v rings"
python dphenome -v stats || die "python dphenome -v stats"
python dphenome -v export || die "python dphenome -v export"

python dape start -s > /dev/null || die "python dape start"

python dape -v export || die "python dape -v export"

python dape -v clear --keep-org || die "python dape -v clear"
sleep 5
python dape -v import kegg.tsv || die "python dape -v import"

rm ductape.db
rm -rf tmp

echo -e $green"Deletion AND insertion mutants"$reset

python dape -v init || die "python dape -v init"
python dape -v import kegg.tsv || die "python dape -v import"
sleep 5
python dape -v add Rm1021 -c red || die "python dape -v add"
python dape -v add-mut -k deletion -c blue -m Rm1021 del || die "python dape -v add-mut"
python dape -v add-mut -k insertion -c green -m Rm1021 add || die "python dape -v add-mut"
python dgenome -v add test/input/Rm1021.faa Rm1021 || die "python dgenome -v add"
python dgenome -v add test/input/del.faa del || die "python dgenome -v add"
python dgenome -v add test/input/add.faa add || die "python dgenome -v add"
python dgenome -v add-ko test/input/ko_Rm1021.tsv || die "python dgenome -v add-ko"
python dgenome -v add-ko test/input/del.tsv || die "python dgenome -v add-ko"
python dgenome -v add-ko test/input/add.tsv || die "python dgenome -v add-ko"
python dgenome start > /dev/null || die "python dgenome start"
python dgenome -v stats || die "python dgenome -v stats"
python dgenome -v export || die "python dgenome -v export"

python dphenome -v add test/input/Rm1021.csv Rm1021 || die "python dphenome -v add"
python dphenome -v add test/input/del.csv del || die "python dphenome -v add"
python dphenome -v add test/input/add.csv add || die "python dphenome -v add"
python dphenome -v zero || die "python dphenome -v zero"
python dphenome start -f > /dev/null || die "python dphenome start"
python dphenome -v purge keep-max || die "python dphenome -v purge"
python dphenome -v restore || die "python dphenome -v restore"
python dphenome plot > /dev/null || die "python dphenome plot"
python dphenome -v rings || die "python dphenome -v rings"
python dphenome -v stats || die "python dphenome -v stats"
python dphenome -v export || die "python dphenome -v export"

python dape start -s > /dev/null || die "python dape start"

python dape -v export || die "python dape -v export"

python dape -v clear --keep-org || die "python dape -v clear"
sleep 5
python dape -v import kegg.tsv || die "python dape -v import"

rm ductape.db
rm -rf tmp

echo -e $green"Pangenome"$reset

python dape -v init || die "python dape -v init"
python dape -v add-multi Rm1021 AK83 AK58 BL225C || die "python dape -v add-multi"
python dape -v import kegg.tsv || die "python dape -v import"
sleep 5

python dgenome -v add-dir test/input/pangenome || die "python dgenome -v add-dir"
python dgenome -v add-ko test/input/pangenome/ko.tab || die "python dgenome -v add-ko"
python dgenome -v add-orth test/input/pangenome/pangenome.tsv || die "python dgenome -v add-orth"
python dgenome start > /dev/null || die "python dgenome start"
python dgenome -v annotate || die "python dgenome -v annotate"
python dgenome -v deannotate || die "python dgenome -v deannotate"
python dgenome -v annotate || die "python dgenome -v annotate"
python dgenome -v stats || die "python dgenome -v stats"
python dgenome -v export || die "python dgenome -v export"

python dphenome -v add-dir test/input/pangenome || die "python dphenome -v add-dir"
python dphenome -v zero || die "python dphenome -v zero"
python dphenome start -f > /dev/null || die "python dphenome start"
python dphenome -v purge keep-max || die "python dphenome -v purge"
python dphenome -v restore PM03B || die "python dphenome -v restore"
python dphenome -v purge keep-min PM03B || die "python dphenome -v purge"
python dphenome plot > /dev/null || die "python dphenome plot"
python dphenome -v rings || die "python dphenome -v rings"
python dphenome -v stats || die "python dphenome -v stats"
python dphenome -v export || die "python dphenome -v export"

python dape start -s > /dev/null || die "python dape start"

python dape -v export || die "python dape -v export"
python dape -v clear --keep-org || die "python dape -v clear"
sleep 5
python dape -v import kegg.tsv || die "python dape -v import"

rm ductape.db

echo -e $green"Custom plates"$reset

python dape -v init || die "dape init"
python dape -v import kegg.tsv || die "python dape -v import"
sleep 5
python dphenome -v import-plates test/input/newplate.tsv || die "dphenome import-plates (good)"
sleep 5
python dphenome -v export || die "dphenome export"
python dphenome -v import-plates test/input/newplate_wrong.tsv && die "dphenome import-plates (wrong)"
python dape -v add Rm1021 -c red || die "dape add"
python dphenome -v add test/input/Rm1021newplate.csv Rm1021 || die "dphenome add"
python dphenome -v add test/input/Rm1021.yml Rm1021 || die "dphenome add (yml)"
python dphenome -v zero || die "dphenome zero"
python dphenome start -f > /dev/null || die "dphenome start"
sleep 5
python dphenome -v purge keep-min || die "dphenome purge"
python dphenome -v restore || die "dphenome restore"
python dphenome plot > /dev/null|| die "dphenome plot"
python dphenome -v rings || die "dphenome rings"
python dphenome -v stats || die "dphenome stats"
python dphenome -v export || die "dphenome export"
python dphenome start -f -r > /dev/null || die "dphenome start"

rm ductape.db

echo -e $green"All tests passed"$reset
