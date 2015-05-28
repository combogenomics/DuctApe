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

/usr/bin/python dape -v init || die "/usr/bin/python dape -v init"
/usr/bin/python dape -v add Rm1021 -c red || die "/usr/bin/python dape -v add"
/usr/bin/python dgenome -v add test/input/Rm1021.faa Rm1021 || die "/usr/bin/python dgenome -v add"
/usr/bin/python dgenome -v add-ko test/input/ko_Rm1021.tsv || die "/usr/bin/python dgenome -v add-ko"
/usr/bin/python dgenome start || die "/usr/bin/python dgenome start"
/usr/bin/python dgenome -v stats || die "/usr/bin/python dgenome -v stats"
/usr/bin/python dgenome -v export || die "/usr/bin/python dgenome -v export"

/usr/bin/python dphenome -v add test/input/Rm1021.csv Rm1021 || die "/usr/bin/python dphenome -v add"
/usr/bin/python dphenome -v zero || die "/usr/bin/python dphenome -v zero"
/usr/bin/python dphenome start -f > /dev/null || die "/usr/bin/python dphenome start"
/usr/bin/python dphenome -v purge keep-max || die "/usr/bin/python dphenome -v purge"
/usr/bin/python dphenome -v restore || die "/usr/bin/python dphenome -v restore"
/usr/bin/python dphenome plot > /dev/null || die "/usr/bin/python dphenome plot"
/usr/bin/python dphenome -v rings || die "/usr/bin/python dphenome -v rings"
/usr/bin/python dphenome -v stats || die "/usr/bin/python dphenome -v stats"
/usr/bin/python dphenome -v export || die "/usr/bin/python dphenome -v export"

/usr/bin/python dape start -s > /dev/null || die "/usr/bin/python dape start"

/usr/bin/python dape -v export || die "/usr/bin/python dape -v export"

/usr/bin/python dape -v clear --keep-org || die "/usr/bin/python dape -v clear"
sleep 5
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"

rm ductape.db
rm -rf tmp

echo -e $green"Deletion AND insertion mutants"$reset

/usr/bin/python dape -v init || die "/usr/bin/python dape -v init"
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"
sleep 5
/usr/bin/python dape -v add Rm1021 -c red || die "/usr/bin/python dape -v add"
/usr/bin/python dape -v add-mut -k deletion -c blue -m Rm1021 del || die "/usr/bin/python dape -v add-mut"
/usr/bin/python dape -v add-mut -k insertion -c green -m Rm1021 add || die "/usr/bin/python dape -v add-mut"
/usr/bin/python dgenome -v add test/input/Rm1021.faa Rm1021 || die "/usr/bin/python dgenome -v add"
/usr/bin/python dgenome -v add test/input/del.faa del || die "/usr/bin/python dgenome -v add"
/usr/bin/python dgenome -v add test/input/add.faa add || die "/usr/bin/python dgenome -v add"
/usr/bin/python dgenome -v add-ko test/input/ko_Rm1021.tsv || die "/usr/bin/python dgenome -v add-ko"
/usr/bin/python dgenome -v add-ko test/input/del.tsv || die "/usr/bin/python dgenome -v add-ko"
/usr/bin/python dgenome -v add-ko test/input/add.tsv || die "/usr/bin/python dgenome -v add-ko"
/usr/bin/python dgenome start > /dev/null || die "/usr/bin/python dgenome start"
/usr/bin/python dgenome -v stats || die "/usr/bin/python dgenome -v stats"
/usr/bin/python dgenome -v export || die "/usr/bin/python dgenome -v export"

/usr/bin/python dphenome -v add test/input/Rm1021.csv Rm1021 || die "/usr/bin/python dphenome -v add"
/usr/bin/python dphenome -v add test/input/del.csv del || die "/usr/bin/python dphenome -v add"
/usr/bin/python dphenome -v add test/input/add.csv add || die "/usr/bin/python dphenome -v add"
/usr/bin/python dphenome -v zero || die "/usr/bin/python dphenome -v zero"
/usr/bin/python dphenome start -f > /dev/null || die "/usr/bin/python dphenome start"
/usr/bin/python dphenome -v purge keep-max || die "/usr/bin/python dphenome -v purge"
/usr/bin/python dphenome -v restore || die "/usr/bin/python dphenome -v restore"
/usr/bin/python dphenome plot > /dev/null || die "/usr/bin/python dphenome plot"
/usr/bin/python dphenome -v rings || die "/usr/bin/python dphenome -v rings"
/usr/bin/python dphenome -v stats || die "/usr/bin/python dphenome -v stats"
/usr/bin/python dphenome -v export || die "/usr/bin/python dphenome -v export"

/usr/bin/python dape start -s > /dev/null || die "/usr/bin/python dape start"

/usr/bin/python dape -v export || die "/usr/bin/python dape -v export"

/usr/bin/python dape -v clear --keep-org || die "/usr/bin/python dape -v clear"
sleep 5
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"

rm ductape.db
rm -rf tmp

echo -e $green"Pangenome"$reset

/usr/bin/python dape -v init || die "/usr/bin/python dape -v init"
/usr/bin/python dape -v add-multi Rm1021 AK83 AK58 BL225C || die "/usr/bin/python dape -v add-multi"
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"
sleep 5

/usr/bin/python dgenome -v add-dir test/input/pangenome || die "/usr/bin/python dgenome -v add-dir"
/usr/bin/python dgenome -v add-ko test/input/pangenome/ko.tab || die "/usr/bin/python dgenome -v add-ko"
/usr/bin/python dgenome -v add-orth test/input/pangenome/pangenome.tsv || die "/usr/bin/python dgenome -v add-orth"
/usr/bin/python dgenome start > /dev/null || die "/usr/bin/python dgenome start"
/usr/bin/python dgenome -v annotate || die "/usr/bin/python dgenome -v annotate"
/usr/bin/python dgenome -v deannotate || die "/usr/bin/python dgenome -v deannotate"
/usr/bin/python dgenome -v annotate || die "/usr/bin/python dgenome -v annotate"
/usr/bin/python dgenome -v stats || die "/usr/bin/python dgenome -v stats"
/usr/bin/python dgenome -v export || die "/usr/bin/python dgenome -v export"

/usr/bin/python dphenome -v add-dir test/input/pangenome || die "/usr/bin/python dphenome -v add-dir"
/usr/bin/python dphenome -v zero || die "/usr/bin/python dphenome -v zero"
/usr/bin/python dphenome start -f > /dev/null || die "/usr/bin/python dphenome start"
/usr/bin/python dphenome -v purge keep-max || die "/usr/bin/python dphenome -v purge"
/usr/bin/python dphenome -v restore PM03B || die "/usr/bin/python dphenome -v restore"
/usr/bin/python dphenome -v purge keep-min PM03B || die "/usr/bin/python dphenome -v purge"
/usr/bin/python dphenome plot > /dev/null || die "/usr/bin/python dphenome plot"
/usr/bin/python dphenome -v rings || die "/usr/bin/python dphenome -v rings"
/usr/bin/python dphenome -v stats || die "/usr/bin/python dphenome -v stats"
/usr/bin/python dphenome -v export || die "/usr/bin/python dphenome -v export"

/usr/bin/python dape start -s > /dev/null || die "/usr/bin/python dape start"

/usr/bin/python dape -v export || die "/usr/bin/python dape -v export"
/usr/bin/python dape -v clear --keep-org || die "/usr/bin/python dape -v clear"
sleep 5
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"

rm ductape.db

echo -e $green"Custom plates"$reset

/usr/bin/python dape -v init || die "dape init"
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"
sleep 5
/usr/bin/python dphenome -v import-plates test/input/newplate.tsv || die "dphenome import-plates (good)"
sleep 5
/usr/bin/python dphenome -v export || die "dphenome export"
/usr/bin/python dphenome -v import-plates test/input/newplate_wrong.tsv && die "dphenome import-plates (wrong)"
/usr/bin/python dape -v add Rm1021 -c red || die "dape add"
/usr/bin/python dphenome -v add test/input/Rm1021newplate.csv Rm1021 || die "dphenome add"
/usr/bin/python dphenome -v add test/input/Rm1021.yml Rm1021 || die "dphenome add (yml)"
/usr/bin/python dphenome -v zero || die "dphenome zero"
/usr/bin/python dphenome start -f > /dev/null || die "dphenome start"
sleep 5
/usr/bin/python dphenome -v purge keep-min || die "dphenome purge"
/usr/bin/python dphenome -v restore || die "dphenome restore"
/usr/bin/python dphenome plot > /dev/null|| die "dphenome plot"
/usr/bin/python dphenome -v rings || die "dphenome rings"
/usr/bin/python dphenome -v stats || die "dphenome stats"
/usr/bin/python dphenome -v export || die "dphenome export"
/usr/bin/python dphenome start -f -r > /dev/null || die "dphenome start"

rm ductape.db

echo -e $green"All tests passed"$reset
