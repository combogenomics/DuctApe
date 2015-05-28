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
/usr/bin/python dphenome trim || die "/usr/bin/python dphenome trim"
/usr/bin/python dphenome start -f > /dev/null || die "/usr/bin/python dphenome start"
/usr/bin/python dphenome -v purge replica -r 2 || die "/usr/bin/python dphenome purge replica"
/usr/bin/python dphenome -v restore -r 2 || die "/usr/bin/python dphenome restore replica"
/usr/bin/python dphenome -v purge keep-max || die "/usr/bin/python dphenome purge"
/usr/bin/python dphenome -v restore PM03B || die "/usr/bin/python dphenome restore"
/usr/bin/python dphenome -v purge keep-min PM03B || die "/usr/bin/python dphenome purge plate"
/usr/bin/python dphenome plot > /dev/null || die "/usr/bin/python dphenome plot"
/usr/bin/python dphenome plot PM01 > /dev/null || die "/usr/bin/python dphenome plot PM01"
/usr/bin/python dphenome plot PM01 H12 > /dev/null || die "/usr/bin/python dphenome plot PM01 H12"
/usr/bin/python dphenome -v rings || die "/usr/bin/python dphenome -v rings"
/usr/bin/python dphenome -v stats || die "/usr/bin/python dphenome -v stats"
/usr/bin/python dphenome -v export || die "/usr/bin/python dphenome -v export"

/usr/bin/python dape start -s > /dev/null || die "/usr/bin/python dape start"

/usr/bin/python dape -v export || die "/usr/bin/python dape -v export"
/usr/bin/python dape -v clear --keep-org || die "/usr/bin/python dape -v clear"
sleep 5
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"

rm ductape.db
