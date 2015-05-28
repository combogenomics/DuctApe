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
python dphenome trim || die "python dphenome trim"
python dphenome start -f > /dev/null || die "python dphenome start"
python dphenome -v purge replica -r 2 || die "python dphenome purge replica"
python dphenome -v restore -r 2 || die "python dphenome restore replica"
python dphenome -v purge keep-max || die "python dphenome purge"
python dphenome -v restore PM03B || die "python dphenome restore"
python dphenome -v purge keep-min PM03B || die "python dphenome purge plate"
python dphenome plot > /dev/null || die "python dphenome plot"
python dphenome plot PM01 > /dev/null || die "python dphenome plot PM01"
python dphenome plot PM01 H12 > /dev/null || die "python dphenome plot PM01 H12"
python dphenome -v rings || die "python dphenome -v rings"
python dphenome -v stats || die "python dphenome -v stats"
python dphenome -v export || die "python dphenome -v export"

python dape start -s > /dev/null || die "python dape start"

python dape -v export || die "python dape -v export"
python dape -v clear --keep-org || die "python dape -v clear"
sleep 5
python dape -v import kegg.tsv || die "python dape -v import"

rm ductape.db
