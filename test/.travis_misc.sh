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

echo -e $green"New CSV format"$reset

python dape init || die "dape init"
python dape add Ecoli -c red || die "dape add"
python dphenome add test/input/new_csv/single.csv Ecoli || die "dphenome add"
python dphenome start -g || die "dphenome start -g"
python dphenome export || die "dphenome export"

rm ductape.db

python dape init || die "dape init"
python dape add Ecoli -c red || die "dape add"
python dphenome add test/input/new_csv/multi.csv Ecoli || die "dphenome add"
python dphenome start -g || die "dphenome start -g"
python dphenome export || die "dphenome export"

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

python dape init || die "dape init"
python dape add Rm1021 -c red || die "dape add"
python dphenome add test/input/Rm1021strangeplate.csv Rm1021 || die "dphenome add (strange)"
python dphenome add test/input/Rm1021strangeplate.yml Rm1021 || die "dphenome add (strange yml)"

rm ductape.db

echo -e $green"All tests passed"$reset
