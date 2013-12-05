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

echo -e $green"Custom plates"$reset

./dape -v init || die "dape init"
./dape -v import kegg.tsv || die "./dape -v import"
sleep 5
./dphenome -v import-plates test/input/newplate.tsv || die "dphenome import-plates (good)"
sleep 5
./dphenome -v export || die "dphenome export"
./dphenome -v import-plates test/input/newplate_wrong.tsv && die "dphenome import-plates (wrong)"
./dape -v add Rm1021 -c red || die "dape add"
./dphenome -v add test/input/Rm1021newplate.csv Rm1021 || die "dphenome add"
./dphenome -v add test/input/Rm1021.yml Rm1021 || die "dphenome add (yml)"
./dphenome -v zero || die "dphenome zero"
./dphenome start -f > /dev/null || die "dphenome start"
sleep 5
./dphenome -v purge keep-min || die "dphenome purge"
./dphenome -v restore || die "dphenome restore"
./dphenome plot > /dev/null|| die "dphenome plot"
./dphenome -v rings || die "dphenome rings"
./dphenome -v stats || die "dphenome stats"
./dphenome -v export || die "dphenome export"
./dphenome start -f -r > /dev/null || die "dphenome start"

rm ductape.db

./dape init || die "dape init"
./dape add Rm1021 -c red || die "dape add"
./dphenome add test/input/Rm1021strangeplate.csv Rm1021 || die "dphenome add (strange)"
./dphenome add test/input/Rm1021strangeplate.yml Rm1021 || die "dphenome add (strange yml)"

rm ductape.db

echo -e $green"All tests passed"$reset
