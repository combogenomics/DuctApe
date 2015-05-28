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

/usr/bin/python dape init || die "dape init"
/usr/bin/python dape add Rm1021 -c red || die "dape add"
/usr/bin/python dphenome add test/input/Rm1021strangeplate.csv Rm1021 || die "dphenome add (strange)"
/usr/bin/python dphenome add test/input/Rm1021strangeplate.yml Rm1021 || die "dphenome add (strange yml)"

rm ductape.db

echo -e $green"All tests passed"$reset
