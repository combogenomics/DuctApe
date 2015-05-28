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
/usr/bin/python dphenome start -g -e || die "/usr/bin/python dphenome start (elbow)"
/usr/bin/python dphenome start -f > /dev/null || die "/usr/bin/python dphenome start"
/usr/bin/python dphenome -v purge keep-max || die "/usr/bin/python dphenome -v purge"
/usr/bin/python dphenome -v restore || die "/usr/bin/python dphenome -v restore"
/usr/bin/python dphenome plot > /dev/null || die "/usr/bin/python dphenome plot"
/usr/bin/python dphenome plot PM01 > /dev/null || die "/usr/bin/python dphenome plot PM01"
/usr/bin/python dphenome plot PM01 H12 > /dev/null || die "/usr/bin/python dphenome plot PM01 H12"
/usr/bin/python dphenome -v rings || die "/usr/bin/python dphenome -v rings"
/usr/bin/python dphenome -v rings -r area || die "/usr/bin/python dphenome -v rings -r area"
/usr/bin/python dphenome -v rings -r areaz && die "/usr/bin/python dphenome -v rings -r area"
/usr/bin/python dphenome -v stats || die "/usr/bin/python dphenome -v stats"
/usr/bin/python dphenome -v export || die "/usr/bin/python dphenome -v export"

/usr/bin/python dape start -s > /dev/null || die "/usr/bin/python dape start"

/usr/bin/python dape -v export || die "/usr/bin/python dape -v export"

/usr/bin/python dape -v clear --keep-org || die "/usr/bin/python dape -v clear"
sleep 5
/usr/bin/python dape -v import kegg.tsv || die "/usr/bin/python dape -v import"

rm ductape.db
rm -rf tmp
