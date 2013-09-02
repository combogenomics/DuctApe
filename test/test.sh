#!/bin/bash

green="\033[1;32m"
red="\033[1;31m"
reset="\033[0m"

die () {
        echo -e $red"############"$reset
	echo -e $red$1$reset
	echo -e $red"Test failed!"$reset
	echo -e $red"############"$reset
	#cleanUp
	exit 1
}

cleanUp () {
	find . ! -name . -prune -type f -not -name '*sh' -not -name "*log" | xargs rm
	rm -rf tmp &> /dev/null
}

echo -e $green"Single organism"$reset

../dape init || die "dape init"
../dape add Rm1021 -c red || die "dape add"
../dgenome add input/Rm1021.faa Rm1021 || die "dgenome add"
../dgenome add-ko input/ko_Rm1021.tsv || die "dgenome add-ko"
../dgenome start || die "dgenome start"
../dgenome stats || die "dgenome stats"
../dgenome export || die "dgenome export"

../dphenome add input/Rm1021.csv Rm1021 || die "dphenome add"
../dphenome zero || die "dphenome zero"
../dphenome start -e -g || die "dphenome start (elbow)"
../dphenome start -f || die "dphenome start"
../dphenome purge keep-max || die "dphenome purge"
../dphenome restore || die "dphenome restore"
../dphenome plot || die "dphenome plot"
../dphenome rings || die "dphenome rings"
../dphenome stats || die "dphenome stats"
../dphenome export || die "dphenome export"

../dape map || die "dape map"
../dape start -p || die "dape start"

../dape export || die "dape export"

cp kegg.tsv input/ &> /dev/null

../dape clear --keep-org || die "dape clear"
../dape import kegg.tsv || die "dape import"
sleep 5

../dape clear --keep-kegg || die "dape clear"

cleanUp

echo -e $green"Deletion AND insertion mutants"$reset

../dape init || die "dape init"
../dape import input/kegg.tsv || die "dape import"
sleep 5
../dape add Rm1021 -c red || die "dape add"
../dape add-mut -k deletion -c blue -m Rm1021 del || die "dape add-mut"
../dape add-mut -k insertion -c green -m Rm1021 add || die "dape add-mut"
../dgenome add input/Rm1021.faa Rm1021 || die "dgenome add"
../dgenome add input/del.faa del || die "dgenome add"
../dgenome add input/add.faa add || die "dgenome add"
../dgenome add-ko input/ko_Rm1021.tsv || die "dgenome add-ko"
../dgenome add-ko input/del.tsv || die "dgenome add-ko"
../dgenome add-ko input/add.tsv || die "dgenome add-ko"
../dgenome start || die "dgenome start"
../dgenome stats || die "dgenome stats"
../dgenome export || die "dgenome export"

../dphenome add input/Rm1021.csv Rm1021 || die "dphenome add"
../dphenome add input/del.csv del || die "dphenome add"
../dphenome add input/add.csv add || die "dphenome add"
../dphenome zero || die "dphenome zero"
../dphenome start -f || die "dphenome start"
../dphenome purge keep-max || die "dphenome purge"
../dphenome restore || die "dphenome restore"
../dphenome plot || die "dphenome plot"
../dphenome rings || die "dphenome rings"
../dphenome stats || die "dphenome stats"
../dphenome export || die "dphenome export"

../dape map || die "dape map"
../dape start -p || die "dape start"

../dape export || die "dape export"

cp kegg.tsv input/ &> /dev/null

../dape clear --keep-org || die "dape clear"
../dape import kegg.tsv || die "dape import"
sleep 5

../dape clear --keep-kegg || die "dape clear"

cleanUp

echo -e $green"Pangenome"$reset

../dape init || die "dape init"
../dape add-multi Rm1021 AK83 AK58 BL225C || die "dape add-multi"
../dape import input/kegg.tsv || die "dape import"
sleep 5

../dgenome add-dir input/pangenome || die "dgenome add-dir"
../dgenome add-ko input/pangenome/ko.tab || die "dgenome add-ko"
../dgenome start || die "dgenome start"
../dgenome annotate || die "dgenome annotate"
../dgenome deannotate || die "dgenome deannotate"
../dgenome annotate || die "dgenome annotate"
../dgenome stats || die "dgenome stats"
../dgenome export || die "dgenome export"

../dphenome add-dir input/pangenome || die "dphenome add-dir"
../dphenome zero || die "dphenome zero"
../dphenome trim || die "dphenome trim"
../dphenome start -f || die "dphenome start"
../dphenome purge replica -r 2 || die "dphenome purge replica"
../dphenome restore -r 2 || die "dphenome restore replica"
../dphenome purge keep-max || die "dphenome purge"
../dphenome restore PM03B || die "dphenome restore"
../dphenome purge keep-min PM03B || die "dphenome purge plate"
../dphenome plot || die "dphenome plot"
../dphenome rings || die "dphenome rings"
../dphenome stats || die "dphenome stats"
../dphenome export || die "dphenome export"

../dape map || die "dape map"
../dape start -p || die "dape start"

../dape export || die "dape export"
../dape clear --keep-org || die "dape clear"
../dape import kegg.tsv || die "dape import"

../dape clear --keep-kegg || die "dape clear"

../dape add-multi Rm1021 AK83 AK58 BL225C || die "dape add-multi"
../dgenome add-dir input/pangenome || die "dgenome add-dir"
../dgenome add-orth pangenome.tsv || die "dgenome add-orth"

cleanUp

echo -e $green"Custom plates"$reset

../dape init || die "dape init"
../dape import input/kegg.tsv || die "dape import"
sleep 5
../dphenome import-plates input/newplate.tsv || die "dphenome import-plates (good)"
../dphenome export || die "dphenome export"
../dphenome import-plates input/newplate_wrong.tsv && die "dphenome import-plates (wrong)"
../dape add Rm1021 -c red || die "dape add"
../dphenome add input/Rm1021newplate.csv Rm1021 || die "dphenome add"
../dphenome add input/Rm1021.yml Rm1021 || die "dphenome add (yml)"
../dphenome zero || die "dphenome zero"
../dphenome start -f || die "dphenome start"
../dphenome purge keep-min || die "dphenome purge"
../dphenome restore || die "dphenome restore"
../dphenome plot || die "dphenome plot"
../dphenome rings || die "dphenome rings"
../dphenome stats || die "dphenome stats"
../dphenome export || die "dphenome export"
../dphenome start -f -r || die "dphenome start"

cleanUp

echo -e $green"All tests passed"$reset
