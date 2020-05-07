#!/bin/bash

## to run: ./Fillhists.sh -o test.root -f "4top2016"

runTime=$(stat -c "%Y" .compile)
modTime=$(ls -tr src/* interface/* | tail -n 1 | xargs stat -c "%Y")

if [ "$modTime" -gt "$runTime" ]; then
    echo Recompiling
    scramv1 b -j 10 1>/dev/null
    if [ $? != 0 ]; then
	echo "Error"
	exit 1
    fi
    touch .compile
    echo
fi

EXE=./Utilities/scripts/makeHistFile.py
# EXE=./Utilities/scripts/makeFakeRates.py
ANA=ThreeLep
YEAR=2016
# INP=TwoLep_Met40
INP=NoSkim_MultYear
# INP=TwoLep_Met25

# # SEL=FourTopCutBasedEl
# SEL=FourTopMVAEl
SEL=MVAStudy

# ANA=FR
# INP=FakeRateEOS
# YEAR=2016
# SEL=FakeRate

# EXE=./Utilities/scripts/makeHistFile.py
# ANA=Efficiency:ThreeLep
# # ANA=ThreeLep
# YEAR=2016
# INP=TwoLep_Met40
# # SEL=FourTopCutBasedEl
# SEL=test


$EXE -a $ANA --input_tier $INP -s $SEL \
     --year $YEAR ${@}

