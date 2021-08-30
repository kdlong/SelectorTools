#!/bin/bash

# Following implementation by N. Smith, Fermilab
# https://gitlab.cern.ch/ncsmith/monoZ/tree/master/selector 

tar xvzf ${tarball}
CMSSW_RELEASE_BASE="${CMSSW_RELEASE_BASE}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd $$CMSSW_RELEASE_BASE
scramv1 runtime -sh
eval `scramv1 runtime -sh`
popd
export LD_LIBRARY_PATH=$$PWD/lib:$$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$$PWD:$$ROOT_INCLUDE_PATH
which python3
echo "CMSSW BASE is at:"
echo $CMSSW_RELEASE_BASE
echo "CMSSW BASE contains:"
ls $CMSSW_RELEASE_BASE
echo "The PATH var:"
echo $$PATH
echo "what env gives:"
env
echo "the SHELL var:"
echo $$SHELL

./Analysis/SelectorTools/Utilities/scripts/makeHistFile.py $$@ || exit $$?

