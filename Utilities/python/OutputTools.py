import ROOT
import datetime
import subprocess
import sys
import logging
import numpy as np
import array

def addMetaInfo(fOut):
    metaInfo = fOut.mkdir("MetaInfo")
    metaInfo.cd()
    time = ROOT.TNamed("datetime", str(datetime.datetime.now()))
    command = ROOT.TNamed("command", getScriptCall())
    try: 
        githash = ROOT.TNamed("githash", gitHash())
        gitdiff = ROOT.TNamed("gitdiff", gitDiff())
        githash.Write()
        gitdiff.Write()
    except subprocess.CalledProcessError:
        logging.warning("Not a git repo, skipping git hash/diff")
        pass

    time.Write()
    command.Write()

def gitHash():
    return subprocess.check_output(['git', 'log', '-1', '--format="%H"'], encoding='UTF-8')

def gitDiff():
    return subprocess.check_output(['git', 'diff',], encoding='UTF-8')

def getScriptCall():
    return ' '.join(sys.argv)

def writeOutputListItem(item, directory):
    if item.ClassName() == "TList":
        d = directory.Get(item.GetName())
        if not d:
            d = directory.mkdir(item.GetName())
            ROOT.SetOwnership(d, False)
        for subItem in item:
            writeOutputListItem(subItem, d)
    elif hasattr(item, 'Write'):
        directory.cd()
        item.Write()
    else:
        logging.warning("Couldn't write output item: %s " % repr(item))
    directory.cd()

def numpy3DHistToRoot(histName, bins, hist):
    rthist = ROOT.TH3D(histName, histName, 
        len(bins[0])-1, array.array('d', bins[0]), len(bins[1])-1, 
        array.array('d', bins[1]), len(bins[2])-1, array.array('d', bins[2]))
    for ix, x in enumerate(hist):
        for iy, y in enumerate(x):
            for iz, val in enumerate(y):
                rthist.SetBinContent(ix+1, iy+1, iz+1, val)
    return rthist
