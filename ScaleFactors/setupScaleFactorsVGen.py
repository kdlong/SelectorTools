#!/usr/bin/env python3
# Example: python3 ScaleFactors/setupScaleFactorsVGen.py --name scetlibCorr3D_Z --hist_names hist mass_y_pT_3D_prefsr_mm -n DYm50_scetlib -d DYm50_minnlo --smooth 35 45 --npOut test.npz /eos/user/k/kelong/HistFiles/ZGen/Scetlib/inclusive_Z_pT.npz /eos/user/k/kelong/HistFiles/ZGen/DY_MiNNLO_3Dreweight.root -a ZGen
# Example Wp: python3 ScaleFactors/setupScaleFactorsVGen.py --name scetlibCorr3D_Wp --hist_names hist mass_y_pT_3D_prefsr_mp -n wpmunu_scetlib -d wpmunu_minnlo --smooth 35 45 --npOut testWp.npz /eos/user/k/kelong/HistFiles/WGen/Scetlib/inclusive_Wp_pT.npz /eos/user/k/kelong/HistFiles/WGen/WToMuNu_MiNNLO_N3LLCorr.root -a WGen
# Example Wm (added to Wp file):  python3 ScaleFactors/setupScaleFactorsVGen.py --name scetlibCorr3D_Wm --hist_names hist mass_y_pT_3D_prefsr_mn -n wmmunu_scetlib -d wmmunu_minnlo --smooth 35 45 --npOut testWm.npz /eos/user/k/kelong/HistFiles/WGen/Scetlib/inclusive_Wm_pT.npz /eos/user/k/kelong/HistFiles/WGen/WToMuNu_MiNNLO_N3LLCorr.root -a WGen --append
import ROOT
import argparse
import os
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import time
from python import UserInput,OutputTools,HistTools,ConfigureJobs,sfhelpers
import array
import uproot
import numpy as np
import logging

def getRootHist(rtfile, dataset, histname, xsec, rebin):
    hist = rtfile.Get("/".join([dataset, histname]))
    sumw = rtfile.Get("/".join([dataset, "sumweights"]))
    if sumw:
        logging.debug("xsec is %0.2f sumw is %0.2f" % (xsec, sumw.Integral()))
        hist.Scale(xsec/sumw.Integral())

    return hist

def getHistUproot(rtfile, dataset, histname, xsec):
    f = uproot.open(rtfile)
    h = f["/".join([dataset, histname])]
    hist, bins = h.numpy()
    bins = bins[0]
    sumw = np.sum(f["/".join([dataset, "sumweights"])])
    hist = hist*xsec/sumw
    logging.debug("xsec is %0.2f sumw is %0.2f" % (xsec, sumw))
    logging.debug("Dataset is %s, integral is %0.2f" % (dataset, np.sum(hist)))
    return hist,bins

def getHistNpz(npzfile, dataset, histname, binsname):
    if binsname == "":
        binsname = histname.replace("hist", "bins")
    f = np.load(npzfile, allow_pickle=True)
    hist = f[histname]
    bins = f[binsname]
    if "scetlib" in dataset and bins[3][0] == 0.25:
        # Assuming these are points
        hist = hist*0.5
        bins[3] = bins[3]-0.25
        # Something weird going on here
    return hist,bins

def getHist(filename, dataset, histname, xsec=1, binsname=""):
    if "root" in filename[-4:]:
        return getHistUproot(filename, dataset, histname, xsec)
    elif "npz" in filename[-3:]:
        return getHistNpz(filename, dataset, histname, binsname)
    else:
        raise ValueError("Invalid file type! Must be either .npz or .root. Found %s" % filename)

def getHists(filename, dataset, histnames, xsec=1, binsname=""):
    hist, bins = getHist(filename, dataset, histnames[0], xsec, binsname)
    if len(histnames) > 1:
        for h in histnames[1:]:
            add,_ = getHist(filename, dataset, h, xsec, binsname)
            hist = np.add(hist, add)
        # Assuming these are from the same channel
        hist /= len(histnames)
    return hist, bins

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, required=True, help="Name of output scale factor")
parser.add_argument("--hist_names", required=True, type=str, nargs=2, help="Name of histograms (numerator denominator)")
parser.add_argument("--append", action='store_true', help="Add to file (by default overwrite file)")
parser.add_argument("-a", "--analysis", type=str, help="Analysis name")
parser.add_argument("-n", "--numerator", required=True, type=str, help="Name of numerator dataset")
parser.add_argument("-d", "--denominator", required=True, type=str, help="Name of denominator dataset")
parser.add_argument("inputFiles", type=str, nargs=2, help="file(s) with histograms. <numerator> <denominator> if there are two")
parser.add_argument("--smooth", type=float, nargs='+', help="<low pt> <high pt>")
parser.add_argument("--npOut", type=str, help="store ratio as npz file")
parser.add_argument("--debug", action='store_true', help="Output debug info")
args = parser.parse_args()

logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO))

output_file = os.environ["CMSSW_BASE"]+'/src/Analysis/SelectorTools/data/%s_scaleFactors.root' % args.analysis
fScales = ROOT.TFile(output_file, 'recreate' if not args.append else 'update')

xsecs  = ConfigureJobs.getListOfFilesWithXSec([args.numerator, args.denominator], 
                                                os.path.expanduser("~/work/"))
numfile = args.inputFiles[0]
denomfile = args.inputFiles[0] if len(args.inputFiles) == 1 else args.inputFiles[1]

histnum, binsn = getHists(numfile, args.numerator, args.hist_names[0].split("+"), xsecs[args.numerator])
histdenom, binsd = getHists(denomfile, args.denominator, args.hist_names[1].split("+"), xsecs[args.denominator])

if len(histnum.shape) == len(histdenom.shape):
    corr = histnum/histdenom
elif len(histnum.shape) == 4 and len(histdenom.shape) == 3:
    corr = histnum/histdenom[np.newaxis,:,:,:]

logging.debug("HistNum is %s" % histnum)
logging.debug("HistDenom is %s" % histdenom)

if args.smooth:
    corr = sfhelpers.smoothingWeights(histnum, histdenom, binsn, *args.smooth)

if args.npOut:
    labels = {args.name : corr,
                "bins" : binsn,
                "numerator" : histnum,
                "denominator" : histdenom,
                "metaData" : np.array([OutputTools.getScriptCall(), OutputTools.gitHash(), OutputTools.gitHash(),])
            }
    np.savez(args.npOut, **labels)

for var in range(histnum.shape[0]):
    varname = args.name+"_var%i" % var

    hist = OutputTools.numpy3DHistToRoot(varname, binsd, corr[var,:,:,:])
    sf = ROOT.ScaleFactor(varname, varname)
    if hist.InheritsFrom("TH3"):
        sf.Set3DHist(hist, 0, 0)
    elif hist.InheritsFrom("TH1"):
        sf.Set1DHist(hist, 0, 0)
    fScales.cd()
    sf.Write()
