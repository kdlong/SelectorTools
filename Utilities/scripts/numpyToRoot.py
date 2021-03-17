import ROOT
import numpy as np
import argparse
import array

parser = argparse.ArgumentParser()
parser.add_argument('--numpyFile', '-f', required=True, type=str, help='.npz file with 3D histogram')
parser.add_argument('--rtfile', '-o', type=str, required=True, help='Output root file')
parser.add_argument('--obs', '-v', choices=['pt', 'eta', 'mass'], default='', help='plot only a 1D observalbe')
parser.add_argument('--histName', '-b', type=str, required=True, help='Name of histogram to write to ROOT file')
parser.add_argument('--altNames', '-a', type=str, nargs='*', help='Write multiple copies of hist with different names')
parser.add_argument('--procName', '-n', type=str, required=True, help='Name of folder in output file (data set name)')
args = parser.parse_args()

infile = np.load(args.numpyFile, allow_pickle=True)
rtfile = ROOT.TFile(args.rtfile, "recreate")
rtfile.mkdir(args.procName)
rtfile.cd(args.procName)

for histname in infile[::2]:
    hist = infile[histname]
    bins = infile[histname.replace("hist", "bins")]

    if args.obs:
        axes = ["mass", "eta", "pt"]
        idx = axes.index(args.obs)
        flatten = tuple([i for i,v in enumerate(axes) if args.obs not in v])
        bins = bins[idx]
        hist = np.sum(hist, axis=flatten)

    rthist = ROOT.TH1D(args.histName, args.histName, len(bins)-1, array.array('d', bins))
    for i, x in enumerate(hist):
        rthist.SetBinContent(i+1, x)
        rthist.SetBinError(i+1, x*1e-8)

    rthist.Write()

    for name in args.altNames:
        alt = rthist.Clone(name)
        alt.Write()
