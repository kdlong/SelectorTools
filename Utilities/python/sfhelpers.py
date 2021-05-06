import ROOT
from scipy import special
import numpy as np

def float2double(hist):
    if hist.ClassName() == 'TH1D' or hist.ClassName() == 'TH2D':
        return hist
    elif hist.ClassName() == 'TH1F':
        new = ROOT.TH1D()
        hist.Copy(new)
    elif hist.ClassName() == 'TH2F':
        new = ROOT.TH2D()
        hist.Copy(new)
    else:
        raise Exception("Bad hist, dummy")
    return new

def invert2DHist(hist):
    new_hist = hist.Clone()
    ROOT.SetOwnership(new_hist, False)
    for x in range(hist.GetNbinsX()+1):
        for y in range(hist.GetNbinsY()+1):
            value = hist.GetBinContent(x, y)
            new_hist.SetBinContent(y, x, value)
    new_hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    new_hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    return new_hist

def setScaleFactorObj(infoVector, scaleFactorsObj, fScales):
    TotalSF = None
    for info in infoVector:
        # filename in 0, histname in 0
        f = ROOT.TFile(info[0])
        hist = float2double(f.Get(info[1]))
        if TotalSF == None:
            TotalSF = hist.Clone()
        else:
            TotalSF.Multiply(hist)
    scaleFactorsObj.Set2DHist(TotalSF)
    fScales.cd()
    scaleFactorsObj.Write()
    
def smoothingWeights(histlow, histhigh, binsl, turnon, turnoff, k=0.2, axishigh=2, axislow=3):
    low = np.digitize(turnon, binsl[axislow], right=True)
    high = np.digitize(turnoff, binsl[axislow], right=True)
    if histlow.shape[axislow] < histhigh.shape[axishigh]:
        temp = np.copy(histhigh)
        # This assumes axis is the Z axis...
        temp[:,:,:,:histlow.shape[axislow]] = histlow
        histlow = temp
    mid = int(low+0.5*(high-low))
    midval = binsl[axislow][mid]
    print(low, high, mid, midval)
    corr = 0.5*(1-special.erf(k*(binsl[axislow][:-1]-midval)))
    weights = (corr*histlow + (1-corr)*histhigh[np.newaxis,:,:,:])/histhigh[np.newaxis,:,:,:]
	# Not sure what's going on in the last bin now
    #weights[:,:,-1] = 1.
    return weights
