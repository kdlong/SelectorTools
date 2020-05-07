#!/usr/bin/env python
import ROOT
import glob
import datetime
from python import UserInput,OutputTools
from python import ConfigureJobs
from python import SelectorTools,HistTools
import sys

ROOT.gROOT.SetBatch(True)


def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--proof", "-p", 
        action='store_true', help="Don't use proof")
    parser.add_argument("--lumi", "-l", type=float,
        default=35.87, help="luminosity value (in fb-1)")
    parser.add_argument("--output_file", "-o", type=str,
        default="", help="Output file name")
    parser.add_argument("--input_tier", type=str,
        default="", help="Selection stage of input files")
    parser.add_argument("--year", type=str,
        default="default", help="Year of Analysis")
    parser.add_argument("--maxEntries", "-m", type=int,
        default=-1, help="Max entries to process")
    return vars(parser.parse_args())

def getHistNames(channels):
    base_hists = [x+y for x in ["passingLoose", "passingTight"] \
            for y in "1DEta", "1DPt", "2D"]
    if len(channels) == 0:
        return base_hists
    return [x+"_"+y for x in base_hists for y in channels]

def addOverflow(hist):
    nxbins = hist.GetNbinsX()
    nybins = hist.GetNbinsY()
    for i in range(1, nybins+1):
        setbin = hist.GetBin(nxbins, i)
        obin = hist.GetBin(nxbins+1, i)
        hist.SetBinContent(setbin, hist.GetBinContent(obin)+hist.GetBinContent(setbin))
    for i in range(1,nxbins):
        setbin = hist.GetBin(i, nybins)
        obin = hist.GetBin(i, nybins+1)
        hist.SetBinContent(setbin, hist.GetBinContent(obin)+hist.GetBinContent(setbin))
    
    return hist

# Turn off overflow for FR hists (> 50 is pretty much all EWK anyway)
def makeCompositeHists(infile, name, members, channels, addRatios=True, overflow=False):
    composite = ROOT.TList()
    composite.SetName(name)
    for directory in [str(i) for i in members.keys()]:
        for histname in getHistNames(channels):
            hist = infile.Get("/".join([directory, histname]))
            if hist:
                sumhist = composite.FindObject(hist.GetName())
                if "data" not in directory and hist.GetEntries() > 0:
                    sumweights_hist = infile.Get("/".join([directory, "sumweights"]))
                    sumweights = sumweights_hist.Integral()
                    hist.Scale(members[directory]*1000*args['lumi']/sumweights)
                if overflow and isinstance(hist, ROOT.TH1):
                    addOverflow(hist)
            else:
                raise RuntimeError("hist %s was not produced for "
                    "dataset %s!" % (histname, directory))
            if not sumhist:
                sumhist = hist.Clone()
                composite.Add(sumhist)
            else:
                sumhist.Add(hist)

    if addRatios:
        ratios = getRatios(composite)
        for ratio in ratios:
            composite.Add(ratio) 
    return composite

def getRatios(hists):
    ratios = []
    for hist in hists:
        if "Tight" not in hist.GetName():
            continue
        ratio = hist.Clone()
        ratio.SetName(hist.GetName().replace("passingTight", "ratio"))
        if not ratio.GetSumw2():
            ratio.Sumw2()
        ratio.Divide(hists.FindObject(hist.GetName().replace("Tight", "Loose")))
        ratios.append(ratio)
    return ratios



runAnalysis = True
channels = ["e", "m"]


manager_path = ConfigureJobs.getManagerPath()
if manager_path not in sys.path:
    sys.path.insert(0, "/".join([manager_path, 
        "AnalysisDatasetManager", "Utilities/python"]))
    
args = getComLineArgs()
proof = 0
if args['proof']:
    ROOT.TProof.Open("workers=12")
    proof = ROOT.gProof
today = datetime.date.today().strftime("%d%b%Y")
fileName = args['output_file']
if not fileName:
    fileName = "data/fakeRate%s-%s.root" % (today, args['selection'])        

# scale factors
fScales = ROOT.TFile('data/scaleFactors.root')
muonIsoSF = fScales.Get('muonIsoSF')
muonIdSF = fScales.Get('muonTightIdSF')
electronTightIdSF = fScales.Get('electronTightIdSF')
electronGsfSF = fScales.Get('electronGsfSF')
pileupSF = fScales.Get('pileupSF')
sf_inputs = [electronTightIdSF, electronGsfSF, muonIsoSF, muonIdSF, pileupSF]

if runAnalysis:
    # Setup and run actual analysis
    selector = SelectorTools.SelectorDriver(args['analysis'], args['selection'], args['input_tier'], args['year'])
    #selector.setNumCores(args['numCores'])
    selector.setOutputfile(fileName)
    selector.setInputs(sf_inputs)
    selector.setMaxEntries(args['maxEntries'])
    selector.setNtupeType("NanoAOD")
    if args['filenames']:
        selector.setDatasets(args['filenames'])
    else:
        selector.setFileList(*args['inputs_from_file'])


    mc = selector.applySelector()

    selector.outputFile().Close()
    #OutputTools.addMetaInfo(ROOT.TFile(fileName, "UPDATE"))


# Create end files
fOut = ROOT.TFile(fileName, "UPDATE")

combinedHist = makeCompositeHists(fOut, "AllData", ConfigureJobs.getListOfFilesWithXSec(args["filenames"], args['analysis']), channels)
for i in combinedHist:
    i.Print()
OutputTools.writeOutputListItem(combinedHist, fOut)
combinedHist.Clear()
# alldata = makeCompositeHists("AllData", ConfigureJobs.getListOfFilesWithXSec(["WZxsec2016data"]), channels)
# OutputTools.writeOutputListItem(alldata, fOut)
# allewk = makeCompositeHists("AllEWK", ConfigureJobs.getListOfFilesWithXSec(
#     ConfigureJobs.getListOfEWKFilenames()), channels, False)
# OutputTools.writeOutputListItem(allewk, fOut)
# allnonprompt = makeCompositeHists("NonpromptMC", ConfigureJobs.getListOfFilesWithXSec(
#     ConfigureJobs.getListOfNonpromptFilenames()), channels)
# OutputTools.writeOutputListItem(allnonprompt, fOut)
# final = HistTools.getDifference(fOut, "DataEWKCorrected", "AllData", "AllEWK", getRatios)
# OutputTools.writeOutputListItem(final, fOut)


