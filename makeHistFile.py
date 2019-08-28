#!/usr/bin/env python
import ROOT
from python import SelectorTools
from python import UserInput
from python import OutputTools
from python import ConfigureJobs
from python import HistTools
import os
import sys
import datetime
channels = ["eeee", "eemm","mmee","mmmm"]
#channels = ["eeee"]
def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--proof", "-p", 
        action='store_true', help="Don't use proof")
    parser.add_argument("--lumi", "-l", type=float,
        default=41.5, help="luminosity value (in fb-1)")
    parser.add_argument("--output_file", "-o", type=str,
        default="", help="Output file name")
    parser.add_argument("--scalefactors_file", "-sf", type=str,
        default="", help="ScaleFactors file name")
    parser.add_argument("--leptonSelections", "-ls", type=str,
        default="TightLeptons", help="Either All Loose or Tight")
    parser.add_argument("--test", action='store_true',
        help="Run test job (no background estimate)")
    parser.add_argument("--output_selection", type=str,
        default="", help="Selection stage of output file "
        "(Same as input if not give)")
    parser.add_argument("-b", "--hist_names", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["all"], help="List of histograms, "
                        "as defined in ZZ4lRun2DatasetManager, separated "
                        "by commas")
    return vars(parser.parse_args())

def getHistExpr(hist_names, selection):
    info = ROOT.TList()
    info.SetName("histinfo")
    for hist_name in hist_names:
        bin_info = ConfigHistTools.getHistBinInfo(manager_path, selection, hist_name)
        if "TH1" in ConfigHistTools.getHistType(manager_path, selection, hist_name):
            bin_expr = "{nbins}, {xmin}, {xmax}".format(**bin_info)
        else:
            bin_expr = "{nbinsx}, {xmin}, {xmax}, {nbinsy}, {ymin}, {ymax}".format(**bin_info)
        info.Add(ROOT.TNamed(hist_name, " $ ".join([hist_name, bin_expr])))
    return info

ROOT.gROOT.SetBatch(True)

args = getComLineArgs()
manager_path = ConfigureJobs.getManagerPath()
sys.path.append("/".join([manager_path, 
    "ZZ4lRun2DatasetManager", "Utilities/python"]))
import ConfigHistTools 

today = datetime.date.today().strftime("%d%b%Y")
tmpFileName = "HistFiles/Hists%s-%s.root" % (today, args['leptonSelections']) if args['output_file'] == "" \
        else args['output_file']
fOut = ROOT.TFile(tmpFileName, "recreate")

fScales = ROOT.TFile(args['scalefactors_file'])
mZZTightFakeRate = fScales.Get("mZZTightFakeRate")
eZZTightFakeRate = fScales.Get("eZZTightFakeRate")
#fakeRateFile=ROOT.TFile('data/fakeRate18Dec2018-ZplusLSkimZZ2018Data.root')
#eZZTightFakeRate=fakeRateFile.Get('DataEWKCorrected/ratioE2D_allE')
#mZZTightFakeRate=fakeRateFile.Get('DataEWKCorrected/ratioMu2D_allMu')
#useSvenjasFRs = False
#useJakobsFRs = False
#if useSvenjasFRs:
#    mZZTightFakeRate = fScales.Get("mZZTightFakeRate_Svenja")
#    eZZTightFakeRate = fScales.Get("eZZTightFakeRate_Svenja")
#elif useJakobsFRs:
#    mZZTightFakeRate = fScales.Get("mZZTightFakeRate_Jakob")
#    eZZTightFakeRate = fScales.Get("eZZTightFakeRate_Jakob")
## For medium muons
##mZZMedFakeRate.SetName("fakeRate_allMu")
if mZZTightFakeRate:
    #print "Yes muon fake rates"
    mZZTightFakeRate.SetName("fakeRate_allMu")
if eZZTightFakeRate:
    #print "Yes electron fake rates"
    eZZTightFakeRate.SetName("fakeRate_allE")
#
#muonIsoSF = fScales.Get('muonIsoSF')
muonMoriond19SF= fScales.Get('muonMoriond19SF')
#Electron Reco SF (POG) https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2018DataRecommendations
electronLowReco19SF = fScales.Get('electronLowReco19SF')
electronReco19SF = fScales.Get('electronReco19SF')
#Electron ID + SIP HZZ Scale factors
electronMoriond19SF = fScales.Get('electronMoriond19SF')
electronMoriond19GapSF = fScales.Get('electronMoriond19GapSF')
pileupSF = fScales.Get('pileupSF')

fr_inputs = [eZZTightFakeRate, mZZTightFakeRate,]
sf_inputs = [electronLowReco19SF,electronReco19SF,electronMoriond19SF, electronMoriond19GapSF,muonMoriond19SF,pileupSF]
#sf_inputs=[pileupSF]
selection = args['selection'].replace("LooseLeptons", "") \
    if args['output_selection'] == "" else args['output_selection'].split("_")[0]
if selection == "":
    selection = "LooseLeptons"
    print "Info: Using BasicZZSelections for hist defintions"
analysis = "/".join([args['analysis'], selection])
hists = ConfigHistTools.getAllHistNames(manager_path, analysis) \
    if "all" in args['hist_names'] else args['hist_names']

hists = [h for h in hists if "unrolled" not in h and "wCR" not in h and h not in  ["YieldByChannel", "CutFlow"]]
#print "hists: ", hists
hist_inputs = [getHistExpr(hists, analysis)]
#print "hist_inputs: ",hist_inputs
tselection = [ROOT.TNamed("selection", args['output_selection'])]
#print "tselection: ",tselection
if args['proof']:
    ROOT.TProof.Open('workers=12')

if "FakeRate" not in args['output_selection'] and not args['test']:
   # background = SelectorTools.applySelector(["ZZ4l2018data"] +
   #     ConfigureJobs.getListOfEWK() + ["wz3lnu-powheg"] +
   #     ConfigureJobs.getListOfNonpromptFilenames(),channels, 
   #         "WZBackgroundSelector", args['selection'], fOut, 
   #         extra_inputs=sf_inputs+fr_inputs+hist_inputs+tselection, 
   #         addSumweights=False,
   #         proof=args['proof'])
    background = SelectorTools.applySelector([args['analysis']+"data"] +
            ConfigureJobs.getListOfEWKFilenames()+ConfigureJobs.getListOfDYFilenames()+ConfigureJobs.getListOfNonpromptFilenames()+ConfigureJobs.getListOfHZZFilenames(),channels, 
            "ZZBackgroundSelector", args['selection'], fOut,args['analysis'],
            extra_inputs=sf_inputs+fr_inputs+hist_inputs+tselection, 
            addSumweights=False,
            proof=args['proof'])
mc = SelectorTools.applySelector(args['filenames'],channels, "ZZSelector", args['selection'], fOut,args['analysis'],
        extra_inputs=sf_inputs+hist_inputs+tselection, 
        addSumweights=True, proof=args['proof'])
if args['test']:
    exit(0)

alldata = HistTools.makeCompositeHists(fOut,"AllData", 
    ConfigureJobs.getListOfFilesWithXSec([args['analysis']+"data"], manager_path), args['lumi'],
    underflow=False, overflow=False)
OutputTools.writeOutputListItem(alldata, fOut)
alldata.Delete()

#nonpromptmc = HistTools.makeCompositeHists(fOut, "NonpromptMC", ConfigureJobs.getListOfFilesWithXSec( 
#    ConfigureJobs.getListOfNonpromptFilenames(), manager_path), args['lumi'],
#    underflow=False, overflow=False)
#nonpromptmc.Delete()
#
#OutputTools.writeOutputListItem(nonpromptmc, fOut)

#allggZZmc= HistTools.makeCompositeHists(fOut,"ggZZ", ConfigureJobs.getListOfFilesWithXSec(
#    ConfigureJobs.getListOfggZZFilenames(), manager_path), args['lumi'],
#    underflow=False, overflow=False)
#OutputTools.writeOutputListItem(allggZZmc, fOut)
#allggZZmc.Delete()
#
#allSignal= HistTools.makeCompositeHists(fOut,"HZZ-signal", ConfigureJobs.getListOfFilesWithXSec(
#    ConfigureJobs.getListOfHZZFilenames(), manager_path), args['lumi'],
#    underflow=False, overflow=False)
#OutputTools.writeOutputListItem(allSignal, fOut)
#allSignal.Delete()

ewkmc = HistTools.makeCompositeHists(fOut,"AllEWK", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWK(), manager_path), args['lumi'],
    underflow=False, overflow=False)
OutputTools.writeOutputListItem(ewkmc, fOut)
ewkmc.Delete()

ewkcorr = HistTools.getDifference(fOut, "DataEWKCorrected", "AllData", "AllEWK")
OutputTools.writeOutputListItem(ewkcorr, fOut)
ewkcorr.Delete()
