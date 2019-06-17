#!/usr/bin/env python
import ROOT
from python import SelectorTools
from python import UserInput
from python import OutputTools
from python import ConfigureJobs
from python import HistTools
import os
import sys

def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--lumi", "-l", type=float,
        default=35.87, help="luminosity value (in fb-1)")
    parser.add_argument("--output_file", "-o", type=str,
        default="test.root", help="Output file name")
    parser.add_argument("--test", action='store_true',
        help="Run test job (no background estimate)")
    parser.add_argument("--uwvv", action='store_true',
        help="Use UWVV format ntuples in stead of NanoAOD")
    parser.add_argument("--noHistConfig", action='store_true',
        help="Don't rely on config file to specify hist info")
    parser.add_argument("-j", "--numCores", type=int, default=1,
        help="Number of cores to use (parallelize by dataset)")
    parser.add_argument("--input_tier", type=str,
        default="", help="Selection stage of input files")
    parser.add_argument("-c", "--channels", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["eee","eem","emm"], help="List of channels"
                        "separated by commas. NOTE: set to Inclusive for NanoAOD")
    parser.add_argument("-b", "--hist_names", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["all"], help="List of histograms, "
                        "as defined in AnalysisDatasetManager, separated "
                        "by commas")
    return vars(parser.parse_args())

ROOT.gROOT.SetBatch(True)

args = getComLineArgs()
manager_path = ConfigureJobs.getManagerPath()
if manager_path not in sys.path:
    sys.path.insert(0, "/".join([manager_path, 
        "AnalysisDatasetManager", "Utilities/python"]))

tmpFileName = args['output_file']
fOut = ROOT.TFile(tmpFileName, "recreate")

addScaleFacs = False
if args['analysis'] == "WZxsec2016" or args['analysis'] == 'Zstudy_2016':
    addScaleFacs = True
addScaleFacs=False
sf_inputs = [ROOT.TParameter(bool)("applyScaleFacs", False)]

if addScaleFacs:
    fScales = ROOT.TFile('data/scaleFactors.root')
    mCBTightFakeRate = fScales.Get("mCBTightFakeRate")
    eCBTightFakeRate = fScales.Get("eCBTightFakeRate")
    useSvenjasFRs = False
    useJakobsFRs = False
    if useSvenjasFRs:
        mCBTightFakeRate = fScales.Get("mCBTightFakeRate_Svenja")
        eCBTightFakeRate = fScales.Get("eCBTightFakeRate_Svenja")
    elif useJakobsFRs:
        mCBTightFakeRate = fScales.Get("mCBTightFakeRate_Jakob")
        eCBTightFakeRate = fScales.Get("eCBTightFakeRate_Jakob")
    # For medium muons
    #mCBMedFakeRate.SetName("fakeRate_allMu")
    if mCBTightFakeRate:
        mCBTightFakeRate.SetName("fakeRate_allMu")
    if eCBTightFakeRate:
        eCBTightFakeRate.SetName("fakeRate_allE")

    muonIsoSF = fScales.Get('muonIsoSF')
    muonIdSF = fScales.Get('muonTightIdSF')
    electronTightIdSF = fScales.Get('electronTightIdSF')
    electronGsfSF = fScales.Get('electronGsfSF')
    pileupSF = fScales.Get('pileupSF')

    #fPrefireEfficiency = ROOT.TFile('data/Map_Jet_L1FinOReff_bxm1_looseJet_JetHT_Run2016B-H.root')
    fPrefireEfficiency = ROOT.TFile('data/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root')
    prefireEff = fPrefireEfficiency.Get('prefireEfficiencyMap')

    fr_inputs = [eCBTightFakeRate, mCBTightFakeRate,]
    sf_inputs = [electronTightIdSF, electronGsfSF, muonIsoSF, muonIdSF, pileupSF, prefireEff]
    sf_inputs.append(ROOT.TParameter(bool)("applyScaleFacs", True))

if args['input_tier'] == '':
    args['input_tier'] = args['selection']
selection = args['selection'].split("_")[0]

if selection == "Inclusive2Jet":
    selection = "Wselection"
    print "Info: Using Wselection for hist defintions"
analysis = "/".join([args['analysis'], selection])
hists, hist_inputs = UserInput.getHistInfo(analysis, args['hist_names'], args['noHistConfig'])


#if "WZxsec2016" in analysis and "FakeRate" not in args['output_selection'] and not args['test']:
#    background = SelectorTools.applySelector(["WZxsec2016data"] +
#        ConfigureJobs.getListOfEWKFilenames() + ["wz3lnu-powheg"] +
#        ConfigureJobs.getListOfNonpromptFilenames(), 
#            "WZBackgroundSelector", args['selection'], fOut, 
#            extra_inputs=sf_inputs+fr_inputs+hist_inputs+tselection, 
#            channels=channels,
#            addSumweights=False,
#            nanoAOD=nanoAOD,
#            parallel=args['parallel'],
#            )

selector = SelectorTools.SelectorDriver(args['analysis'], args['selection'], args['input_tier'])
selector.setOutputfile(fOut.GetName())
selector.setInputs(sf_inputs+hist_inputs)
selector.setNtupeType("UWVV" if args['uwvv'] else "NanoAOD")
if args['uwvv']:
    selector.setChannels(args['channels'])
selector.setNumCores(args['numCores'])

mc = selector.applySelector(args['filenames'])

if args['test']:
    fOut.Close()
    sys.exit(0)

alldata = HistTools.makeCompositeHists(fOut,"AllData", 
    ConfigureJobs.getListOfFilesWithXSec(["WZxsec2016data"], manager_path), args['lumi'],
    underflow=False, overflow=False)
OutputTools.writeOutputListItem(alldata, fOut)
alldata.Delete()

nonpromptmc = HistTools.makeCompositeHists(fOut, "NonpromptMC", ConfigureJobs.getListOfFilesWithXSec( 
    ConfigureJobs.getListOfNonpromptFilenames(), manager_path), args['lumi'],
    underflow=False, overflow=False)
nonpromptmc.Delete()

OutputTools.writeOutputListItem(nonpromptmc, fOut)
ewkmc = HistTools.makeCompositeHists(fOut,"AllEWK", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWKFilenames(), manager_path), args['lumi'],
    underflow=False, overflow=False)
OutputTools.writeOutputListItem(ewkmc, fOut)
ewkmc.Delete()

ewkcorr = HistTools.getDifference(fOut, "DataEWKCorrected", "AllData", "AllEWK")
OutputTools.writeOutputListItem(ewkcorr, fOut)
ewkcorr.Delete()
