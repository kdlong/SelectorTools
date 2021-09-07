#!/usr/bin/env python3
from python import ConfigureJobs,CombineCardTools,UserInput
import sys
import ROOT
import logging
import array
import argparse
import os

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("--debug", action='store_true',
    help="Print debug info")
parser.add_argument("--mc2hes", action='store_true',
    help="Convert MC errors to hessian")
parser.add_argument("-c", "--central", type=str, default="wlnu_jetbinned_nlo_cp5",
    help="Sample to use as central value")
parser.add_argument("--files", type=lambda x: [i.strip() for i in x.split(",")], 
    default=[], help="Samples to add to output file")
parser.add_argument("-d", "--data", type=str, default="wlnu_nlo",
    help="Sample to use as dummy data")
parser.add_argument("-a", "--append", type=str, default="",
    help="Append to output folder name")
parser.add_argument("-b", "--fitvar", type=str, default="ptWmet",
    help="Variable to use in the fit")
parser.add_argument("-f", "--input_file", type=str, required=True,
    help="Input hist file")
parser.add_argument("--channels", type=lambda x: x.split(','), 
    default=["mp","mn"], help="List of channels (separated by comma)")
parser.add_argument("-l", "--lumi", type=float, 
    default=35.9*0.7, help="lumi")
parser.add_argument("--pdf", type=str, default="nnpdf31", help="PDF to store",
    choices=["nnpdf31", "cmsw1", "cmsw2", "cmsw3", "cms4", "ct18", "ct18z", "mmht", "all"])
parser.add_argument("--outFolder", type=str, default="/data/shared/{user}",
    help="Output folder")
parser.add_argument("--storePdfCenValues", action='store_true', 
    help="Store all central PDF variations")
parser.add_argument("--splitPtV", action='store_true', 
    help="Don't split scale uncertainties by pt(V)")
parser.add_argument("--allHessianVars", action='store_true', 
    help="store all hessian variations")
parser.add_argument("--addEff", action='store_true', 
    help="add dummy efficiency uncertainties")
parser.add_argument("--theoryOnly", action='store_true', 
    help="Only add theory uncertainties")
parser.add_argument("--scetlibUnc", action='store_true', 
    help="Use SCETlib scale uncertainties")
parser.add_argument("-r", "--rebin", 
                    type=str, default=None, help="Rebin array: "
                    "values (bin edges) separated by commas.")
args = parser.parse_args()

logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO))

cardtool = CombineCardTools.CombineCardTools()
cardtool.setCorrelateScaleUnc(False)

manager_path = ConfigureJobs.getManagerPath() 
sys.path.append("/".join([manager_path, "AnalysisDatasetManager",
    "Utilities"]))

from configTools import ConfigHistFactory
config_factory = ConfigHistFactory.ConfigHistFactory(
    "%s/AnalysisDatasetManager" % manager_path,
    "WGen/NanoAOD",
)

#plot_groups = ["wlnu_lo", "wlnu_lo_cp5", "wlnu_nlo", "wlnu_jetbinned_nlo", "wlnu_jetbinned_nlo_cp5", ]
plot_groups = args.files if args.files else ["wpmunu_minnlo_nnlopslike_photos", "wpmunu_nnlops_photos", "wpmunu_nnlops_nlow"]
plotGroupsMap = {name : config_factory.getPlotGroupMembers(name) for name in plot_groups}

xsecs  = ConfigureJobs.getListOfFilesWithXSec([f for files in plotGroupsMap.values() for f in files])

if args.rebin and "unrolled" not in args.fitvar:
    if ":" in args.rebin:
        bins = UserInput.getRebin(args.rebin)
        args.rebin = array.array('d', bins)
    elif "," in args.rebin:
        args.rebin = array.array('d', [float(i.strip()) for i in args.rebin.split(",")])
    else:
        args.rebin = int(args.rebin)
    cardtool.setRebin(args.rebin)

cardtool.setFitVariable(args.fitvar)
if "unrolled" in args.fitvar:
    #cardtool.setUnrolled([-2.5+0.2*i for i in range(0,26)], range(26, 56, 1))
    cardtool.setUnrolled([-2.4+0.2*i for i in range(0,25)], range(26, 57, 1))
    #cardtool.setUnrolled([-2.5+0.5*i for i in range(0,11)], range(26, 56, 3))
cardtool.setProcesses(plotGroupsMap)
cardtool.setChannels(args.channels)
print("Channels are", args.channels)
cardtool.setCrosSectionMap(xsecs)

variations = [] if args.theoryOnly else ["CMS_scale_m", "CMS_res_m"] 
cardtool.setVariations(variations)
normVariations = [] if args.theoryOnly else ["mWBWShift100MeV", "mWBWShift50MeV"] 
cardtool.setNormalizedVariations(normVariations)

folder_name = "_".join([args.fitvar,args.pdf]) 
if args.append != "":
    folder_name += "_%s" % args.append

basefolder = args.outFolder.format(user=os.getlogin())
cardtool.setOutputFolder(basefolder+"/CombineStudies/WGen/%s" % folder_name)

cardtool.setLumi(args.lumi)
cardtool.setInputFile(args.input_file)
cardtool.setOutputFile("WGenCombineInput.root")
if "mp" in args.channels and "mn" in args.channels:
    cardtool.setCombineChannels({"m" : ["mp", "mn"]})
    args.channels = args.channels+["m"]
cardtool.setRemoveZeros(False)
cardtool.setAddOverflow(False)

ptbins = [0,3,5,7,9,12,15,20,27,40,100]
ptbinPairs = [(x,y) for x,y in zip(ptbins[:-1], ptbins[1:])]

firstPdfIdx = 10
pdfIdxMap = {
        "nnpdf31" : {
            "name" : "NNPDF31",
            "unc" : "pdf_hessian",
            "cenidx" : firstPdfIdx,
            "nsets" : 103,
        },
        "cmsw1" : {
            "name" : "NNPDF31CMSW1",
            "unc" : "pdf_hessian",
            "cenidx" : firstPdfIdx+1*args.storePdfCenValues,
            "nsets" : 103,
        },
        "cmsw2" : {
            "name" : "NNPDF31CMSW2",
            "unc" : "pdf_hessian",
            "cenidx" : firstPdfIdx+2*args.storePdfCenValues,
            "nsets" : 103,
        },
        "cmsw3" : {
            "name" : "NNPDF31CMSW3",
            "unc" : "pdf_hessian",
            "cenidx" : firstPdfIdx+3*args.storePdfCenValues,
            "nsets" : 103,
        },
        "cmsw4" : {
            "name" : "NNPDF31CMSW4",
            "unc" : "pdf_hessian",
            "cenidx" : firstPdfIdx+4*args.storePdfCenValues,
            "nsets" : 103,
        },
        "cmsw4" : {
            "name" : "NNPDF30",
            "unc" : "pdf_hessian",
            "cenidx" : firstPdfIdx+5*args.storePdfCenValues,
            "nsets" : 103,
        },
        "ct18" : {
            "name" : "CT18",
            "cenidx" : firstPdfIdx+6*args.storePdfCenValues,
            "unc" : "pdf_assymhessian",
            "nsets" : 61,
        },
        "ct18z" : {
            "name" : "CT18Z",
            "cenidx" : firstPdfIdx+6*args.storePdfCenValues+61,
            "unc" : "pdf_assymhessian",
            "nsets" : 61,
        },
        "mmht" : {
            "name" : "MMHT",
            "cenidx" : firstPdfIdx+7*args.storePdfCenValues,
            "unc" : "pdf_assymhessian",
            "nsets" : 51,
        },
        "hera" : {
            "name" : "HERA",
            "cenidx" : firstPdfIdx+8*args.storePdfCenValues,
            "unc" : "pdf_assymhessian",
            "nsets" : 51,
        },
}
for process in plot_groups:
    if "matrix" in process:
        cardtool.setVariations(variations+["QCDscale_"+process])
    #Turn this back on when the theory uncertainties are added
    if "minnlo" in process:
        # If using the NanoGen
        #cardtool.addTheoryVar(process, 'scale', range(1, 10) , exclude=[6, 8], central=0)
        cardtool.addTheoryVar(process, 'scale', range(1,19)[::2], exclude=[2, 6], central=4)
        cardtool.setScaleVarGroups(process, [(1,7), (3,5), (0,8)])
        # NNPDF3.0 scale unc
        # cardtool.addTheoryVar(process, 'scale', range(10, 19), exclude=[6, 8], central=0, specName="NNPDF30")
        if args.pdf != "none" and "all" not in args.pdf:
            info = pdfIdxMap[args.pdf]
            firstAlphaIdx = info["cenidx"]+info["nsets"]-2
            indices = range(info["cenidx"], firstAlphaIdx)
            alphaIndices = range(firstAlphaIdx, firstAlphaIdx+2)
            cardtool.addTheoryVar(process, info["unc"], indices, central=0, specName=info["name"])
            cardtool.addTheoryVar(process, 'other', alphaIndices, central=0, specName=info["name"]+"_alphas")
        if args.storePdfCenValues:
            for label, pdfInfo in pdfIdxMap.items():
                cardtool.addTheoryVar(process, 'other', [pdfInfo["cenidx"]], central=0, specName=pdfInfo["name"]+"Cen")

        cenMassIdx = firstPdfIdx+8*args.storePdfCenValues+pdfIdxMap[args.pdf]["nsets"]+11-1
        print("central is", cenMassIdx)
        massVars = lambda i: [cenMassIdx+i, cenMassIdx-i]
        cardtool.addTheoryVar(process, 'other', massVars(0), exclude=[], central=0, specName="massShift0MeV")
        cardtool.addTheoryVar(process, 'other', massVars(1), exclude=[], central=0, specName="massShift10MeV")
        cardtool.addTheoryVar(process, 'other', massVars(2), exclude=[], central=0, specName="massShift20MeV")
        cardtool.addTheoryVar(process, 'other', massVars(3), exclude=[], central=0, specName="massShift30MeV")
        cardtool.addTheoryVar(process, 'other', massVars(5), exclude=[], central=0, specName="massShift50MeV")
        # Width weights are broken for now
        # width = (18+890+21+3) if not isAltTh else (18+nsets+21+3)
        # cardtool.addTheoryVar(process, 'other', [width, width], exclude=[], central=0, specName="width2043")
        if "N3LLCorr" in process:
            resumIdx = cenMassIdx+11
            cardtool.addTheoryVar(process, 'resumscale', range(resumIdx, resumIdx+45), central=0)
            cardtool.addTheoryVar(process, 'resumscaleDFO', range(resumIdx, resumIdx+3), central=0)
            cardtool.addTheoryVar(process, 'resumscaleDLambda', range(resumIdx+3, resumIdx+5), central=-1)
            cardtool.addTheoryVar(process, 'resumscaleDMatch', range(resumIdx+5, resumIdx+9), central=-1)
            cardtool.addTheoryVar(process, 'resumscaleDResum', range(resumIdx+9, resumIdx+45), central=-1)
    elif process not in ["nonprompt", "data"]:
        cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[3, 7], central=4)
        if args.pdf != "none":
            cardtool.addTheoryVar(process, 'pdf_mc' if "cp5" not in process else "pdf_hessian", range(10,111), central=0)

    if args.splitPtV:
        for pair in ptbinPairs:
            varName = 'ptV%ito%i' % pair
            varName = varName.replace("100", "Inf")
            cardtool.addScaleBasedVar(process, varName) 
    if process in args.central.split(",") and args.addEff:
        cardtool.addPerBinVariation(process, "CMS_eff_m", 0.0025, False)

    cardtool.loadHistsForProcess(process, expandedTheory=args.allHessianVars)
    cardtool.writeProcessHistsToOutput(process)

scriptdir = os.path.dirname(os.path.realpath(__file__)) 
path = "/".join(scriptdir.split("/")[:-2]+["Templates", "CombineCards", "VGen"])

nnu = 2
if args.splitPtV:
    nnu += cardtool.addCustomizeCard(path+"/Customize/PtV_template.txt")
if not args.theoryOnly:
    nnu += cardtool.addCustomizeCard(path+"/Customize/muscale_template.txt")
    cardtool.addCardGroup("CMS_scale_m group = CMS_scale_m")
    cardtool.addCardGroup("CMS_res_m group = CMS_res_m")

if args.pdf != "none":
    if args.allHessianVars:
        # TODO: More sets of course
        pdflabel = "NNPDF" if "nnpdf" in args.pdf else "CT18"
        npdfs = cardtool.addCustomizeCard(f"{path}/Customize/pdf{pdflabel}_template.txt")
        nnu += npdfs
        cardtool.addCardGroup("pdf group = %s" % " ".join(["pdf%i" % i for i in range(1,npdfs+1)]))
    else:
        nnu += cardtool.addCustomizeCard(path+"/Customize/pdf_template.txt")

cen = args.central
if not args.scetlibUnc:
    nnu += cardtool.addCustomizeCard(path+"/Customize/scale_template.txt")
    cardtool.addCardGroup(f"QCDscale group = QCDscale_muR_{cen} QCDscale_muF_{cen} QCDscale_muRmuF_{cen}")
else:
    nnu += cardtool.addCustomizeCard(path+"/Customize/scetlibscale_template.txt")
    cardtool.addCardGroup(f"QCDscale group = resumscaleDLambda_{cen} resumscaleDFO_{cen} resumscaleDMatch_{cen} resumscaleDResum_{cen}")

cardtool.addCardGroup("massnoi noiGroup = massShift100MeV")

nuissance_map = {"mn" : nnu, "mp" : nnu, "m" : nnu}
for i, chan in enumerate(args.channels):
    data = args.data if "," not in args.data else args.data.split(",")[i]
    central = args.central if "," not in args.central else args.data.split(",")[i]
    cardtool.setTemplateFileName("%s/WGen_template_{channel}.txt" % path)
    logging.info("Writting cards for channel %s" % chan)
    cardtool.writeCards(chan, nuissance_map[chan], 
        extraArgs={"data_name" : data, 
            "w_sample" : central, "w_yield" : "yield:%s" % central, 
            "pdfName" : pdfIdxMap[args.pdf]["name"]+"Hes"
        }
    )
