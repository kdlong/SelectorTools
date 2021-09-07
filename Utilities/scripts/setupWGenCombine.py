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
parser.add_argument("--pdfs", type=str, default="NNPDF31", help="List of PDFs to store")
parser.add_argument("--outFolder", type=str, default="/data/shared/{user}",
    help="Output folder")
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

folder_name = "_".join([args.fitvar,args.append]) if args.append != "" else args.fitvar

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
        #isAltTh = "lhe" in args.fitvar or "prefsr" in args.fitvar
        #isAltTh = "allth" not in process
        if args.pdfs != "none":
            cenidx=19 # 9+9 scale vars
            # NNPDF3.1
            nsets = 103
            if "nnpdf31" in args.pdfs.lower():
                cardtool.addTheoryVar(process, 'pdf_hessian', range(cenidx, cenidx+nsets-2), central=0, specName="NNPDF31")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="NNPDF31_alphas")
                # NNPDF31_nnlo_as_0118_CMSW1_hessian_100; LHAPDFID = 325700, AltSet5
            if "cmsw1" in args.pdfs.lower():
                # Offset by 1 central value (NNPDF3.1)
                cenidx += 1
                cardtool.addTheoryVar(process, 'pdf_hessian', range(cenidx, cenidx+nsets-2), central=0, specName="CMSW1")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="CMSW1_alphas")
                # NNPDF31_nnlo_as_0118_CMSW2_hessian_100; LHAPDFID = 325900, AltSet6
            if "cmsw2" in args.pdfs.lower():
                cenidx += 2
                cardtool.addTheoryVar(process, 'pdf_hessian', range(cenidx, cenidx+nsets-2), central=0, specName="CMSW2")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="CMSW2_alphas")
                # NNPDF31_nnlo_as_0118_CMSW3_hessian_100; LHAPDFID = 326100, AltSet7
            if "cmsw3" in args.pdfs.lower():
                cenidx += 3
                cardtool.addTheoryVar(process, 'pdf_hessian', range(cenidx, cenidx+nsets-2), central=0, specName="CMSW3")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="CMSW3_alphas")
                # NNPDF31_nnlo_as_0118_CMSW3_hessian_100; LHAPDFID = 326300, AltSet8
            if "cmsw4" in args.pdfs.lower():
                cenidx += 4
                cardtool.addTheoryVar(process, 'pdf_hessian', range(cenidx, cenidx+nsets-2), central=0, specName="CMSW4")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="CMSW4_alphas")
                # NNPDF30_nnlo_as_0118_CMSW3_hessian_100; LHAPDFID = 303200, AltSet9
            if "nnpdf30" in args.pdfs.lower():
                cenidx += 5
                cardtool.addTheoryVar(process, 'pdf_hessian', range(cenidx, cenidx+nsets-2), central=0, specName="NNPDF30")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="NNPDF30")
                # CT18, LHAPDF ID = 14000
            if args.pdfs.lower() == "ct18":
                nsets=61
                cenidx += 6
                cardtool.addTheoryVar(process, 'pdf_assymhessian', range(cenidx, cenidx+nsets-2), central=0, specName="CT18")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="CT18")
                # CT18Z, LHAPDF ID = 14000
            if "ct18z" in args.pdfs.lower():
                nsets=61
                cenidx = cenidx+6+61
                cardtool.addTheoryVar(process, 'pdf_assymhessian', range(cenidx, cenidx+nsets-2), central=0, specName="CT18Z")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="CT18Z")
                # MMHT
            if "mmht" in args.pdfs.lower():
                cenidx += 7
                cardtool.addTheoryVar(process, 'pdf_assymhessian', range(cenidx, cenidx+nsets-2), central=0, specName="MMHT")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="MMHT")
                # HERA20_EIG, LHAPDF=61200
            if "hera" in args.pdfs.lower():
                cenidx += 8
                nsets = 30
                cardtool.addTheoryVar(process, 'pdf_assymhessian', range(cenix, cenidx+nsets-2), central=0, specName="HERA2")
                cardtool.addTheoryVar(process, 'other', range(cenidx+nsets-2, cenidx+nsets), central=0, specName="HERA2")
        isAltTh = True
        cenMassIdx = 919 if not isAltTh else 18+nsets+11
        massVars = lambda i: [cenMassIdx+i, cenMassIdx-i]
        cardtool.addTheoryVar(process, 'other', massVars(0), exclude=[], central=0, specName="massShift0MeV")
        cardtool.addTheoryVar(process, 'other', massVars(1), exclude=[], central=0, specName="massShift10MeV")
        cardtool.addTheoryVar(process, 'other', massVars(2), exclude=[], central=0, specName="massShift20MeV")
        cardtool.addTheoryVar(process, 'other', massVars(3), exclude=[], central=0, specName="massShift30MeV")
        cardtool.addTheoryVar(process, 'other', massVars(5), exclude=[], central=0, specName="massShift50MeV")
        cardtool.addTheoryVar(process, 'other', massVars(10), exclude=[], central=0, specName="massShift100MeV")
        # This is broken for now
        # width = (18+890+21+3) if not isAltTh else (18+nsets+21+3)
        # cardtool.addTheoryVar(process, 'other', [width, width], exclude=[], central=0, specName="width2043")
        if "N3LLCorr" in process:
            resumIdx = cenMassIdx+11
            cardtool.addTheoryVar(process, 'resumscale', range(resumIdx, resumIdx+45), central=0)
            cardtool.addTheoryVar(process, 'resumscaleDFO', range(resumIdx, resumIdx+3), central=0)
            cardtool.addTheoryVar(process, 'resumscaleDLambda', range(resumIdx+3, resumIdx+5), central=-1)
            cardtool.addTheoryVar(process, 'resumscaleDMatch', range(resumIdx+5, resumIdx+9), central=-1)
            cardtool.addTheoryVar(process, 'resumscaleDResum', range(resumIdx+9, resumIdx+45), central=-1)

    elif "nnlops" in process:
        cardtool.addTheoryVar(process, 'scale', range(10, 19), exclude=[15, 17], central=0)
        cardtool.setScaleVarGroups(process, [(3,6), (1,2), (4,8)])
        if args.pdfs != "none":
            # NNPDF3.1
            cardtool.addTheoryVar(process, 'pdf_hessian', range(885, 986), central=0, specName="NNPDF31")
    elif process not in ["nonprompt", "data"]:
        cardtool.addTheoryVar(process, 'scale', range(1, 10), exclude=[3, 7], central=4)
        if args.pdfs != "none":
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

if args.pdfs != "none":
    if args.allHessianVars:
        # TODO: More sets of course
        pdflabel = "NNPDF" if "nnpdf" in args.pdfs else "CT18"
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
        }
    )
