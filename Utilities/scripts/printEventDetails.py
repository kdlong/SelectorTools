import ROOT
import os
from python import ConfigureJobs,ApplySelection
from python.prettytable import PrettyTable
import argparse
import datetime

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def getEventSelectionExpr(path, comparison, channel):
    events = []
    
    file_name = "%s/%s_%s.txt" % (path, comparison, channel)
    print file_name
    with open(file_name, "r") as events_file:
        for line in events_file.readlines():
            event_num = line.strip().split(":")
            if len(event_num) == 3:
                events.append("(run == {0} && lumi == {1} && evt == {2})".format(*event_num))
            elif len(event_num) == 1:
                events.append("evt == {0}".format(*event_num))
    return events

parser = argparse.ArgumentParser() 
parser.add_argument("-f", "--filelist", 
                    type=lambda x : [i.strip() for i in x.split(',')],
                    required=True, help="List of input file names "
                    "to be processed (separated by commas)")
parser.add_argument("-s", "--selection", required=True)
parser.add_argument("--output_selection", required=False, default="")
parser.add_argument("-m", "--cut_string", required=False, type=str,
                    default="")
parser.add_argument("-c", "--channels", required=False, type=str,
                    default="eee,eem,emm,mmm")
parser.add_argument("-o", "--output_file", required=False, type=str,
                    default=".")
args = parser.parse_args()

isfile = any(os.path.isfile(name) or os.path.exists(name.rstrip("/*")) 
                for name in args.filelist)
filelist = ConfigureJobs.getListOfFiles(args.filelist, input_tier=args.selection) if \
    not isfile else args.filelist
print filelist
states = [x.strip() for x in args.channels.split(",")]
file_paths = []
output_dir = os.path.dirname(args.output_file)
for name in filelist:
    if not isfile:
        try:
            label = name
            file_paths.append((label, ConfigureJobs.getInputFilesPath(name, 
                args.selection, "WZxsec2016")))
        except ValueError as e:
            print e
            continue
    else:
        label = filelist[0].split(".")[0].split("/")[-1]
        file_paths.append((label, name))
    try:
        os.makedirs(output_dir + "/" + label)
    except OSError as e:
        print e
        pass

for state in states:
    first = True
    for file_path in file_paths:
        state = state.strip()
        chain = ROOT.TChain("%s/ntuple" % state)
        if "store" in file_path[1][:6]:
            path = "root://cmsxrootd.hep.wisc.edu/" + file_path[1]
        else:
            path = file_path[1]
        chain.Add(path)
        print "Path is", path
        ApplySelection.setAliases(chain, state, "Cuts/WZxsec2016/aliases.json")
        run_expr = "run:lumi:evt"
        trig_expr = "singleEPass:singleMuPass:singleIsoMuPass:singleESingleMuPass:singleMuSingleEPass:doubleMuPass:doubleMuDZPass:doubleEPass:(singleESingleMuPass || singleMuSingleEPass || doubleMuPass || doubleMuDZPass || doubleEPass || singleEPass || singleMuPass || singleIsoMuPass)"
        filter_expr = "Flag_duplicateMuonsPass:Flag_BadPFMuonFilterPass:Flag_badMuonsPass:Flag_goodVerticesPass"
        all_filters = "Flag_BadChargedCandidateFilterPass:Flag_HBHENoiseFilterPass:Flag_HBHENoiseIsoFilterPass:Flag_BadPFMuonFilterPass:Flag_EcalDeadCellTriggerPrimitiveFilterPass:Flag_goodVerticesPass:Flag_globalTightHalo2016FilterPass:Flag_eeBadScFilterPass:Flag_duplicateMuonsPass:Flag_badMuonsPass"
        filter_expr = all_filters
        filter_expr += ":metFiltersData" if "data_" in file_path[0] else ":metFiltersMC"
        lepid_expr = "Zlep1IsLoose:Zlep1IsTight:Zlep1_Pt:Zlep1_Eta:Zlep1_Phi:Zlep1_PVDXY:Zlep1_PVDZ:" \
                            "Zlep2IsLoose:Zlep2IsTight:Zlep2_Pt:Zlep2_Eta:Zlep2_Phi:Zlep2_PVDXY:Zlep2_PVDZ:" \
                            "WlepIsLoose:WlepIsTight:Wlep_Pt:Wlep_Eta:Wlep_Phi:Wlep_PVDXY:Wlep_PVDZ"
        other_kinematics = "Mass:ZMass:type1_pfMETEt:Max$(jetCSVv2)"
        jet_vars = "mjj:dEtajj:jetPt[0]:jetEta[0]:jetPt[1]:jetEta[1]:zep3l"
        veto_expr = "nCBVIDHLTSafeElec:nWZMediumMuon:nCBVIDTightElec:nWZTightMuon"
        scan_expr = ":".join([run_expr,other_kinematics,lepid_expr,veto_expr,jet_vars,filter_expr,trig_expr])

        outfile_name = "/".join([output_dir, file_path[0], args.output_file.split("/")[-1].replace("chan", state)])
        events = getEventSelectionExpr("/eos/user/k/kelong/WZAnalysisData/SyncWithJakob/Differences_2018Apr/DataLooseControl/",
                #"KennethNotJakob",
                "MPnotWisc",
        #        "WiscnotMP",
                state
        )
        print events
        for i, eventlist in enumerate(chunks(events, 40)):
            if i > 0:
                outfile_name = outfile_name.replace(".txt", "%i.txt" % i)
            evtstring = " || ".join(eventlist)
            print chain, chain.GetEntries()
            print "Output file", outfile_name
            chain.GetPlayer().SetScanRedirect(True)
            chain.GetPlayer().SetScanFileName(outfile_name)
            chain.Scan(scan_expr, evtstring,"colsize=30")
            with open(args.output_file.replace("chan", state), "w" if first else "a") as allfiles:
                if not os.path.isfile(outfile_name):
                    continue
                with open(outfile_name) as subfile:
                    lines = subfile.readlines()
                    if not first:
                        lines = lines[3:-1]
                    allfiles.write("\n".join(lines))
                    first = False
