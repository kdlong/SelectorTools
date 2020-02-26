from python import HistTools, OutputTools, ConfigureJobs
import ROOT
import array
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", type=str,
    default="test.root", help="Output file name")
args = parser.parse_args()

saveToFile = True
adddRjj = True
input_file_name = args.input_file
input_file = ROOT.TFile(input_file_name, "update" if saveToFile else "read")

variations = [i for x in ["CMS_scale_j", "CMS_res_j", \
    "CMS_eff_m", "CMS_scale_m", "CMS_eff_e", "CMS_scale_e", "CMS_pileup", "CMS_scale_unclEnergy"] for i in [x+"Up", x+"Down"]]
jeVariations = [i for x in ["CMS_scale_j", "CMS_res_j"] for i in [x+"Up", x+"Down"]]

transformed_mjj_etajj_hists = HistTools.getTransformedHists(input_file, 
        ConfigureJobs.getListOfFiles(ConfigureJobs.getListOfEWKFilenames() + \
            ConfigureJobs.getListOfNonpromptFilenames() + \
            ['wzjj-vbfnlo-sf', 'wzjj-vbfnlo-of', ] + \
                ['wz3lnu-mg5amcnlo','wz3lnu-powheg', 'zz4l-mg5amcnlo'] + \
                ['AllData', 'WZxsec2016data', 'DataEWKCorrected'], 
            input_tier='Wselection'),
        ["mjj_etajj_2D_%s" % c for c in ConfigureJobs.getChannels()] + \
        ["mjj_etajj_2D_Fakes_%s" % c for c in ConfigureJobs.getChannels()] + \
        ["mjj_etajj_2D_%s_%s" % (var, c) for c in ConfigureJobs.getChannels()
                for var in variations] + \
        ["mjj_etajj_2D_%s_Fakes_%s" % (var, c) for c in ConfigureJobs.getChannels()
                for var in variations],
        HistTools.makeUnrolledHist, ConfigureJobs.get2DBinning()
)
transformed_mjj_dRjj_hists = []
if adddRjj:
    transformed_mjj_dRjj_hists = HistTools.getTransformedHists(input_file, 
            ConfigureJobs.getListOfFiles(ConfigureJobs.getListOfEWKFilenames() + \
                ['wzjj-vbfnlo-sf', 'wzjj-vbfnlo-of', ] + \
                    ['wz3lnu-mg5amcnlo','wz3lnu-powheg', 'zz4l-mg5amcnlo'] + \
                    ['AllData', 'WZxsec2016data', 'DataEWKCorrected'], 
                input_tier='Wselection'),
            ["mjj_dRjj_2D_%s" % c for c in ConfigureJobs.getChannels()] + \
            ["mjj_dRjj_2D_Fakes_%s" % c for c in ConfigureJobs.getChannels()] + \
            ["mjj_dRjj_2D_%s_%s" % (var, c) for c in ConfigureJobs.getChannels()
                    for var in variations] + \
            ["mjj_dRjj_2D_%s_Fakes_%s" % (var, c) for c in ConfigureJobs.getChannels()
                    for var in variations],
            HistTools.makeUnrolledHist, ConfigureJobs.get2DBinning("mjj", "dRjj")
    )


transformed_hists = transformed_mjj_etajj_hists+transformed_mjj_dRjj_hists
addControlRegion=True
unrolled_hists_wcontrol = []
if addControlRegion:
    for folder in transformed_hists:
        new_folder = folder.Clone()
        unrolled_hists_wcontrol.append(new_folder)
        for h in folder:
            append = h.GetName().split("_")
            append = append[append.index("unrolled")+1:]
            hist_name = "_".join(["backgroundControlYield"] + append)
            control_hist = input_file.Get(folder.GetName() + "/" +hist_name)
            if not control_hist:
                raise RuntimeError("failed to find hist %s for folder %s" % (hist_name, folder))
            hist = HistTools.addControlRegionToFitHist(control_hist, h)
            new_folder.Add(hist)
    transformed_hists += unrolled_hists_wcontrol

#transformed_mjj_mtwz_hists = HistTools.getTransformedHists(input_file, 
#        ConfigureJobs.getListOfFiles(ConfigureJobs.getListOfEWKFilenames() + \
#            ['wzjj-vbfnlo-sf', 'wzjj-vbfnlo-of', ] + \
#                ['wz3lnu-mg5amcnlo','wz3lnu-powheg', 'zz4l-mg5amcnlo'] + \
#                ['AllData', 'WZxsec2016data', 'DataEWKCorrected'], 
#             'Wselection'),
#        ["mjj_mtwz_2D_%s" % c for c in ConfigureJobs.getChannels()] + \
#            ["mjj_mtwz_2D_Fakes_%s" % c for c in ConfigureJobs.getChannels()],
#        HistTools.makeUnrolledHist, [
#            array.array('d', [500, 1000,1500, 2000, 2500]),
#            [0, 150, 300, 450]
#        ]
#)

if saveToFile:
    for output in transformed_hists:
        OutputTools.writeOutputListItem(output, input_file)
