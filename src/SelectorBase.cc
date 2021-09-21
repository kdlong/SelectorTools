#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include <boost/algorithm/string.hpp>
#include <TStyle.h>
#include <regex>
#include "TParameter.h"

void SelectorBase::Begin(TTree * /*tree*/)
{
    TString option = GetOption();
}

void SelectorBase::SlaveBegin(TTree * /*tree*/)
{
    if (GetInputList() != nullptr) {
        TParameter<bool>* applyScaleFactors = (TParameter<bool>*) GetInputList()->FindObject("applyScaleFacs");
        if (applyScaleFactors != nullptr && applyScaleFactors->GetVal()) {
           SetScaleFactors();
        }
    }
}

void SelectorBase::Init(TTree *tree)
{
    if (!tree) return;
    fChain = tree;
    
    TString option = GetOption();

    if (GetInputList() != nullptr) {
	TNamed* ntupleType = (TNamed *) GetInputList()->FindObject("ntupleType");
    TNamed* name = (TNamed *) GetInputList()->FindObject("name");
    TNamed* chan = (TNamed *) GetInputList()->FindObject("channel");
    TNamed* selection = (TNamed *) GetInputList()->FindObject("selection");
	TNamed* year = (TNamed *) GetInputList()->FindObject("year");

	if (ntupleType != nullptr) {
	    std::string ntupleName = ntupleType->GetTitle();
	    if (ntupleName == "NanoAOD")
		    ntupleType_ = NanoAOD;
	    else if (ntupleName  == "UWVV")
		ntupleType_ = UWVV;
	    else if (ntupleName  == "Bacon")
		ntupleType_ = Bacon;
	    else
		    throw std::invalid_argument("Unsupported ntuple type!");
	}
	else {
	    std::cerr << "INFO: Assuming NanoAOD ntuples" << std::endl;
	    ntupleType_ = NanoAOD;
	}

        if (name != nullptr) {
            name_ = name->GetTitle();
        }
        else {
            name_ = GetNameFromFile();
        }
        if (name_.empty()){
            std::cerr << "INFO: Using default name \"Unknown\" for file" << std::endl;
            name_ = "Unknown";
        }

	if(year != nullptr) {
	    year_ = yearMap_[year->GetTitle()];
	}
	
	if (chan != nullptr) {
	    channelName_ = chan->GetTitle();
	}
	else if (ntupleType_ == UWVV)
            channelName_ = fChain->GetTree()->GetDirectory()->GetName();
        if (selection != nullptr) {
            selectionName_ = selection->GetTitle();
        }
    }
    
    if (selectionMap_.find(selectionName_) != selectionMap_.end()) {
	    selection_ = selectionMap_[selectionName_];
    }
    else
        throw std::invalid_argument(selectionName_ + " is not a valid selection!");
    
    isMC_ = false;
    if (name_.find("data") == std::string::npos){
        isMC_ = true;
    }
    if (doSystematics_ && isMC_ && !isNonprompt_)
        variations_.insert(systematics_.begin(), systematics_.end());

    if (channelMap_.find(channelName_) != channelMap_.end())
        channel_ = channelMap_[channelName_];
    else {
        std::string message = "Invalid channel choice! ";
        message += "Choice was " + channelName_ + "\n";
        message += "Valid choices are: ";
        for (const auto& chan : channelMap_)
            message += chan.first + ", " ;
        throw std::invalid_argument(message);
    }

    // only make the directory iff class isn't being run as a slave class /////
    TNamed* isSlaveClass = (TNamed *) GetInputList()->FindObject("isSlaveClass");
    if(isSlaveClass != nullptr) return;
    
    makeOutputDirs();
    SetBranches();
}

void SelectorBase::addSubprocesses(std::vector<std::string> processes) {
    subprocesses_ = processes;
}

void SelectorBase::setSubprocesses(std::string process) {
    currentHistDir_ = dynamic_cast<TList*>(fOutput->FindObject(process.c_str()));
    if (currentHistDir_ == nullptr) {
        currentHistDir_ = new TList();
        currentHistDir_->SetName(process.c_str());
        fOutput->Add(currentHistDir_);
    }
        //throw std::invalid_argument(process + " is not a valid subprocess for selector!");
}

void SelectorBase::makeOutputDirs() {
    std::vector<std::string> allProcesses = {name_};
    allProcesses.insert(allProcesses.begin(), subprocesses_.begin(), subprocesses_.end());
    
    for (auto& name : allProcesses) { 
        currentHistDir_ = dynamic_cast<TList*>(fOutput->FindObject(name.c_str()));
        
        if ( currentHistDir_ == nullptr ) {
            currentHistDir_ = new TList();
            currentHistDir_->SetName(name.c_str());
            fOutput->Add(currentHistDir_);
            if (name == name_)
                SetupNewDirectory();
        }
    }
    //TODO: Determine why this fails and if it's really needed
    //size_t existingObjectPtrsSize = allObjects_.size();
    // Watch for something that I hope never happens,
    //if ( existingObjectPtrsSize > 0 && allObjects_.size() != existingObjectPtrsSize ) {
    //    throw std::invalid_argument(Form("SelectorBase: Size of allObjects has changed!: %lu to %lu", existingObjectPtrsSize, allObjects_.size()));
    //}
    //UpdateDirectory();
}

void SelectorBase::SetScaleFactors() {
    std::invalid_argument("No scale factors defined for selector!");
}

void SelectorBase::SetBranches() {
    if (ntupleType_ == UWVV)
        SetBranchesUWVV();
    else if (ntupleType_ == NanoAOD)
        SetBranchesNanoAOD();
    else if (ntupleType_ == Bacon)
        SetBranchesBacon();
    else 
        throw std::invalid_argument("Undefined ntuple type!");
        
}

void SelectorBase::LoadBranches(Long64_t entry, SystPair variation) {
    if (ntupleType_ == NanoAOD)
        LoadBranchesNanoAOD(entry, variation);
    else if (ntupleType_ == Bacon)
        LoadBranchesBacon(entry, variation);
    else {
        LoadBranchesUWVV(entry, variation);
    }
}

Bool_t SelectorBase::Process(Long64_t entry)
{
    for (const auto& variation : variations_) {
        LoadBranches(entry, variation);
        FillHistograms(entry, variation);
    }
    return kTRUE;
}

Bool_t SelectorBase::Notify()
{
    return kTRUE;
}

float SelectorBase::GetPrefiringEfficiencyWeight(
        std::vector<float>* jetPt, std::vector<float>* jetEta) {
    float prefire_weight = 1;
    for (size_t i = 0; i < jetPt->size(); i++) {
        float jPt = jetPt->at(i);
        float jEta = std::abs(jetEta->at(i));
        prefire_weight *= (1-prefireEff_->GetEfficiency(prefireEff_->FindFixBin(jEta, jPt)));
    }
    return prefire_weight;
}

void SelectorBase::Terminate()
{
}
    
void SelectorBase::SlaveTerminate()
{
}
void SelectorBase::UpdateDirectory()
{
  for(TNamed** objPtrPtr : allObjects_) {
    if ( *objPtrPtr == nullptr ) std::invalid_argument("SelectorBase: Call to UpdateObject but existing pointer is null");
    *objPtrPtr = (TNamed *) currentHistDir_->FindObject((*objPtrPtr)->GetName());
    if ( *objPtrPtr == nullptr ) std::invalid_argument("SelectorBase: Call to UpdateObject but current directory has no instance");
  }
}

template<typename T>
void SelectorBase::InitializeHistMap(std::vector<std::string>& labels, std::unordered_map<HistLabel, T*>& histMap) {
    for (auto& label : labels) {
        if (channel_ != Inclusive) {
            HistLabel histlabel = {label, channel_, Central};
            histMap[histlabel] = {};
        }
        else {
            for (auto& chan : allChannels_) {
                HistLabel histlabel = {label, chan.first, Central};
                histMap[histlabel] = {};
            }
        }
    }
}

void SelectorBase::InitializeHistogramsFromConfig() {
    TList* histInfo = (TList *) GetInputList()->FindObject("histinfo");
    if (histInfo == nullptr ) 
        throw std::domain_error("Can't initialize histograms without passing histogram information to TSelector");

    InitializeHistMap(hists1D_,histMap1D_);
    InitializeHistMap(hists2D_,histMap2D_);
    InitializeHistMap(hists3D_,histMap3D_);
    InitializeHistMap(weighthists1D_, weighthistMap1D_);
    InitializeHistMap(weighthists2D_, weighthistMap2D_);

    std::vector<std::string> tempSystHistNames;
    for (auto& syst : systematics_) {
        for (auto hist : systHists_) {
            tempSystHistNames.push_back(hist + "_" + syst.second);
        }
        for (auto hist : systHists2D_) {
            tempSystHistNames.push_back(hist + "_" + syst.second);
        }
    }

    for (auto && entry : *histInfo) {  
        TNamed* currentHistInfo = dynamic_cast<TNamed*>(entry);
        const char* name = currentHistInfo->GetName();
        std::vector<std::string> histData = ReadHistDataFromConfig(currentHistInfo->GetTitle());
        
        std::vector<ChannelPair> channels = {{channel_, channelName_}};
        if (channel_ == Inclusive) {
            channels = allChannels_;
        }

        for (const auto& chan : channels) {
            HistLabel centralLabel = {name, chan.first, Central};
            if (histMap1D_.find(centralLabel) != histMap1D_.end() || histMap2D_.find(centralLabel) != histMap2D_.end()
                                                                        || histMap3D_.find(centralLabel) != histMap3D_.end()) { 
                InitializeHistogramFromConfig(name, chan, histData);
            }
            //Don't print out a ton of annoying errors if it's a syst hist
            else {
                if (std::find(tempSystHistNames.begin(), tempSystHistNames.end(), name) == tempSystHistNames.end())
                    break;
                size_t idx = centralLabel.name.find_last_of("_");
                HistLabel tmplabel = {idx != std::string::npos ? centralLabel.name.substr(0, idx) : "", chan.first, Central};
                if (idx == std::string::npos || (histMap1D_.find(tmplabel) == histMap1D_.end() && histMap2D_.find(tmplabel) == histMap2D_.end()))
                    std::cerr << "Skipping invalid histogram '" << name << "'" << std::endl;
            }
        }
    }

    for (auto& subprocess : subprocesses_) {
        setSubprocesses(subprocess);
        auto& subprocessMap = subprocessHistMaps1D_[subprocess];
        subprocessMap = {};
        for (auto& hist : histMap1D_) {
            subprocessMap[hist.first] = {};
            if (hist.second == nullptr)
                continue;
            AddObject<TH1D>(subprocessMap[hist.first], const_cast<TH1D&>(*hist.second));
            //subprocessWeightHistMaps1D_[subporcess] = weighthistMap1D_;
        }
    }

    currentHistDir_ = dynamic_cast<TList*>(fOutput->FindObject(name_.c_str()));
}

void SelectorBase::InitializeHistogramFromConfig(std::string name, ChannelPair channel, std::vector<std::string>& histData) {
    if (histData.size() != 4 && histData.size() != 7 && histData.size() != 10) {
        std::cerr << "Malformed data string for histogram '" << name
                    << ".' Must have form: 'Title; (optional info) $ nbins, xmin, xmax'"
                    << "\n   OR form: 'Title; (optional info) $ nbins, xmin, xmax nbinsy ymin ymax'"
                    << std::endl;
        exit(1);
    }
    
    int nbins = std::stoi(histData.at(1));
    float xmin = std::stof(histData.at(2));
    float xmax = std::stof(histData.at(3));
    
    bool is1D = histData.size() == 4;
    int nbinsy = is1D ? 1 : std::stoi(histData[4]);
    float ymin = is1D ? 1. : std::stof(histData[5]);
    float ymax = is1D ? 1. : std::stof(histData[6]);

    bool is3D = histData.size() == 10;
    int nbinsz = !is3D ? 1 : std::stoi(histData[7]);
    float zmin = !is3D ? 1. : std::stof(histData[8]);
    float zmax = !is3D ? 1. : std::stof(histData[9]);

    for (auto& variation : variations_) {
        std::string histName = getHistName(name, variation.second, channel.second);
        HistLabel histlabel = {name, channel.first, variation.first};

        if (variation.first == Central) {
            if (is1D) {
                AddObject<TH1D>(histMap1D_[histlabel], histName.c_str(), histData.at(0).c_str(),
                        nbins, xmin, xmax);
            }
            // For now, only allow 3D without variations
            else if (is3D) {
                // Ths is obviously a hack, but it's too much trouble to support this generally
                if (nbins == 4 && xmin == 5 && xmax == 1000) {
                    float xrange[nbins+1]  = {5., 75., 80, 85., 1000,};
                    float yrange[nbinsy+1];
                    float zrange[nbinsz+1];
                    for (int i = 0; i <= nbinsy; i++) {
                        yrange[i] = ymin + (ymax - ymin)/nbinsy*i;
                    }
                    for (int i = 0; i <= nbinsz; i++) {
                        zrange[i] = zmin + (zmax - zmin)/nbinsz*i;
                    }

                    AddObject<TH3D>(histMap3D_[histlabel], histName.c_str(), histData.at(0).c_str(), 
                            nbins, xrange, nbinsy, yrange, nbinsz, zrange);
                }
                else
                    AddObject<TH3D>(histMap3D_[histlabel], histName.c_str(), histData.at(0).c_str(), 
                            nbins, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
            }
            else {
                AddObject<TH2D>(histMap2D_[histlabel], histName.c_str(), histData.at(0).c_str(), 
                        nbins, xmin, xmax, nbinsy, ymin, ymax);
            }
        }
        else if ((is1D && std::find(systHists_.begin(), systHists_.end(), name) != systHists_.end()) ||
                    (!is1D && !is3D && std::find(systHists2D_.begin(), systHists2D_.end(), name) != systHists2D_.end()) ||
                    (is3D && std::find(systHists3D_.begin(), systHists3D_.end(), name) != systHists3D_.end())) {
            if (is1D) {
                histMap1D_[histlabel] = {};
                AddObject<TH1D>(histMap1D_[histlabel], histName.c_str(), 
                    histData.at(0).c_str(), nbins, xmin, xmax);
            }
            else if (is3D) {
                histMap3D_[histlabel] = {};
                // Ths is obviously a hack, but it's too much trouble to support this generally
                if (((nbins == 4 && xmin == 5) || (nbins == 5 && xmin == 50)) && xmax == 1000) {
                    std::vector<float> xrange = {50., 75., 85., 90., 95., 1000.};
                    if (nbins == 4)
                        xrange = {5., 75., 80, 85., 1000,};
                    float yrange[nbinsy+1];
                    float zrange[nbinsz+1];
                    for (int i = 0; i <= nbinsy; i++) {
                        yrange[i] = ymin + (ymax - ymin)/nbinsy*i;
                    }
                    for (int i = 0; i <= nbinsz; i++) {
                        zrange[i] = zmin + (zmax - zmin)/nbinsz*i;
                    }

                    AddObject<TH3D>(histMap3D_[histlabel], histName.c_str(), histData.at(0).c_str(), 
                            nbins, xrange.data(), nbinsy, yrange, nbinsz, zrange);
                }
                else
                    AddObject<TH3D>(histMap3D_[histlabel], histName.c_str(), histData.at(0).c_str(), 
                            nbins, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
            }
            else {
                histMap2D_[histlabel] = {};
                AddObject<TH2D>(histMap2D_[histlabel], histName.c_str(), 
                    histData.at(0).c_str(), nbins, xmin, xmax, nbinsy, ymin, ymax);
            }
        }
        bool isPtVvar = (variation.first == ptV0to3 || variation.first == ptV0to3_lhe ||
                variation.first == ptV3to5 || variation.first == ptV3to5_lhe ||
                variation.first == ptV5to7 || variation.first == ptV5to7_lhe ||
                variation.first == ptV7to9 || variation.first == ptV7to9_lhe ||
                variation.first == ptV9to12 || variation.first == ptV9to12_lhe ||
                variation.first == ptV12to15 || variation.first == ptV12to15_lhe ||
                variation.first == ptV15to20 || variation.first == ptV15to20_lhe ||
                variation.first == ptV20to27 || variation.first == ptV20to27_lhe ||
                variation.first == ptV27to40 || variation.first == ptV27to40_lhe ||
                variation.first == ptV40toInf || variation.first == ptV40toInf_lhe);

        if (std::find(theoryVarSysts_.begin(), theoryVarSysts_.end(), variation.first) != theoryVarSysts_.end()) {
            // Should change this to support SCETlib weights at some point
            int nWeights = isPtVvar ? 10 : nTheoryWeights_;

            size_t pos = variation.first == Central ? name.size() : (name.size() + variation.second.size()+1);
            std::string weighthistName = histName.insert(pos, "_lheWeights");
            if (is1D && std::find(weighthists1D_.begin(), weighthists1D_.end(), name) != weighthists1D_.end()) {
                weighthistMap1D_[histlabel] = {};
                AddObject<TH2D>(weighthistMap1D_[histlabel], 
                    weighthistName.c_str(), histData[0].c_str(),
                    nbins, xmin, xmax, nWeights, 0, nWeights);
            }
            // 3D weight hists must be subset of 2D hists!
            else if (std::find(weighthists2D_.begin(), weighthists2D_.end(), name) != weighthists2D_.end()) {
                weighthistMap2D_[histlabel] = {};
                AddObject<TH3D>(weighthistMap2D_[histlabel], 
                    weighthistName.c_str(), histData[0].c_str(),
                    nbins, xmin, xmax, nbinsy, ymin, ymax, nWeights, 0, nWeights);

            }
        }
    }
}

std::vector<std::string> SelectorBase::ReadHistDataFromConfig(std::string histDataString) {
    std::vector<std::string> histData;
    boost::split(histData, histDataString, boost::is_any_of("$"));
    std::vector<std::string> binInfo;
    if (histData.size() != 2)
        return {};
    
    boost::split(binInfo, histData[1], boost::is_any_of(","));
   
    histData.pop_back();
    for (const auto& x : binInfo) {
        histData.push_back(x);
    }
    
    return histData;
}

void SelectorBase::SetupNewDirectory()
{
}

std::string SelectorBase::concatenateNames(const std::string& baseName, std::string& toAppend) {
    return concatenateNames(baseName.c_str(), toAppend);
}

std::string SelectorBase::concatenateNames(const std::string& baseName, const char* toAppend) {
    std::string app = toAppend;
    return concatenateNames(baseName.c_str(), toAppend);
}

std::string SelectorBase::concatenateNames(const char* baseName, std::string& toAppend) {
    if (toAppend.empty())
        return baseName;
    const std::string delimit = "_";
    std::string name = baseName;
    name.reserve(name.size()+delimit.size()+toAppend.size());
    name.append(delimit);
    name.append(toAppend);
    return name;
}

std::string SelectorBase::getHistName(std::string histName, std::string variationName) {
    return getHistName(histName, variationName, channelName_);
}

std::string SelectorBase::getHistName(std::string histName, std::string variationName, std::string channel) {
    histName.append("_");
    if (!variationName.empty()) {
        histName.append(variationName);
        histName.append("_");
    }
    histName.append(channel);
    return histName;
}
