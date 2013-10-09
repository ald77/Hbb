/* Example of how to use the cutflow class*/
#include "cutflow.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TChain.h"


typedef unsigned int uint;
using namespace std;

int main(int argc, char* argv[]) {
	
	// Option 1: print the cutflow from vectors of root files
	vector<TFile*> ttbar, qcd, wt, diboson, znn, SM, data;

	ttbar.push_back(new TFile("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_SyncSkim.root","read"));
	ttbar.push_back(new TFile("raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_SyncSkim.root","read"));
    ttbar.push_back(new TFile("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_SyncSkim.root","read"));
    ttbar.push_back(new TFile("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71_SyncSkim.root","read"));
    ttbar.push_back(new TFile("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71_SyncSkim.root","read"));
    ttbar.push_back(new TFile("raw_plots_and_values/TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71_SyncSkim.root","read"));
    
	qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1897_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1898_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1899_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1900_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1901_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1902_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1903_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1904_v71_SyncSkim.root","read"));
    qcd.push_back(new TFile("raw_plots_and_values/QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1905_v71_SyncSkim.root","read"));

	
	wt.push_back(new TFile("raw_plots_and_values/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71_SyncSkim.root","read"));
	wt.push_back(new TFile("raw_plots_and_values/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71_SyncSkim.root","read"));
    wt.push_back(new TFile("raw_plots_and_values/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859_v71_SyncSkim.root","read"));
  
    znn.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71_SyncSkim.root","read"));
    znn.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71_SyncSkim.root","read"));
    znn.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71_SyncSkim.root","read"));
    znn.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71_SyncSkim.root","read"));
    znn.push_back(new TFile("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71_SyncSkim.root","read"));
    
    diboson.push_back(new TFile("raw_plots_and_values/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71_SyncSkim.root","read"));
    diboson.push_back(new TFile("raw_plots_and_values/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71_SyncSkim.root","read"));
    diboson.push_back(new TFile("raw_plots_and_values/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71_SyncSkim.root","read"));
    
    for (uint f = 0; f < 20; f++) { // fill the SM vector, save a few lines of code
    	if (f < qcd.size()) SM.push_back(qcd[f]);
    	if (f < ttbar.size()) SM.push_back(ttbar[f]);
    	if (f < wt.size()) SM.push_back(wt[f]);
    	if (f < znn.size()) SM.push_back(znn[f]);
    	if (f < diboson.size()) SM.push_back(ttbar[f]);
    }
    
    data.push_back(new TFile("raw_plots_and_values/MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71_SyncSkim.root","read"));
    data.push_back(new TFile("raw_plots_and_values/MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71_SyncSkim.root","read"));
    data.push_back(new TFile("raw_plots_and_values/MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71_SyncSkim.root","read"));
    data.push_back(new TFile("raw_plots_and_values/MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71_SyncSkim.root","read"));
    data.push_back(new TFile("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71_SyncSkim.root","read"));
    data.push_back(new TFile("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71_SyncSkim.root","read"));    
    // buggy! diboson.push_back(new TFile("raw_plots_and_values/ZH","read"));	
    
    /*Cutflow* ttbar_cutflow = new Cutflow(ttbar);
    ttbar_cutflow->PrintCSV();
    Cutflow* qcd_cutflow = new Cutflow(qcd);
    qcd_cutflow->PrintCSV();
    Cutflow* wt_cutflow = new Cutflow(wt);
    wt_cutflow->PrintCSV();
    Cutflow* znn_cutflow = new Cutflow(znn);
    znn_cutflow->PrintCSV();
    Cutflow* diboson_cutflow = new Cutflow(diboson);
    diboson_cutflow->PrintCSV();
    Cutflow* SM_cutflow = new Cutflow(SM);
    SM_cutflow->PrintCSV();*/
    
    // Option 2: use a list of files specified by the 'i' flag
    char in_filename[1024]="";
	char c(' ');
	bool fileList(false);
	while((c=getopt(argc, argv, "i:"))!=-1){
	switch(c){
	case 'i':
	  sprintf(in_filename, "%s", optarg);
	  break;
	default:
	  break;
	}
	}
	if (strlen(in_filename)>0) {
		// check to see if input is a text file or just a path to a root file
		char* txt = strstr(in_filename,".txt");
		char* list = strstr(in_filename,".list");
		cout << txt << " " << list << endl;
		if (txt||list) fileList = true;
		cout << fileList << endl;
		cout << in_filename << endl;
		Cutflow* in_cutflow = new Cutflow(in_filename,fileList);
		in_cutflow->Print();  
	}
    
	return 0;
}
