#! /bin/bash

./scripts/make_reduced_tree.exe -c -i WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859_v71 &> logs/make_reduced_tree_1859.log &

./scripts/make_reduced_tree.exe -c -i TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71 &> logs/make_reduced_tree_1883.log &
./scripts/make_reduced_tree.exe -c -i TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71 &> logs/make_reduced_tree_1880.log &
./scripts/make_reduced_tree.exe -c -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71 &> logs/make_reduced_tree_1884.log &
./scripts/make_reduced_tree.exe -c -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext2-v1_AODSIM_UCSB1959_v71 &> logs/make_reduced_tree_1959.log &
./scripts/make_reduced_tree.exe -c -i TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext1-v1_AODSIM_UCSB1962_v71 &> logs/make_reduced_tree_1962.log &
./scripts/make_reduced_tree.exe -c -i TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1850_v71 &> logs/make_reduced_tree_1962.log &
./scripts/make_reduced_tree.exe -c -i TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V7C_FSIM-v2_AODSIM_UCSB1976_v71 &> logs/make_reduced_tree_1962.log &
./scripts/make_reduced_tree.exe -c -i TTJets_WToBC_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1966_v71 &> logs/make_reduced_tree_1966.log &

./scripts/make_reduced_tree.exe -c -i TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71 &> logs/make_reduced_tree_1857.log &
./scripts/make_reduced_tree.exe -c -i TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71 &> logs/make_reduced_tree_1856.log &

./scripts/make_reduced_tree.exe -c -i ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71 &> logs/make_reduced_tree_1887.log &
./scripts/make_reduced_tree.exe -c -i ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71 &> logs/make_reduced_tree_1889.log &
./scripts/make_reduced_tree.exe -c -i ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71 &> logs/make_reduced_tree_1888.log &
./scripts/make_reduced_tree.exe -c -i ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71 &> logs/make_reduced_tree_1891.log &
./scripts/make_reduced_tree.exe -c -i ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71 &> logs/make_reduced_tree_1890.log &

./scripts/make_reduced_tree.exe -c -i ZH_ZToBB_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1868_v71 &> logs/make_reduced_tree_1868.log &
./scripts/make_reduced_tree.exe -c -i WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71 &> logs/make_reduced_tree_1858.log &

./scripts/make_reduced_tree.exe -c -i WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71 &> logs/make_reduced_tree_1874.log &
./scripts/make_reduced_tree.exe -c -i WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1875_v71 &> logs/make_reduced_tree_1875.log &
./scripts/make_reduced_tree.exe -c -i ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71 &> logs/make_reduced_tree_1876.log &

./scripts/make_reduced_tree.exe -c -i TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71 &> logs/make_reduced_tree_1855.log &

#./scripts/make_reduced_tree.exe -c -i QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1897_v71 &> logs/make_reduced_tree_1897.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1898_v71 &> logs/make_reduced_tree_1898.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1899_v71 &> logs/make_reduced_tree_1899.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1900_v71 &> logs/make_reduced_tree_1900.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1901_v71 &> logs/make_reduced_tree_1901.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1902_v71 &> logs/make_reduced_tree_1902.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1903_v71 &> logs/make_reduced_tree_1903.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1904_v71 &> logs/make_reduced_tree_1904.log &
#./scripts/make_reduced_tree.exe -c -i QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1905_v71 &> logs/make_reduced_tree_1905.log &

./scripts/make_reduced_tree.exe -c -i BJets_HT-250To500_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1893_v71 &> logs/make_reduced_tree_1893.log &
./scripts/make_reduced_tree.exe -c -i BJets_HT-500To1000_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1894_v71 &> logs/make_reduced_tree_1894.log &
./scripts/make_reduced_tree.exe -c -i BJets_HT-1000ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1895_v71 &> logs/make_reduced_tree_1895.log &

./scripts/make_reduced_tree.exe -c -i W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71 &> logs/make_reduced_tree_1877.log &
./scripts/make_reduced_tree.exe -c -i W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71 &> logs/make_reduced_tree_1878.log &
./scripts/make_reduced_tree.exe -c -i W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71 &> logs/make_reduced_tree_1879.log &

./scripts/make_reduced_tree.exe -c -i T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71 &> logs/make_reduced_tree_1860.log &
./scripts/make_reduced_tree.exe -c -i T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71 &> logs/make_reduced_tree_1861.log &
./scripts/make_reduced_tree.exe -c -i T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71 &> logs/make_reduced_tree_1862.log &
./scripts/make_reduced_tree.exe -c -i Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71 &> logs/make_reduced_tree_1864.log &
./scripts/make_reduced_tree.exe -c -i Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71 &> logs/make_reduced_tree_1865.log &
./scripts/make_reduced_tree.exe -c -i Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71 &> logs/make_reduced_tree_1866.log &

./scripts/make_reduced_tree.exe -c -i MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71 &> logs/make_reduced_tree_1852.log &
./scripts/make_reduced_tree.exe -c -i MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71 &> logs/make_reduced_tree_1853.log &
./scripts/make_reduced_tree.exe -c -i MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71 &> logs/make_reduced_tree_1854.log &
./scripts/make_reduced_tree.exe -c -i MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71 &> logs/make_reduced_tree_1867.log &
./scripts/make_reduced_tree.exe -c -i MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71 &> logs/make_reduced_tree_1869.log &
./scripts/make_reduced_tree.exe -c -i MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71 &> logs/make_reduced_tree_1870.log &

# for trigger eff
#./scripts/make_reduced_tree.exe -c -i 

wait

echo "Done"
exit 0;
