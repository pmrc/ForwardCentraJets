// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 09 Jul 2013
//
// run_analysis()
// compiles the analysis routines and runs them to do the analysis

void run_piotr_analysis()
{
  gROOT -> Reset();
  gROOT->SetStyle("Plain");

 //------------ FWLite libraries ------
  gSystem->Load("libFWCoreFWLite.so");
  gSystem->Load("libKKousourQCDAnalysis.so");
  gSystem->Load("../RooUnfold-1.1.1/libRooUnfold.so");
  AutoLibraryLoader::enable();

//verbose level
  bool detail = true;
  bool show_results = true;
  bool show_steps = true;

//enable testing
  bool test = false; //if set to true will get only some entries from histograms and save outputs with extension _test

  bool read_mc_ntuples_gen = false;
  bool read_mc_ntuples_det = false;
  bool read_data_ntuples = true;

//directory configuration
  string output_dir = "../piotr_output/";

// histogram directories
  string histograms_dir = output_dir + "histograms/";
  string raw_data_dir = histograms_dir + "raw_xsec_data/";

//ntuple directories
  string ntuple_storage_path = "root://eoscms//eos/cms/store/user/";
  string ntuple_mc_path = ntuple_storage_path + "ilknur/7TeVMC/CMSSW_4_2_4/2010_MC_with_LowPU/GlobalTag_V16_RemoteGlidein/";
  string ntuple_pythia6_path = ntuple_mc_path + "pythia6/TuneZ2star_pythia6_START42_V16/";
  string ntuple_pythia8_path = ntuple_mc_path + "pythia8/Tune4C_pythia8_START42_V16/";
  string ntuple_data_path = ntuple_storage_path + "ilknur/7TeVdata/CMSSW_4_2_4/2010_DATA_with_QCD_NTuples_with_pTCut_5GeV/";
  string ntuple_trigeff_path = ntuple_storage_path + "ilknur/Trigger_Eff_Ntuple/";

//MC Ntuples locations and luminosities

   //Pythia6_TuneZ2star
   int n_files_p6_z2 = 5;
   string mc_p6_z2[5];
   double lumi_p6_z2[5];
   mc_p6_z2[0] = ntuple_pythia6_path + "QCD_Pt_10to25_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_10to25_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[0] = 0.06727 * 3.651e9/5205128.0;
   mc_p6_z2[1] = ntuple_pythia6_path + "QCD_Pt_25to40_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_25to40_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[1] = 0.27707 * 1.059e8/5194564.0;
   mc_p6_z2[2] = ntuple_pythia6_path + "QCD_Pt_40to80_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_40to80_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[2] = 0.25234 * 1.771e7/5011067.0;
   mc_p6_z2[3] = ntuple_pythia6_path + "QCD_Pt_80to150_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_80to150_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[3] = 0.19052 * 8.756e5/4876756.0;
   mc_p6_z2[4] = "../../../../../../../../work/c/cipriano/public/pythia6_150toinf/QCD_Pt_150toInf_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[4] = 0.14464 * 4.76e4/5710726.0;


   //Pythia8_Tune4C
   int n_files_p8_4c = 4;
   string mc_p8_4c[5];
   double lumi_p8_4c[5];
   mc_p8_4c[0] = ntuple_pythia8_path + "QCD_Pt_10to25_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_10to25_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[0] = 0.07463 * 4.334e9/9695621.0;
   mc_p8_4c[1] = ntuple_pythia8_path + "QCD_Pt_25to40_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_25to40_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[1] = 0.29023 * 1.226e8/9949451.0;
   mc_p8_4c[2] = ntuple_pythia8_path + "QCD_Pt_40to80_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_40to80_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[2] = 0.27368 * 2.021e7/9658287.0;
   mc_p8_4c[3] = ntuple_pythia8_path + "QCD_Pt_80to150_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_80to150_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[3] = 0.21046 * 9.855e5/9117322.0;
   mc_p8_4c[4] = "../../../../../../../../work/p/pgunnell/public/PYTHIA8PEDRO/JetTreePYTHIA8.root";
   lumi_p8_4c[4] = 0.16428 * 5.283e4/9761171.0;


//JetMETTau_2010A
   string jetmettau[3];
   double lumi_jetmettau[3];
   jetmettau[0] = ntuple_data_path + "JetMETTau_Run2010A_Apr21ReReco_v1_AOD/JetMETTau_Run2010A_Apr21ReReco_v1_AOD_ProcessedTree_data.root"; //new
   jetmettau[1] = jetmettau[0];
   jetmettau[2] = jetmettau[0];
   int n_files_jetmettau = 3;
   lumi_jetmettau[0] = 0.013799; //cipriano
   lumi_jetmettau[1] = 0.117223; //cipriano
   lumi_jetmettau[2] = 0.278789; //cipriano
   string out_jetmettau_allvertex = raw_data_dir + "xsec_JetMETTau2010A_allvertex.root";
   string out_jetmettau_allvertex_25gev = raw_data_dir + "xsec_JetMETTau2010A_allvertex_25GeVcut.root";
   string out_jetmettau_allvertex_fwd20gev = raw_data_dir + "xsec_JetMETTau2010A_allvertex_20GeVFwDcut.root";

   //JetMET_2010A
   string jetmet[3];
   double lumi_jetmet[3];
   jetmet[0] = ntuple_data_path + "JetMET_Run2010A_Apr21ReReco_v1_AOD/JETMET_AOD_21Apr_MERGE_ALL_ROOT_ProcessedTree_data.root";
   jetmet[1] = jetmet[0];
   jetmet[2] = jetmet[0];
   int n_files_jetmet = 3;
   //lumi_jetmet[0] = 2.896; // /pb cipriano
   //lumi_jetmet[0] = 0.013799; // walter
   //lumi_jetmet[1] = 0.117223; // walter
   //lumi_jetmet[2] = 0.278789; // walter
   lumi_jetmet[0] = 0.009645; // cipriano
   lumi_jetmet[1] = 0.192895; // cipriano
   lumi_jetmet[2] = 2.869; // cipriano
   string out_jetmet_allvertex = raw_data_dir + "xsec_JetMET2010A_allvertex.root";
   string out_jetmet_allvertex_25gev = raw_data_dir + "xsec_JetMET2010A_allvertex_25GeVcut.root";
   string out_jetmet_allvertex_fwd20gev = raw_data_dir + "xsec_JetMET2010A_allvertex_20GeVFwDcut.root";

//compile the different routines
if (show_steps) { cout << "Compiling code..."<<endl; }
  gROOT -> ProcessLine(".L create_directories.C++");
  gROOT -> ProcessLine(".L common_methods.h++");
  if (read_mc_ntuples_gen || read_mc_ntuples_det || read_data_ntuples) { gROOT -> ProcessLine(".L piotr_analysis.C++"); }


  create_directories("", output_dir, "root");

//read the data Ntuples
if (read_data_ntuples)
{
if (show_steps) { cout << "Reading Data Ntuples..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, raw_data_dir, "root");
if (show_steps) { cout << "JetMETTau"<<endl; }
//piotr_analysis(jetmettau, out_jetmettau_1vertex, lumi_jetmettau, n_files_jetmettau, "DATA", "1vertex", "", "", detail, test);
//piotr_analysis(jetmettau, out_jetmettau_allvertex, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
//piotr_analysis(jetmettau, out_jetmettau_allvertex_25gev, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
//piotr_analysis(jetmettau, out_jetmettau_allvertex_fwd20gev, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
//piotr_analysis(jetmettau, out_jetmettau_up, lumi_jetmettau, n_files_jetmettau, "DATA", "up", "", "", detail, test);
//piotr_analysis(jetmettau, out_jetmettau_down, lumi_jetmettau, n_files_jetmettau, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
//piotr_analysis(jetmet, out_jetmet_1vertex, lumi_jetmet, n_files_jetmet, "DATA", "1vertex", "", "", detail, test);
//piotr_analysis(jetmet, out_jetmet_allvertex, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
//piotr_analysis(jetmet, out_jetmet_allvertex_25gev, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
piotr_analysis(jetmet, out_jetmet_allvertex_fwd20gev, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
//piotr_analysis(jetmet, out_jetmet_up, lumi_jetmet, n_files_jetmet, "DATA", "up", "", "", detail, test);
//piotr_analysis(jetmet, out_jetmet_down, lumi_jetmet, n_files_jetmet, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "Jet"<<endl; }
//piotr_analysis(jet, out_jet_1vertex, lumi_jet, n_files_jet, "DATA", "1vertex", "", "", detail, test);
//piotr_analysis(jet, out_jet_allvertex, lumi_jet, n_files_jet, "DATA", "allvertex", "", "", detail, test);
//piotr_analysis(jet, out_jet_up, lumi_jet, n_files_jet, "DATA", "up", "", "", detail, test);
//piotr_analysis(jet, out_jet_down, lumi_jet, n_files_jet, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "All the Data Ntuples were sucessfully read!"<<endl; }
}


}
