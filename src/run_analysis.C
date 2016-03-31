// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 15 Jun 2015
//
// run_analysis()
// compiles the analysis routines and runs them to do the analysis

void run_analysis()
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

//steps to be done
  bool get_trigger_efficiencies = false; //not used anymore
  bool compute_trigger_efficiency = false; //not used anymore
  bool compute_trigger_turn_on = false;
  bool compute_trigger_correction_v1 = false; //not used anymore
  bool check_trigger_correction_v1 = false; //not used anymore
  bool compute_trigger_correction_v2 = false; //not used anymore
  bool check_trigger_correction_v2 = false; //not used anymore
  bool apply_trigger_correction = false; //not used anymore
  bool apply_delta_phi_factor = false; // not needed or recomended
  bool plot_final_trigger_efficiency = false;
  bool counting_combination_statistics = false; //need to have root output files set
  bool get_vertex_distribution_v0 = false;
  bool compute_vertex_weights_v1 = false;
  bool get_vertex_distribution_v1 = false;
  bool compute_vertex_weights_v2 = false;
  bool get_vertex_distribution_v2 = false; 
  bool compute_vertex_weights_v3 = false;
  bool get_vertex_distribution_v3 = false;
  bool compute_vertex_weights_v4 = false;
  bool get_vertex_distribution_v4 = false;   // here
  bool compute_vertex_weights_v5 = false;
  bool plot_status_vertex_reweight = false;
  bool compute_resolution = false;
  bool plot_resolution = false;
  bool get_pileup_normalization = false;
  bool read_mc_ntuples_gen = true; //needs p6 for unfolding
  bool read_mc_ntuples_det = true; //needs p6 for unfolding
  bool read_data_ntuples = false; //done
  bool check_cross_section = false; //done
  bool normalize_mc = false; //done
  bool plot_control_dist = false; //done
  bool compute_corrections = false; //done
  bool get_mc_matched_events = false; //done
  bool compute_psba = false; // later if needed: compile latex and check
  bool apply_corrections = false; //might need to run it
  bool create_unfolding_response = false; //done
  bool check_response_matrix = false; //might be done
  bool unfold = false; //here
  bool compare_unfolding_results = false; //here
  bool merge_data = false; //suspeita
  bool compute_model_uncertainty = false;
  bool compute_jes_uncertainty = false;
  bool compute_jer_uncertainty = false;
  bool compute_corr_uncertainty = false;
  bool estimate_combination_systematic = false; //need to change and run
  bool run_stat_vs_unf_unc = false; //need to change and run
  bool compute_uncorr_uncertainty = false;
  bool compute_total_uncertainty = false;
  bool apply_uncertainty = false;
  bool create_ratio = false;  //
  bool do_final_plots = false; //need to run

//directory configuration
  string output_dir = "../output/";

// histogram directories
  string histograms_dir = output_dir + "histograms/";
  string trigger_syst_dir = histograms_dir + "trigger_systematic/";
  string trigger_dir = histograms_dir + "trigger_eff/";
  string trigger_corr_dir = histograms_dir + "trigger_correction/";
  string eff_dir = histograms_dir + "efficiency/";
  string vertex_dir = histograms_dir + "vertex/";
  string vertex_weights_dir = histograms_dir + "vertex_weights/";
  string resolution_dir = histograms_dir + "resolution/";
  string pileup_norm_dir = histograms_dir + "pileup_norm/";
  string mc_det_dir = histograms_dir + "xsec_mc_det/";
  string mc_gen_dir = histograms_dir + "xsec_mc_gen/";
  string mc_norm_dir = histograms_dir + "normalized_mc/";
  string raw_data_dir = histograms_dir + "raw_xsec_data/";
  string corrections_dir = histograms_dir + "corrections/";
  string closure_dir = histograms_dir + "closure_test/";
  string matched_dir = histograms_dir + "mc_matched/";
  string corrected_dir = histograms_dir + "corrected_data/";
  string model_unc_dir = histograms_dir + "model_uncertainty/";
  string jes_unc_dir = histograms_dir + "jes_uncertainty/";
  string jer_unc_dir = histograms_dir + "jer_uncertainty/";
  string corr_unc_dir = histograms_dir + "corr_uncertainty/";
  string uncorr_unc_dir = histograms_dir + "uncorr_uncertainty/";
  string total_unc_dir = histograms_dir + "total_uncertainty/";
  string mc_pred_dir = histograms_dir + "mc_prediction/";
  string ratio_dir = histograms_dir + "ratio_to_data/";
  string response_dir = histograms_dir + "unfolding_response/";
  string unfolded_dir = histograms_dir + "unfolded/";


//plot directories
  string trigger_syst_plots = output_dir + "trigger_systematic/";
  string trigger_plots = output_dir + "trigger_eff/";
  string trigger_corr_plots = output_dir + "../../../../../scratch0/output/trigger_correction/";
  string vertex_weights_plots = output_dir + "vertex_weights/";
  string resolution_plots = output_dir + "../../../../../scratch0/output/resolution/";
  string mc_norm_plots = output_dir + "../../../../../scratch0/output/normalized_mc/";
  string control_plots = output_dir + "control_dist/";
  string corrections_plots = output_dir + "corrections/";
  string sel_corrections_plots = output_dir + "sel_corrections/";
  string closure_plots = output_dir + "closure/";
  string corrected_plots = output_dir + "corrected_data/";
  string psba_plots = output_dir + "psba/";
  string merged_data_plots = output_dir + "merged_data/";
  string merged_uncorrected_data_plots = output_dir + "merged_uncorrected_data/";
  string model_unc_plots = output_dir + "model_uncertainty/";
  string jes_unc_plots = output_dir + "jes_uncertainty/";
  string jer_unc_plots = output_dir + "jer_uncertainty/";
  string corr_unc_plots = output_dir + "corr_uncertainty/";
  string uncorr_unc_plots = output_dir + "uncorr_uncertainty/";
  string stat_vs_unf_unc_plots = output_dir + "stat_vs_unf_unc/";
  string total_unc_plots = output_dir + "total_uncertainty/";
  string data_unc_plots = output_dir + "data_uncertainty/";
  string ratio_plots = output_dir + "ratios_to_data/";
  string correction_ratio_plots = output_dir + "correction_ratios/";
  string final_plots = output_dir + "final_xsec/";
  string check_response_dir = output_dir + "check_response/";
  string unfolding_plots = output_dir + "../../../../../scratch0/output/unfolding/";
  string compare_unfolding_plots = output_dir + "compare_unfolding/";

//ntuple directories
  string ntuple_storage_path = "root://eoscms//eos/cms/store/user/";
  string ntuple_mc_path = ntuple_storage_path + "ilknur/7TeVMC/CMSSW_4_2_4/2010_MC_with_LowPU/GlobalTag_V16_RemoteGlidein/";
  string ntuple_pythia6_path = ntuple_mc_path + "pythia6/TuneZ2star_pythia6_START42_V16/";
  string ntuple_pythia8_path = ntuple_mc_path + "pythia8/Tune4C_pythia8_START42_V16/";
  string ntuple_data_path = ntuple_storage_path + "ilknur/7TeVdata/CMSSW_4_2_4/2010_DATA_with_QCD_NTuples_with_pTCut_5GeV/";
  string ntuple_trigeff_path = ntuple_storage_path + "ilknur/Trigger_Eff_Ntuple/";

  string ntuple_pythia6_v16_path =  "root://eoscms//eos/cms/store/group/phys_diffraction/walterjr/MCOfficial2010LowPU_NewJEC/";

//temp vars
  string prefix;


//MC Ntuples locations and luminosities

   //Pythia6_TuneZ2Star_v16
   int n_files_p6_z2 = 5;
   string mc_p6_z2[5];
   double lumi_p6_z2[5];
   mc_p6_z2[0] = ntuple_pythia6_v16_path + "QCD_Ht-50to100_TuneZ2star_7TeV_madgraph_pythia6_V16.root";
   lumi_p6_z2[0] = 1.0 * 1.741e8/11111394.0;
   mc_p6_z2[1] = ntuple_pythia6_v16_path + "QCD_Ht-100to250_TuneZ2star_7TeV_madgraph_pythia6_V16.root";
   lumi_p6_z2[1] = 1.0 * 4.069e6/10353295.0;
   mc_p6_z2[2] = ntuple_pythia6_v16_path + "QCD_Ht-250to500_TuneZ2star_7TeV_madgraph_pythia6_V16.root";
   lumi_p6_z2[2] = 1.0 * 1.978e5/5602739.0;
   mc_p6_z2[3] = ntuple_pythia6_v16_path + "QCD_Ht-500to1000_TuneZ2star_7TeV_madgraph_pythia6_V16.root";
   lumi_p6_z2[3] = 1.0 * 5.803e3/4651372.0;
   mc_p6_z2[4] = ntuple_pythia6_v16_path + "QCD_Ht-1000toInf_TuneZ2star_7TeV_madgraph_pythia6_V16.root";
   lumi_p6_z2[4] = 1.0 * 1.237e2/3478280.0;


   //Pythia6_TuneZ2star
//   int n_files_p6_z2 = 5;
//   string mc_p6_z2[5];
//   double lumi_p6_z2[5];
//   mc_p6_z2[0] = ntuple_pythia6_path + "QCD_Pt_10to25_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_10to25_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
//   lumi_p6_z2[0] = 0.06727 * 3.651e9/5205128.0;
//   mc_p6_z2[1] = ntuple_pythia6_path + "QCD_Pt_25to40_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_25to40_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
//   lumi_p6_z2[1] = 0.27707 * 1.059e8/5194564.0;
//   mc_p6_z2[2] = ntuple_pythia6_path + "QCD_Pt_40to80_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_40to80_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
//   lumi_p6_z2[2] = 0.25234 * 1.771e7/5011067.0;
//   mc_p6_z2[3] = ntuple_pythia6_path + "QCD_Pt_80to150_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_80to150_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
//   lumi_p6_z2[3] = 0.19052 * 8.756e5/4876756.0;
//   mc_p6_z2[4] = "../../../../../../../../work/c/cipriano/public/pythia6_150toinf/QCD_Pt_150toInf_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
//   lumi_p6_z2[4] = 0.14464 * 4.76e4/5710726.0;
   string vertex_p6_z2_allvertex = vertex_dir + "v0/vertex_Pythia6_TuneZ2star_allvertex.root";
   string vertex_p6_z2_allvertex_jetmettau_v1 = vertex_dir + "v1/vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v1.root";
   string vertex_p6_z2_allvertex_jetmet_v1 = vertex_dir + "v1/vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v1.root";
   string vertex_p6_z2_allvertex_jet_v1 = vertex_dir + "v1/vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v1.root";
   string vertex_p6_z2_allvertex_jetmettau_v2 = vertex_dir + "v2/vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v2.root";
   string vertex_p6_z2_allvertex_jetmet_v2 = vertex_dir + "v2/vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v2.root";
   string vertex_p6_z2_allvertex_jet_v2 = vertex_dir + "v2/vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v2.root";
   string vertex_p6_z2_allvertex_jetmettau_v3 = vertex_dir + "v3/vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v3.root";
   string vertex_p6_z2_allvertex_jetmet_v3 = vertex_dir + "v3/vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v3.root";
   string vertex_p6_z2_allvertex_jet_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v3.root";
   string vertex_p6_z2_allvertex_jetmettau_v4 = vertex_dir + "v4/vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v4.root";
   string vertex_p6_z2_allvertex_jetmet_v4 = vertex_dir + "v4/vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v4.root";
   string vertex_p6_z2_allvertex_jet_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v4.root";
   string vertex_weights_p6_z2_jetmettau_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v1.root";
   string vertex_weights_p6_z2_jetmet_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v1.root";
   string vertex_weights_p6_z2_jet_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v1.root";
   string vertex_weights_p6_z2_jetmettau_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v2.root";
   string vertex_weights_p6_z2_jetmet_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v2.root";
   string vertex_weights_p6_z2_jet_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v2.root";
   string vertex_weights_p6_z2_jetmettau_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v3.root";
   string vertex_weights_p6_z2_jetmet_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v3.root";
   string vertex_weights_p6_z2_jet_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v3.root";
   string vertex_weights_p6_z2_jetmettau_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v4.root";
   string vertex_weights_p6_z2_jetmet_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v4.root";
   string vertex_weights_p6_z2_jet_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v4.root";
   string vertex_weights_p6_z2_jetmettau_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v5.root";
   string vertex_weights_p6_z2_jetmet_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v5.root";
   string vertex_weights_p6_z2_jet_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v5.root";
   string vertex_p6_z2_allvertex_jetmettau_notnorm_v1 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_notnorm_v1.root";
   string vertex_p6_z2_allvertex_jetmet_notnorm_v1 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_notnorm_v1.root";
   string vertex_p6_z2_allvertex_jet_notnorm_v1 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_notnorm_v1.root";
   string vertex_p6_z2_allvertex_jetmettau_notnorm_v2 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_notnorm_v2.root";
   string vertex_p6_z2_allvertex_jetmet_notnorm_v2 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_notnorm_v2.root";
   string vertex_p6_z2_allvertex_jet_notnorm_v2 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_notnorm_v2.root";
   string vertex_p6_z2_allvertex_jetmettau_notnorm_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_notnorm_v3.root";
   string vertex_p6_z2_allvertex_jetmet_notnorm_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_notnorm_v3.root";
   string vertex_p6_z2_allvertex_jet_notnorm_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_notnorm_v3.root";
   string vertex_p6_z2_allvertex_jetmettau_notnorm_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_notnorm_v4.root";
   string vertex_p6_z2_allvertex_jetmet_notnorm_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_notnorm_v4.root";
   string vertex_p6_z2_allvertex_jet_notnorm_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_notnorm_v4.root";
   string vertex_weights_p6_z2_jetmettau_notnorm_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_notnorm_v1.root";
   string vertex_weights_p6_z2_jetmet_notnorm_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_notnorm_v1.root";
   string vertex_weights_p6_z2_jet_notnorm_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_notnorm_v1.root";
   string vertex_weights_p6_z2_jetmettau_notnorm_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_notnorm_v2.root";
   string vertex_weights_p6_z2_jetmet_notnorm_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_notnorm_v2.root";
   string vertex_weights_p6_z2_jet_notnorm_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_notnorm_v2.root";
   string vertex_weights_p6_z2_jetmettau_notnorm_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_notnorm_v3.root";
   string vertex_weights_p6_z2_jetmet_notnorm_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_notnorm_v3.root";
   string vertex_weights_p6_z2_jet_notnorm_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_notnorm_v3.root";
   string vertex_weights_p6_z2_jetmettau_notnorm_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_notnorm_v4.root";
   string vertex_weights_p6_z2_jetmet_notnorm_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_notnorm_v4.root";
   string vertex_weights_p6_z2_jet_notnorm_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_notnorm_v4.root";
   string vertex_weights_p6_z2_jetmettau_notnorm_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_notnorm_v5.root";
   string vertex_weights_p6_z2_jetmet_notnorm_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_notnorm_v5.root";
   string vertex_weights_p6_z2_jet_notnorm_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_notnorm_v5.root";
   string resolution_p6_z2_jetmettau = resolution_dir + "resolution_Pythia6_TuneZ2star_JetMETTau_2010A.root";
   string resolution_p6_z2_jetmet = resolution_dir + "resolution_Pythia6_TuneZ2star_JetMET_2010A.root";
   string resolution_p6_z2_jet = resolution_dir + "resolution_Pythia6_TuneZ2star_Jet_2010B.root";
   string norm_p6_z2 = pileup_norm_dir + "norm_Pythia6_TuneZ2star.root";
   string norm_p6_z2_jetmettau = pileup_norm_dir + "norm_Pythia6_TuneZ2star_JetMETTau_2010A.root";
   string norm_p6_z2_jetmet = pileup_norm_dir + "norm_Pythia6_TuneZ2star_JetMET_2010A.root";
   string norm_p6_z2_jet = pileup_norm_dir + "norm_Pythia6_TuneZ2star_Jet_2010B.root";
   string norm_p6_z2_notnorm_jetmettau = pileup_norm_dir + "norm_Pythia6_TuneZ2star_notnorm_JetMETTau_2010A.root";
   string norm_p6_z2_notnorm_jetmet = pileup_norm_dir + "norm_Pythia6_TuneZ2star_notnorm_JetMET_2010A.root";
   string norm_p6_z2_notnorm_jet = pileup_norm_dir + "norm_Pythia6_TuneZ2star_notnorm_Jet_2010B.root";
   string out_p6_z2_gen_nopileup = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup.root";
   string out_p6_z2_gen_allvertex = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex.root";
   string out_p6_z2_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup_JetMETTau_2010A.root"; 
   string out_p6_z2_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex_JetMETTau_2010A.root";
   string out_p6_z2_gen_nopileup_jetmet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup_JetMET_2010A.root"; 
   string out_p6_z2_gen_allvertex_jetmet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex_JetMET_2010A.root";
   string out_p6_z2_gen_nopileup_jet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup_Jet_2010B.root"; 
   string out_p6_z2_gen_allvertex_jet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex_Jet_2010B.root";
   string out_p6_z2_notnorm_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Pythia6_TuneZ2star_notnorm_gen_nopileup_JetMETTau_2010A.root"; 
   string out_p6_z2_notnorm_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Pythia6_TuneZ2star_notnorm_gen_allvertex_JetMETTau_2010A.root";
   string out_p6_z2_notnorm_gen_nopileup_jetmet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_notnorm_gen_nopileup_JetMET_2010A.root"; 
   string out_p6_z2_notnorm_gen_allvertex_jetmet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_notnorm_gen_allvertex_JetMET_2010A.root";
   string out_p6_z2_notnorm_gen_nopileup_jet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_notnorm_gen_nopileup_Jet_2010B.root"; 
   string out_p6_z2_notnorm_gen_allvertex_jet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_notnorm_gen_allvertex_Jet_2010B.root";
   string out_p6_z2_gen_norm = output_dir + "histograms/normalized_mc/xsec_p6_z2_gen.root";
   string out_p6_z2_gen_jetmettau_norm = mc_norm_dir + "xsec_p6_z2_JetMETTau_2010A_gen.root";
   string out_p6_z2_gen_jetmet_norm = mc_norm_dir + "xsec_p6_z2_JetMET_2010A_gen.root";
   string out_p6_z2_gen_jet_norm = mc_norm_dir + "xsec_p6_z2_Jet_2010B_gen.root";
   string out_p6_z2_det_1vertex = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex.root"; 
   string out_p6_z2_det_allvertex = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex.root";
   string out_p6_z2_det_allvertex_up = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_up.root";
   string out_p6_z2_det_allvertex_down = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_down.root";
   string out_p6_z2_det_1vertex_jetmettau = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex_JetMETTau_2010A.root"; 
   string out_p6_z2_det_allvertex_jetmettau = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMETTau_2010A.root";
   string out_p6_z2_det_allvertex_jetmettau_up = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMETTau_2010A_up.root";
   string out_p6_z2_det_allvertex_jetmettau_down = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMETTau_2010A_down.root";
   string out_p6_z2_det_1vertex_jetmet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex_JetMET_2010A.root"; 
   string out_p6_z2_det_allvertex_jetmet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMET_2010A.root";
   string out_p6_z2_det_allvertex_jetmet_up = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMET_2010A_up.root";
   string out_p6_z2_det_allvertex_jetmet_down = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMET_2010A_down.root";
   string out_p6_z2_det_1vertex_jet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex_Jet_2010B.root"; 
   string out_p6_z2_det_allvertex_jet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_Jet_2010B.root";
   string out_p6_z2_det_allvertex_jet_up = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_Jet_2010B_up.root";
   string out_p6_z2_det_allvertex_jet_down = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_Jet_2010B_down.root";
   string out_p6_z2_det_norm = mc_norm_dir + "xsec_p6_z2_det.root";
   string out_p6_z2_det_jetmettau_norm = mc_norm_dir + "xsec_p6_z2_JetMETTau_2010A_det.root";
   string out_p6_z2_det_jetmet_norm = mc_norm_dir + "xsec_p6_z2_JetMET_2010A_det.root";
   string out_p6_z2_det_jet_norm = mc_norm_dir + "xsec_p6_z2_Jet_2010B_det.root";
   string corr_pileup_p6_z2 = corrections_dir + "correction_pileup_p6_z2.root";
   string corr_pileup_p6_z2_jetmettau = corrections_dir + "correction_pileup_p6_z2_JetMETTau_2010A.root";
   string corr_pileup_p6_z2_jetmet = corrections_dir + "correction_pileup_p6_z2_JetMET_2010A.root";
   string corr_pileup_p6_z2_jet = corrections_dir + "correction_pileup_p6_z2_Jet_2010B.root";
   string corr_detector_p6_z2 = corrections_dir + "correction_detector_p6_z2.root";
   string corr_detector_p6_z2_jetmettau = corrections_dir + "correction_detector_p6_z2_JetMETTau_2010A.root";
   string corr_detector_p6_z2_jetmet = corrections_dir + "correction_detector_p6_z2_JetMET_2010A.root";
   string corr_detector_p6_z2_jet = corrections_dir + "correction_detector_p6_z2_Jet_2010B.root";
   string corr_final_p6_z2 = corrections_dir + "correction_final_p6_z2.root";
   string corr_final_p6_z2_jetmettau = corrections_dir + "correction_final_p6_z2_JetMETTau_2010A.root";
   string corr_final_p6_z2_jetmet = corrections_dir + "correction_final_p6_z2_JetMET_2010A.root";
   string corr_final_p6_z2_jet = corrections_dir + "correction_final_p6_z2_Jet_2010B.root";
   string corrected_p6_z2_jetmettau = closure_dir + "closure_final_p6_z2_JetMETTau_2010A.root";
   string corrected_p6_z2_jetmet = closure_dir + "closure_final_p6_z2_JetMET_2010A.root";
   string corrected_p6_z2_jet = closure_dir + "closure_final_p6_z2_Jet_2010B.root";
   string closure_ratio_p6_z2_jetmettau = closure_dir + "closure_ratio_p6_z2_JetMETTau_2010A.root";
   string closure_ratio_p6_z2_jetmet = closure_dir + "closure_ratio_p6_z2_JetMET_2010A.root";
   string closure_ratio_p6_z2_jet = closure_dir + "closure_ratio_p6_z2_Jet_2010B.root";
   string matched_p6_z2 = matched_dir + "matched_p6_z2.root";
   string matched_p6_z2_jetmettau = matched_dir + "matched_p6_z2_JetMETTau_2010A.root";
   string matched_p6_z2_jetmet = matched_dir + "matched_p6_z2_JetMET_2010A.root";
   string matched_p6_z2_jet = matched_dir + "matched_p6_z2_Jet_2010B.root";
   string out_p6_z2_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_jetmettau_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_jetmet_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_jet_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_bbb_data = unfolded_dir + "Data_unfolded_BinByBin_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_bbb_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_bbb_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_bbb_data_jetmettau = unfolded_dir + "Data_unfolded_BinByBin_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_bbb_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_bbb_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_bbb_data_jetmet = unfolded_dir + "Data_unfolded_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_bbb_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_bbb_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_bbb_data_jet = unfolded_dir + "Data_unfolded_BinByBin_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_bbb_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_bbb_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_tunfold_data = unfolded_dir + "Data_unfolded_TUnfold_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p8_4c = unfolded_dir + "Pythia8_Tune4C_TUnfold_unfolded_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_tunfold_data_jetmettau = unfolded_dir + "Data_unfolded_TUnfold_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_data_jetmet = unfolded_dir + "Data_unfolded_TUnfold_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_tunfold_data_jet = unfolded_dir + "Data_unfolded_TUnfold_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_tunfold_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_svd_data = unfolded_dir + "Data_unfolded_SVD_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_svd_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_svd_p8_4c = unfolded_dir + "Pythia8_Tune4C_SVD_unfolded_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_svd_data_jetmettau = unfolded_dir + "Data_unfolded_SVD_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_svd_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_svd_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_svd_data_jetmet = unfolded_dir + "Data_unfolded_SVD_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_svd_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_svd_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_svd_data_jet = unfolded_dir + "Data_unfolded_SVD_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_svd_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_svd_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_bayes_data = unfolded_dir + "Data_unfolded_Bayes_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_bayes_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_bayes_p8_4c = unfolded_dir + "Pythia8_Tune4C_Bayes_unfolded_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_bayes_data_jetmettau = unfolded_dir + "Data_unfolded_Bayes_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_bayes_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_bayes_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_bayes_data_jetmet = unfolded_dir + "Data_unfolded_Bayes_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_bayes_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_bayes_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_bayes_data_jet = unfolded_dir + "Data_unfolded_Bayes_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_bayes_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_bayes_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";


   //Pythia8_Tune4C
   int n_files_p8_4c = 1;
   string mc_p8_4c[1];
   double lumi_p8_4c[1];
   mc_p8_4c[0] = ntuple_pythia6_v16_path + "QCD_Pt-15to3000_Tune4C_Flat_7TeV_pythia8_V16.root";
   lumi_p8_4c[0] = 1.0 * 2.84e10/4989964.0;

   // int n_files_p8_4c = 4;
   // string mc_p8_4c[5];
   // double lumi_p8_4c[5];
   // mc_p8_4c[0] = ntuple_pythia8_path + "QCD_Pt_10to25_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_10to25_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   // lumi_p8_4c[0] = 0.07463 * 4.334e9/9695621.0;
   // mc_p8_4c[1] = ntuple_pythia8_path + "QCD_Pt_25to40_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_25to40_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   // lumi_p8_4c[1] = 0.29023 * 1.226e8/9949451.0;
   // mc_p8_4c[2] = ntuple_pythia8_path + "QCD_Pt_40to80_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_40to80_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   // lumi_p8_4c[2] = 0.27368 * 2.021e7/9658287.0;
   // mc_p8_4c[3] = ntuple_pythia8_path + "QCD_Pt_80to150_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_80to150_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   // lumi_p8_4c[3] = 0.21046 * 9.855e5/9117322.0;
   // mc_p8_4c[4] = "../../../../../../../../work/p/pgunnell/public/PYTHIA8PEDRO/JetTreePYTHIA8.root";
   // lumi_p8_4c[4] = 0.16428 * 5.283e4/9761171.0;

   string vertex_p8_4c_allvertex = vertex_dir + "vertex_Pythia8_Tune4C_allvertex.root";
   string vertex_p8_4c_allvertex_jetmettau_v1 = vertex_dir + "v1/vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v1.root";
   string vertex_p8_4c_allvertex_jetmet_v1 = vertex_dir + "v1/vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v1.root";      
   string vertex_p8_4c_allvertex_jet_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v1.root";
   string vertex_p8_4c_allvertex_jetmettau_v2 = vertex_dir + "v2/vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v2.root";
   string vertex_p8_4c_allvertex_jetmet_v2 = vertex_dir + "v2/vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v2.root";      
   string vertex_p8_4c_allvertex_jet_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v2.root";
   string vertex_p8_4c_allvertex_jetmettau_v3 = vertex_dir + "v3/vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v3.root";
   string vertex_p8_4c_allvertex_jetmet_v3 = vertex_dir + "v3/vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v3.root";      
   string vertex_p8_4c_allvertex_jet_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v3.root";
   string vertex_p8_4c_allvertex_jetmettau_v4 = vertex_dir + "v4/vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v4.root";
   string vertex_p8_4c_allvertex_jetmet_v4 = vertex_dir + "v4/vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v4.root";      
   string vertex_p8_4c_allvertex_jet_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v4.root";
   string vertex_weights_p8_4c_jetmettau_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v1.root";
   string vertex_weights_p8_4c_jetmet_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v1.root";
   string vertex_weights_p8_4c_jet_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v1.root";
   string vertex_weights_p8_4c_jetmettau_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v2.root";
   string vertex_weights_p8_4c_jetmet_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v2.root"; 
   string vertex_weights_p8_4c_jet_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v2.root";
   string vertex_weights_p8_4c_jetmettau_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v3.root";
   string vertex_weights_p8_4c_jetmet_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v3.root"; 
   string vertex_weights_p8_4c_jet_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v3.root";
   string vertex_weights_p8_4c_jetmettau_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v4.root";
   string vertex_weights_p8_4c_jetmet_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v4.root";
   string vertex_weights_p8_4c_jet_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v4.root";
   string vertex_weights_p8_4c_jetmettau_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v5.root";
   string vertex_weights_p8_4c_jetmet_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v5.root";
   string vertex_weights_p8_4c_jet_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v5.root";
   string vertex_p8_4c_allvertex_jetmettau_notnorm_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_notnorm_v1.root";
   string vertex_p8_4c_allvertex_jetmet_notnorm_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_notnorm_v1.root";      
   string vertex_p8_4c_allvertex_jet_notnorm_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_notnorm_v1.root";
   string vertex_p8_4c_allvertex_jetmettau_notnorm_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_notnorm_v2.root";
   string vertex_p8_4c_allvertex_jetmet_notnorm_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_notnorm_v2.root";      
   string vertex_p8_4c_allvertex_jet_notnorm_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_notnorm_v2.root";
   string vertex_p8_4c_allvertex_jetmettau_notnorm_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_notnorm_v3.root";
   string vertex_p8_4c_allvertex_jetmet_notnorm_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_notnorm_v3.root";      
   string vertex_p8_4c_allvertex_jet_notnorm_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_notnorm_v3.root";
   string vertex_p8_4c_allvertex_jetmettau_notnorm_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_notnorm_v4.root";
   string vertex_p8_4c_allvertex_jetmet_notnorm_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_notnorm_v4.root";      
   string vertex_p8_4c_allvertex_jet_notnorm_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_notnorm_v4.root";
   string vertex_weights_p8_4c_jetmettau_notnorm_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_notnorm_v1.root";
   string vertex_weights_p8_4c_jetmet_notnorm_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_notnorm_v1.root";
   string vertex_weights_p8_4c_jet_notnorm_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_notnorm_v1.root";
   string vertex_weights_p8_4c_jetmettau_notnorm_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_notnorm_v2.root";
   string vertex_weights_p8_4c_jetmet_notnorm_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_notnorm_v2.root"; 
   string vertex_weights_p8_4c_jet_notnorm_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_notnorm_v2.root";
   string vertex_weights_p8_4c_jetmettau_notnorm_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_notnorm_v3.root";
   string vertex_weights_p8_4c_jetmet_notnorm_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_notnorm_v3.root"; 
   string vertex_weights_p8_4c_jet_notnorm_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_notnorm_v3.root";
   string vertex_weights_p8_4c_jetmettau_notnorm_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_notnorm_v4.root";
   string vertex_weights_p8_4c_jetmet_notnorm_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_notnorm_v4.root";
   string vertex_weights_p8_4c_jet_notnorm_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_notnorm_v4.root";
   string vertex_weights_p8_4c_jetmettau_notnorm_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_notnorm_v5.root";
   string vertex_weights_p8_4c_jetmet_notnorm_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_notnorm_v5.root";
   string vertex_weights_p8_4c_jet_notnorm_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_notnorm_v5.root";
   string resolution_p8_4c_jetmettau = resolution_dir + "resolution_Pythia8_Tune4C_JetMETTau_2010A.root";
   string resolution_p8_4c_jetmet = resolution_dir + "resolution_Pythia8_Tune4C_JetMET_2010A.root";
   string resolution_p8_4c_jet = resolution_dir + "resolution_Pythia8_Tune4C_Jet_2010B.root";
   string norm_p8_4c = pileup_norm_dir + "norm_Pythia8_Tune4C.root";
   string norm_p8_4c_jetmettau = pileup_norm_dir + "norm_Pythia8_Tune4C_JetMETTau_2010A.root";
   string norm_p8_4c_jetmet = pileup_norm_dir + "norm_Pythia8_Tune4C_JetMET_2010A.root";
   string norm_p8_4c_jet = pileup_norm_dir + "norm_Pythia8_Tune4C_Jet_2010B.root";
   string norm_p8_4c_notnorm_jetmettau = pileup_norm_dir + "norm_Pythia8_Tune4C_notnorm_JetMETTau_2010A.root";
   string norm_p8_4c_notnorm_jetmet = pileup_norm_dir + "norm_Pythia8_Tune4C_notnorm_JetMET_2010A.root";
   string norm_p8_4c_notnorm_jet = pileup_norm_dir + "norm_Pythia8_Tune4C_notnorm_Jet_2010B.root";
   string out_p8_4c_gen_nopileup = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup.root";
   string out_p8_4c_gen_allvertex = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex.root";
   string out_p8_4c_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup_JetMETTau_2010A.root";
   string out_p8_4c_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex_JetMETTau_2010A.root";
   string out_p8_4c_gen_nopileup_jetmet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup_JetMET_2010A.root";
   string out_p8_4c_gen_allvertex_jetmet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex_JetMET_2010A.root";
   string out_p8_4c_gen_nopileup_jet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup_Jet_2010B.root";
   string out_p8_4c_gen_allvertex_jet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex_Jet_2010B.root";
   string out_p8_4c_notnorm_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Pythia8_Tune4C_notnorm_gen_nopileup_JetMETTau_2010A.root";
   string out_p8_4c_notnorm_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Pythia8_Tune4C_notnorm_gen_allvertex_JetMETTau_2010A.root";
   string out_p8_4c_notnorm_gen_nopileup_jetmet = mc_gen_dir + "xsec_Pythia8_Tune4C_notnorm_gen_nopileup_JetMET_2010A.root";
   string out_p8_4c_notnorm_gen_allvertex_jetmet = mc_gen_dir + "xsec_Pythia8_Tune4C_notnorm_gen_allvertex_JetMET_2010A.root";
   string out_p8_4c_notnorm_gen_nopileup_jet = mc_gen_dir + "xsec_Pythia8_Tune4C_notnorm_gen_nopileup_Jet_2010B.root";
   string out_p8_4c_notnorm_gen_allvertex_jet = mc_gen_dir + "xsec_Pythia8_Tune4C_notnorm_gen_allvertex_Jet_2010B.root";
   string out_p8_4c_gen_norm = mc_norm_dir + "xsec_p8_4c_gen.root";
   string out_p8_4c_gen_jetmettau_norm = mc_norm_dir + "xsec_p8_4c_JetMETTau_2010A_gen.root";
   string out_p8_4c_gen_jetmet_norm = mc_norm_dir + "xsec_p8_4c_JetMET_2010A_gen.root";
   string out_p8_4c_gen_jet_norm = mc_norm_dir + "xsec_p8_4c_Jet_2010B_gen.root";
   string out_p8_4c_det_1vertex = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex.root";
   string out_p8_4c_det_allvertex = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex.root";
   string out_p8_4c_det_1vertex_jetmettau = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex_JetMETTau_2010A.root";
   string out_p8_4c_det_allvertex_jetmettau = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMETTau_2010A.root";
   string out_p8_4c_det_allvertex_jetmettau_up = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMETTau_2010A_up.root";
   string out_p8_4c_det_allvertex_jetmettau_down = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMETTau_2010A_down.root";
   string out_p8_4c_det_1vertex_jetmet = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex_JetMET_2010A.root";
   string out_p8_4c_det_allvertex_jetmet = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMET_2010A.root";
   string out_p8_4c_det_allvertex_jetmet_up = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMET_2010A_up.root";
   string out_p8_4c_det_allvertex_jetmet_down = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMET_2010A_down.root";
   string out_p8_4c_det_1vertex_jet = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex_Jet_2010B.root";
   string out_p8_4c_det_allvertex_jet = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_Jet_2010B.root";
   string out_p8_4c_det_allvertex_jet_up = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_Jet_2010B_up.root";
   string out_p8_4c_det_allvertex_jet_down = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_Jet_2010B_down.root";
   string out_p8_4c_det_norm = mc_norm_dir + "xsec_p8_4c_det.root";
   string out_p8_4c_det_jetmettau_norm = mc_norm_dir + "xsec_p8_4c_JetMETTau_2010A_det.root";
   string out_p8_4c_det_jetmet_norm = mc_norm_dir + "xsec_p8_4c_JetMET_2010A_det.root";
   string out_p8_4c_det_jet_norm = mc_norm_dir + "xsec_p8_4c_Jet_2010B_det.root";
   string corr_pileup_p8_4c = corrections_dir + "correction_pileup_p8_4c.root";
   string corr_pileup_p8_4c_jetmettau = corrections_dir + "correction_pileup_p8_4c_JetMETTau_2010A.root";
   string corr_pileup_p8_4c_jetmet = corrections_dir + "correction_pileup_p8_4c_JetMET_2010A.root";
   string corr_pileup_p8_4c_jet = corrections_dir + "correction_pileup_p8_4c_Jet_2010B.root";
   string corr_detector_p8_4c = corrections_dir + "correction_detector_p8_4c.root";
   string corr_detector_p8_4c_jetmettau = corrections_dir + "correction_detector_p8_4c_JetMETTau_2010A.root";
   string corr_detector_p8_4c_jetmet = corrections_dir + "correction_detector_p8_4c_JetMET_2010A.root";
   string corr_detector_p8_4c_jet = corrections_dir + "correction_detector_p8_4c_Jet_2010B.root";
   string corr_final_p8_4c = corrections_dir + "correction_final_p8_4c.root";
   string corr_final_p8_4c_jetmettau = corrections_dir + "correction_final_p8_4c_JetMETTau_2010A.root";
   string corr_final_p8_4c_jetmet = corrections_dir + "correction_final_p8_4c_JetMET_2010A.root";
   string corr_final_p8_4c_jet = corrections_dir + "correction_final_p8_4c_Jet_2010B.root";
   string corrected_p8_4c_jetmettau = closure_dir + "closure_final_p8_4c_JetMETTau_2010A.root";
   string corrected_p8_4c_jetmet = closure_dir + "closure_final_p8_4c_JetMET_2010A.root";
   string corrected_p8_4c_jet = closure_dir + "closure_final_p8_4c_Jet_2010B.root";
   string closure_ratio_p8_4c_jetmettau = closure_dir + "closure_ratio_p8_4c_JetMETTau_2010A.root";
   string closure_ratio_p8_4c_jetmet = closure_dir + "closure_ratio_p8_4c_JetMET_2010A.root";
   string closure_ratio_p8_4c_jet = closure_dir + "closure_ratio_p8_4c_Jet_2010B.root";
   string matched_p8_4c = matched_dir + "matched_p8_4c.root";
   string matched_p8_4c_jetmettau = matched_dir + "matched_p8_4c_JetMETTau_2010A.root";
   string matched_p8_4c_jetmet = matched_dir + "matched_p8_4c_JetMET_2010A.root";
   string matched_p8_4c_jet = matched_dir + "matched_p8_4c_Jet_2010B.root";
   string out_p8_4c_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_jetmettau_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_jetmet_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_jet_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_bbb_data = unfolded_dir + "Data_unfolded_BinByBin_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bbb_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bbb_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bbb_data_jetmettau = unfolded_dir + "Data_unfolded_BinByBin_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_BinByBin_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_data_jetmet = unfolded_dir + "Data_unfolded_BinByBin_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_BinByBin_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bbb_data_jet = unfolded_dir + "Data_unfolded_BinByBin_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_bbb_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_bbb_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_tunfold_data = unfolded_dir + "Data_unfolded_TUnfold_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_tunfold_data_jetmettau = unfolded_dir + "Data_unfolded_TUnfold_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";

//method test
   string out_p8_4c_unfolded_tunfold1_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold1_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold2_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold2_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold3_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold3_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold4_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold4_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";

   string out_p8_4c_unfolded_tunfold_data_jetmet = unfolded_dir + "Data_unfolded_TUnfold_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia8_Tune4C_JetMET_2010A_allvertex.root";

//method test
   string out_p8_4c_unfolded_tunfold1_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold1_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold2_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold2_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold3_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold3_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_tunfold4_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold4_Pythia8_Tune4C_JetMET_2010A_allvertex.root";

   string out_p8_4c_unfolded_tunfold_data_jet = unfolded_dir + "Data_unfolded_TUnfold_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_tunfold_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_svd_data = unfolded_dir + "Data_unfolded_SVD_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_svd_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_svd_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_svd_data_jetmettau = unfolded_dir + "Data_unfolded_SVD_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";

//method test
   string out_p8_4c_unfolded_svd1_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD1_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd2_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD2_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd3_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD3_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd4_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD4_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";


   string out_p8_4c_unfolded_svd_data_jetmet = unfolded_dir + "Data_unfolded_SVD_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia8_Tune4C_JetMET_2010A_allvertex.root";

//method test
   string out_p8_4c_unfolded_svd1_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD1_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd2_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD2_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd3_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD3_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_svd4_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD4_Pythia8_Tune4C_JetMET_2010A_allvertex.root";


   string out_p8_4c_unfolded_svd_data_jet = unfolded_dir + "Data_unfolded_SVD_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_svd_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_svd_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_bayes_data = unfolded_dir + "Data_unfolded_Bayes_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes_hw_eec3 = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes_hw_23 = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Pythia_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia8_Tune4C_allvertex.root";

//method test
   string out_p8_4c_unfolded_bayes1_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes1_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes2_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes2_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes3_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes3_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_bayes4_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes4_Pythia8_Tune4C_allvertex.root";

   string out_p8_4c_unfolded_bayes_data_jetmettau = unfolded_dir + "Data_unfolded_Bayes_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";

//method test
   string out_p8_4c_unfolded_bayes1_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes1_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes2_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes2_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes3_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes3_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes4_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes4_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";

   string out_p8_4c_unfolded_bayes_data_jetmet = unfolded_dir + "Data_unfolded_Bayes_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Pythia8_Tune4c_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Pythia8_Tun4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia8_Tune4C_JetMET_2010A_allvertex.root";

//method test
   string out_p8_4c_unfolded_bayes1_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes1_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes2_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes2_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes3_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes3_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_bayes4_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes4_Pythia8_Tune4C_JetMET_2010A_allvertex.root";

   string out_p8_4c_unfolded_bayes_data_jet = unfolded_dir + "Data_unfolded_Bayes_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_bayes_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_bayes_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Pythia8_Tune4C_Jet_2010B_allvertex.root";


// Herwig pp EEC3
   string mc_hw_eec3[1];
   double lumi_hw_eec3[1];
   int n_files_hw_eec3 = 1;
   mc_hw_eec3[0] = ntuple_pythia6_v16_path + "QCD_Pt-15to1000_EEC3_Flat_7TeV_herwigpp_V16.root";
   lumi_hw_eec3[0] = 1.0 * 1.69532e10/10801533.0;

   string vertex_hw_eec3_allvertex = vertex_dir + "v0/vertex_Herwigpp_TuneEEC3_allvertex.root";
   string vertex_hw_eec3_allvertex_jetmettau_v1 = vertex_dir + "v1/vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_v1.root";
   string vertex_hw_eec3_allvertex_jetmet_v1 = vertex_dir + "v1/vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_v1.root";      
   string vertex_hw_eec3_allvertex_jet_v1 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_v1.root";
   string vertex_hw_eec3_allvertex_jetmettau_v2 = vertex_dir + "v2/vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_v2.root";
   string vertex_hw_eec3_allvertex_jetmet_v2 = vertex_dir + "v2/vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_v2.root";      
   string vertex_hw_eec3_allvertex_jet_v2 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_v2.root";
   string vertex_hw_eec3_allvertex_jetmettau_v3 = vertex_dir + "v3/vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_v3.root";
   string vertex_hw_eec3_allvertex_jetmet_v3 = vertex_dir + "v3/vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_v3.root";      
   string vertex_hw_eec3_allvertex_jet_v3 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_v3.root";
   string vertex_hw_eec3_allvertex_jetmettau_v4 = vertex_dir + "v4/vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_v4.root";
   string vertex_hw_eec3_allvertex_jetmet_v4 = vertex_dir + "v4/vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_v4.root";      
   string vertex_hw_eec3_allvertex_jet_v4 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_v4.root";
   string vertex_weights_hw_eec3_jetmettau_v1 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_v1.root";
   string vertex_weights_hw_eec3_jetmet_v1 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_v1.root";
   string vertex_weights_hw_eec3_jet_v1 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_v1.root";
   string vertex_weights_hw_eec3_jetmettau_v2 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_v2.root";
   string vertex_weights_hw_eec3_jetmet_v2 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_v2.root"; 
   string vertex_weights_hw_eec3_jet_v2 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_v2.root";
   string vertex_weights_hw_eec3_jetmettau_v3 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_v3.root";
   string vertex_weights_hw_eec3_jetmet_v3 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_v3.root"; 
   string vertex_weights_hw_eec3_jet_v3 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_v3.root";
   string vertex_weights_hw_eec3_jetmettau_v4 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_v4.root";
   string vertex_weights_hw_eec3_jetmet_v4 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_v4.root";
   string vertex_weights_hw_eec3_jet_v4 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_v4.root";
   string vertex_weights_hw_eec3_jetmettau_v5 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_v5.root";
   string vertex_weights_hw_eec3_jetmet_v5 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_v5.root";
   string vertex_weights_hw_eec3_jet_v5 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_v5.root";
   string vertex_hw_eec3_allvertex_jetmettau_notnorm_v1 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_notnorm_v1.root";
   string vertex_hw_eec3_allvertex_jetmet_notnorm_v1 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_notnorm_v1.root";      
   string vertex_hw_eec3_allvertex_jet_notnorm_v1 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_notnorm_v1.root";
   string vertex_hw_eec3_allvertex_jetmettau_notnorm_v2 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_notnorm_v2.root";
   string vertex_hw_eec3_allvertex_jetmet_notnorm_v2 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_notnorm_v2.root";      
   string vertex_hw_eec3_allvertex_jet_notnorm_v2 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_notnorm_v2.root";
   string vertex_hw_eec3_allvertex_jetmettau_notnorm_v3 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_notnorm_v3.root";
   string vertex_hw_eec3_allvertex_jetmet_notnorm_v3 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_notnorm_v3.root";      
   string vertex_hw_eec3_allvertex_jet_notnorm_v3 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_notnorm_v3.root";
   string vertex_hw_eec3_allvertex_jetmettau_notnorm_v4 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMETTau_2010A_notnorm_v4.root";
   string vertex_hw_eec3_allvertex_jetmet_notnorm_v4 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_JetMET_2010A_notnorm_v4.root";      
   string vertex_hw_eec3_allvertex_jet_notnorm_v4 = vertex_dir + "vertex_Herwigpp_TuneEEC3_allvertex_Jet_2010B_notnorm_v4.root";
   string vertex_weights_hw_eec3_jetmettau_notnorm_v1 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_notnorm_v1.root";
   string vertex_weights_hw_eec3_jetmet_notnorm_v1 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_notnorm_v1.root";
   string vertex_weights_hw_eec3_jet_notnorm_v1 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_notnorm_v1.root";
   string vertex_weights_hw_eec3_jetmettau_notnorm_v2 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_notnorm_v2.root";
   string vertex_weights_hw_eec3_jetmet_notnorm_v2 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_notnorm_v2.root"; 
   string vertex_weights_hw_eec3_jet_notnorm_v2 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_notnorm_v2.root";
   string vertex_weights_hw_eec3_jetmettau_notnorm_v3 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_notnorm_v3.root";
   string vertex_weights_hw_eec3_jetmet_notnorm_v3 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_notnorm_v3.root"; 
   string vertex_weights_hw_eec3_jet_notnorm_v3 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_notnorm_v3.root";
   string vertex_weights_hw_eec3_jetmettau_notnorm_v4 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_notnorm_v4.root";
   string vertex_weights_hw_eec3_jetmet_notnorm_v4 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_notnorm_v4.root";
   string vertex_weights_hw_eec3_jet_notnorm_v4 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_notnorm_v4.root";
   string vertex_weights_hw_eec3_jetmettau_notnorm_v5 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMETTau2010A_notnorm_v5.root";
   string vertex_weights_hw_eec3_jetmet_notnorm_v5 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_JetMET2010A_notnorm_v5.root";
   string vertex_weights_hw_eec3_jet_notnorm_v5 = vertex_weights_dir + "vertex_Herwigpp_TuneEEC3_Jet2010B_notnorm_v5.root";
   string resolution_hw_eec3_jetmettau = resolution_dir + "resolution_Herwigpp_TuneEEC3_JetMETTau_2010A.root";
   string resolution_hw_eec3_jetmet = resolution_dir + "resolution_Herwigpp_TuneEEC3_JetMET_2010A.root";
   string resolution_hw_eec3_jet = resolution_dir + "resolution_Herwigpp_TuneEEC3_Jet_2010B.root";
   string norm_hw_eec3 = pileup_norm_dir + "norm_Herwigpp_TuneEEC3.root";
   string norm_hw_eec3_jetmettau = pileup_norm_dir + "norm_Herwigpp_TuneEEC3_JetMETTau_2010A.root";
   string norm_hw_eec3_jetmet = pileup_norm_dir + "norm_Herwigpp_TuneEEC3_JetMET_2010A.root";
   string norm_hw_eec3_jet = pileup_norm_dir + "norm_Herwigpp_TuneEEC3_Jet_2010B.root";
   string norm_hw_eec3_notnorm_jetmettau = pileup_norm_dir + "norm_Herwigpp_TuneEEC3_notnorm_JetMETTau_2010A.root";
   string norm_hw_eec3_notnorm_jetmet = pileup_norm_dir + "norm_Herwigpp_TuneEEC3_notnorm_JetMET_2010A.root";
   string norm_hw_eec3_notnorm_jet = pileup_norm_dir + "norm_Herwigpp_TuneEEC3_notnorm_Jet_2010B.root";

   string out_hw_eec3_gen_nopileup = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_nopileup.root";
   string out_hw_eec3_gen_allvertex = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_allvertex.root";
   string out_hw_eec3_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_nopileup_JetMETTau_2010A.root"; 
   string out_hw_eec3_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_allvertex_JetMETTau_2010A.root";
   string out_hw_eec3_gen_nopileup_jetmet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_nopileup_JetMET_2010A.root"; 
   string out_hw_eec3_gen_allvertex_jetmet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_allvertex_JetMET_2010A.root";
   string out_hw_eec3_gen_nopileup_jet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_nopileup_Jet_2010B.root"; 
   string out_hw_eec3_gen_allvertex_jet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_gen_allvertex_Jet_2010B.root";
   string out_hw_eec3_notnorm_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_notnorm_gen_nopileup_JetMETTau_2010A.root"; 
   string out_hw_eec3_notnorm_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_notnorm_gen_allvertex_JetMETTau_2010A.root";
   string out_hw_eec3_notnorm_gen_nopileup_jetmet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_notnorm_gen_nopileup_JetMET_2010A.root"; 
   string out_hw_eec3_notnorm_gen_allvertex_jetmet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_notnorm_gen_allvertex_JetMET_2010A.root";
   string out_hw_eec3_notnorm_gen_nopileup_jet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_notnorm_gen_nopileup_Jet_2010B.root"; 
   string out_hw_eec3_notnorm_gen_allvertex_jet = mc_gen_dir + "xsec_Herwigpp_TuneEEC3_notnorm_gen_allvertex_Jet_2010B.root";
   string out_hw_eec3_gen_norm = output_dir + "histograms/normalized_mc/xsec_hw_eec3_gen.root";
   string out_hw_eec3_gen_jetmettau_norm = mc_norm_dir + "xsec_hw_eec3_JetMETTau_2010A_gen.root";
   string out_hw_eec3_gen_jetmet_norm = mc_norm_dir + "xsec_hw_eec3_JetMET_2010A_gen.root";
   string out_hw_eec3_gen_jet_norm = mc_norm_dir + "xsec_hw_eec3_Jet_2010B_gen.root";
   string out_hw_eec3_det_1vertex = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_1vertex.root"; 
   string out_hw_eec3_det_allvertex = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex.root";
   string out_hw_eec3_det_allvertex_up = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_up.root";
   string out_hw_eec3_det_allvertex_down = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_down.root";
   string out_hw_eec3_det_1vertex_jetmettau = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_1vertex_JetMETTau_2010A.root"; 
   string out_hw_eec3_det_allvertex_jetmettau = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_JetMETTau_2010A.root";
   string out_hw_eec3_det_allvertex_jetmettau_up = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_JetMETTau_2010A_up.root";
   string out_hw_eec3_det_allvertex_jetmettau_down = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_JetMETTau_2010A_down.root";
   string out_hw_eec3_det_1vertex_jetmet = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_1vertex_JetMET_2010A.root"; 
   string out_hw_eec3_det_allvertex_jetmet = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_JetMET_2010A.root";
   string out_hw_eec3_det_allvertex_jetmet_up = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_JetMET_2010A_up.root";
   string out_hw_eec3_det_allvertex_jetmet_down = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_JetMET_2010A_down.root";
   string out_hw_eec3_det_1vertex_jet = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_1vertex_Jet_2010B.root"; 
   string out_hw_eec3_det_allvertex_jet = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_Jet_2010B.root";
   string out_hw_eec3_det_allvertex_jet_up = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_Jet_2010B_up.root";
   string out_hw_eec3_det_allvertex_jet_down = mc_det_dir + "xsec_Herwigpp_TuneEEC3_det_allvertex_Jet_2010B_down.root";
   string out_hw_eec3_det_norm = mc_norm_dir + "xsec_hw_eec3_det.root";
   string out_hw_eec3_det_jetmettau_norm = mc_norm_dir + "xsec_hw_eec3_JetMETTau_2010A_det.root";
   string out_hw_eec3_det_jetmet_norm = mc_norm_dir + "xsec_hw_eec3_JetMET_2010A_det.root";
   string out_hw_eec3_det_jet_norm = mc_norm_dir + "xsec_hw_eec3_Jet_2010B_det.root";

   string corr_pileup_hw_eec3 = corrections_dir + "correction_pileup_hw_eec3.root";
   string corr_pileup_hw_eec3_jetmettau = corrections_dir + "correction_pileup_hw_eec3_JetMETTau_2010A.root";
   string corr_pileup_hw_eec3_jetmet = corrections_dir + "correction_pileup_hw_eec3_JetMET_2010A.root";
   string corr_pileup_hw_eec3_jet = corrections_dir + "correction_pileup_hw_eec3_Jet_2010B.root";
   string corr_detector_hw_eec3 = corrections_dir + "correction_detector_hw_eec3.root";
   string corr_detector_hw_eec3_jetmettau = corrections_dir + "correction_detector_hw_eec3_JetMETTau_2010A.root";
   string corr_detector_hw_eec3_jetmet = corrections_dir + "correction_detector_hw_eec3_JetMET_2010A.root";
   string corr_detector_hw_eec3_jet = corrections_dir + "correction_detector_hw_eec3_Jet_2010B.root";
   string corr_final_hw_eec3 = corrections_dir + "correction_final_hw_eec3.root";
   string corr_final_hw_eec3_jetmettau = corrections_dir + "correction_final_hw_eec3_JetMETTau_2010A.root";
   string corr_final_hw_eec3_jetmet = corrections_dir + "correction_final_hw_eec3_JetMET_2010A.root";
   string corr_final_hw_eec3_jet = corrections_dir + "correction_final_hw_eec3_Jet_2010B.root";
   string corrected_hw_eec3_jetmettau = closure_dir + "closure_final_hw_eec3_JetMETTau_2010A.root";
   string corrected_hw_eec3_jetmet = closure_dir + "closure_final_hw_eec3_JetMET_2010A.root";
   string corrected_hw_eec3_jet = closure_dir + "closure_final_hw_eec3_Jet_2010B.root";
   string closure_ratio_hw_eec3_jetmettau = closure_dir + "closure_ratio_hw_eec3_JetMETTau_2010A.root";
   string closure_ratio_hw_eec3_jetmet = closure_dir + "closure_ratio_hw_eec3_JetMET_2010A.root";
   string closure_ratio_hw_eec3_jet = closure_dir + "closure_ratio_hw_eec3_Jet_2010B.root";
   string matched_hw_eec3 = matched_dir + "matched_hw_eec3.root";
   string matched_hw_eec3_jetmettau = matched_dir + "matched_hw_eec3_JetMETTau_2010A.root";
   string matched_hw_eec3_jetmet = matched_dir + "matched_hw_eec3_JetMET_2010A.root";
   string matched_hw_eec3_jet = matched_dir + "matched_hw_eec3_Jet_2010B.root";

   string out_hw_eec3_unf_allvertex = response_dir + "unfresp_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_jetmettau_unf_allvertex = response_dir + "unfresp_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_jetmet_unf_allvertex = response_dir + "unfresp_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_jet_unf_allvertex = response_dir + "unfresp_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_bbb_data = unfolded_dir + "Data_unfolded_BinByBin_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_bbb_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_bbb_hw_eec3 = unfolded_dir + "Herwig_TuneEEC3_unfolded_BinByBin_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_bbb_data_jetmettau = unfolded_dir + "Data_unfolded_BinByBin_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_BinByBin_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_data_jetmet = unfolded_dir + "Data_unfolded_BinByBin_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_BinByBin_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bbb_data_jet = unfolded_dir + "Data_unfolded_BinByBin_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_bbb_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_bbb_hw_eec3_jet = unfolded_dir + "Herwig_TuneEEC3_unfolded_BinByBin_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_data = unfolded_dir + "Data_unfolded_TUnfold_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_hw_eec3 = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_data_jetmettau = unfolded_dir + "Data_unfolded_TUnfold_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";

//method test
   string out_hw_eec3_unfolded_tunfold1_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold1_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold2_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold2_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold3_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold3_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold4_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold4_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";

   string out_hw_eec3_unfolded_tunfold_data_jetmet = unfolded_dir + "Data_unfolded_TUnfold_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";

//method test
   string out_hw_eec3_unfolded_tunfold1_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold1_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold2_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold2_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold3_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold3_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_tunfold4_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold4_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";

   string out_hw_eec3_unfolded_tunfold_data_jet = unfolded_dir + "Data_unfolded_TUnfold_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_tunfold_hw_eec3_jet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_svd_data = unfolded_dir + "Data_unfolded_SVD_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_svd_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_svd_hw_eec3 = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_svd_data_jetmettau = unfolded_dir + "Data_unfolded_SVD_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_p8_4c_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_hw_23_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";

//method test
   string out_hw_eec3_unfolded_svd1_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD1_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd2_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD2_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd3_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD3_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd4_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD4_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";

   string out_hw_eec3_unfolded_svd_data_jetmet = unfolded_dir + "Data_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_data_jetmet = unfolded_dir + "Data_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";

//method test
   string out_hw_eec3_unfolded_svd1_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD1_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd2_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD2_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd3_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD3_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_svd4_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD4_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";

   string out_hw_eec3_unfolded_bayes_data = unfolded_dir + "Data_unfolded_Bayes_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_bayes_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_bayes_hw_eec3 = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Herwig_TuneEEC3_allvertex.root";
   string out_hw_eec3_unfolded_bayes_data_jetmettau = unfolded_dir + "Data_unfolded_Bayes_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";

//method test
   string out_hw_eec3_unfolded_bayes1_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes1_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes2_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes2_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes3_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes3_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes4_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes4_Herwig_TuneEEC3_JetMETTau_2010A_allvertex.root";


   string out_hw_eec3_unfolded_bayes_data_jetmet = unfolded_dir + "Data_unfolded_Bayes_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";

//method test
   string out_hw_eec3_unfolded_bayes1_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes1_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes2_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes2_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes3_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes3_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";
   string out_hw_eec3_unfolded_bayes4_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes4_Herwig_TuneEEC3_JetMET_2010A_allvertex.root";

   string out_hw_eec3_unfolded_bayes_data_jet = unfolded_dir + "Data_unfolded_Bayes_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_bayes_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_TuneEEC3_Jet_2010B_allvertex.root";
   string out_hw_eec3_unfolded_bayes_hw_eec3_jet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Herwig_TuneEEC3_Jet_2010B_allvertex.root";


// Herwig pp 23
   string mc_hw_23[1];
   double lumi_hw_23[1];
   int n_files_hw_23 = 1;
   mc_hw_23[0] = ntuple_pythia6_v16_path + "QCD_Pt-15to3000_Tune23_Flat_7TeV_herwigpp_V16.root";
   lumi_hw_23[0] = 1.0 * 2.31e10/9966899.0;

   string vertex_hw_23_allvertex = vertex_dir + "v0/vertex_Herwigpp_Tune23_allvertex.root";
   string vertex_hw_23_allvertex_jetmettau_v1 = vertex_dir + "v1/vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_v1.root";
   string vertex_hw_23_allvertex_jetmet_v1 = vertex_dir + "v1/vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_v1.root";      
   string vertex_hw_23_allvertex_jet_v1 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_v1.root";
   string vertex_hw_23_allvertex_jetmettau_v2 = vertex_dir + "v2/vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_v2.root";
   string vertex_hw_23_allvertex_jetmet_v2 = vertex_dir + "v2/vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_v2.root";      
   string vertex_hw_23_allvertex_jet_v2 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_v2.root";
   string vertex_hw_23_allvertex_jetmettau_v3 = vertex_dir + "v3/vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_v3.root";
   string vertex_hw_23_allvertex_jetmet_v3 = vertex_dir + "v3/vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_v3.root";      
   string vertex_hw_23_allvertex_jet_v3 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_v3.root";
   string vertex_hw_23_allvertex_jetmettau_v4 = vertex_dir + "v4/vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_v4.root";
   string vertex_hw_23_allvertex_jetmet_v4 = vertex_dir + "v4/vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_v4.root";      
   string vertex_hw_23_allvertex_jet_v4 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_v4.root";
   string vertex_weights_hw_23_jetmettau_v1 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_v1.root";
   string vertex_weights_hw_23_jetmet_v1 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_v1.root";
   string vertex_weights_hw_23_jet_v1 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_v1.root";
   string vertex_weights_hw_23_jetmettau_v2 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_v2.root";
   string vertex_weights_hw_23_jetmet_v2 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_v2.root"; 
   string vertex_weights_hw_23_jet_v2 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_v2.root";
   string vertex_weights_hw_23_jetmettau_v3 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_v3.root";
   string vertex_weights_hw_23_jetmet_v3 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_v3.root"; 
   string vertex_weights_hw_23_jet_v3 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_v3.root";
   string vertex_weights_hw_23_jetmettau_v4 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_v4.root";
   string vertex_weights_hw_23_jetmet_v4 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_v4.root";
   string vertex_weights_hw_23_jet_v4 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_v4.root";
   string vertex_weights_hw_23_jetmettau_v5 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_v5.root";
   string vertex_weights_hw_23_jetmet_v5 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_v5.root";
   string vertex_weights_hw_23_jet_v5 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_v5.root";
   string vertex_hw_23_allvertex_jetmettau_notnorm_v1 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_notnorm_v1.root";
   string vertex_hw_23_allvertex_jetmet_notnorm_v1 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_notnorm_v1.root";      
   string vertex_hw_23_allvertex_jet_notnorm_v1 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_notnorm_v1.root";
   string vertex_hw_23_allvertex_jetmettau_notnorm_v2 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_notnorm_v2.root";
   string vertex_hw_23_allvertex_jetmet_notnorm_v2 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_notnorm_v2.root";      
   string vertex_hw_23_allvertex_jet_notnorm_v2 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_notnorm_v2.root";
   string vertex_hw_23_allvertex_jetmettau_notnorm_v3 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_notnorm_v3.root";
   string vertex_hw_23_allvertex_jetmet_notnorm_v3 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_notnorm_v3.root";      
   string vertex_hw_23_allvertex_jet_notnorm_v3 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_notnorm_v3.root";
   string vertex_hw_23_allvertex_jetmettau_notnorm_v4 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMETTau_2010A_notnorm_v4.root";
   string vertex_hw_23_allvertex_jetmet_notnorm_v4 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_JetMET_2010A_notnorm_v4.root";      
   string vertex_hw_23_allvertex_jet_notnorm_v4 = vertex_dir + "vertex_Herwigpp_Tune23_allvertex_Jet_2010B_notnorm_v4.root";
   string vertex_weights_hw_23_jetmettau_notnorm_v1 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_notnorm_v1.root";
   string vertex_weights_hw_23_jetmet_notnorm_v1 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_notnorm_v1.root";
   string vertex_weights_hw_23_jet_notnorm_v1 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_notnorm_v1.root";
   string vertex_weights_hw_23_jetmettau_notnorm_v2 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_notnorm_v2.root";
   string vertex_weights_hw_23_jetmet_notnorm_v2 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_notnorm_v2.root"; 
   string vertex_weights_hw_23_jet_notnorm_v2 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_notnorm_v2.root";
   string vertex_weights_hw_23_jetmettau_notnorm_v3 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_notnorm_v3.root";
   string vertex_weights_hw_23_jetmet_notnorm_v3 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_notnorm_v3.root"; 
   string vertex_weights_hw_23_jet_notnorm_v3 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_notnorm_v3.root";
   string vertex_weights_hw_23_jetmettau_notnorm_v4 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_notnorm_v4.root";
   string vertex_weights_hw_23_jetmet_notnorm_v4 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_notnorm_v4.root";
   string vertex_weights_hw_23_jet_notnorm_v4 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_notnorm_v4.root";
   string vertex_weights_hw_23_jetmettau_notnorm_v5 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMETTau2010A_notnorm_v5.root";
   string vertex_weights_hw_23_jetmet_notnorm_v5 = vertex_weights_dir + "vertex_Herwigpp_Tune23_JetMET2010A_notnorm_v5.root";
   string vertex_weights_hw_23_jet_notnorm_v5 = vertex_weights_dir + "vertex_Herwigpp_Tune23_Jet2010B_notnorm_v5.root";
   string resolution_hw_23_jetmettau = resolution_dir + "resolution_Herwigpp_Tune23_JetMETTau_2010A.root";
   string resolution_hw_23_jetmet = resolution_dir + "resolution_Herwigpp_Tune23_JetMET_2010A.root";
   string resolution_hw_23_jet = resolution_dir + "resolution_Herwigpp_Tune23_Jet_2010B.root";
   string norm_hw_23 = pileup_norm_dir + "norm_Herwigpp_Tune23.root";
   string norm_hw_23_jetmettau = pileup_norm_dir + "norm_Herwigpp_Tune23_JetMETTau_2010A.root";
   string norm_hw_23_jetmet = pileup_norm_dir + "norm_Herwigpp_Tune23_JetMET_2010A.root";
   string norm_hw_23_jet = pileup_norm_dir + "norm_Herwigpp_Tune23_Jet_2010B.root";
   string norm_hw_23_notnorm_jetmettau = pileup_norm_dir + "norm_Herwigpp_Tune23_notnorm_JetMETTau_2010A.root";
   string norm_hw_23_notnorm_jetmet = pileup_norm_dir + "norm_Herwigpp_Tune23_notnorm_JetMET_2010A.root";
   string norm_hw_23_notnorm_jet = pileup_norm_dir + "norm_Herwigpp_Tune23_notnorm_Jet_2010B.root";

   string out_hw_23_gen_nopileup = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_nopileup.root";
   string out_hw_23_gen_allvertex = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_allvertex.root";
   string out_hw_23_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_nopileup_JetMETTau_2010A.root"; 
   string out_hw_23_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_allvertex_JetMETTau_2010A.root";
   string out_hw_23_gen_nopileup_jetmet = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_nopileup_JetMET_2010A.root"; 
   string out_hw_23_gen_allvertex_jetmet = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_allvertex_JetMET_2010A.root";
   string out_hw_23_gen_nopileup_jet = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_nopileup_Jet_2010B.root"; 
   string out_hw_23_gen_allvertex_jet = mc_gen_dir + "xsec_Herwigpp_Tune23_gen_allvertex_Jet_2010B.root";
   string out_hw_23_notnorm_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Herwigpp_Tune23_notnorm_gen_nopileup_JetMETTau_2010A.root"; 
   string out_hw_23_notnorm_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Herwigpp_Tune23_notnorm_gen_allvertex_JetMETTau_2010A.root";
   string out_hw_23_notnorm_gen_nopileup_jetmet = mc_gen_dir + "xsec_Herwigpp_Tune23_notnorm_gen_nopileup_JetMET_2010A.root"; 
   string out_hw_23_notnorm_gen_allvertex_jetmet = mc_gen_dir + "xsec_Herwigpp_Tune23_notnorm_gen_allvertex_JetMET_2010A.root";
   string out_hw_23_notnorm_gen_nopileup_jet = mc_gen_dir + "xsec_Herwigpp_Tune23_notnorm_gen_nopileup_Jet_2010B.root"; 
   string out_hw_23_notnorm_gen_allvertex_jet = mc_gen_dir + "xsec_Herwigpp_Tune23_notnorm_gen_allvertex_Jet_2010B.root";
   string out_hw_23_gen_norm = output_dir + "histograms/normalized_mc/xsec_hw_23_gen.root";
   string out_hw_23_gen_jetmettau_norm = mc_norm_dir + "xsec_hw_23_JetMETTau_2010A_gen.root";
   string out_hw_23_gen_jetmet_norm = mc_norm_dir + "xsec_hw_23_JetMET_2010A_gen.root";
   string out_hw_23_gen_jet_norm = mc_norm_dir + "xsec_hw_23_Jet_2010B_gen.root";
   string out_hw_23_det_1vertex = mc_det_dir + "xsec_Herwigpp_Tune23_det_1vertex.root"; 
   string out_hw_23_det_allvertex = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex.root";
   string out_hw_23_det_allvertex_up = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_up.root";
   string out_hw_23_det_allvertex_down = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_down.root";
   string out_hw_23_det_1vertex_jetmettau = mc_det_dir + "xsec_Herwigpp_Tune23_det_1vertex_JetMETTau_2010A.root"; 
   string out_hw_23_det_allvertex_jetmettau = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_JetMETTau_2010A.root";
   string out_hw_23_det_allvertex_jetmettau_up = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_JetMETTau_2010A_up.root";
   string out_hw_23_det_allvertex_jetmettau_down = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_JetMETTau_2010A_down.root";
   string out_hw_23_det_1vertex_jetmet = mc_det_dir + "xsec_Herwigpp_Tune23_det_1vertex_JetMET_2010A.root"; 
   string out_hw_23_det_allvertex_jetmet = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_JetMET_2010A.root";
   string out_hw_23_det_allvertex_jetmet_up = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_JetMET_2010A_up.root";
   string out_hw_23_det_allvertex_jetmet_down = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_JetMET_2010A_down.root";
   string out_hw_23_det_1vertex_jet = mc_det_dir + "xsec_Herwigpp_Tune23_det_1vertex_Jet_2010B.root"; 
   string out_hw_23_det_allvertex_jet = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_Jet_2010B.root";
   string out_hw_23_det_allvertex_jet_up = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_Jet_2010B_up.root";
   string out_hw_23_det_allvertex_jet_down = mc_det_dir + "xsec_Herwigpp_Tune23_det_allvertex_Jet_2010B_down.root";
   string out_hw_23_det_norm = mc_norm_dir + "xsec_hw_23_det.root";
   string out_hw_23_det_jetmettau_norm = mc_norm_dir + "xsec_hw_23_JetMETTau_2010A_det.root";
   string out_hw_23_det_jetmet_norm = mc_norm_dir + "xsec_hw_23_JetMET_2010A_det.root";
   string out_hw_23_det_jet_norm = mc_norm_dir + "xsec_hw_23_Jet_2010B_det.root";


   string corr_pileup_hw_23 = corrections_dir + "correction_pileup_hw_23.root";
   string corr_pileup_hw_23_jetmettau = corrections_dir + "correction_pileup_hw_23_JetMETTau_2010A.root";
   string corr_pileup_hw_23_jetmet = corrections_dir + "correction_pileup_hw_23_JetMET_2010A.root";
   string corr_pileup_hw_23_jet = corrections_dir + "correction_pileup_hw_23_Jet_2010B.root";
   string corr_detector_hw_23 = corrections_dir + "correction_detector_hw_23.root";
   string corr_detector_hw_23_jetmettau = corrections_dir + "correction_detector_hw_23_JetMETTau_2010A.root";
   string corr_detector_hw_23_jetmet = corrections_dir + "correction_detector_hw_23_JetMET_2010A.root";
   string corr_detector_hw_23_jet = corrections_dir + "correction_detector_hw_23_Jet_2010B.root";
   string corr_final_hw_23 = corrections_dir + "correction_final_hw_23.root";
   string corr_final_hw_23_jetmettau = corrections_dir + "correction_final_hw_23_JetMETTau_2010A.root";
   string corr_final_hw_23_jetmet = corrections_dir + "correction_final_hw_23_JetMET_2010A.root";
   string corr_final_hw_23_jet = corrections_dir + "correction_final_hw_23_Jet_2010B.root";
   string corrected_hw_23_jetmettau = closure_dir + "closure_final_hw_23_JetMETTau_2010A.root";
   string corrected_hw_23_jetmet = closure_dir + "closure_final_hw_23_JetMET_2010A.root";
   string corrected_hw_23_jet = closure_dir + "closure_final_hw_23_Jet_2010B.root";
   string closure_ratio_hw_23_jetmettau = closure_dir + "closure_ratio_hw_23_JetMETTau_2010A.root";
   string closure_ratio_hw_23_jetmet = closure_dir + "closure_ratio_hw_23_JetMET_2010A.root";
   string closure_ratio_hw_23_jet = closure_dir + "closure_ratio_hw_23_Jet_2010B.root";
   string matched_hw_23 = matched_dir + "matched_hw_23.root";
   string matched_hw_23_jetmettau = matched_dir + "matched_hw_23_JetMETTau_2010A.root";
   string matched_hw_23_jetmet = matched_dir + "matched_hw_23_JetMET_2010A.root";
   string matched_hw_23_jet = matched_dir + "matched_hw_23_Jet_2010B.root";

   string out_hw_23_unf_allvertex = response_dir + "unfresp_Herwig_Tune23_allvertex.root";
   string out_hw_23_jetmettau_unf_allvertex = response_dir + "unfresp_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_jetmet_unf_allvertex = response_dir + "unfresp_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_jet_unf_allvertex = response_dir + "unfresp_Herwig_Tune23_Jet_2010B_allvertex.root";
   string out_hw_23_unfolded_bbb_data = unfolded_dir + "Data_unfolded_BinByBin_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_bbb_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_bbb_hw_23 = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_bbb_data_jetmettau = unfolded_dir + "Data_unfolded_BinByBin_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_hw_eec3_jetmettau = unfolded_dir + "Harwig_TuneEEC3_TuneZ2star_unfolded_BinByBin_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_data_jetmet = unfolded_dir + "Data_unfolded_BinByBin_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_BinByBin_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_TuneZ2star_unfolded_BinByBin_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bbb_data_jet = unfolded_dir + "Data_unfolded_BinByBin_Herwig_Tune23_Jet_2010B_allvertex.root";
   string out_hw_23_unfolded_bbb_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_BinByBin_Herwig_Tune23_Jet_2010B_allvertex.root";
   string out_hw_23_unfolded_bbb_hw_23_jet = unfolded_dir + "Herwig_Tune23_unfolded_BinByBin_Herwig_Tune23_Jet_2010B_allvertex.root";
   string out_hw_23_unfolded_tunfold_data = unfolded_dir + "Data_unfolded_TUnfold_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_tunfold_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_tunfold_hw_23 = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_tunfold_data_jetmettau = unfolded_dir + "Data_unfolded_TUnfold_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Herwig_Tune23_JetMETTau_2010A_allvertex.root";

//method test
   string out_hw_23_unfolded_tunfold1_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold1_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold2_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold2_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold3_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold3_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold4_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold4_Herwig_Tune23_JetMETTau_2010A_allvertex.root";

   string out_hw_23_unfolded_tunfold_data_jetmet = unfolded_dir + "Data_unfolded_TUnfold_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_TUnfold_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_TUnfold_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold_data_jet = unfolded_dir + "Data_unfolded_TUnfold_Herwig_Tune23_JetMET_2010B_allvertex.root";
   string out_hw_23_unfolded_tunfold_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_TUnfold_Herwig_Tune23_JetMET_2010B_allvertex.root";
   string out_hw_23_unfolded_tunfold_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold_Herwig_Tune23_JetMET_2010A_allvertex.root";

//method test
   string out_hw_23_unfolded_tunfold1_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold1_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold2_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold2_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold3_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold3_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_tunfold4_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_TUnfold4_Herwig_Tune23_JetMET_2010A_allvertex.root";


   string out_hw_23_unfolded_svd_data = unfolded_dir + "Data_unfolded_SVD_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_svd_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_svd_hw_23 = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_svd_data_jetmettau = unfolded_dir + "Data_unfolded_SVD_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Herwig_Tune23_JetMETTau_2010A_allvertex.root";

//method test
   string out_hw_23_unfolded_svd1_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_SVD1_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd2_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_SVD2_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd3_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_SVD3_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_svd4_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_SVD4_Herwig_Tune23_JetMETTau_2010A_allvertex.root";

   string out_hw_23_unfolded_svd_data_jetmet = unfolded_dir + "Data_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_data_jetmet = unfolded_dir + "Data_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD_Herwig_Tune23_JetMET_2010A_allvertex.root";

//method test
   string out_hw_23_unfolded_svd1_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD1_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd2_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD2_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd3_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD3_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_svd4_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_SVD4_Herwig_Tune23_JetMET_2010A_allvertex.root";

   string out_hw_23_unfolded_bayes_data = unfolded_dir + "Data_unfolded_Bayes_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_bayes_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_bayes_hw_23 = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Herwig_Tune23_allvertex.root";
   string out_hw_23_unfolded_bayes_data_jetmettau = unfolded_dir + "Data_unfolded_Bayes_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_hw_eec3_jetmettau = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Herwig_Tune23_JetMETTau_2010A_allvertex.root";

//method test
   string out_hw_23_unfolded_bayes1_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes1_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes2_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes2_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes3_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes3_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes4_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes4_Herwig_Tune23_JetMETTau_2010A_allvertex.root";

   string out_hw_23_unfolded_bayes_hw_23_jetmettau = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Herwig_Tune23_JetMETTau_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_data_jetmet = unfolded_dir + "Data_unfolded_Bayes_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Bayes_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_hw_eec3_jetmet = unfolded_dir + "Herwig_TuneEEC3_unfolded_Bayes_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Herwig_Tune23_JetMET_2010A_allvertex.root";

//method test
   string out_hw_23_unfolded_bayes1_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes1_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes2_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes2_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes3_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes3_Herwig_Tune23_JetMET_2010A_allvertex.root";
   string out_hw_23_unfolded_bayes4_hw_23_jetmet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes4_Herwig_Tune23_JetMET_2010A_allvertex.root";

   string out_hw_23_unfolded_bayes_data_jet = unfolded_dir + "Data_unfolded_Bayes_Herwig_Tune23_Jet_2010B_allvertex.root";
   string out_hw_23_unfolded_bayes_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Bayes_Herwig_Tune23_Jet_2010B_allvertex.root";
   string out_hw_23_unfolded_bayes_hw_23_jet = unfolded_dir + "Herwig_Tune23_unfolded_Bayes_Herwig_Tune23_Jet_2010B_allvertex.root";
//Data Ntuples locations and luminosities

//JetMETTauMonitor_2010A
   string jetmettaumon[1];
   double lumi_jetmettaumon[2];
   //jetmettaumon[0] = ntuple_trigeff_path + "JetTreeJETMETTAUMONITOR.root";
   jetmettaumon[0] = "../../../../../../../../work/c/cipriano/public/JetTreeJETMETTAUMONITOR.root";
   string systematics_trigger_jetmettaumon = trigger_syst_dir + "systematic_trigger_JetMETTauMonitor_2010A.root";
   string out_jetmettaumon_emu = trigger_dir + "events_JetMETTauMonitor2010A_emulated.root";
   string correction_jetmettaumon_emu_v1 = trigger_corr_dir + "correction_JetMETTauMonitor2010A_emulated_v1.txt";
   string correction_sj_jetmettaumon_emu_v1 = trigger_corr_dir + "correction_sj_JetMETTauMonitor2010A_emulated_v1.txt";
   string check_jetmettaumon_emu_v1 = trigger_dir + "check_JetMETTauMonitor2010A_emulated_v1.root";
   string check_sj_jetmettaumon_emu_v1 = trigger_dir + "check_sj_JetMETTauMonitor2010A_emulated_v1.root";
   string correction_jetmettaumon_emu_v2 = trigger_corr_dir + "correction_JetMETTauMonitor2010A_emulated_v2.txt";
   string correction_sj_jetmettaumon_emu_v2 = trigger_corr_dir + "correction_sj_JetMETTauMonitor2010A_emulated_v2.txt";
   string check_jetmettaumon_emu_v2 = trigger_dir + "check_JetMETTauMonitor2010A_emulated_v2.root";
   string check_sj_jetmettaumon_emu_v2 = trigger_dir + "check_sj_JetMETTauMonitor2010A_emulated_v2.root";
   int n_files_jetmettaumon = 1;
   lumi_jetmettaumon[0] = 1.0; // /pb cipriano

//JetMETTau_2010A
   string jetmettau[3];
   double lumi_jetmettau[3];
   //jetmettau[0] = ntuple_data_path + "JetMETTau_Run2010A_Apr21ReReco_v1_AOD/JetMETTau_Run2010A_Apr21ReReco_v1_AOD_ProcessedTree_data.root"; //new
   jetmettau[0] = "/nfs/dust/cms/user/haeverm/QCDntuple_JetMETTau_Run2010A-Apr21ReReco-v1_AOD_NewJEC-Jan2015_BTag_v1/JetMerged.root";
   jetmettau[1] = jetmettau[0];
   jetmettau[2] = jetmettau[0];
   string systematics_trigger_jetmettau = trigger_syst_dir + "systematic_trigger_JetMETTau_2010A.root";
   string vertex_jetmettau_allvertex = vertex_dir + "vertex_JetMETTau2010A_allvertex.root";
   string out_jetmettau_1vertex = raw_data_dir + "xsec_JetMETTau2010A_1vertex.root";
   string out_jetmettau_allvertex = raw_data_dir + "xsec_JetMETTau2010A_allvertex.root";
   string out_jetmettau_up = raw_data_dir + "xsec_JetMETTau2010A_up.root";
   string out_jetmettau_down = raw_data_dir + "xsec_JetMETTau2010A_down.root";
   string out_jetmettau_emu = trigger_dir + "events_JetMETTau2010A_emulated.root";
   string correction_jetmettau_emu_v1 = trigger_corr_dir + "correction_JetMETTau2010A_emulated_v1.txt";
   string correction_sj_jetmettau_emu_v1 = trigger_corr_dir + "correction_sj_JetMETTau2010A_emulated_v1.txt";
   string check_jetmettau_emu_v1 = trigger_dir + "check_JetMETTau2010A_emulated_v1.root";
   string check_sj_jetmettau_emu_v1 = trigger_dir + "check_sj_JetMETTau2010A_emulated_v1.root";
   string correction_jetmettau_emu_v2 = trigger_corr_dir + "correction_JetMETTau2010A_emulated_v2.txt";
   string correction_sj_jetmettau_emu_v2 = trigger_corr_dir + "correction_sj_JetMETTau2010A_emulated_v2.txt";
   string check_jetmettau_emu_v2 = trigger_dir + "check_JetMETTau2010A_emulated_v2.root";
   string check_sj_jetmettau_emu_v2 = trigger_dir + "check_sj_JetMETTau2010A_emulated_v2.root";
   string jes_up_unc_jetmettau = jes_unc_dir + "jes_up_unc_JetMETTau_2010A.root";
   string jes_down_unc_jetmettau = jes_unc_dir + "jes_down_unc_JetMETTau_2010A.root";
   string jer_up_unc_jetmettau = jer_unc_dir + "jer_up_unc_JetMETTau_2010A.root";
   string jer_down_unc_jetmettau = jer_unc_dir + "jer_down_unc_JetMETTau_2010A.root";
   int n_files_jetmettau = 3;
   //lumi_jetmettau[0] = 0.284161; // /pb cipriano
   //lumi_jetmettau[0] = 0.009645; //walter
   //lumi_jetmettau[1] = 0.192895; //walter
   //lumi_jetmettau[2] = 2.869; //walter
   lumi_jetmettau[0] = 0.013799; //cipriano
   lumi_jetmettau[1] = 0.117223; //cipriano
   lumi_jetmettau[2] = 0.278789; //cipriano
   
   //JetMET_2010A
   string jetmet[3];
   double lumi_jetmet[3];
   jetmet[0] = "/nfs/dust/cms/user/haeverm/QCDntuple_JetMET_Run2010A-Apr21ReReco-v1_AOD_NewJEC-Jan2015_BTag_v1/JetMerged.root";
   jetmet[1] = jetmet[0];
   jetmet[2] = jetmet[0];
   string systematics_trigger_jetmet = trigger_syst_dir + "systematic_trigger_JetMET_2010A.root";
   string vertex_jetmet_allvertex = vertex_dir + "vertex_JetMET2010A_allvertex.root";
   string out_jetmet_1vertex = raw_data_dir + "xsec_JetMET2010A_1vertex.root";
   string out_jetmet_allvertex = raw_data_dir + "xsec_JetMET2010A_allvertex.root";
   string out_jetmet_up = raw_data_dir + "xsec_JetMET2010A_up.root";
   string out_jetmet_down = raw_data_dir + "xsec_JetMET2010A_down.root";
   string out_jetmet_emu = trigger_dir + "events_JetMET2010A_emulated.root";
   string correction_jetmet_emu_v1 = trigger_corr_dir + "correction_JetMET2010A_emulated_v1.txt";
   string correction_sj_jetmet_emu_v1 = trigger_corr_dir + "correction_sj_JetMET2010A_emulated_v1.txt";
   string check_jetmet_emu_v1 = trigger_dir + "check_JetMET2010A_emulated_v1.root";
   string check_sj_jetmet_emu_v1 = trigger_dir + "check_sj_JetMET2010A_emulated_v1.root";
   string correction_jetmet_emu_v2 = trigger_corr_dir + "correction_JetMET2010A_emulated_v2.txt";
   string correction_sj_jetmet_emu_v2 = trigger_corr_dir + "correction_sj_JetMET2010A_emulated_v2.txt";
   string check_jetmet_emu_v2 = trigger_dir + "check_JetMET2010A_emulated_v2.root";
   string check_sj_jetmet_emu_v2 = trigger_dir + "check_sj_JetMET2010A_emulated_v2.root";
   string jes_up_unc_jetmet = jes_unc_dir + "jes_up_unc_JetMET_2010A.root";
   string jes_down_unc_jetmet = jes_unc_dir + "jes_down_unc_JetMET_2010A.root";
   string jer_up_unc_jetmet = jer_unc_dir + "jer_up_unc_JetMET_2010A.root";
   string jer_down_unc_jetmet = jer_unc_dir + "jer_down_unc_JetMET_2010A.root";
   int n_files_jetmet = 3;
   //lumi_jetmet[0] = 2.896; // /pb cipriano
   //lumi_jetmet[0] = 0.013799; // walter
   //lumi_jetmet[1] = 0.117223; // walter
   //lumi_jetmet[2] = 0.278789; // walter
   lumi_jetmet[0] = 0.009645; // cipriano
   lumi_jetmet[1] = 0.192895; // cipriano
   lumi_jetmet[2] = 2.869; // cipriano
   
   //Jet_2010B
   string jet[3];
   double lumi_jet[3];
   jet[0] = ntuple_data_path + "Jet_Run2010B_Apr21ReReco_v1_AOD_New/JET_2010_AOD_21APR_MERGE_ALL_ROOT_ProcessedTree_data.root";
   jet[1] = jet[0];
   jet[2] = jet[0];   
   string systematics_trigger_jet = trigger_syst_dir + "systematic_trigger_Jet2010B.root";
   string vertex_jet_allvertex = vertex_dir + "vertex_Jet2010B_allvertex.root";
   string out_jet_1vertex = raw_data_dir + "xsec_Jet2010B_1vertex.root";
   string out_jet_allvertex = raw_data_dir + "xsec_Jet2010B_allvertex.root";
   string out_jet_up = raw_data_dir + "xsec_Jet2010B_up.root";
   string out_jet_down = raw_data_dir + "xsec_Jet2010B_down.root";
   string out_jet_emu = trigger_dir + "events_Jet2010B_emulated.root";
   string correction_jet_emu_v1 = trigger_corr_dir + "correction_Jet2010B_emulated_v1.txt";
   string correction_sj_jet_emu_v1 = trigger_corr_dir + "correction_sj_Jet2010B_emulated_v1.txt";
   string check_jet_emu_v1 = trigger_dir + "check_Jet2010B_emulated_v1.root";
   string check_sj_jet_emu_v1 = trigger_dir + "check_sj_Jet2010B_emulated_v1.root";
   string correction_jet_emu_v2 = trigger_corr_dir + "correction_Jet2010B_emulated_v2.txt";
   string correction_sj_jet_emu_v2 = trigger_corr_dir + "correction_sj_Jet2010B_emulated_v2.txt";
   string check_jet_emu_v2 = trigger_dir + "check_Jet2010B_emulated_v2.root";
   string check_sj_jet_emu_v2 = trigger_dir + "check_sj_Jet2010B_emulated_v2.root";
   string jes_up_unc_jet = jes_unc_dir + "jes_up_unc_Jet_2010B.root";
   string jes_down_unc_jet = jes_unc_dir + "jes_down_unc_Jet_2010B.root";
   string jer_up_unc_jet = jer_unc_dir + "jer_up_unc_Jet_2010B.root";
   string jer_down_unc_jet = jer_unc_dir + "jer_down_unc_Jet_2010B.root";
   int n_files_jet = 3;
   //lumi_jet[0] = 32.072; // /pb cipriano
   //lumi_jet[0] = 0.001855; //walter
   //lumi_jet[1] = 0.026783; //walter
   //lumi_jet[2] = 0.239874; //walter
   lumi_jet[0] = 0.001855;
   lumi_jet[1] = 0.026783;
   lumi_jet[2] = 0.239874;

//histogram locations
   string corrected_jetmettau = corrected_dir + "corrected_xsec_data_JetMETTau_2010A.root";
   string corrected_jetmet = corrected_dir + "corrected_xsec_data_JetMET_2010A.root";
   string corrected_jet = corrected_dir + "corrected_xsec_data_Jet_2010B.root";

   string merged_uncorrected_data = histograms_dir + "xsec_uncorrected_data_2010.root";
   string merged_data = histograms_dir + "xsec_data_2010.root";
   string merged_unfolded_data = histograms_dir + "xsec_unfolded_data_2010.root";
   string merged_uncorr_unfolded_data = histograms_dir + "xsec_uncorr_unfolded_data_2010.root";
   string ratio_merged_data = ratio_dir + "xsec_data_2010.root";
   string ratio_merged_unfolded_data = ratio_dir + "xsec_unfolded_data_2010.root";
   string data_unc = histograms_dir + "xsec_data_2010_unc.root";
   string ratio_data_unc = ratio_dir + "xsec_data_2010_unc.root";
   string merged_trig_comb_syst = trigger_syst_dir + "trig_comb_syst_data_2010.root";

   string model_unc_jetmettau = model_unc_dir + "model_uncertainty_JetMETTau_2010A.root";
   string model_unc_jetmet = model_unc_dir + "model_uncertainty_JetMET_2010A.root";
   string model_unc_jet = model_unc_dir + "model_uncertainty_Jet_2010B.root";
   string model_unc_merged = model_unc_dir + "model_uncertainty.root";
   string jes_unc_up_merged = jes_unc_dir + "jes_uncertainty_up.root";
   string jes_unc_down_merged = jes_unc_dir + "jes_uncertainty_down.root";
   string jer_unc_p6_jetmettau = jer_unc_dir + "jer_uncertainty_p6_jetmettau.root";
   string jer_unc_p8_jetmettau = jer_unc_dir + "jer_uncertainty_p8_jetmettau.root";
   string jer_unc_p6_jetmet = jer_unc_dir + "jer_uncertainty_p6_jetmet.root";
   string jer_unc_p8_jetmet = jer_unc_dir + "jer_uncertainty_p8_jetmet.root";
   string jer_unc_p6_jet = jer_unc_dir + "jer_uncertainty_p6_jet.root";
   string jer_unc_p8_jet = jer_unc_dir + "jer_uncertainty_p8_jet.root";
   string jer_unc_p6_merged = jer_unc_dir + "jer_uncertainty_p6_merged.root";
   string jer_unc_p8_merged = jer_unc_dir + "jer_uncertainty_p8_merged.root";
   string jer_unc_merged = jer_unc_dir + "jer_uncertainty_merged.root";
   string corr_unc_up_merged = corr_unc_dir + "corr_uncertainty_up.root";
   string corr_unc_down_merged = corr_unc_dir + "corr_uncertainty_down.root";
   string uncorr_unc_merged = uncorr_unc_dir + "uncorr_uncertainty.root";
   string total_unc_up = total_unc_dir + "total_uncertainty_up.root";
   string total_unc_down = total_unc_dir + "total_uncertainty_down.root";
   string total_unc = total_unc_dir + "total_uncertainty.root";

//monte carlo predictions locations
   string albert_public_dir = "../../../../../../../../work/k/knutsson/public/mcforpedro/v2/1M/";
   string mc_dir = "../mc/1M/";
   string pedro_generated_dir = "../output/histograms/mc_generated/";

   string mc_pred_herwig6 = mc_dir + "herwig6_pt30.root";
   string mc_pred_herwigpp = mc_dir + "herwigpp_pt30.root";
   string mc_pred_py6_ambt1 = mc_dir + "py6_ambt1_pt30.root";
   string mc_pred_py6_p11 = mc_dir + "py6_p11_pt30.root";
   string mc_pred_py6_z2_nompi = mc_dir + "py6_z2_nompi_pt30.root";
   string mc_pred_py6_z2 = mc_dir + "py6_z2_pt30.root";
   string mc_pred_py8_4c = mc_dir + "py8_4c_pt30.root";
   string mc_pred_py6_z2_nompi2 = pedro_generated_dir + "7tev_z2star_nompi_2.96984e08.root";
   string mc_pred_py6_z2_nompi_nohad = pedro_generated_dir + "7tev_z2star_nompi_nohad_2.964422e08.root";
   string mc_pred_py6_z2_nompi_nohad_nofsr = pedro_generated_dir + "7tev_z2star_nompi_nohad_nofsr_2.95915e08.root";
   string mc_pred_py6_z2_nompi_nohad_nofsr_noisr = pedro_generated_dir + "7tev_z2star_nompi_nohad_nofsr_noisr_2.96715e08.root";

   string mc_ratio_herwig6 = ratio_dir + "herwig6_pt30.root";
   string mc_ratio_herwigpp = ratio_dir + "herwigpp_pt30.root";
   string mc_ratio_py6_ambt1 = ratio_dir + "pythia6_ambt1_pt30.root";
   string mc_ratio_py6_p11 = ratio_dir + "pythia6_p11_pt30.root";
   string mc_ratio_py6_z2_nompi = ratio_dir + "pythia6_z2starnompi_pt30.root";
   string mc_ratio_py6_z2 = ratio_dir + "pythia6_z2star_pt30.root";
   string mc_ratio_py8_4c = ratio_dir + "pythia8_4c_pt30.root";
   string mc_ratio_py6_z2_off = ratio_dir + "pythia6_z2star_off.root";
   string mc_ratio_py8_4c_off = ratio_dir + "pythia8_4c_off.root";
   string mc_ratio_py6_z2_nompi2 = ratio_dir + "pythia6_z2_nompi2.root";
   string mc_ratio_py6_z2_nompi_nohad = ratio_dir + "pythia6_z2_nompi_nohad.root";
   string mc_ratio_py6_z2_nompi_nohad_nofsr = ratio_dir + "pythia6_z2_nompi_nohad_nofsr.root";
   string mc_ratio_py6_z2_nompi_nohad_nofsr_noisr = ratio_dir + "pythia6_z2_nompi_nohad_nofsr_noisr.root";

   string mc_ratio_py6_z2_jetmettau = ratio_dir + "pythia6_z2star_jetmettau.root";
   string mc_ratio_py6_z2_jetmet = ratio_dir + "pythia6_z2star_jetmet.root";
   string mc_ratio_py6_z2_jet = ratio_dir + "pythia6_z2star_jet.root";
   string mc_ratio_py8_4c_jetmettau = ratio_dir + "pythia8_4c_jetmemttau.root";
   string mc_ratio_py8_4c_jetmet = ratio_dir + "pythia8_4c_jetmet.root";
   string mc_ratio_py8_4c_jet = ratio_dir + "pythia8_4c_jet.root";

//predictions list
   string final_mc_list[20];
   string final_label_list[20];
   string final_prefix_list[20];

   final_mc_list[0] = out_p6_z2_gen_allvertex;
   final_mc_list[1] = out_p8_4c_gen_allvertex;
   final_mc_list[2] = mc_pred_herwig6;
   final_mc_list[3] = mc_pred_herwigpp;
   final_mc_list[4] = mc_pred_py6_p11;
   final_mc_list[5] = mc_pred_py6_z2_nompi;
   final_mc_list[6] = mc_pred_py6_z2;
   final_mc_list[7] = mc_pred_py6_ambt1;
   final_mc_list[8] = mc_pred_py8_4c;
   final_mc_list[9] = mc_pred_py6_z2_nompi_nohad;
   final_mc_list[10] = mc_pred_py6_z2_nompi_nohad_nofsr;
   final_mc_list[11] = mc_pred_py6_z2_nompi_nohad_nofsr_noisr;
   final_mc_list[12] = "";
   final_mc_list[13] = "";
   final_mc_list[14] = "";
   final_mc_list[15] = "";
   final_mc_list[16] = "";
   final_mc_list[17] = "";
   final_mc_list[18] = "";
   final_mc_list[19] = "";

   final_label_list[0] = "Pythia 6 - Z2* Tune (Official)";
   final_label_list[1] = "Pythia 8 - 4C Tune (Official)";
   final_label_list[2] = "Herwig 6";
   final_label_list[3] = "Herwig ++";
   final_label_list[4] = "Pythia 6 - P11 Tune";
   final_label_list[5] = "Pythia 6 - Z2* Tune (No MPI)";
   final_label_list[6] = "Pythia 6 - Z2* Tune";
   final_label_list[7] = "Pythia 6 - AMBT1 Tune";
   final_label_list[8] = "Pythia 8 - 4C";
   final_label_list[9] = "Same - No MPI, No Hadronization";
   final_label_list[10] = "Same - No MPI, No Hadronization, No FSR";
   final_label_list[11] = "Same - No MPI, No Hadronization, No FSR, No ISR";
   final_label_list[12] = "";
   final_label_list[13] = "";
   final_label_list[14] = "";
   final_label_list[15] = "";
   final_label_list[16] = "";
   final_label_list[17] = "";
   final_label_list[18] = "";
   final_label_list[19] = "";


   final_prefix_list[0] = "ak5Gen_";
   final_prefix_list[1] = "ak5Gen_";
   final_prefix_list[2] = "";
   final_prefix_list[3] = "";
   final_prefix_list[4] = "";
   final_prefix_list[5] = "";
   final_prefix_list[6] = "";
   final_prefix_list[7] = "";
   final_prefix_list[8] = "";
   final_prefix_list[9] = "";
   final_prefix_list[10] = "";
   final_prefix_list[11] = "";
   final_prefix_list[12] = "";
   final_prefix_list[13] = "";
   final_prefix_list[14] = "";
   final_prefix_list[15] = "";
   final_prefix_list[16] = "";
   final_prefix_list[17] = "";
   final_prefix_list[18] = "";
   final_prefix_list[19] = "";

//rations list
   string ratio_mc_list[20];
   string ratio_label_list[20];
   string ratio_prefix_list[20];

   ratio_mc_list[0] = mc_ratio_py6_z2_off;
   ratio_mc_list[1] = mc_ratio_py8_4c_off;
   ratio_mc_list[2] = mc_ratio_herwig6;
   ratio_mc_list[3] = mc_ratio_herwigpp;
   ratio_mc_list[4] = mc_ratio_py6_p11;
   ratio_mc_list[5] = mc_ratio_py6_z2_nompi;
   ratio_mc_list[6] = mc_ratio_py6_z2;
   ratio_mc_list[7] = mc_ratio_py6_ambt1;
   ratio_mc_list[8] = mc_ratio_py8_4c;
   ratio_mc_list[9] = mc_ratio_py6_z2_nompi_nohad;
   ratio_mc_list[12] = mc_ratio_py6_z2_nompi_nohad_nofsr;
   ratio_mc_list[11] = mc_ratio_py6_z2_nompi_nohad_nofsr_noisr;
   ratio_mc_list[12] = "";
   ratio_mc_list[13] = "";
   ratio_mc_list[14] = "";
   ratio_mc_list[15] = "";
   ratio_mc_list[16] = "";
   ratio_mc_list[17] = "";
   ratio_mc_list[18] = "";
   ratio_mc_list[19] = "";

   ratio_label_list[0] = "Pythia 6 - Z2* Tune (Official)";
   ratio_label_list[1] = "Pythia 8 - 4C Tune (Official)";
   ratio_label_list[2] = "Herwig 6";
   ratio_label_list[3] = "Herwig ++";
   ratio_label_list[4] = "Pythia 6 - P11 Tune";
   ratio_label_list[5] = "Pythia 6 - Z2* Tune (No MPI)";
   ratio_label_list[6] = "Pythia 6 - Z2* Tune";
   ratio_label_list[7] = "Pythia 6 - AMBT1 Tune";
   ratio_label_list[8] = "Pythia 8 - 4C";
   ratio_label_list[9] = "Same - No MPI, No Hadronization";
   ratio_label_list[10] = "Same - No MPI, No Hadronization, No FSR";
   ratio_label_list[11] = "Same - No MPI, No Hadronization, No FSR, No ISR";
   ratio_label_list[12] = "";
   ratio_label_list[13] = "";
   ratio_label_list[14] = "";
   ratio_label_list[15] = "";
   ratio_label_list[16] = "";
   ratio_label_list[17] = "";
   ratio_label_list[18] = "";
   ratio_label_list[19] = "";

   ratio_prefix_list[0] = "ak5Gen_";
   ratio_prefix_list[1] = "ak5Gen_";
   ratio_prefix_list[2] = "";
   ratio_prefix_list[3] = "";
   ratio_prefix_list[4] = "";
   ratio_prefix_list[5] = "";
   ratio_prefix_list[6] = "";
   ratio_prefix_list[7] = "";
   ratio_prefix_list[8] = "";
   ratio_prefix_list[9] = "";
   ratio_prefix_list[10] = "";
   ratio_prefix_list[11] = "";
   ratio_prefix_list[12] = "";
   ratio_prefix_list[13] = "";
   ratio_prefix_list[14] = "";
   ratio_prefix_list[15] = "";
   ratio_prefix_list[16] = "";
   ratio_prefix_list[17] = "";
   ratio_prefix_list[18] = "";
   ratio_prefix_list[19] = "";

//compile the different routines
if (show_steps) { cout << "Compiling code..."<<endl; }
  gROOT -> ProcessLine(".L create_directories.C++");
  gROOT -> ProcessLine(".L common_methods.h++");
  if (estimate_combination_systematic) { gROOT -> ProcessLine(".L estimate_combination_systematic.C++"); }
  if (get_trigger_efficiencies) { gROOT -> ProcessLine(".L get_triggered_events.C++"); }
  if (compute_trigger_efficiency) { gROOT -> ProcessLine(".L compute_trigger_efficiency.C++"); }
  if (compute_trigger_turn_on || check_trigger_correction_v1 || check_trigger_correction_v2 || apply_trigger_correction || plot_final_trigger_efficiency)
  {
  gROOT -> ProcessLine(".L compute_trigger_turn_on.C++");
  }
  if (counting_combination_statistics) { gROOT -> ProcessLine(".L counting_combination_statistics.C++"); }
  if (compute_trigger_correction_v1 || compute_trigger_correction_v2) { gROOT -> ProcessLine(".L compute_trigger_correction.C++"); }
  if (get_vertex_distribution_v0 || get_vertex_distribution_v1 || get_vertex_distribution_v2 || get_vertex_distribution_v3 || get_vertex_distribution_v4)
  {
  gROOT -> ProcessLine(".L get_vertex_distribution.C++");
  }
  if (compute_vertex_weights_v1 || compute_vertex_weights_v2 || compute_vertex_weights_v3 || compute_vertex_weights_v4 || compute_vertex_weights_v5)
  {
  gROOT -> ProcessLine(".L compute_vertex_weights.C++");
  gROOT -> ProcessLine(".L plot_status_vertex_reweight.C++");
  }
  if (plot_status_vertex_reweight) { gROOT -> ProcessLine(".L plot_evolution_vertex_reweight.C++"); }
  if (compute_resolution) { gROOT -> ProcessLine(".L compute_resolution.C++"); }
  if (plot_resolution) { gROOT -> ProcessLine(".L plot_resolution.C++"); }
  if (get_pileup_normalization) { gROOT -> ProcessLine(".L get_pileup_normalization.C++"); }
  if (read_mc_ntuples_gen || read_mc_ntuples_det || read_data_ntuples) { gROOT -> ProcessLine(".L read_ntuple.C++"); }
  if (check_cross_section) { gROOT -> ProcessLine(".L check_cross_section.C++"); }
  if (normalize_mc) { gROOT -> ProcessLine(".L normalize_mc.C++"); }
  
  if (plot_control_dist)
  {
  gROOT -> ProcessLine(".L plot_mc_gen_level.C++");
  gROOT -> ProcessLine(".L control_plots.C++");
  }
  if (compute_corrections)   { gROOT -> ProcessLine(".L compute_correction.C++"); }
  if (get_mc_matched_events)   { gROOT -> ProcessLine(".L get_matched_events.C++"); }
  if (compute_psba)   { gROOT -> ProcessLine(".L compute_psba.C++"); }
  if (apply_corrections) { gROOT -> ProcessLine(".L apply_correction.C++"); }
  if (create_unfolding_response) { gROOT -> ProcessLine(".L create_unfolding_response.C++"); }
  if (check_response_matrix) { gROOT -> ProcessLine(".L check_response_matrix.C++"); }
  if (unfold) { gROOT -> ProcessLine(".L unfolding.C++"); }
  if (compare_unfolding_results) { gROOT -> ProcessLine(".L compare_unfolding_results.C++"); }
  if (merge_data)
  {
  gROOT -> ProcessLine(".L merge_data.C++");
  gROOT -> ProcessLine(".L merge_unfolded_data.C++");
  }
  if (compute_model_uncertainty) { gROOT -> ProcessLine(".L estimate_model_uncertainty.C++"); }
  if (compute_jes_uncertainty) { gROOT -> ProcessLine(".L estimate_jes_uncertainty.C++"); }
  if (compute_jer_uncertainty) { gROOT -> ProcessLine(".L estimate_jer_uncertainty.C++"); }
  if (compute_corr_uncertainty) { gROOT -> ProcessLine(".L compute_corr_uncertainty.C++"); }
  if (compute_model_uncertainty || compute_jes_uncertainty || compute_jer_uncertainty) { gROOT -> ProcessLine(".L merge_uncertainties.C++"); }
  if (run_stat_vs_unf_unc) { gROOT -> ProcessLine(".L stat_vs_unf_unc.C++"); }
  if (compute_uncorr_uncertainty) { gROOT -> ProcessLine(".L compute_uncorr_uncertainty.C++"); }
  if (compute_total_uncertainty) { gROOT -> ProcessLine(".L total_uncertainty.C++"); }
  if (apply_uncertainty) { gROOT -> ProcessLine(".L apply_uncertainty.C++"); }
  if (create_ratio) { gROOT -> ProcessLine(".L create_ratios.C++"); }
  if (do_final_plots) { gROOT -> ProcessLine(".L final_plots.C++"); }



//get the trigger efficiencies
if (get_trigger_efficiencies)
{
if (show_steps) { cout << "Get Trigger Efficiencies..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//get_triggered_events(jetmettaumon, out_jetmettaumon_none, n_files_jetmettaumon, "NONE", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U, n_files_jetmettaumon, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U, n_files_jetmettaumon, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet15U, n_files_jetmettaumon, "HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jetmettaumon, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jetmettaumon, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet30U, n_files_jetmettaumon, "HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jetmettaumon, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jetmettaumon, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U, n_files_jetmettaumon, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet50U, n_files_jetmettaumon, "HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jetmettaumon, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jetmettaumon, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U, n_files_jetmettaumon, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_DiJetAve15U, n_files_jetmettaumon, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_DiJetAve30U, n_files_jetmettaumon, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_FwdJet20U, n_files_jetmettaumon, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_DoubleJet15U_ForwardBackward, n_files_jetmettaumon, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_all, n_files_jetmettaumon, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "JetMETTau"<<endl; }
//get_triggered_events(jetmettau, out_jetmettau_none, n_files_jetmettau, "NONE", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U, n_files_jetmettau, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U, n_files_jetmettau, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet15U, n_files_jetmettau, "HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jetmettau, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jetmettau, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet30U, n_files_jetmettau, "HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jetmettau, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jetmettau, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet15U_AND_HLT_Jet30U, n_files_jetmettau, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet50U, n_files_jetmettau, "HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jetmettau, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jetmettau, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet15U_AND_HLT_Jet50U, n_files_jetmettau, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_DiJetAve15U, n_files_jetmettau, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_DiJetAve30U, n_files_jetmettau, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_FwdJet20U, n_files_jetmettau, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_DoubleJet15U_ForwardBackward, n_files_jetmettau, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_all, n_files_jetmettau, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
//get_triggered_events(jetmet, out_jetmet_none, n_files_jetmet, "NONE", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U, n_files_jetmet, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U, n_files_jetmet, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_Jet15U, n_files_jetmet, "HLT_Jet15U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jetmet, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jetmet, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_Jet30U, n_files_jetmet, "HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jetmet, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jetmet, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_Jet15U_AND_HLT_Jet30U, n_files_jetmet, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_Jet50U, n_files_jetmet, "HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jetmet, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jetmet, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_Jet15U_AND_HLT_Jet50U, n_files_jetmet, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_DiJetAve15U, n_files_jetmet, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_DiJetAve30U, n_files_jetmet, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_FwdJet20U, n_files_jetmet, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_DoubleJet15U_ForwardBackward, n_files_jetmet, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_all, n_files_jetmet, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "Jet"<<endl; }
//get_triggered_events(jet, out_jet_none, n_files_jet, "NONE", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet6U, n_files_jet, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet10U, n_files_jet, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet15U, n_files_jet, "HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jet, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jet, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet30U, n_files_jet, "HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jet, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jet, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet15U_AND_HLT_Jet30U, n_files_jet, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet50U, n_files_jet, "HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jet, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jet, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_Jet15U_AND_HLT_Jet50U, n_files_jet, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_DiJetAve15U, n_files_jet, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_DiJetAve30U, n_files_jet, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_FwdJet20U, n_files_jet, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_DoubleJet15U_ForwardBackward, n_files_jet, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_all, n_files_jet, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "All triggered were sucessfully extracted!"<<endl; }
}


//compute the trigger efficiencies
if (compute_trigger_efficiency)
{
if (show_steps) { cout << "Compute Trigger Efficiencies..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, eff_dir, "root");
create_directories(output_dir, trigger_plots, "plots");

if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_L1Jet6U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_L1Jet10U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet15U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_Jet15U, "JetMETTauMonitor_2010A_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jetmettaumon_HLT_L1Jet6U, false, eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jetmettaumon_HLT_L1Jet10U, false, eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet30U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jetmettaumon_HLT_L1Jet6U, false, eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jetmettaumon_HLT_L1Jet10U, false, eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U, true, out_jetmettaumon_HLT_Jet15U, false, eff_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet50U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jetmettaumon_HLT_L1Jet6U, false, eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jetmettaumon_HLT_L1Jet10U, false, eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U, true, out_jetmettaumon_HLT_Jet15U, false, eff_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_DiJetAve15U, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_DiJetAve15U, "JetMETTauMonitor_2010A_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_DiJetAve30U, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_DiJetAve30U, "JetMETTauMonitor_2010A_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_FwdJet20U, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_FwdJet20U, "JetMETTauMonitor_2010A_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_DoubleJet15U_ForwardBackward, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_DoubleJet15U_ForwardBackward, "JetMETTauMonitor_2010A_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_all, false, out_jetmettaumon_none, false, eff_jetmettaumon_all, "JetMETTauMonitor_2010A_all_", trigger_plots, detail);

if (show_steps) { cout << "JetMETTau"<<endl; }
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U, false, out_jetmettau_none, false, eff_jetmettau_HLT_L1Jet6U, "JetMETTau_2010A_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U, false, out_jetmettau_none, false, eff_jetmettau_HLT_L1Jet10U, "JetMETTau_2010A_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet15U, false, out_jetmettau_none, false, eff_jetmettau_HLT_Jet15U, "JetMETTau_2010A_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jetmettau_HLT_L1Jet6U, false, eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U, "JetMETTau_2010A_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jetmettau_HLT_L1Jet10U, false, eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U, "JetMETTau_2010A_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet30U, true, out_jetmettau_none, false, eff_jetmettau_HLT_Jet30U, "JetMETTau_2010A_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jetmettau_HLT_L1Jet6U, false, eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U, "JetMETTau_2010A_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jetmettau_HLT_L1Jet10U, false, eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U, "JetMETTau_2010A_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet15U_AND_HLT_Jet30U, true, out_jetmettau_HLT_Jet15U, false, eff_jetmettau_HLT_Jet15U_AND_HLT_Jet30U, "JetMETTau_2010A_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet50U, true, out_jetmettau_none, false, eff_jetmettau_HLT_Jet50U, "JetMETTau_2010A_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jetmettau_HLT_L1Jet6U, false, eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U, "JetMETTau_2010A_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jetmettau_HLT_L1Jet10U, false, eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U, "JetMETTau_2010A_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet15U_AND_HLT_Jet50U, true, out_jetmettau_HLT_Jet15U, false, eff_jetmettau_HLT_Jet15U_AND_HLT_Jet50U, "JetMETTau_2010A_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_DiJetAve15U, false, out_jetmettau_none, false, eff_jetmettau_HLT_DiJetAve15U, "JetMETTau_2010A_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_DiJetAve30U, false, out_jetmettau_none, false, eff_jetmettau_HLT_DiJetAve30U, "JetMETTau_2010A_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_FwdJet20U, false, out_jetmettau_none, false, eff_jetmettau_HLT_FwdJet20U, "JetMETTau_2010A_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_DoubleJet15U_ForwardBackward, false, out_jetmettau_none, false, eff_jetmettau_HLT_DoubleJet15U_ForwardBackward, "JetMETTau_2010A_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_all, false, out_jetmettau_none, false, eff_jetmettau_all, "JetMETTau_2010A_all_", trigger_plots, detail);

if (show_steps) { cout << "JetMET"<<endl; }
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U, false, out_jetmet_none, false, eff_jetmet_HLT_L1Jet6U, "JetMET_2010A_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U, false, out_jetmet_none, false, eff_jetmet_HLT_L1Jet10U, "JetMET_2010A_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet15U, false, out_jetmet_none, false, eff_jetmet_HLT_Jet15U, "JetMET_2010A_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jetmet_HLT_L1Jet6U, false, eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U, "JetMET_2010A_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jetmet_HLT_L1Jet10U, false, eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U, "JetMET_2010A_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet30U, true, out_jetmet_none, false, eff_jetmet_HLT_Jet30U, "JetMET_2010A_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jetmet_HLT_L1Jet6U, false, eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U, "JetMET_2010A_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jetmet_HLT_L1Jet10U, false, eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U, "JetMET_2010A_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet15U_AND_HLT_Jet30U, true, out_jetmet_HLT_Jet15U, false, eff_jetmet_HLT_Jet15U_AND_HLT_Jet30U, "JetMET_2010A_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet50U, true, out_jetmet_none, false, eff_jetmet_HLT_Jet50U, "JetMET_2010A_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jetmet_HLT_L1Jet6U, false, eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U, "JetMET_2010A_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jetmet_HLT_L1Jet10U, false, eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U, "JetMET_2010A_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet15U_AND_HLT_Jet50U, true, out_jetmet_HLT_Jet15U, false, eff_jetmet_HLT_Jet15U_AND_HLT_Jet50U, "JetMET_2010A_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_DiJetAve15U, false, out_jetmet_none, false, eff_jetmet_HLT_DiJetAve15U, "JetMET_2010A_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_DiJetAve30U, false, out_jetmet_none, false, eff_jetmet_HLT_DiJetAve30U, "JetMET_2010A_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_FwdJet20U, false, out_jetmet_none, false, eff_jetmet_HLT_FwdJet20U, "JetMET_2010A_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_DoubleJet15U_ForwardBackward, false, out_jetmet_none, false, eff_jetmet_HLT_DoubleJet15U_ForwardBackward, "JetMET_2010A_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_all, false, out_jetmet_none, false, eff_jetmet_all, "JetMET_2010A_all_", trigger_plots, detail);

if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_efficiency(out_jet_HLT_L1Jet6U, false, out_jet_none, false, eff_jet_HLT_L1Jet6U, "Jet_2010B_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U, false, out_jet_none, false, eff_jet_HLT_L1Jet10U, "Jet_2010B_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet15U, false, out_jet_none, false, eff_jet_HLT_Jet15U, "Jet_2010B_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jet_HLT_L1Jet6U, false, eff_jet_HLT_L1Jet6U_AND_HLT_Jet15U, "Jet_2010B_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jet_HLT_L1Jet10U, false, eff_jet_HLT_L1Jet10U_AND_HLT_Jet15U, "Jet_2010B_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet30U, true, out_jet_none, false, eff_jet_HLT_Jet30U, "Jet_2010B_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jet_HLT_L1Jet6U, false, eff_jet_HLT_L1Jet6U_AND_HLT_Jet30U, "Jet_2010B_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jet_HLT_L1Jet10U, false, eff_jet_HLT_L1Jet10U_AND_HLT_Jet30U, "Jet_2010B_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet15U_AND_HLT_Jet30U, true, out_jet_HLT_Jet15U, false, eff_jet_HLT_Jet15U_AND_HLT_Jet30U, "Jet_2010B_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet50U, true, out_jet_none, false, eff_jet_HLT_Jet50U, "Jet_2010B_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jet_HLT_L1Jet6U, false, eff_jet_HLT_L1Jet6U_AND_HLT_Jet50U, "Jet_2010B_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jet_HLT_L1Jet10U, false, eff_jet_HLT_L1Jet10U_AND_HLT_Jet50U, "Jet_2010B_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet15U_AND_HLT_Jet50U, true, out_jet_HLT_Jet15U, false, eff_jet_HLT_Jet15U_AND_HLT_Jet50U, "Jet_2010B_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_DiJetAve15U, false, out_jet_none, false, eff_jet_HLT_DiJetAve15U, "Jet_2010B_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_DiJetAve30U, false, out_jet_none, false, eff_jet_HLT_DiJetAve30U, "Jet_2010B_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_FwdJet20U, false, out_jet_none, false, eff_jet_HLT_FwdJet20U, "Jet_2010B_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_DoubleJet15U_ForwardBackward, false, out_jet_none, false, eff_jet_HLT_DoubleJet15U_ForwardBackward, "Jet_2010B_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_all, false, out_jet_none, false, eff_jet_HLT_DiJetAve30U, "Jet_2010B_all_", trigger_plots, detail);
if (show_steps) { cout << "All trigger efficiencies were sucessfully computed!"<<endl; }
}


//get the trigger turn on
if (compute_trigger_turn_on)
{
if (show_steps) { cout << "Get Trigger turn on..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//compute_trigger_turn_on(jetmettaumon, out_jetmettaumon_emu, n_files_jetmettaumon, "", "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//compute_trigger_turn_on(jetmettau, out_jetmettau_emu, 1, "", "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//compute_trigger_turn_on(jetmet, out_jetmet_emu, 1, "", "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_turn_on(jet, out_jet_emu, 1, "", "allvertex", detail, test);
}


//get the trigger correction v1
if (compute_trigger_correction_v1)
{
if (show_steps) { cout << "Compute Trigger Correction v1..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_corr_dir, "root");
create_directories(output_dir, trigger_corr_plots, "plots");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }

prefix = "JetMETTauMonitor_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jetmettaumon_emu, correction_jetmettaumon_emu_v1, trigger_corr_plots, prefix, "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jetmettaumon_emu, correction_sj_jetmettaumon_emu_v1, trigger_corr_plots, prefix, "_v1", true, detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }

prefix = "JetMETTau_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jetmettau_emu, correction_jetmettau_emu_v1, trigger_corr_plots, prefix, "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jetmettau_emu, correction_sj_jetmettau_emu_v1, trigger_corr_plots, prefix, "_v1", true, detail, test);

if (show_steps) { cout << "JetMET"<<endl; }

prefix = "JetMET_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jetmet_emu, correction_jetmet_emu_v1, trigger_corr_plots, prefix, "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jetmet_emu, correction_sj_jetmet_emu_v1, trigger_corr_plots, prefix, "_v1", true, detail, test);

if (show_steps) { cout << "Jet"<<endl; }

prefix = "Jet_2010B";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jet_emu, correction_jet_emu_v1, trigger_corr_plots, "Jet_2010B", "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jet_emu, correction_sj_jet_emu_v1, trigger_corr_plots, "Jet_2010B", "_v1", true, detail, test);

}

//check the trigger correction v1
if (check_trigger_correction_v1)
{
if (show_steps) { cout << "Check Trigger Correction v1..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
compute_trigger_turn_on(jetmettaumon, check_jetmettaumon_emu_v1, n_files_jetmettaumon, correction_jetmettaumon_emu_v1, "allvertex", detail, test);
compute_trigger_turn_on(jetmettaumon, check_sj_jetmettaumon_emu_v1, n_files_jetmettaumon, correction_sj_jetmettaumon_emu_v1, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//compute_trigger_turn_on(jetmettau, check_jetmettau_emu_v1, n_files_jetmettau, correction_jetmettau_emu_v1, "allvertex", detail, test);
//compute_trigger_turn_on(jetmettau, check_sj_jetmettau_emu_v1, n_files_jetmettau, correction_sj_jetmettau_emu_v1, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//compute_trigger_turn_on(jetmet, check_jetmet_emu_v1, n_files_jetmet, correction_jetmet_emu_v1, "allvertex", detail, test);
//compute_trigger_turn_on(jetmet, check_sj_jetmet_emu_v1, n_files_jetmet, correction_sj_jetmet_emu_v1, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
//compute_trigger_turn_on(jet, check_jet_emu_v1, n_files_jet, correction_jet_emu_v1, "allvertex", detail, test);
//compute_trigger_turn_on(jet, check_sj_jet_emu_v1, n_files_jet, correction_sj_jet_emu_v1, "allvertex", detail, test);
}


//get the trigger correction v2
if (compute_trigger_correction_v2)
{
if (show_steps) { cout << "Compute Trigger Correction v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_corr_dir, "root");
create_directories(output_dir, trigger_corr_plots, "plots");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }

prefix = "JetMETTauMonitor_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(check_jetmettaumon_emu_v1, correction_jetmettaumon_emu_v2, trigger_corr_plots, prefix, "_v2", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(check_sj_jetmettaumon_emu_v1, correction_sj_jetmettaumon_emu_v2, trigger_corr_plots, prefix, "_v2", true, detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
compute_trigger_correction(check_jetmettau_emu_v1, correction_jetmettau_emu_v2, trigger_corr_plots, "JetMETTau_2010A", "_v2", false, detail, test);
compute_trigger_correction(check_jetmettau_emu_v1, correction_sj_jetmettau_emu_v2, trigger_corr_plots, "JetMETTau_2010A", "_v2", true, detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
compute_trigger_correction(check_jetmet_emu_v1, correction_jetmet_emu_v2, trigger_corr_plots, "JetMET_2010A", "_v2", false, detail, test);
compute_trigger_correction(check_jetmet_emu_v1, correction_sj_jetmet_emu_v2, trigger_corr_plots, "JetMET_2010A", "_v2", true, detail, test);
if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_correction(check_jet_emu_v1, correction_jet_emu_v2, trigger_corr_plots, "Jet_2010B", "_v2", false, detail, test);
compute_trigger_correction(check_jet_emu_v1, correction_sj_jet_emu_v2, trigger_corr_plots, "Jet_2010B", "_v2", true, detail, test);
}


//check the trigger correction v2
if (check_trigger_correction_v2)
{
if (show_steps) { cout << "Check Trigger Correction v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
///compute_trigger_turn_on(jetmettaumon, check_jetmettaumon_emu_v2, n_files_jetmettaumon, correction_jetmettaumon_emu_v2, "allvertex", detail, test);
///compute_trigger_turn_on(jetmettaumon, check_sj_jetmettaumon_emu_v2, n_files_jetmettaumon, correction_sj_jetmettaumon_emu_v2, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//compute_trigger_turn_on(jetmettau, check_jetmettau_emu_v2, n_files_jetmettau, correction_jetmettau_emu_v2, "allvertex", detail, test);
//compute_trigger_turn_on(jetmettau, check_sj_jetmettau_emu_v2, n_files_jetmettau, correction_sj_jetmettau_emu_v2, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//compute_trigger_turn_on(jetmet, check_jetmet_emu_v2, n_files_jetmet, correction_jetmet_emu_v2, "allvertex", detail, test);
//compute_trigger_turn_on(jetmet, check_sj_jetmet_emu_v2, n_files_jetmet, correction_sj_jetmet_emu_v2, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
//compute_trigger_turn_on(jet, check_jet_emu_v2, n_files_jet, correction_jet_emu_v2, "allvertex", detail, test);
//compute_trigger_turn_on(jet, check_sj_jet_emu_v2, n_files_jet, correction_sj_jet_emu_v2, "allvertex", detail, test);
}


//apply the trigger correction
if (apply_trigger_correction)
{
if (show_steps) { cout << "Apply Trigger Correction..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
///compute_trigger_turn_on(jetmettaumon, check_jetmettaumon_emu_v1, n_files_jetmettaumon, out_jetmettaumon_emu, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
///compute_trigger_turn_on(jetmettau, check_jetmettau_emu_v1, n_files_jetmettau, out_jetmettau_emu, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
///compute_trigger_turn_on(jetmet, check_jetmet_emu_v1, n_files_jetmet, out_jetmet_emu, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_turn_on(jet, check_jet_emu_v1, n_files_jet, out_jet_emu, "allvertex", detail, test);
}


// plot final trigger efficiency
if (plot_final_trigger_efficiency)
{
if (show_steps) { cout << "Plot Final Trigger Efficiencies..."<<endl; }
create_directories(output_dir, trigger_plots, "plots");

plot_final_trigger_efficiency(out_jetmettau_emu, out_jetmet_emu, out_jet_emu, trigger_plots, detail, test);

plot_final_trigger_efficiency2(out_jetmettau_emu, out_jetmet_emu, out_jet_emu, trigger_plots, detail, test);

if (show_steps) { cout << "Done!"<<endl; }
}


//counting combination statistics
if (counting_combination_statistics)
{
if (show_steps) { cout << "Counting Combination Statistics..."<<endl; }
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//counting_combination_statistics(jetmettaumon, "", n_files_jetmettaumon, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//counting_combination_statistics(jetmettau, "", n_files_jetmettau, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//counting_combination_statistics(jetmet, "", n_files_jetmet, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
counting_combination_statistics(jet, "", n_files_jet, "allvertex", detail, test);
}

//get the vertex distribution
if (get_vertex_distribution_v0)
{
if (show_steps) { cout << "Get the vertex distribution..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", "", "", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex, lumi_p6_z2, 1, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Herwig pp EEC3"<<endl; }
get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex, lumi_hw_eec3, 1, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Herwig pp 23"<<endl; }
//get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex, lumi_mc_hw_23, 1, "MC_DET", "allvertex", "", "", detail, test); //done
if (show_steps) { cout << "DATA"<<endl; }
if (show_steps) { cout << "JetMETTAU 2010A"<<endl; }
//get_vertex_distribution(jetmettau, vertex_jetmettau_allvertex, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
if (show_steps) { cout << "JetMET 2010A"<<endl; }
//get_vertex_distribution(jetmet, vertex_jetmet_allvertex, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Jet 2010B"<<endl; }
//get_vertex_distribution(jet, vertex_jet_allvertex, lumi_jet, n_files_jet, "DATA", "allvertex", "", "", detail, test);



if (show_steps) { cout << "The vertex distribution was sucessfully computed!"<<endl; }
}


// compute the vertex weights v1
if (compute_vertex_weights_v1)
{
if (show_steps) { cout << "Compute the vertex weights..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
///compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmettau_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jetmettau_v1_", false, detail);
///compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmet_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jetmet_v1_", false, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jet_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jet_v1_", false, detail);
///plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v1, vertex_weights_p6_z2_jetmet_v1, "_v1", vertex_weights_plots, "p6_z2_v1_", detail);

//compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmettau_notnorm_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jetmettau_notnerm_v1_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmet_notnorm_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jetmet_notnorm_v1_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jet_notnorm_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jet_notnorm_v1_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_notnorm_v1, vertex_weights_p6_z2_jetmet_notnorm_v1, vertex_weights_p6_z2_jet_notnorm_v1, "_v1", vertex_weights_plots, "p6_z2_notnorm_v1_", detail);

if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
///compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmettau_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jetmettau_v1_", false, detail);
///compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmet_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jetmet_v1_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jet_v1_", false, detail);
///plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v1, vertex_weights_p8_4c_jetmet_v1, "_v1", vertex_weights_plots, "p8_4c_v1_", detail);

//compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmettau_notnorm_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jetmettau_notnorm_v1_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmet_notnorm_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jetmet_notnorm_v1_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_notnorm_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jet_notnorm_v1_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_notnorm_v1, vertex_weights_p8_4c_jetmet_notnorm_v1, vertex_weights_p8_4c_jet_notnorm_v1, "_v1", vertex_weights_plots, "p8_4c_notnorm_v1_", detail);


if (show_steps) { cout << "Herwig pp TuneEEC3"<<endl; }
compute_vertex_weights(vertex_hw_eec3_allvertex, "Herwig pp TuneEEC3", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_hw_eec3_allvertex, vertex_weights_hw_eec3_jetmettau_v1, "_v1", "", "", vertex_weights_plots, "hw_eec3_jetmettau_v1_", true, detail);
compute_vertex_weights(vertex_hw_eec3_allvertex, "Herwig pp TuneEEC3", vertex_jetmet_allvertex, "JetMET 2010A", vertex_hw_eec3_allvertex, vertex_weights_hw_eec3_jetmet_v1, "_v1", "", "", vertex_weights_plots, "hw_eec3_jetmet_v1_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jet_v1_", false, detail);
plot_status_vertex_reweight(vertex_weights_hw_eec3_jetmettau_v1, vertex_weights_hw_eec3_jetmet_v1, "_v1", vertex_weights_plots, "hw_eec3_v1_", detail);


if (show_steps) { cout << "Herwig pp Tune23"<<endl; }
///compute_vertex_weights(vertex_hw_23_allvertex, "Herwig pp Tune23", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_hw_23_allvertex, vertex_weights_hw_23_jetmettau_v1, "_v1", "", "", vertex_weights_plots, "hw_23_jetmettau_v1_", false, detail);
///compute_vertex_weights(vertex_hw_23_allvertex, "Herwig pp Tune23", vertex_jetmet_allvertex, "JetMET 2010A", vertex_hw_23_allvertex, vertex_weights_hw_23_jetmet_v1, "_v1", "", "", vertex_weights_plots, "hw_23_jetmet_v1_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jet_v1_", false, detail);
///plot_status_vertex_reweight(vertex_weights_hw_23_jetmettau_v1, vertex_weights_hw_23_jetmet_v1, "_v1", vertex_weights_plots, "hw_23_v1_", detail);

if (show_steps) { cout << "All the vertex weights were sucessfully read!"<<endl; }
}


//get the vertex distribution v1
if (get_vertex_distribution_v1)
{
if (show_steps) { cout << "Get the vertex distribution v1..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v1, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_notnorm_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_notnorm_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v1, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_notnorm_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_notnorm_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_notnorm_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_notnorm_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v1, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_notnorm_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_notnorm_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v1, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_notnorm_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_notnorm_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v1, "_v1", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_notnorm_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_notnorm_v1, "_v1", detail, test);

if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmettau_v1, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmettau_v1, "_v1", detail, test);
if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMET_2010A"<<endl; }
get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmet_v1, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmet_v1, "_v1", detail, test);


if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmettau_v1, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmettau_v1, "_v1", detail, test);
if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMET_2010A"<<endl; }
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmet_v1, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmet_v1, "_v1", detail, test);

if (show_steps) { cout << "The vertex distribution v1 was sucessfully computed!"<<endl; }
}

// compute the vertex weights v2
if (compute_vertex_weights_v2)
{
if (show_steps) { cout << "Compute the vertex weights v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
///compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v1, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmettau_v2, "_v2", vertex_weights_p6_z2_jetmettau_v1,  "_v1",vertex_weights_plots, "p6_z2_jetmettau_v2_", false, detail);
///compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v1, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmet_v2, "_v2", vertex_weights_p6_z2_jetmet_v1, "_v1", vertex_weights_plots, "p6_z2_jetmet_v2_", false, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jet_v1, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jet_v2, "_v2", vertex_weights_p6_z2_jet_v1, "_v1", vertex_weights_plots, "p6_z2_jet_v2_", true, detail);
///plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v2, vertex_weights_p6_z2_jetmet_v2, "_v2", vertex_weights_plots, "p6_z2_v2_", detail);

//compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_notnorm_v1, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmettau_notnorm_v2, "_v2", vertex_weights_p6_z2_jetmettau_notnorm_v1,  "_v1",vertex_weights_plots, "p6_z2_jetmettau_notnorm_v2_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_notnorm_v1, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmet_notnorm_v2, "_v2", vertex_weights_p6_z2_jetmet_notnorm_v1, "_v1", vertex_weights_plots, "p6_z2_jetmet_notnorm_v2_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jet_notnorm_v1, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jet_notnorm_v2, "_v2", vertex_weights_p6_z2_jet_notnorm_v1, "_v1", vertex_weights_plots, "p6_z2_jet_notnorm_v2_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_notnorm_v2, vertex_weights_p6_z2_jetmet_notnorm_v2, vertex_weights_p6_z2_jet_notnorm_v2, "_v2", vertex_weights_plots, "p6_z2_notnorm_v2_", detail);

if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
///compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v1, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmettau_v2, "_v2", vertex_weights_p8_4c_jetmettau_v1, "_v1", vertex_weights_plots, "p8_4c_jetmettau_v2_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v1, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmet_v2, "_v2", vertex_weights_p8_4c_jetmet_v1, "_v1", vertex_weights_plots, "p8_4c_jetmet_v2_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_v1, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_p8_4c_jet_v1, "_v1", vertex_weights_plots, "p8_4c_jet_v2_", false, detail);
///plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v2, vertex_weights_p8_4c_jetmet_v2, "_v2", vertex_weights_plots, "p8_4c_v2_", detail);

//compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_notnorm_v1, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmettau_notnorm_v2, "_v2", vertex_weights_p8_4c_jetmettau_notnorm_v1, "_v1", vertex_weights_plots, "p8_4c_jetmettau_notnorm_v2_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_notnorm_v1, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmet_notnorm_v2, "_v2", vertex_weights_p8_4c_jetmet_notnorm_v1, "_v1", vertex_weights_plots, "p8_4c_jetmet_notnorm_v2_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_notnorm_v1, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_notnorm_v2, "_v2", vertex_weights_p8_4c_jet_notnorm_v1, "_v1", vertex_weights_plots, "p8_4c_jet_notnorm_v2_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_notnorm_v2, vertex_weights_p8_4c_jetmet_notnorm_v2, vertex_weights_p8_4c_jet_notnorm_v2, "_v2", vertex_weights_plots, "p8_4c_notnorm_v2_", detail);


if (show_steps) { cout << "Herwig pp Tune EEC3 reweighted to data"<<endl; }
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmettau_v1, "Herwig pp - Tune EEC3", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_hw_eec3_allvertex, vertex_weights_hw_eec3_jetmettau_v2, "_v2", vertex_weights_hw_eec3_jetmettau_v1, "_v1", vertex_weights_plots, "hw_eec3_jetmettau_v2_", false, detail);
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmet_v1, "Herwig pp - Tune EEC3", vertex_jetmet_allvertex, "JetMET 2010A", vertex_hw_eec3_allvertex, vertex_weights_hw_eec3_jetmet_v2, "_v2", vertex_weights_hw_eec3_jetmet_v1, "_v1", vertex_weights_plots, "hw_eec3_jetmet_v2_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_v1, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_p8_4c_jet_v1, "_v1", vertex_weights_plots, "p8_4c_jet_v2_", false, detail);
plot_status_vertex_reweight(vertex_weights_hw_eec3_jetmettau_v2, vertex_weights_hw_eec3_jetmet_v2, "_v2", vertex_weights_plots, "hw_eec3_v2_", detail);

if (show_steps) { cout << "Herwig pp Tune 23 reweighted to data"<<endl; }
compute_vertex_weights(vertex_hw_23_allvertex_jetmettau_v1, "Herwig pp - Tune 23", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_hw_23_allvertex, vertex_weights_hw_23_jetmettau_v2, "_v2", vertex_weights_hw_23_jetmettau_v1, "_v1", vertex_weights_plots, "hw_23_jetmettau_v2_", false, detail);
compute_vertex_weights(vertex_hw_23_allvertex_jetmet_v1, "Herwig pp - Tune 23", vertex_jetmet_allvertex, "JetMET 2010A", vertex_hw_23_allvertex, vertex_weights_hw_23_jetmet_v2, "_v2", vertex_weights_hw_23_jetmet_v1, "_v1", vertex_weights_plots, "hw_23_jetmet_v2_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_v1, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_p8_4c_jet_v1, "_v1", vertex_weights_plots, "p8_4c_jet_v2_", false, detail);
///plot_status_vertex_reweight(vertex_weights_hw_23_jetmettau_v2, vertex_weights_hw_23_jetmet_v2, "_v2", vertex_weights_plots, "hw_23_v2_", detail);


if (show_steps) { cout << "All the vertex weights v2 were sucessfully read!"<<endl; }
}


//get the vertex distribution v2
if (get_vertex_distribution_v2)
{
if (show_steps) { cout << "Get the vertex distribution v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v2, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_notnorm_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_notnorm_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v2, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_notnorm_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_notnorm_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_notnorm_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_notnorm_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v2, "_v2", detail, test);
///get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v2, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_notnorm_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_notnorm_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v2, "_v2", detail, test);
////get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v2, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_notnorm_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_notnorm_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v2, "_v2", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_notnorm_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_notnorm_v2, "_v2", detail, test);

if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmettau_v2, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmettau_v2, "_v2", detail, test);
if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMET_2010A"<<endl; }
////get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmet_v2, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmet_v2, "_v2", detail, test);


if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmettau_v2, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmettau_v2, "_v2", detail, test);
if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMET_2010A"<<endl; }
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmet_v2, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmet_v2, "_v2", detail, test);

if (show_steps) { cout << "The vertex distribution v2 was sucessfully computed!"<<endl; }
}


// compute the vertex weights
if (compute_vertex_weights_v3)
{
if (show_steps) { cout << "Compute the vertex weights v3..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
///compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v2, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmettau_v3, "_v3", vertex_weights_p6_z2_jetmettau_v2,  "_v2",vertex_weights_plots, "p6_z2_jetmettau_v3_", false, detail);
////compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v2, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmet_v3, "_v3", vertex_weights_p6_z2_jetmet_v2, "_v2", vertex_weights_plots, "p6_z2_jetmet_v3_", false, detail);
///compute_vertex_weights(vertex_p6_z2_allvertex_jet_v2, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jet_v3, "_v3", vertex_weights_p6_z2_jet_v2, "_v2", vertex_weights_plots, "p6_z2_jet_v3_", false, detail);
////plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v3, vertex_weights_p6_z2_jetmet_v3, "_v3", vertex_weights_plots, "p6_z2_v3_", detail);

//compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_notnorm_v2, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmettau_notnorm_v3, "_v3", vertex_weights_p6_z2_jetmettau_notnorm_v2,  "_v2",vertex_weights_plots, "p6_z2_jetmettau_notnorm_v3_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_notnorm_v2, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jetmet_notnorm_v3, "_v3", vertex_weights_p6_z2_jetmet_notnorm_v2, "_v2", vertex_weights_plots, "p6_z2_jetmet_notnorm_v3_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jet_notnorm_v2, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_p6_z2_allvertex, vertex_weights_p6_z2_jet_notnorm_v3, "_v3", vertex_weights_p6_z2_jet_notnorm_v2, "_v2", vertex_weights_plots, "p6_z2_jet_notnorm_v3_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_notnorm_v3, vertex_weights_p6_z2_jetmet_notnorm_v3, vertex_weights_p6_z2_jet_notnorm_v3, "_v3", vertex_weights_plots, "p6_z2_notnorm_v3_", detail);

if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
////compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v2, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmettau_v3, "_v3", vertex_weights_p8_4c_jetmettau_v2, "_v2", vertex_weights_plots, "p8_4c_jetmettau_v3_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v2, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmet_v3, "_v3", vertex_weights_p8_4c_jetmet_v2, "_v2", vertex_weights_plots, "p8_4c_jetmet_v3_", false, detail);
///compute_vertex_weights(vertex_p8_4c_allvertex_jet_v2, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_plots, "p8_4c_jet_v3_", false, detail);
////plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v3, vertex_weights_p8_4c_jetmet_v3, "_v3", vertex_weights_plots, "p8_4c_v3_", detail);

//compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_notnorm_v2, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmettau_notnorm_v3, "_v3", vertex_weights_p8_4c_jetmettau_notnorm_v2, "_v2", vertex_weights_plots, "p8_4c_jetmettau_notnorm_v3_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_notnorm_v2, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jetmet_notnorm_v3, "_v3", vertex_weights_p8_4c_jetmet_notnorm_v2, "_v2", vertex_weights_plots, "p8_4c_jetmet_notnorm_v3_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_notnorm_v2, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_notnorm_v3, "_v3", vertex_weights_p8_4c_jet_notnorm_v2, "_v2", vertex_weights_plots, "p8_4c_jet_notnorm_v3_", false, detail);
//plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_notnorm_v3, vertex_weights_p8_4c_jetmet_notnorm_v3, vertex_weights_p8_4c_jet_notnorm_v3, "_v3", vertex_weights_plots, "p8_4c_notnorm_v3_", detail);


if (show_steps) { cout << "Herwig pp Tune 23 reweighted to data"<<endl; }
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmettau_v2, "Herwig pp Tune EEC3", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_hw_eec3_allvertex, vertex_weights_hw_eec3_jetmettau_v3, "_v3", vertex_weights_hw_eec3_jetmettau_v2, "_v2", vertex_weights_plots, "hw_eec3_jetmettau_v3_", false, detail);
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmet_v2, "Herwig pp Tune EEC3", vertex_jetmet_allvertex, "JetMET 2010A", vertex_hw_eec3_allvertex, vertex_weights_hw_eec3_jetmet_v3, "_v3", vertex_weights_hw_eec3_jetmet_v2, "_v2", vertex_weights_plots, "hw_eec3_jetmet_v3_", false, detail);
///compute_vertex_weights(vertex_p8_4c_allvertex_jet_v2, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_plots, "p8_4c_jet_v3_", false, detail);
plot_status_vertex_reweight(vertex_weights_hw_eec3_jetmettau_v3, vertex_weights_hw_eec3_jetmet_v3, "_v3", vertex_weights_plots, "hw_eec3_v3_", detail);


if (show_steps) { cout << "Herwig pp Tune 23 reweighted to data"<<endl; }
compute_vertex_weights(vertex_hw_23_allvertex_jetmettau_v2, "Herwig pp Tune 23", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_hw_23_allvertex, vertex_weights_hw_23_jetmettau_v3, "_v3", vertex_weights_hw_23_jetmettau_v2, "_v2", vertex_weights_plots, "hw_23_jetmettau_v3_", false, detail);
compute_vertex_weights(vertex_hw_23_allvertex_jetmet_v2, "Herwig pp Tune 23", vertex_jetmet_allvertex, "JetMET 2010A", vertex_hw_23_allvertex, vertex_weights_hw_23_jetmet_v3, "_v3", vertex_weights_hw_23_jetmet_v2, "_v2", vertex_weights_plots, "hw_23_jetmet_v3_", false, detail);
///compute_vertex_weights(vertex_p8_4c_allvertex_jet_v2, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_p8_4c_allvertex, vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_plots, "p8_4c_jet_v3_", false, detail);
plot_status_vertex_reweight(vertex_weights_hw_23_jetmettau_v3, vertex_weights_hw_23_jetmet_v3, "_v3", vertex_weights_plots, "hw_23_v3_", detail);


if (show_steps) { cout << "All the new vertex weights v3 were sucessfully read!"<<endl; }
}

//get the vertex distribution v3
if (get_vertex_distribution_v3)
{
if (show_steps) { cout << "Get the vertex distribution v3..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v3, "_v3", detail, test);
/////get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v3, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_notnorm_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_notnorm_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v3, "_v3", detail, test);
////get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v3, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_notnorm_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_notnorm_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_notnorm_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_notnorm_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v3, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_notnorm_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_notnorm_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v3, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_notnorm_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_notnorm_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v3, "_v3", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_notnorm_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_notnorm_v3, "_v3", detail, test);

if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
////get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmettau_v3, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmettau_v3, "_v3", detail, test);
if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMET_2010A"<<endl; }
get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmet_v3, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmet_v3, "_v3", detail, test);


if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmettau_v3, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmettau_v3, "_v3", detail, test);
if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMET_2010A"<<endl; }
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmet_v3, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmet_v3, "_v3", detail, test);
if (show_steps) { cout << "The vertex distribution v3 was sucessfully computed!"<<endl; }
}


// compute the vertex weights v4
if (compute_vertex_weights_v4)
{
if (show_steps) { cout << "Compute the new vertex weights v4..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
////compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v3, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_p6_z2_jetmettau_v4, "_v4", vertex_weights_p6_z2_jetmettau_v3,  "_v3",vertex_weights_plots, "p6_z2_jetmettau_v4_", false, detail);
////compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v3, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_p6_z2_jetmet_v4, "_v4", vertex_weights_p6_z2_jetmet_v3, "_v3", vertex_weights_plots, "p6_z2_jetmet_v4_", false, detail);
////compute_vertex_weights(vertex_p6_z2_allvertex_jet_v3, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p6_z2_jet_v4, "_v4", vertex_weights_p6_z2_jet_v3, "_v3", vertex_weights_plots, "p6_z2_jet_v4_", false, detail);
////plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v4, vertex_weights_p6_z2_jetmet_v4, "_v4", vertex_weights_plots, "p6_z2_v4_", detail);

//compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_notnorm_v3, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_notnorm_v4, "_v4", vertex_weights_p6_z2_jetmettau_notnorm_v3,  "_v3", vertex_weights_plots, "p6_z2_jetmettau_notnorm_v4_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_notnorm_v3, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_notnorm_v4, "_v4", vertex_weights_p6_z2_jetmet_notnorm_v3, "_v3", vertex_weights_plots, "p6_z2_jetmet_notnorm_v4_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jet_notnorm_v3, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p6_z2_jet_notnorm_v4, "_v4", vertex_weights_p6_z2_jet_notnorm_v3, "_v3", vertex_weights_plots, "p6_z2_jet_notnorm_v4_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_notnorm_v4, vertex_weights_p6_z2_jetmet_notnorm_v4, vertex_weights_p6_z2_jet_notnorm_v4, "_v4", vertex_weights_plots, "p6_z2_notnorm_v4_", detail);

if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
///compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v3, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_p8_4c_jetmettau_v4, "_v4", vertex_weights_p8_4c_jetmettau_v3, "_v3", vertex_weights_plots, "p8_4c_jetmettau_v4_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v3, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_p8_4c_jetmet_v4, "_v4", vertex_weights_p8_4c_jetmet_v3, "_v3", vertex_weights_plots, "p8_4c_jetmet_v4_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_v3, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_plots, "p8_4c_jet_v4_", false, detail);
////plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v4, vertex_weights_p8_4c_jetmet_v4, "_v4", vertex_weights_plots, "p8_4c_v4_", detail);

//compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_notnorm_v3, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_notnorm_v4, "_v4", vertex_weights_p8_4c_jetmettau_notnorm_v3, "_v3", vertex_weights_plots, "p8_4c_jetmettau_notnorm_v4_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_notnorm_v3, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_notnorm_v4, "_v4", vertex_weights_p8_4c_jetmet_notnorm_v3, "_v3", vertex_weights_plots, "p8_4c_jetmet_notnorm_v4_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_notnorm_v3, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_notnorm_v4, "_v4", vertex_weights_p8_4c_jet_notnorm_v3, "_v3", vertex_weights_plots, "p8_4c_jet_notnorm_v4_", false, detail);
//plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_notnorm_v4, vertex_weights_p8_4c_jetmet_notnorm_v4, vertex_weights_p8_4c_jet_notnorm_v4, "_v4", vertex_weights_plots, "p8_4c_notnorm_v4_", detail);

if (show_steps) { cout << "Herwig pp Tune EEC3 reweighted to data"<<endl; }
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmettau_v3, "Herwig pp Tune EEC3", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_hw_eec3_jetmettau_v4, "_v4", vertex_weights_hw_eec3_jetmettau_v3, "_v3", vertex_weights_plots, "hw_eec3_jetmettau_v4_", false, detail);
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmet_v3, "Herwig pp Tune EEC3", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_hw_eec3_jetmet_v4, "_v4", vertex_weights_hw_eec3_jetmet_v3, "_v3", vertex_weights_plots, "hw_eec3_jetmet_v4_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_v3, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_plots, "p8_4c_jet_v4_", false, detail);
plot_status_vertex_reweight(vertex_weights_hw_eec3_jetmettau_v4, vertex_weights_hw_eec3_jetmet_v4, "_v4", vertex_weights_plots, "hw_eec3_v4_", detail);



if (show_steps) { cout << "Herwig pp Tune 23 reweighted to data"<<endl; }
/////compute_vertex_weights(vertex_hw_23_allvertex_jetmettau_v3, "Herwig pp Tune 23", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_hw_23_jetmettau_v4, "_v4", vertex_weights_hw_23_jetmettau_v3, "_v3", vertex_weights_plots, "hw_23_jetmettau_v4_", false, detail);
/////compute_vertex_weights(vertex_hw_23_allvertex_jetmet_v3, "Herwig pp Tune 23", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_hw_23_jetmet_v4, "_v4", vertex_weights_hw_23_jetmet_v3, "_v3", vertex_weights_plots, "hw_23_jetmet_v4_", false, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_v3, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_plots, "p8_4c_jet_v4_", false, detail);
/////plot_status_vertex_reweight(vertex_weights_hw_23_jetmettau_v4, vertex_weights_hw_23_jetmet_v4, "_v4", vertex_weights_plots, "hw_23_v4_", detail);


if (show_steps) { cout << "All the new vertex weights v4 were sucessfully read!"<<endl; }
}


//get the vertex distribution v4
if (get_vertex_distribution_v4)
{
if (show_steps) { cout << "Get the vertex distribution v4..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v4, "_v4", detail, test);
/////get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v4, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v4, "_v4", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_notnorm_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_notnorm_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v4, "_v4", detail, test);
/////get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v4, lumi_p6_z2, 1, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v4, "_v4", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_notnorm_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_notnorm_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v4, "_v4", detail, test);
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_notnorm_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_notnorm_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v4, "_v4", detail, test);
///get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v4, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v4, "_v4", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_notnorm_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_notnorm_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v4, "_v4", detail, test);
////get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v4, lumi_p8_4c, 1, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v4, "_v4", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_notnorm_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_notnorm_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v4, "_v4", detail, test);
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_notnorm_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_notnorm_v4, "_v4", detail, test);

if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v4, "_v4", detail, test);
get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmettau_v4, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmettau_v4, "_v4", detail, test);
if (show_steps) { cout << "Herwig pp Tune EEC3 with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v4, "_v4", detail, test);
get_vertex_distribution(mc_hw_eec3, vertex_hw_eec3_allvertex_jetmet_v4, lumi_hw_eec3, 1, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmet_v4, "_v4", detail, test);



if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v4, "_v4", detail, test);
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmettau_v4, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmettau_v4, "_v4", detail, test);
if (show_steps) { cout << "Herwig pp Tune 23 with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v4, "_v4", detail, test);
///get_vertex_distribution(mc_hw_23, vertex_hw_23_allvertex_jetmet_v4, lumi_hw_23, 1, "MC_DET", "allvertex", vertex_weights_hw_23_jetmet_v4, "_v4", detail, test);


if (show_steps) { cout << "The vertex distribution v4 was sucessfully computed!"<<endl; }
}

// compute the vertex weights v5
if (compute_vertex_weights_v5)
{
if (show_steps) { cout << "Compute the new vertex weights v5..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
///compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v4, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_p6_z2_jetmettau_v5, "_v5", vertex_weights_p6_z2_jetmettau_v4,  "_v4",vertex_weights_plots, "p6_z2_jetmettau_v5_", false, detail);
///compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v4, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_p6_z2_jetmet_v5, "_v5", vertex_weights_p6_z2_jetmet_v4, "_v4", vertex_weights_plots, "p6_z2_jetmet_v5_", false, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jet_v4, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p6_z2_jet_v5, "_v5", vertex_weights_p6_z2_jet_v4, "_v4", vertex_weights_plots, "p6_z2_jet_v5_", false, detail);
////plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v5, vertex_weights_p6_z2_jetmet_v5, "_v5", vertex_weights_plots, "p6_z2_v5_", detail);

//compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_notnorm_v4, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_notnorm_v5, "_v5", vertex_weights_p6_z2_jetmettau_notnorm_v4,  "_v4",vertex_weights_plots, "p6_z2_jetmettau_notnorm_v5_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_notnorm_v4, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_notnorm_v5, "_v5", vertex_weights_p6_z2_jetmet_notnorm_v4, "_v4", vertex_weights_plots, "p6_z2_jetmet_notnorm_v5_", true, detail);
//compute_vertex_weights(vertex_p6_z2_allvertex_jet_notnorm_v4, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p6_z2_jet_notnorm_v5, "_v5", vertex_weights_p6_z2_jet_notnorm_v4, "_v4", vertex_weights_plots, "p6_z2_jet_notnorm_v5_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_notnorm_v5, vertex_weights_p6_z2_jetmet_notnorm_v5, vertex_weights_p6_z2_jet_notnorm_v5, "_v5", vertex_weights_plots, "p6_z2_notnorm_v5_", detail);

if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
///compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v4, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_p8_4c_jetmettau_v5, "_v5", vertex_weights_p8_4c_jetmettau_v4, "_v4", vertex_weights_plots, "p8_4c_jetmettau_v5_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v4, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_p8_4c_jetmet_v5, "_v5", vertex_weights_p8_4c_jetmet_v4, "_v4", vertex_weights_plots, "p8_4c_jetmet_v5_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jet_v4, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p8_4c_jet_v5, "_v5", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_plots, "p8_4c_jet_v5_", false, detail);
///plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v5, vertex_weights_p8_4c_jetmet_v5, "_v5", vertex_weights_plots, "p8_4c_v5_", detail);

//compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_notnorm_v4, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_notnorm_v5, "_v5", vertex_weights_p8_4c_jetmettau_notnorm_v4, "_v4", vertex_weights_plots, "p8_4c_jetmettau_notnorm_v5_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_notnorm_v4, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_notnorm_v5, "_v5", vertex_weights_p8_4c_jetmet_notnorm_v4, "_v4", vertex_weights_plots, "p8_4c_jetmet_notnorm_v5_", true, detail);
//compute_vertex_weights(vertex_p8_4c_allvertex_jet_notnorm_v4, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_notnorm_v5, "_v5", vertex_weights_p8_4c_jet_notnorm_v4, "_v4", vertex_weights_plots, "p8_4c_jet_notnorm_v5_", true, detail);
//plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_notnorm_v5, vertex_weights_p8_4c_jetmet_notnorm_v5, vertex_weights_p8_4c_jet_notnorm_v5, "_v5", vertex_weights_plots, "p8_4c_notnorm_v5_", detail);


if (show_steps) { cout << "Herwig pp Tune EEC3 reweighted to data"<<endl; }
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmettau_v4, "Herwig pp Tune EEC3", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_hw_eec3_jetmettau_v5, "_v5", vertex_weights_hw_eec3_jetmettau_v4, "_v4", vertex_weights_plots, "hw_eec3_jetmettau_v5_", false, detail);
compute_vertex_weights(vertex_hw_eec3_allvertex_jetmet_v4, "Herwig pp Tune EEC3", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_hw_eec3_jetmet_v5, "_v5", vertex_weights_hw_eec3_jetmet_v4, "_v4", vertex_weights_plots, "hw_eec3_jetmet_v5_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jet_v4, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p8_4c_jet_v5, "_v5", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_plots, "p8_4c_jet_v5_", false, detail);
plot_status_vertex_reweight(vertex_weights_hw_eec3_jetmettau_v5, vertex_weights_hw_eec3_jetmet_v5, "_v5", vertex_weights_plots, "hw_eec3_v5_", detail);

if (show_steps) { cout << "Herwig pp Tune 23 reweighted to data"<<endl; }
///compute_vertex_weights(vertex_hw_23_allvertex_jetmettau_v4, "Herwig pp Tune 23", vertex_jetmettau_allvertex, "JetMETTau 2010A", "", vertex_weights_hw_23_jetmettau_v5, "_v5", vertex_weights_hw_23_jetmettau_v4, "_v4", vertex_weights_plots, "hw_23_jetmettau_v5_", false, detail);
///compute_vertex_weights(vertex_hw_23_allvertex_jetmet_v4, "Herwig pp Tune 23", vertex_jetmet_allvertex, "JetMET 2010A", "", vertex_weights_hw_23_jetmet_v5, "_v5", vertex_weights_hw_23_jetmet_v4, "_v4", vertex_weights_plots, "hw_23_jetmet_v5_", false, detail);
////compute_vertex_weights(vertex_p8_4c_allvertex_jet_v4, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", "", vertex_weights_p8_4c_jet_v5, "_v5", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_plots, "p8_4c_jet_v5_", false, detail);
////plot_status_vertex_reweight(vertex_weights_hw_23_jetmettau_v5, vertex_weights_hw_23_jetmet_v5, "_v5", vertex_weights_plots, "hw_23_v5_", detail);


if (show_steps) { cout << "All the new vertex weights v5 were sucessfully read!"<<endl; }
}


//plots the status of the vertex reweighting process
if (plot_status_vertex_reweight)
{
if (show_steps) { cout << "Plot the Status of the Vertex Reweight..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighting evolution..."<<endl; }
////plot_evolution_vertex_reweight(vertex_weights_p6_z2_jetmettau_v1, vertex_weights_p6_z2_jetmettau_v2, vertex_weights_p6_z2_jetmettau_v3, vertex_weights_p6_z2_jetmettau_v4, vertex_weights_p6_z2_jetmettau_v5, vertex_weights_plots, "p6_z2_jetmettau_", detail);
///plot_evolution_vertex_reweight(vertex_weights_p6_z2_jetmet_v1, vertex_weights_p6_z2_jetmet_v2, vertex_weights_p6_z2_jetmet_v3, vertex_weights_p6_z2_jetmet_v4, vertex_weights_p6_z2_jetmet_v5, vertex_weights_plots, "p6_z2_jetmet_", detail);
///plot_evolution_vertex_reweight(vertex_weights_p6_z2_jet_v1, vertex_weights_p6_z2_jet_v2, vertex_weights_p6_z2_jet_v3, vertex_weights_p6_z2_jet_v4, vertex_weights_p6_z2_jet_v5, vertex_weights_plots, "p6_z2_jet_", detail);

//plot_evolution_vertex_reweight(vertex_weights_p6_z2_jetmettau_notnorm_v1, vertex_weights_p6_z2_jetmettau_notnorm_v2, vertex_weights_p6_z2_jetmettau_notnorm_v3, vertex_weights_p6_z2_jetmettau_notnorm_v4, vertex_weights_p6_z2_jetmettau_notnorm_v5, vertex_weights_plots, "p6_z2_jetmettau_notnorm_", detail);
//plot_evolution_vertex_reweight(vertex_weights_p6_z2_jetmet_notnorm_v1, vertex_weights_p6_z2_jetmet_notnorm_v2, vertex_weights_p6_z2_jetmet_notnorm_v3, vertex_weights_p6_z2_jetmet_notnorm_v4, vertex_weights_p6_z2_jetmet_notnorm_v5, vertex_weights_plots, "p6_z2_jetmet_notnorm_", detail);
//plot_evolution_vertex_reweight(vertex_weights_p6_z2_jet_notnorm_v1, vertex_weights_p6_z2_jet_notnorm_v2, vertex_weights_p6_z2_jet_notnorm_v3, vertex_weights_p6_z2_jet_notnorm_v4, vertex_weights_p6_z2_jet_notnorm_v5, vertex_weights_plots, "p6_z2_jet_notnorm_", detail);

if (show_steps) { cout << "Pythia8 4CTune reweighting evolution..."<<endl; }
////plot_evolution_vertex_reweight(vertex_weights_p8_4c_jetmettau_v1, vertex_weights_p8_4c_jetmettau_v2, vertex_weights_p8_4c_jetmettau_v3, vertex_weights_p8_4c_jetmettau_v4, vertex_weights_p8_4c_jetmettau_v5, vertex_weights_plots, "p8_4c_jetmettau_", detail);
///plot_evolution_vertex_reweight(vertex_weights_p8_4c_jetmet_v1, vertex_weights_p8_4c_jetmet_v2, vertex_weights_p8_4c_jetmet_v3, vertex_weights_p8_4c_jetmet_v4, vertex_weights_p8_4c_jetmet_v5, vertex_weights_plots, "p8_4c_jetmet_",detail);
///plot_evolution_vertex_reweight(vertex_weights_p8_4c_jet_v1, vertex_weights_p8_4c_jet_v2, vertex_weights_p8_4c_jet_v3, vertex_weights_p8_4c_jet_v4, vertex_weights_p8_4c_jet_v5, vertex_weights_plots, "p8_4c_jet_", detail);

//plot_evolution_vertex_reweight(vertex_weights_p8_4c_jetmettau_notnorm_v1, vertex_weights_p8_4c_jetmettau_notnorm_v2, vertex_weights_p8_4c_jetmettau_notnorm_v3, vertex_weights_p8_4c_jetmettau_notnorm_v4, vertex_weights_p8_4c_jetmettau_notnorm_v5, vertex_weights_plots, "p8_4c_jetmettau_notnorm_", detail);
//plot_evolution_vertex_reweight(vertex_weights_p8_4c_jetmet_notnorm_v1, vertex_weights_p8_4c_jetmet_notnorm_v2, vertex_weights_p8_4c_jetmet_notnorm_v3, vertex_weights_p8_4c_jetmet_notnorm_v4, vertex_weights_p8_4c_jetmet_notnorm_v5, vertex_weights_plots, "p8_4c_jetmet_notnorm_",detail);
//plot_evolution_vertex_reweight(vertex_weights_p8_4c_jet_notnorm_v1, vertex_weights_p8_4c_jet_notnorm_v2, vertex_weights_p8_4c_jet_notnorm_v3, vertex_weights_p8_4c_jet_notnorm_v4, vertex_weights_p8_4c_jet_notnorm_v5, vertex_weights_plots, "p8_4c_jet_notnorm_", detail);

if (show_steps) { cout << "Herwig pp Tune EEC3 reweighting evolution..."<<endl; }
plot_evolution_vertex_reweight(vertex_weights_hw_eec3_jetmettau_v1, vertex_weights_hw_eec3_jetmettau_v2, vertex_weights_hw_eec3_jetmettau_v3, vertex_weights_hw_eec3_jetmettau_v4, vertex_weights_hw_eec3_jetmettau_v5, vertex_weights_plots, "hw_eec3_jetmettau_", detail);
plot_evolution_vertex_reweight(vertex_weights_hw_eec3_jetmet_v1, vertex_weights_hw_eec3_jetmet_v2, vertex_weights_hw_eec3_jetmet_v3, vertex_weights_hw_eec3_jetmet_v4, vertex_weights_hw_eec3_jetmet_v5, vertex_weights_plots, "hw_eec3_jetmet_",detail);


if (show_steps) { cout << "Herwig pp Tune 23 reweighting evolution..."<<endl; }
////plot_evolution_vertex_reweight(vertex_weights_hw_23_jetmettau_v1, vertex_weights_hw_23_jetmettau_v2, vertex_weights_hw_23_jetmettau_v3, vertex_weights_hw_23_jetmettau_v4, vertex_weights_hw_23_jetmettau_v5, vertex_weights_plots, "hw_23_jetmettau_", detail);
////plot_evolution_vertex_reweight(vertex_weights_hw_23_jetmet_v1, vertex_weights_hw_23_jetmet_v2, vertex_weights_hw_23_jetmet_v3, vertex_weights_hw_23_jetmet_v4, vertex_weights_hw_23_jetmet_v5, vertex_weights_plots, "hw_23_jetmet_",detail);


if (show_steps) { cout << "The status was sucessfully plotted!"<<endl; }
}


if (compute_resolution)
{
if (show_steps) { cout << "Compute Resolution..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, resolution_dir, "root");
if (show_steps) { cout << "Compute Resolution for Pythia 6 - Tune Z2Star..."<<endl; }
//compute_resolution(mc_p6_z2, resolution_p6_z2_jetmettau, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//compute_resolution(mc_p6_z2, resolution_p6_z2_jetmet, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//compute_resolution(mc_p6_z2, resolution_p6_z2_jet, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "Compute Resolution for Pythia 8 - Tune 4C..."<<endl; }
///compute_resolution(mc_p8_4c, resolution_p8_4c_jetmettau, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///compute_resolution(mc_p8_4c, resolution_p8_4c_jetmet, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
//compute_resolution(mc_p8_4c, resolution_p8_4c_jet, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);

if (show_steps) { cout << "Compute Resolution for Herwig pp - Tune EEC3..."<<endl; }
compute_resolution(mc_hw_eec3, resolution_hw_eec3_jetmettau, lumi_hw_eec3, n_files_hw_eec3, "allvertex", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);
//compute_resolution(mc_hw_eec3, resolution_hw_eec3_jetmet, lumi_hw_eec3, n_files_hw_eec3, "allvertex", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);

if (show_steps) { cout << "Compute Resolution for Herwig pp - Tune 23..."<<endl; }
//compute_resolution(mc_hw_23, resolution_hw_23_jetmettau, lumi_hw_23, n_files_hw_23, "allvertex", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);
///compute_resolution(mc_hw_23, resolution_hw_23_jetmet, lumi_hw_23, n_files_hw_23, "allvertex", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);


if (show_steps) { cout << "Resolution Sucessfully Computed!"<<endl; }
}


if (plot_resolution)
{
if (show_steps) { cout << "Plotting Resolution..."<<endl; }
create_directories(output_dir, resolution_plots, "plots");

//plot_resolution(resolution_p6_z2_jetmettau, resolution_p6_z2_jetmet, resolution_p8_4c_jetmettau, resolution_p8_4c_jetmet, resolution_plots, "resolution_", true, detail);
plot_resolution(resolution_hw_eec3_jetmettau, resolution_hw_eec3_jetmet, resolution_p8_4c_jetmettau, resolution_p8_4c_jetmet, resolution_plots, "resolution_", true, detail);

if (show_steps) { cout << "Resolution Sucessfully Plotted!"<<endl; }
}


if (get_pileup_normalization)
{
if (show_steps) { cout << "Getting the Pileup Normalization..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, pileup_norm_dir, "root");
if (show_steps) { cout << "For Pythia 6 - Tune Z2Star..."<<endl; }
//get_pileup_normalization(mc_p6_z2, norm_p6_z2, lumi_p6_z2, n_files_p6_z2, "", "", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_jetmettau, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_jetmet, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_jet, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jet_v5, "_v5", detail, test);

//get_pileup_normalization(mc_p6_z2, norm_p6_z2_notnorm_jetmettau, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jetmettau_notnorm_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_notnorm_jetmet, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jetmet_notnorm_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_notnorm_jet, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jet_notnorm_v5, "_v5", detail, test);

if (show_steps) { cout << "For Pythia 8 - Tune 4C..."<<endl; }
///get_pileup_normalization(mc_p8_4c, norm_p8_4c, lumi_p8_4c, n_files_p8_4c, "", "", detail, test);
///get_pileup_normalization(mc_p8_4c, norm_p8_4c_jetmettau, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///get_pileup_normalization(mc_p8_4c, norm_p8_4c_jetmet, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
///get_pileup_normalization(mc_p8_4c, norm_p8_4c_jet, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jet_v5, "_v5", detail, test);

if (show_steps) { cout << "For Herwig pp - EEC3..."<<endl; }
///get_pileup_normalization(mc_hw_eec3, norm_hw_eec3, lumi_hw_eec3, n_files_hw_eec3, "", "", detail, test);
///get_pileup_normalization(mc_hw_eec3, norm_hw_eec3_jetmettau, lumi_hw_eec3, n_files_hw_eec3, vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);
///get_pileup_normalization(mc_hw_eec3, norm_hw_eec3_jetmet, lumi_hw_eec3, n_files_hw_eec3, vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);

if (show_steps) { cout << "For Herwig pp - 23..."<<endl; }
///get_pileup_normalization(mc_hw_23, norm_hw_23, lumi_hw_23, n_files_hw_23, "", "", detail, test);
///get_pileup_normalization(mc_hw_23, norm_hw_23_jetmettau, lumi_hw_23, n_files_hw_23, vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);
///get_pileup_normalization(mc_hw_23, norm_hw_23_jetmet, lumi_hw_23, n_files_hw_23, vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);

//get_pileup_normalization(mc_p8_4c, norm_p8_4c_notnorm_jetmettau, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jetmettau_notnorm_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p8_4c, norm_p8_4c_notnorm_jetmet, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jetmet_notnorm_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p8_4c, norm_p8_4c_notnorm_jet, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jet_notnorm_v5, "_v5", detail, test);
if (show_steps) { cout << "The Pileup Normalization was been sucessfully obtained!"<<endl; }
if (show_steps) { cout << "Showing Pileup Normalization ..."<<endl; }
show_pileup_normalization(norm_hw_eec3, norm_hw_eec3_jetmettau, norm_hw_eec3_jetmet, norm_p8_4c, norm_p8_4c_jetmettau, norm_p8_4c_jetmet, norm_hw_23, norm_hw_23_jetmettau, norm_hw_23_jetmet,detail);
///show_pileup_normalization(norm_p6_z2, norm_p6_z2_notnorm_jetmettau, norm_p6_z2_notnorm_jetmet, norm_p6_z2_notnorm_jet, norm_p8_4c, norm_p8_4c_notnorm_jetmettau, norm_p8_4c_notnorm_jetmet, norm_p8_4c_notnorm_jet, detail);
if (show_steps) { cout << "The Pileup Normalization was been sucessfully outputed!"<<endl; }
}

//read the MC Ntuples on the Generator Level
if (read_mc_ntuples_gen)
{
if (show_steps) { cout << "Reading MC Ntuples on Generator Level..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_gen_dir, "root");
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", "", "", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_notnorm_gen_nopileup_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jetmettau_notnorm_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_notnorm_gen_allvertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jetmettau_notnorm_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_notnorm_gen_nopileup_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jetmet_notnorm_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_notnorm_gen_allvertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jetmet_notnorm_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup_jet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_notnorm_gen_nopileup_jet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jet_notnorm_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_notnorm_gen_allvertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jet_notnorm_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
///read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", "", "", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_notnorm_gen_nopileup_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jetmettau_notnorm_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_notnorm_gen_allvertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jetmettau_notnorm_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
///read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_notnorm_gen_nopileup_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jetmet_notnorm_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_notnorm_gen_allvertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jetmet_notnorm_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
////read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup_jet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_notnorm_gen_nopileup_jet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jet_notnorm_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_notnorm_gen_allvertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jet_notnorm_v5, "_v5", detail, test);



if (show_steps) { cout << "Herwig pp EEC3 Tune"<<endl; }
///read_ntuple(mc_hw_eec3, out_hw_eec3_gen_nopileup, lumi_hw_eec3, n_files_hw_eec3, "MC_GEN", "nopileup", "", "", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_gen_allvertex, lumi_hw_eec3, n_files_hw_eec3, "MC_GEN", "allvertex", "", "", detail, test);

if (show_steps) { cout << "Herwig pp EEC3 Tune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///read_ntuple(mc_hw_eec3, out_hw_eec3_gen_nopileup_jetmettau, lumi_hw_eec3, n_files_hw_eec3, "MC_GEN", "nopileup", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_hw_eec3, out_hw_eec3_gen_allvertex_jetmettau, lumi_hw_eec3, n_files_hw_eec3, "MC_GEN", "allvertex", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);

if (show_steps) { cout << "Herwig pp EEC3 Tune with Vertex Reweight to match JetMET_2010A"<<endl; }
///read_ntuple(mc_hw_eec3, out_hw_eec3_gen_nopileup_jetmet, lumi_hw_eec3, n_files_hw_eec3, "MC_GEN", "nopileup", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_gen_allvertex_jetmet, lumi_hw_eec3, n_files_hw_eec3, "MC_GEN", "allvertex", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);

if (show_steps) { cout << "Herwig pp 23 Tune"<<endl; }
///read_ntuple(mc_hw_23, out_hw_23_gen_nopileup, lumi_hw_23, n_files_hw_23, "MC_GEN", "nopileup", "", "", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_gen_allvertex, lumi_hw_23, n_files_hw_23, "MC_GEN", "allvertex", "", "", detail, test);

if (show_steps) { cout << "Herwig pp 23 Tune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///read_ntuple(mc_hw_23, out_hw_23_gen_nopileup_jetmettau, lumi_hw_23, n_files_hw_23, "MC_GEN", "nopileup", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_hw_23, out_hw_23_gen_allvertex_jetmettau, lumi_hw_23, n_files_hw_23, "MC_GEN", "allvertex", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);

if (show_steps) { cout << "Herwig pp EEC3 Tune with Vertex Reweight to match JetMET_2010A"<<endl; }
///read_ntuple(mc_hw_23, out_hw_23_gen_nopileup_jetmet, lumi_hw_23, n_files_hw_23, "MC_GEN", "nopileup", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_gen_allvertex_jetmet, lumi_hw_23, n_files_hw_23, "MC_GEN", "allvertex", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);


if (show_steps) { cout << "All the MC Ntuples were sucessfully read!"<<endl; }
}


//read the MC Ntuples on the Detector Level
if (read_mc_ntuples_det)
{
if (show_steps) { cout << "Reading MC Ntuples on Detector Level..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_det_dir, "root");
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", "", "", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmettau_up, lumi_p6_z2, n_files_p6_z2, "MC_DET", "up", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmettau_down, lumi_p6_z2, n_files_p6_z2, "MC_DET", "down", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmet_up, lumi_p6_z2, n_files_p6_z2, "MC_DET", "up", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmet_down, lumi_p6_z2, n_files_p6_z2, "MC_DET", "down", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jet_up, lumi_p6_z2, n_files_p6_z2, "MC_DET", "up", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jet_down, lumi_p6_z2, n_files_p6_z2, "MC_DET", "down", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
///read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex", "", "", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex",  vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmettau_up, lumi_p8_4c, n_files_p8_4c, "MC_DET", "up", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmettau_down, lumi_p8_4c, n_files_p8_4c, "MC_DET", "down", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
///read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex",  vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
////read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmet_up, lumi_p8_4c, n_files_p8_4c, "MC_DET", "up", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmet_down, lumi_p8_4c, n_files_p8_4c, "MC_DET", "down", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex",  vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jet_up, lumi_p8_4c, n_files_p8_4c, "MC_DET", "up", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jet_down, lumi_p8_4c, n_files_p8_4c, "MC_DET", "down", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);


if (show_steps) { cout << "Herwig pp - EEC3 Tune"<<endl; }
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_1vertex, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "1vertex", "", "", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "allvertex", "", "", detail, test);

if (show_steps) { cout << "Herwig pp EEC3 Tune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_1vertex_jetmettau, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "1vertex",  vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex_jetmettau, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex_jetmettau_up, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "up", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex_jetmettau_down, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "down", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);


if (show_steps) { cout << "Herwig pp EEC3 Tune with Vertex Reweight to match JetMET_2010A"<<endl; }
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_1vertex_jetmet, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "1vertex",  vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex_jetmet, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "allvertex", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex_jetmet_up, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "up", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_eec3, out_hw_eec3_det_allvertex_jetmet_down, lumi_hw_eec3, n_files_hw_eec3, "MC_DET", "down", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);


if (show_steps) { cout << "Herwig pp - 23 Tune"<<endl; }
///read_ntuple(mc_hw_23, out_hw_23_det_1vertex, lumi_hw_23, n_files_hw_23, "MC_DET", "1vertex", "", "", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex, lumi_hw_23, n_files_hw_23, "MC_DET", "allvertex", "", "", detail, test);

if (show_steps) { cout << "Herwig pp 23 Tune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
///read_ntuple(mc_hw_23, out_hw_23_det_1vertex_jetmettau, lumi_hw_23, n_files_hw_23, "MC_DET", "1vertex",  vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex_jetmettau, lumi_hw_23, n_files_hw_23, "MC_DET", "allvertex", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex_jetmettau_up, lumi_hw_23, n_files_hw_23, "MC_DET", "up", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex_jetmettau_down, lumi_hw_23, n_files_hw_23, "MC_DET", "down", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);


if (show_steps) { cout << "Herwig pp 23 Tune with Vertex Reweight to match JetMET_2010A"<<endl; }
///read_ntuple(mc_hw_23, out_hw_23_det_1vertex_jetmet, lumi_hw_23, n_files_hw_23, "MC_DET", "1vertex",  vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex_jetmet, lumi_hw_23, n_files_hw_23, "MC_DET", "allvertex", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex_jetmet_up, lumi_hw_23, n_files_hw_23, "MC_DET", "up", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);
///read_ntuple(mc_hw_23, out_hw_23_det_allvertex_jetmet_down, lumi_hw_23, n_files_hw_23, "MC_DET", "down", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);


if (show_steps) { cout << "All the MC Ntuples were sucessfully read!"<<endl; }
}


//read the data Ntuples
if (read_data_ntuples)
{
if (show_steps) { cout << "Reading Data Ntuples..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, raw_data_dir, "root");
if (show_steps) { cout << "JetMETTau"<<endl; }
//read_ntuple(jetmettau, out_jetmettau_1vertex, lumi_jetmettau, n_files_jetmettau, "DATA", "1vertex", "", "", detail, test);
//read_ntuple(jetmettau, out_jetmettau_allvertex, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
//read_ntuple(jetmettau, out_jetmettau_up, lumi_jetmettau, n_files_jetmettau, "DATA", "up", "", "", detail, test);
//read_ntuple(jetmettau, out_jetmettau_down, lumi_jetmettau, n_files_jetmettau, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
//read_ntuple(jetmet, out_jetmet_1vertex, lumi_jetmet, n_files_jetmet, "DATA", "1vertex", "", "", detail, test);
//read_ntuple(jetmet, out_jetmet_allvertex, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
//read_ntuple(jetmet, out_jetmet_up, lumi_jetmet, n_files_jetmet, "DATA", "up", "", "", detail, test);
//read_ntuple(jetmet, out_jetmet_down, lumi_jetmet, n_files_jetmet, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "Jet"<<endl; }
//read_ntuple(jet, out_jet_1vertex, lumi_jet, n_files_jet, "DATA", "1vertex", "", "", detail, test);
//read_ntuple(jet, out_jet_allvertex, lumi_jet, n_files_jet, "DATA", "allvertex", "", "", detail, test);
//read_ntuple(jet, out_jet_up, lumi_jet, n_files_jet, "DATA", "up", "", "", detail, test);
//read_ntuple(jet, out_jet_down, lumi_jet, n_files_jet, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "All the Data Ntuples were sucessfully read!"<<endl; }
}


//normalize the different mc dataset
if (check_cross_section)
{
if (show_steps) { cout << "Check Cross Section..."<<endl; }

check_cross_section(out_jetmettau_allvertex, "JetMETTau_2010A", false, out_jetmet_allvertex, "JetMET_2010A", false, out_p8_4c_det_allvertex_jetmettau, "Pythia 8 - 4C", false, true);

check_cross_section(out_hw_eec3_det_allvertex_jetmettau, "Herwig - EEC3", false, out_hw_23_det_allvertex_jetmettau, "Herwig - 23", false, out_p8_4c_det_allvertex_jetmettau, "Pythia 8 - 4C", false, true);

check_cross_section(out_hw_eec3_gen_allvertex_jetmettau, "Herwig - EEC3", true, out_hw_23_gen_allvertex_jetmettau, "Herwig - 23", true, out_p8_4c_gen_allvertex_jetmettau, "Pythia 8 - 4C", true, true);
}


//normalize the different mc dataset
if (normalize_mc)
{
if (show_steps) { cout << "Normalize MC..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_norm_dir, "root");
create_directories(output_dir, mc_norm_plots, "plots");


if (show_steps) { cout << "Generator Level"<<endl; }
if (show_steps) { cout << "Pythia 6 - Z2* Tune"<<endl; }
//normalize_mc(out_p6_z2_gen_allvertex, out_p6_z2_gen_nopileup, mc_norm_plots, out_p6_z2_gen_norm , "pythia6_z2_", "ak5Gen_", detail,show_results);
///normalize_mc(out_p6_z2_gen_allvertex_jetmettau, out_p6_z2_gen_nopileup_jetmettau, mc_norm_plots, out_p6_z2_gen_jetmettau_norm , "pythia6_z2_JetMETTau_2010A_", "ak5Gen_", detail,show_results);
///normalize_mc(out_p6_z2_gen_allvertex_jetmet, out_p6_z2_gen_nopileup_jetmet, mc_norm_plots, out_p6_z2_gen_jetmet_norm , "pythia6_z2_JetMET_2010A_", "ak5Gen_", detail,show_results);
//normalize_mc(out_p6_z2_gen_allvertex_jet, out_p6_z2_gen_nopileup_jetmet, mc_norm_plots, out_p6_z2_gen_jet_norm , "pythia6_z2_Jet_2010B_", "ak5Gen_", detail,show_results);

if (show_steps) { cout << "Pythia 8 - 4C Tune"<<endl; }
normalize_mc(out_p8_4c_gen_allvertex, out_p8_4c_gen_nopileup, mc_norm_plots, out_p8_4c_gen_norm , "pythia8_4c_", "ak5Gen_", detail,show_results);
/normalize_mc(out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_gen_nopileup_jetmettau, mc_norm_plots, out_p8_4c_gen_jetmettau_norm , "pythia8_4c_JetMETTau_2010A", "ak5Gen_", detail,show_results);
normalize_mc(out_p8_4c_gen_allvertex_jetmet, out_p8_4c_gen_nopileup_jetmet, mc_norm_plots, out_p8_4c_gen_jetmet_norm , "pythia8_4c_JetMET_2010A", "ak5Gen_", detail,show_results);
//normalize_mc(out_p8_4c_gen_allvertex_jet, out_p8_4c_gen_nopileup_jetmet, mc_norm_plots, out_p8_4c_gen_jet_norm , "pythia8_4c_Jet_2010B_", "ak5Gen_", detail,show_results);

if (show_steps) { cout << "Herwig pp - EEC3 Tune"<<endl; }
normalize_mc(out_hw_eec3_gen_allvertex, out_hw_eec3_gen_nopileup, mc_norm_plots, out_hw_eec3_gen_norm , "herwig_eec3_", "ak5Gen_", detail,show_results);
normalize_mc(out_hw_eec3_gen_allvertex_jetmettau, out_hw_eec3_gen_nopileup_jetmettau, mc_norm_plots, out_hw_eec3_gen_jetmettau_norm , "herwig_eec3_JetMETTau_2010A", "ak5Gen_", detail,show_results);
normalize_mc(out_hw_eec3_gen_allvertex_jetmet, out_hw_eec3_gen_nopileup_jetmet, mc_norm_plots, out_hw_eec3_gen_jetmet_norm , "herwig_eec3_JetMET_2010A", "ak5Gen_", detail,show_results);


if (show_steps) { cout << "Herwig pp - 23 Tune"<<endl; }
normalize_mc(out_hw_23_gen_allvertex, out_hw_23_gen_nopileup, mc_norm_plots, out_hw_23_gen_norm , "herwig_23_", "ak5Gen_", detail,show_results);
normalize_mc(out_hw_23_gen_allvertex_jetmettau, out_hw_23_gen_nopileup_jetmettau, mc_norm_plots, out_hw_23_gen_jetmettau_norm , "herwig_23_JetMETTau_2010A", "ak5Gen_", detail,show_results);
///normalize_mc(out_hw_23_gen_allvertex_jetmet, out_hw_23_gen_nopileup_jetmet, mc_norm_plots, out_hw_23_gen_jetmet_norm , "herwig_23_JetMET_2010A", "ak5Gen_", detail,show_results);


if (show_steps) { cout << "Detector Level"<<endl; }
if (show_steps) { cout << "Pythia 6 - Z2* Tune"<<endl; }
///normalize_mc(out_p6_z2_det_allvertex, out_p6_z2_det_1vertex, mc_norm_plots, out_p6_z2_det_norm , "pythia6_z2_", "ak5PF_", detail,show_results);
///normalize_mc(out_p6_z2_det_allvertex_jetmettau, out_p6_z2_det_1vertex_jetmettau, mc_norm_plots, out_p6_z2_det_jetmettau_norm , "pythia6_z2_JetMETTau_2010A_", "ak5PF_", detail,show_results);
///normalize_mc(out_p6_z2_det_allvertex_jetmet, out_p6_z2_det_1vertex_jetmet, mc_norm_plots, out_p6_z2_det_jetmet_norm , "pythia6_z2_JetMET_2010A_", "ak5PF_", detail,show_results);
//normalize_mc(out_p6_z2_det_allvertex_jet, out_p6_z2_det_1vertex_jet, mc_norm_plots, out_p6_z2_det_jet_norm , "pythia6_z2_Jet_2010B_", "ak5PF_", detail,show_results);

if (show_steps) { cout << "Pythia 8 - 4C Tune"<<endl; }
//normalize_mc(out_p8_4c_det_allvertex, out_p8_4c_det_1vertex, mc_norm_plots, out_p8_4c_det_norm , "pythia8_4c_", "ak5PF_", detail,show_results);
///normalize_mc(out_p8_4c_det_allvertex_jetmettau, out_p8_4c_det_1vertex_jetmettau, mc_norm_plots, out_p8_4c_det_jetmettau_norm , "pythia8_4c_JetMETTau_2010A_", "ak5PF_", detail,show_results);
///normalize_mc(out_p8_4c_det_allvertex_jetmet, out_p8_4c_det_1vertex_jetmet, mc_norm_plots, out_p8_4c_det_jetmet_norm , "pythia8_4c_JetMET_2010A_", "ak5PF_", detail,show_results);
///normalize_mc(out_p8_4c_det_allvertex_jet, out_p8_4c_det_1vertex_jet, mc_norm_plots, out_p8_4c_det_jet_norm , "pythia8_4c_Jet_2010B_", "ak5PF_", detail,show_results);

if (show_steps) { cout << "Herwig pp - EEC3 Tune"<<endl; }
///normalize_mc(out_hw_eec3_det_allvertex, out_hw_eec3_det_1vertex, mc_norm_plots, out_hw_eec3_det_norm , "herwig_eec3_", "ak5PF_", detail,show_results);
///normalize_mc(out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_det_1vertex_jetmettau, mc_norm_plots, out_hw_eec3_det_jetmettau_norm , "herwig_eec3_JetMETTau_2010A_", "ak5PF_", detail,show_results);
///normalize_mc(out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_det_1vertex_jetmet, mc_norm_plots, out_hw_eec3_det_jetmet_norm , "herwig_eec3_JetMET_2010A_", "ak5PF_", detail,show_results);

if (show_steps) { cout << "Herwig pp - 23 Tune"<<endl; }
///normalize_mc(out_hw_23_det_allvertex, out_hw_23_det_1vertex, mc_norm_plots, out_hw_23_det_norm , "herwig_23_", "ak5PF_", detail,show_results);
///normalize_mc(out_hw_23_det_allvertex_jetmettau, out_hw_23_det_1vertex_jetmettau, mc_norm_plots, out_hw_23_det_jetmettau_norm , "herwig_23_JetMETTau_2010A_", "ak5PF_", detail,show_results);
///normalize_mc(out_hw_23_det_allvertex_jetmet, out_hw_23_det_1vertex_jetmet, mc_norm_plots, out_hw_23_det_jetmet_norm , "herwig_23_JetMET_2010A_", "ak5PF_", detail,show_results);

if (show_steps) { cout << "The MC datasets were sucessfully normalized!"<<endl; }
}

//compute correction
if (compute_corrections)
{
if (show_steps) { cout << "Computing Corrections..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, corrections_dir, "root");
create_directories(output_dir, corrections_plots, "plots");
create_directories(output_dir, sel_corrections_plots, "plots");

if (show_steps) { cout << "Computing the Pileup Correction - No reweighting..."<<endl; }
//compute_correction(out_p6_z2_gen_nopileup, "ak5Gen_", "Pythia 6 Tune Z2* - Generator Level with no pileup", out_p6_z2_gen_allvertex, "ak5Gen_", "Pythia 6 Tune Z2* - Generator Level with all vertex", corrections_plots, corr_pileup_p6_z2, "pileup_pythia6_z2_", "Pileup Correction with Pythia 6 - Tune Z2*", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup, "ak5Gen_", "Pythia 8 Tune 4C - Generator Level with no pileup", out_p8_4c_gen_allvertex, "ak5Gen_", "Pythia 8 Tune 4C - Generator Level with all vertex", corrections_plots, corr_pileup_p8_4c, "pileup_pythia8_4c_", "Pileup Correction with Pythia 8 - Tune 4C", false, detail, show_results);
compute_correction(out_hw_eec3_gen_nopileup, "ak5Gen_", "Herwig Tune EEC3 - Generator Level with no pileup", out_hw_eec3_gen_allvertex, "ak5Gen_", "Herwig Tune EEC3 - Generator Level with all vertex", corrections_plots, corr_pileup_hw_eec3, "pileup_herwig_eec3_", "Pileup Correction with Herwig Tune EEC3", false, detail, show_results);
compute_correction(out_hw_23_gen_nopileup, "ak5Gen_", "Herwig Tune 23 - Generator Level with no pileup", out_hw_23_gen_allvertex, "ak5Gen_", "Herwig Tune 23 - Generator Level with all vertex", corrections_plots, corr_pileup_hw_23, "pileup_herwig_23_", "Pileup Correction with Herwig Tune 23", false, detail, show_results);

if (show_steps) { cout << "Computing the Pileup Correction - JetMETTau..."<<endl; }
//compute_correction(out_p6_z2_gen_nopileup_jetmettau, "ak5Gen_", "Pythia 6 Tune Z2* -  Gen. Level no pileup - JetMETTau_2010A", out_p6_z2_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level all vertex - JetMETTau_2010A", corrections_plots, corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pileup Correction with Pythia 6 - Tune Z2* - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup_jetmettau, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level no pileup - JetMETTau_2010A", out_p8_4c_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level all vertex - JetMETTau_2010A", corrections_plots, corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pileup Correction with Pythia 8 - Tune 4C - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_hw_eec3_gen_nopileup_jetmettau, "ak5Gen_", "Herwig Tune EEC3 -  Gen. Level no pileup - JetMETTau_2010A", out_hw_eec3_gen_allvertex_jetmettau, "ak5Gen_", "Herwig Tune EEC3 - Gen. Level all vertex - JetMETTau_2010A", corrections_plots, corr_pileup_hw_eec3_jetmettau, "pileup_herwig_eec3_JetMETTau_2010A_", "Pileup Correction with Herwig Tune EEC3 - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_hw_23_gen_nopileup_jetmettau, "ak5Gen_", "Herwig Tune 23 -  Gen. Level no pileup - JetMETTau_2010A", out_hw_23_gen_allvertex_jetmettau, "ak5Gen_", "Herwig Tune 23 - Gen. Level all vertex - JetMETTau_2010A", corrections_plots, corr_pileup_hw_23_jetmettau, "pileup_herwig_23_JetMETTau_2010A_", "Pileup Correction with Herwig Tune 23 - JetMETTau_2010A", false, detail, show_results);

if (show_steps) { cout << "Computing the Pileup Correction - JetMET..."<<endl; }
//compute_correction(out_p6_z2_gen_nopileup_jetmet, "ak5Gen_", "Pythia 6 Tune Z2* -  Gen. Level no pileup JetMET_2010A", out_p6_z2_gen_allvertex_jetmet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level all vertex - JetMET_2010A", corrections_plots, corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pileup Correction with Pythia 6 - Tune Z2* - JetMET_2010A", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup_jetmet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level no pileup JetMET_2010A", out_p8_4c_gen_allvertex_jetmet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level all vertex - JetMET_2010A", corrections_plots, corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pileup Correction with Pythia 8 - Tune 4C - JetMET_2010A", false, detail, show_results);
compute_correction(out_hw_eec3_gen_nopileup_jetmet, "ak5Gen_", "Herwig Tune EEC3 - Gen. Level no pileup JetMET_2010A", out_hw_eec3_gen_allvertex_jetmet, "ak5Gen_", "Herwig Tune EEC3 - Gen. Level all vertex - JetMET_2010A", corrections_plots, corr_pileup_hw_eec3_jetmet, "pileup_herwig_eec3_JetMET_2010A_", "Pileup Correction with Herwig - Tune EEC3 - JetMET_2010A", false, detail, show_results);
compute_correction(out_hw_23_gen_nopileup_jetmet, "ak5Gen_", "Herwig Tune 23 - Gen. Level no pileup JetMET_2010A", out_hw_23_gen_allvertex_jetmet, "ak5Gen_", "Herwig Tune 23 - Gen. Level all vertex - JetMET_2010A", corrections_plots, corr_pileup_hw_23_jetmet, "pileup_herwig_23_JetMET_2010A_", "Pileup Correction with Herwig - Tune 23 - JetMET_2010A", false, detail, show_results);

if (show_steps) { cout << "Computing the Pileup Correction - Jet..."<<endl; }
//compute_correction(out_p6_z2_gen_nopileup_jet, "ak5Gen_", "Pythia 6 Tune Z2* -  Gen. Level no pileup Jet_2010B", out_p6_z2_gen_allvertex_jet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level all vertex - Jet_2010B", corrections_plots, corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pileup Correction with Pythia 6 - Tune Z2* - Jet_2010B", false, detail, show_results);
//compute_correction(out_p8_4c_gen_nopileup_jet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level no pileup Jet_2010B", out_p8_4c_gen_allvertex_jet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level all vertex - Jet_2010B", corrections_plots, corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pileup Correction with Pythia 8 - Tune 4C - Jet_2010B", false, detail, show_results);

if (show_steps) { cout << "Computing the Detector Correction - No Reweighting..."<<endl; }
//compute_correction(out_p6_z2_det_allvertex, "ak5PF_", "Pythia 6 Tune Z2* - Detector Level", out_p6_z2_gen_allvertex, "ak5Gen_", "Pythia 6 Tune Z2* - Generator Level", corrections_plots, corr_detector_p6_z2, "detector_pythia6_z2_", "Detector Correction with Pythia 6 - Tune Z2*", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex, "ak5PF_", "Pythia 8  - Tune 4C - Detector Level", out_p8_4c_gen_allvertex, "ak5Gen_", "Pythia 8 Tune 4C - Generator Level", corrections_plots, corr_detector_p8_4c, "detector_pythia8_4c_", "Detector Correction with Pythia 8 - Tune 4C", false, detail, show_results);
compute_correction(out_hw_eec3_det_allvertex, "ak5PF_", "Herwig Tune EEC3 - Detector Level", out_hw_eec3_gen_allvertex, "ak5Gen_", "Herwig Tune EEC3 - Generator Level", corrections_plots, corr_detector_hw_eec3, "detector_herwig_eec3_", "Detector Correction with Herwig - Tune EEC3", false, detail, show_results);
compute_correction(out_hw_23_det_allvertex, "ak5PF_", "Herwig Tune 23 - Detector Level", out_hw_23_gen_allvertex, "ak5Gen_", "Herwig Tune 23 - Generator Level", corrections_plots, corr_detector_hw_23, "detector_herwig_23_", "Detector Correction with Herwig - Tune 23", false, detail, show_results);

if (show_steps) { cout << "Computing the Detector Correction - JetMETTau..."<<endl; }
//compute_correction(out_p6_z2_det_allvertex_jetmettau, "ak5PF_", "Pythia 6 Tune Z2* - Det. Level JetMETTau_2010A", out_p6_z2_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level JetMETTau_2010A", corrections_plots, corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Detector Correction with Pythia 6 - Tune Z2* - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex_jetmettau, "ak5PF_", "Pythia 8 - Tune 4C - Det. Level JetMETTau_2010A", out_p8_4c_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level JetMETTau_2010A", corrections_plots, corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Detector Correction with Pythia 8 - Tune 4C - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_hw_eec3_det_allvertex_jetmettau, "ak5PF_", "Herwig - Tune EEC3 - Det. Level JetMETTau_2010A", out_hw_eec3_gen_allvertex_jetmettau, "ak5Gen_", "Herwig - Tune EEC3 - Gen. Level JetMETTau_2010A", corrections_plots, corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", "Detector Correction with Herwig - Tune EEC3 - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_hw_23_det_allvertex_jetmettau, "ak5PF_", "Herwig - Tune 23 - Det. Level JetMETTau_2010A", out_hw_23_gen_allvertex_jetmettau, "ak5Gen_", "Herwig - Tune 23 - Gen. Level JetMETTau_2010A", corrections_plots, corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", "Detector Correction with Herwig - Tune 23 - JetMETTau_2010A", false, detail, show_results);

if (show_steps) { cout << "Computing the Detector Correction - JetMET..."<<endl; }
//compute_correction(out_p6_z2_det_allvertex_jetmet, "ak5PF_", "Pythia 6 Tune Z2* - Det. Level JetMET_2010A", out_p6_z2_gen_allvertex_jetmet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level JetMET_2010A", corrections_plots, corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Detector Correction with Pythia 6 - Tune Z2* - JetMET_2010A", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex_jetmet, "ak5PF_", "Pythia 8 - Tune 4C - Det. Level JetMET_2010A", out_p8_4c_gen_allvertex_jetmet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level JetMET_2010A", corrections_plots, corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Detector Correction with Pythia 8 - Tune 4C - JetMET_2010A", false, detail, show_results);
compute_correction(out_hw_eec3_det_allvertex_jetmet, "ak5PF_", "Herwig - Tune EEC3 - Det. Level JetMET_2010A", out_hw_eec3_gen_allvertex_jetmet, "ak5Gen_", "Herwig - Tune EEC3 - Gen. Level JetMET_2010A", corrections_plots, corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", "Detector Correction with Herwig - Tune EEC3 - JetMET_2010A", false, detail, show_results);
compute_correction(out_hw_23_det_allvertex_jetmet, "ak5PF_", "Herwig - Tune 23 - Det. Level JetMET_2010A", out_hw_23_gen_allvertex_jetmet, "ak5Gen_", "Herwig - Tune 23 - Gen. Level JetMET_2010A", corrections_plots, corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", "Detector Correction with Herwig - Tune 23 - JetMET_2010A", false, detail, show_results);


if (show_steps) { cout << "Computing the Detector Correction - Jet..."<<endl; }
//compute_correction(out_p6_z2_det_allvertex_jet, "ak5PF_", "Pythia 6 Tune Z2* - Det. Level Jet_2010B", out_p6_z2_gen_allvertex_jet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level Jet_2010B", corrections_plots, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Detector Correction with Pythia 6 - Tune Z2* - Jet_2010B", false, detail, show_results);
//compute_correction(out_p8_4c_det_allvertex_jet, "ak5PF_", "Pythia 8  - Tune 4C - Det. Level Jet_2010B", out_p8_4c_gen_allvertex_jet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level Jet_2010B", corrections_plots, corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Detector Correction with Pythia 8 - Tune 4C - Jet_2010B", false, detail, show_results);

if (show_steps) { cout << "Computing the Overall Correction - No Rewieghting..."<<endl; }
//compute_correction(corr_detector_p6_z2, "detector_pythia6_z2_", "Pythia 6 Tune Z2* - Detector Correction", corr_pileup_p6_z2, "pileup_pythia6_z2_", "Pythia 6 Tune Z2* - Pileup Correction", corrections_plots, corr_final_p6_z2, "final_pythia6_z2_", "Final Correction with Pythia 6 - Tune Z2*", true, detail, show_results);
compute_correction(corr_detector_p8_4c, "detector_pythia8_4c_", "Pythia 8 - Tune 4C - Detector Correction", corr_pileup_p8_4c, "pileup_pythia8_4c_", "Pythia 8 Tune 4C - Pileup Correction", corrections_plots, corr_final_p8_4c, "final_pythia8_4c_", "Final Correction with Pythia 8 - Tune 4C", true, detail, show_results);
compute_correction(corr_detector_hw_eec3, "detector_herwig_eec3_", "Herwig - Tune EEC3 - Detector Correction", corr_pileup_hw_eec3, "pileup_herwig_eec3_", "Herwig Tune EEC3 - Pileup Correction", corrections_plots, corr_final_hw_eec3, "final_herwig_eec3_", "Final Correction with Herwig Tune EEC3", true, detail, show_results);
compute_correction(corr_detector_hw_23, "detector_herwig_23_", "Herwig - Tune 23 - Detector Correction", corr_pileup_hw_23, "pileup_herwig_23_", "Herwig Tune 23 - Pileup Correction", corrections_plots, corr_final_hw_23, "final_herwig_23_", "Final Correction with Herwig Tune 23", true, detail, show_results);


if (show_steps) { cout << "Computing the Overall Correction - JetMETTau..."<<endl; }
//compute_correction(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Pythia 6 Tune Z2* - Detector Correction - JetMETTau_2010A", corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pythia 6 Tune Z2* - Pileup Correction - JetMETTau_2010A", corrections_plots, corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "Final Correction with Pythia 6 - Tune Z2* - JetMETTau_2010A", true, detail, show_results);
compute_correction(corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Pythia 8 Tune 4C - Detector Correction - JetMETTau_2010A", corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pythia 8 Tune 4C - Pileup Correction - JetMETTau_2010A", corrections_plots, corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "Final Correction with Pythia 8 - Tune 4C - JetMETTau_2010A", true, detail, show_results);
compute_correction(corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", "Herwig Tune EEC3 - Detector Correction - JetMETTau_2010A", corr_pileup_hw_eec3_jetmettau, "pileup_herwig_eec3_JetMETTau_2010A_", "Herwig Tune EEC3 - Pileup Correction - JetMETTau_2010A", corrections_plots, corr_final_hw_eec3_jetmettau, "final_herwig_eec3_JetMETTau_2010A_", "Final Correction with Herwig - Tune EEC3 - JetMETTau_2010A", true, detail, show_results);
compute_correction(corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", "Herwig Tune 23 - Detector Correction - JetMETTau_2010A", corr_pileup_hw_23_jetmettau, "pileup_herwig_23_JetMETTau_2010A_", "Herwig Tune 23 - Pileup Correction - JetMETTau_2010A", corrections_plots, corr_final_hw_23_jetmettau, "final_herwig_23_JetMETTau_2010A_", "Final Correction with Herwig - Tune 23 - JetMETTau_2010A", true, detail, show_results);

if (show_steps) { cout << "Computing the Overall Correction - JetMET..."<<endl; }
//compute_correction(corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Pythia 6 Tune Z2* - Detector Correction - JetMET_2010A", corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pythia 6 Tune Z2* - Pileup Correction - JetMET_2010A", corrections_plots, corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "Final Correction with Pythia 6 - Tune Z2* - JetMET_2010A", true, detail, show_results);
compute_correction(corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Pythia 8 Tune 4C - Detector Correction - JetMET_2010A", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pythia 8 Tune 4C - Pileup Correction - JetMET_2010A", corrections_plots, corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "Final Correction with Pythia 8 - Tune 4C - JetMET_2010A", true, detail, show_results);
compute_correction(corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", "Herwig Tune EEC3 - Detector Correction - JetMET_2010A", corr_pileup_hw_eec3_jetmet, "pileup_herwig_eec3_JetMET_2010A_", "Herwig Tune EEC3 - Pileup Correction - JetMET_2010A", corrections_plots, corr_final_hw_eec3_jetmet, "final_herwig_eec3_JetMET_2010A_", "Final Correction with Herwig - Tune EEC3 - JetMET_2010A", true, detail, show_results);
compute_correction(corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", "Herwig Tune 23 - Detector Correction - JetMET_2010A", corr_pileup_hw_23_jetmet, "pileup_herwig_23_JetMET_2010A_", "Herwig Tune 23 - Pileup Correction - JetMET_2010A", corrections_plots, corr_final_hw_23_jetmet, "final_herwig_23_JetMET_2010A_", "Final Correction with Herwig - Tune 23 - JetMET_2010A", true, detail, show_results);


if (show_steps) { cout << "Computing the Overall Correction - Jet..."<<endl; }
//compute_correction(corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Pythia 6 Tune Z2* - Detector Correction - Jet_2010B", corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pythia 6 Tune Z2* - Pileup Correction - Jet_2010B", corrections_plots, corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Final Correction with Pythia 6 - Tune Z2* - Jet_2010B", true, detail, show_results);
//compute_correction(corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Pythia 8 Tune 4C - Detector Correction - Jet_2010B", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pythia 8 Tune 4C - Pileup Correction - Jet_2010B", corrections_plots, corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Final Correction with Pythia 8 - Tune 4C - Jet_2010B", true, detail, show_results);

if (show_steps) { cout << "Ploting Final Results..."<<endl; }
if (show_steps) { cout << "Ploting Final Results - Pythia 6..."<<endl; }
//display_final_corrections(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Detector Correction", corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pileup Correction", corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "Final Correction", corrections_plots, "corrections_pythia6_z2_JetMETTau_2010A_", detail);
//display_final_corrections(corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Detector Correction", corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pileup Correction", corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "Final Correction", corrections_plots, "corrections_pythia6_z2_JetMET_2010A_", detail);
//display_final_corrections(corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Detector Correction", corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pileup Correction", corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Final Correction", corrections_plots, "corrections_pythia6_z2_Jet_2010B_", detail);

//display_final_corrections(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "JetMETTau_2010A", corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "JetMET_2010A", corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia6_z2_detector_", detail);
//display_final_corrections(corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "JetMETTau_2010A", corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "JetMET_2010A", corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia6_z2_pileup_", detail);
//display_final_corrections(corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "JetMETTau_2010A", corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "JetMET_2010A", corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia6_z2_final_", detail);

if (show_steps) { cout << "Ploting Final Results - Pythia 8..."<<endl; }
display_final_corrections(corr_detector_p8_4c, "detector_pythia8_4c_", "Detector Correction", corr_pileup_p8_4c, "pileup_pythia8_4c_", "Pileup Correction", corr_final_p8_4c, "final_pythia8_4c_", "Final Correction", sel_corrections_plots, "corrections_pythia8_4c_", detail);
display_final_corrections(corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Detector Correction", corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pileup Correction", corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "Final Correction", sel_corrections_plots, "corrections_pythia8_4c_JetMETTau_2010A_", detail);
display_final_corrections(corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Detector Correction", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pileup Correction", corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "Final Correction", sel_corrections_plots, "corrections_pythia8_4c_JetMET_2010A_", detail);
//display_final_corrections(corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Detector Correction", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pileup Correction", corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Final Correction", corrections_plots, "corrections_pythia8_4c_Jet_2010B_", detail);

display_final_corrections(corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "JetMETTau_2010A", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "JetMET_2010A", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Jet_2010B", sel_corrections_plots, "corrections_pythia8_4c_detector_", detail);
display_final_corrections(corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMETTau_2010A_", "JetMETTau_2010A", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "JetMET_2010A", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Jet_2010B", sel_corrections_plots, "corrections_pythia8_4c_pileup_", detail);
display_final_corrections(corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "JetMETTau_2010A", corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "JetMET_2010A", corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Jet_2010B", sel_corrections_plots, "corrections_pythia8_4c_final_", detail);


if (show_steps) { cout << "Ploting Final Results - Herwig EEC3..."<<endl; }
display_final_corrections(corr_detector_hw_eec3, "detector_herwig_eec3_", "Detector Correction", corr_pileup_hw_eec3, "pileup_herwig_eec3_", "Pileup Correction", corr_final_hw_eec3, "final_herwig_eec3_", "Final Correction", sel_corrections_plots, "corrections_herwig_eec3_", detail);
display_final_corrections(corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", "Detector Correction", corr_pileup_hw_eec3_jetmettau, "pileup_herwig_eec3_JetMETTau_2010A_", "Pileup Correction", corr_final_hw_eec3_jetmettau, "final_herwig_eec3_JetMETTau_2010A_", "Final Correction", sel_corrections_plots, "corrections_herwig_eec3_JetMETTau_2010A_", detail);
display_final_corrections(corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", "Detector Correction", corr_pileup_hw_eec3_jetmet, "pileup_herwig_eec3_JetMET_2010A_", "Pileup Correction", corr_final_hw_eec3_jetmet, "final_herwig_eec3_JetMET_2010A_", "Final Correction", sel_corrections_plots, "corrections_herwig_eec3_JetMET_2010A_", detail);

display_final_corrections(corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", "JetMETTau_2010A", corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", "JetMET_2010A", "", "", "", sel_corrections_plots, "corrections_herwig_eec3_detector_", detail);
display_final_corrections(corr_pileup_hw_eec3_jetmet, "pileup_herwig_eec3_JetMETTau_2010A_", "JetMETTau_2010A", corr_pileup_hw_eec3_jetmet, "pileup_herwig_eec3_JetMET_2010A_", "JetMET_2010A", "", "", "", sel_corrections_plots, "corrections_herwig_eec3_pileup_", detail);
display_final_corrections(corr_final_hw_eec3_jetmettau, "final_herwig_eec3_JetMETTau_2010A_", "JetMETTau_2010A", corr_final_hw_eec3_jetmet, "final_herwig_eec3_JetMET_2010A_", "JetMET_2010A", "", "", "", sel_corrections_plots, "corrections_herwig_eec3_final_", detail);

if (show_steps) { cout << "Ploting Final Results - Herwig 23..."<<endl; }
display_final_corrections(corr_detector_hw_23, "detector_herwig_23_", "Detector Correction", corr_pileup_hw_23, "pileup_herwig_23_", "Pileup Correction", corr_final_hw_23, "final_herwig_23_", "Final Correction", sel_corrections_plots, "corrections_herwig_23_", detail);
display_final_corrections(corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", "Detector Correction", corr_pileup_hw_23_jetmettau, "pileup_herwig_23_JetMETTau_2010A_", "Pileup Correction", corr_final_hw_23_jetmettau, "final_herwig_23_JetMETTau_2010A_", "Final Correction", sel_corrections_plots, "corrections_herwig_23_JetMETTau_2010A_", detail);
display_final_corrections(corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", "Detector Correction", corr_pileup_hw_23_jetmet, "pileup_herwig_23_JetMET_2010A_", "Pileup Correction", corr_final_hw_23_jetmet, "final_herwig_23_JetMET_2010A_", "Final Correction", sel_corrections_plots, "corrections_herwig_23_JetMET_2010A_", detail);

display_final_corrections(corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", "JetMETTau_2010A", corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", "JetMET_2010A", "", "", "", sel_corrections_plots, "corrections_herwig_23_detector_", detail);
display_final_corrections(corr_pileup_hw_23_jetmettau, "pileup_herwig_23_JetMETTau_2010A_", "JetMETTau_2010A", corr_pileup_hw_23_jetmet, "pileup_herwig_23_JetMET_2010A_", "JetMET_2010A", "", "", "", sel_corrections_plots, "corrections_herwig_23_pileup_", detail);
display_final_corrections(corr_final_hw_23_jetmettau, "final_herwig_23_JetMETTau_2010A_", "JetMETTau_2010A", corr_final_hw_23_jetmet, "final_herwig_23_JetMET_2010A_", "JetMET_2010A", "", "", "", sel_corrections_plots, "corrections_herwig_23_final_", detail);



if (show_steps) { cout << "Ploting Final Results - Detector..."<<endl; }

///display_final_corrections(corr_detector_hw_eec3_jetmettau, "detector_hw_eec3_JetMETTau_2010A_", "Herwig - Tune EEC3", corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Pythia8 - 4C", corr_detector_hw_23_jetmettau, "detector_hw_23_JetMETTau_2010A_", "Herwig - Tune 23", sel_corrections_plots, "corrections_detector_JetMETTau_2010A_", detail);
///display_final_corrections(corr_detector_hw_eec3_jetmet, "detector_hw_eec3_JetMET_2010A_", "Herwig - EEC3", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Pythia8 - 4C", corr_detector_hw_23_jetmettau, "detector_hw_23_JetMETTau_2010A_", "Herwig - Tune 23", sel_corrections_plots, "corrections_detector_JetMET_2010_", detail);
//display_final_corrections(corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Pythia6 - Z2*", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_detector_Jet_2010B_", detail);

if (show_steps) { cout << "Ploting Final Results - Pileup..."<<endl; }
///display_final_corrections(corr_pileup_hw_eec3_jetmettau, "pileup_hw_eec3_JetMETTau_2010A_", "Herwig - Tune EEC3", corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pythia8 - 4C", "", "", "", sel_corrections_plots, "corrections_pileup_JetMETTau_2010A_", detail);
///display_final_corrections(corr_pileup_hw_eec3_jetmet, "pileup_hw_eec3_JetMET_2010A_", "Herwig - EEC3", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pythia8 - 4C", corr_pileup_hw_23_jetmettau, "pileup_hw_23_JetMETTau_2010A_", "Herwig - Tune 23", sel_corrections_plots, "corrections_pileup_JetMET_2010_", detail);
//display_final_corrections(corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pythia6 - Z2*", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_pileup_Jet_2010B_", detail);

if (show_steps) { cout << "Ploting Final Results - Final..."<<endl; }
///display_final_corrections(corr_final_hw_eec3_jetmettau, "final_hw_eec3_JetMETTau_2010A_", "Herwig - Tune EEC3", corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "Pythia8 - 4C", corr_final_hw_23_jetmettau, "final_hw_23_JetMETTau_2010A_", "Herwig - Tune 23", sel_corrections_plots, "corrections_final_JetMETTau_2010A_", detail);
///display_final_corrections(corr_final_hw_eec3_jetmet, "final_hw_eec3_JetMET_2010A_", "Herwig - Tune EEC3", corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "Pythia8 - 4C", corr_final_hw_23_jetmet, "final_hw_23_JetMET_2010A_", "Herwig - Tune 23", sel_corrections_plots, "corrections_final_JetMET_2010_", detail);
//display_final_corrections(corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Pythia6 - Z2*", corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_final_Jet_2010B_", detail);
if (show_steps) { cout << "The correction factors were sucessfully computed!"<<endl; }
}


//get mc matched events
if (get_mc_matched_events)
{
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, matched_dir, "root");
if (show_steps) { cout << "Get MC Matched Events..."<<endl; }
if (show_steps) { cout << "Get MC Matched Events Pythia 6 - Z2*..."<<endl; }
///get_matched_events(mc_p6_z2, matched_p6_z2, lumi_p6_z2, n_files_p6_z2, "allvertex", "", "", detail, test);
//get_matched_events(mc_p6_z2, matched_p6_z2_jetmettau, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//get_matched_events(mc_p6_z2, matched_p6_z2_jetmet, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
///get_matched_events(mc_p6_z2, matched_p6_z2_jet, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);

if (show_steps) { cout << "Get MC Matched Events Pythia 8 - 4C..."<<endl; }
///get_matched_events(mc_p8_4c, matched_p8_4c, lumi_p8_4c, n_files_p8_4c, "allvertex", "", "", detail, test);
///get_matched_events(mc_p8_4c, matched_p8_4c_jetmettau, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
///get_matched_events(mc_p8_4c, matched_p8_4c_jetmet, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
///get_matched_events(mc_p8_4c, matched_p8_4c_jet, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);

if (show_steps) { cout << "Get MC Matched Events Herwig - EEC3..."<<endl; }
///get_matched_events(mc_hw_eec3, matched_hw_eec3, lumi_hw_eec3, n_files_hw_eec3, "allvertex", "", "", detail, test);
///get_matched_events(mc_hw_eec3, matched_hw_eec3_jetmettau, lumi_hw_eec3, n_files_hw_eec3, "allvertex", vertex_weights_hw_eec3_jetmettau_v5, "_v5", detail, test);
///get_matched_events(mc_hw_eec3, matched_hw_eec3_jetmet, lumi_hw_eec3, n_files_hw_eec3, "allvertex", vertex_weights_hw_eec3_jetmet_v5, "_v5", detail, test);

if (show_steps) { cout << "Get MC Matched Events Herwig - 23..."<<endl; }
///get_matched_events(mc_hw_23, matched_hw_23, lumi_hw_23, n_files_hw_23, "allvertex", "", "", detail, test);
///get_matched_events(mc_hw_23, matched_hw_23_jetmettau, lumi_hw_23, n_files_hw_23, "allvertex", vertex_weights_hw_23_jetmettau_v5, "_v5", detail, test);
get_matched_events(mc_hw_23, matched_hw_23_jetmet, lumi_hw_23, n_files_hw_23, "allvertex", vertex_weights_hw_23_jetmet_v5, "_v5", detail, test);

if (show_steps) { cout << "All MC Matched Events computed sucessfully!"<<endl; }
}


//compute psba
if (compute_psba)
{
create_directories(output_dir, psba_plots, "plots");
if (show_steps) { cout << "Compute PSBA..."<<endl; }
if (show_steps) { cout << "Compute PSBA - Pythia 6 - Tune Z2*..."<<endl; }
//compute_psba(matched_p6_z2, corr_detector_p6_z2, "detector_pythia6_z2_", psba_plots, "pythia6_z2_", detail, test);
//compute_psba(matched_p6_z2_jetmettau, corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", psba_plots, "pythia6_z2_jetmettau_", detail, test);
//compute_psba(matched_p6_z2_jetmet, corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", psba_plots, "pythia6_z2_jetmet_", detail, test);
//compute_psba(matched_p6_z2_jet, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", psba_plots, "pythia6_z2_jet_", detail, test);
if (show_steps) { cout << "Compute PSBA - Pythia 8 - Tune 4C..."<<endl; }
///compute_psba(matched_p8_4c, corr_detector_p8_4c, "detector_pythia8_4c_", psba_plots, "pythia8_4c_", detail, test);
compute_psba(matched_p8_4c_jetmettau, corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", psba_plots, "pythia8_4c_jetmettau_", detail, test);
//compute_psba(matched_p8_4c_jetmet, corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", psba_plots, "pythia8_4c_jetmet_", detail, test);
//compute_psba(matched_p8_4c_jet, corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", psba_plots, "pythia8_4c_jet_", detail, test);
if (show_steps) { cout << "Compute PSBA - Herwig pp - Tune EEC3..."<<endl; }
//compute_psba(matched_hw_eec3, corr_detector_hw_eec3, "detector_herwig_eec3_", psba_plots, "herwig_eec3_", detail, test);
//compute_psba(matched_hw_eec3_jetmettau, corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", psba_plots, "herwig_eec3_jetmettau_", detail, test);
//compute_psba(matched_hw_eec3_jetmet, corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", psba_plots, "herwig_eec3_jetmet_", detail, test);
if (show_steps) { cout << "Compute PSBA - Herwig pp - Tune 23..."<<endl; }
//compute_psba(matched_hw_23, corr_detector_hw_23, "detector_herwig_23_", psba_plots, "herwig_23_", detail, test);
//compute_psba(matched_hw_23_jetmettau, corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", psba_plots, "herwig_23_jetmettau_", detail, test);
//compute_psba(matched_hw_23_jetmet, corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", psba_plots, "herwig_23_jetmet_", detail, test);
if (show_steps) { cout << "PSBA Computed sucessfully!"<<endl; }
}

//apply correction
if (apply_corrections)
{
if (show_steps) { cout << "Apply Corrections..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, corrected_dir, "root");
create_directories(output_dir, closure_dir, "root");
create_directories(output_dir, corrected_plots, "plots");
//offical jetmettau
//apply_correction(out_jetmettau_allvertex, corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", corrected_jetmettau, corrected_plots, "corrected_jetmettau_", detail, test);

//individual jetmettau p6
apply_correction(out_jetmettau_allvertex, corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", corrected_p6_z2_jetmettau, corrected_plots, "corrected_p6_z2_jetmettau_", detail, test);

//individual jetmettau p8
apply_correction(out_jetmettau_allvertex, corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", corrected_p8_4c_jetmettau, corrected_plots, "corrected_p8_4c_jetmettau_", detail, test);

//individual jetmettau hw eec3
apply_correction(out_jetmettau_allvertex, corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", corr_detector_hw_eec3_jetmettau, "detector_herwig_eec3_JetMETTau_2010A_", corrected_hw_eec3_jetmettau, corrected_plots, "corrected_herwig_eec3_jetmettau_", detail, test);

//individual jetmettau hw 23
apply_correction(out_jetmettau_allvertex, corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", corr_detector_hw_23_jetmettau, "detector_herwig_23_JetMETTau_2010A_", corrected_hw_23_jetmettau, corrected_plots, "corrected_herwig_23_jetmettau_", detail, test);

//official jetmet
//apply_correction(out_jetmet_allvertex, corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", corrected_jetmet, corrected_plots, "corrected_jetmet_", detail, test);

//individual jetmet p6
apply_correction(out_jetmet_allvertex, corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", corrected_p6_z2_jetmet, corrected_plots, "corrected_p6_z2_jetmet_", detail, test);

//individual jetmet p8
apply_correction(out_jetmet_allvertex, corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", corrected_p8_4c_jetmet, corrected_plots, "corrected_p8_4c_jetmet_", detail, test);

//individual jetmet hw eec3
apply_correction(out_jetmet_allvertex, corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", corr_detector_hw_eec3_jetmet, "detector_herwig_eec3_JetMET_2010A_", corrected_hw_eec3_jetmet, corrected_plots, "corrected_herwig_eec3_jetmet_", detail, test);

//individual jetmet hw 23
apply_correction(out_jetmet_allvertex, corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", corr_detector_hw_23_jetmet, "detector_herwig_23_JetMET_2010A_", corrected_hw_23_jetmet, corrected_plots, "corrected_herwig_23_jetmet_", detail, test);

//apply_correction(out_jet_allvertex, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", corrected_jet, corrected_plots, "corrected_jet_", detail, test);

//apply_correction(out_p6_z2_det_allvertex_jet, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", corrected_p6_z2_jet, corrected_plots, "corrected_p6_z2_jet_", detail, test);

//apply_correction(out_p8_4c_det_allvertex_jet, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", corrected_p6_z2_jet, corrected_plots, "corrected_p8_4c_jetmet_", detail, test);

if (show_steps) { cout << "The correction factors were sucessfully applied!"<<endl; }
}


//plot uncorrected control distributions
if (plot_control_dist)
{
if (show_steps) { cout << "Ploting Control Distributions..."<<endl; }
create_directories(output_dir, control_plots, "plots");
if (show_steps) { cout << "Ploting MC on Generator Level..."<<endl; }
///plot_mc_gen_level(out_p8_4c_gen_allvertex, out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmet, out_p8_4c_gen_allvertex_jet, control_plots, detail, show_results);


if (show_steps) { cout << "Ploting Detector Level and Data..."<<endl; }

//control_plots(out_p6_z2_gen_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_det_allvertex, out_p8_4c_det_allvertex, out_jetmet_allvertex, control_plots, "control_rawvsjetmet_", detail, show_results);

//control_plots(out_p6_z2_gen_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_det_allvertex, out_p8_4c_det_allvertex, merged_data, control_plots, "control_", detail, show_results);

control_plots(out_hw_eec3_gen_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_eec3_det_allvertex_jetmettau, out_p8_4c_det_allvertex_jetmettau, out_jetmettau_allvertex, control_plots, "control_JetMETTau_2010A_", detail, show_results);

control_plots(out_hw_eec3_gen_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_eec3_det_allvertex_jetmet, out_p8_4c_det_allvertex_jetmet, out_jetmet_allvertex, control_plots, "control_JetMET_2010A_", detail, show_results);

//control_plots(out_p6_z2_gen_allvertex_jet, out_p8_4c_gen_allvertex_jet, out_p6_z2_det_allvertex_jet, out_p8_4c_det_allvertex_jet, out_jet_allvertex, control_plots, "control_Jet_2010B_", detail, show_results);
if (show_steps) { cout << "The Control Plots have been sucessfully generated!"<<endl; }
}


//creates the unfolding response matrix
if (create_unfolding_response)
{
if (show_steps) { cout << "Creating the RooUnfold Response MAtrixes..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, response_dir, "root");
if (show_steps) { cout << "MC without vertex reweighting..."<<endl; }
///create_unfolding_response(mc_p6_z2, out_p6_z2_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", "",  "", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", "",  "", detail, test);
if (show_steps) { cout << "MC with vertex reweighting as JetMETTau_2010A..."<<endl; }
//create_unfolding_response(mc_p6_z2, out_p6_z2_jetmettau_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmettau_v5,  "_v5", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_jetmettau_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmettau_v5,  "_v5", detail, test);
///create_unfolding_response(mc_hw_eec3, out_hw_eec3_jetmettau_unf_allvertex, lumi_hw_eec3, n_files_hw_eec3, "allvertex", vertex_weights_hw_eec3_jetmettau_v5,  "_v5", detail, test);
///create_unfolding_response(mc_hw_23, out_hw_23_jetmettau_unf_allvertex, lumi_hw_23, n_files_hw_23, "allvertex", vertex_weights_hw_23_jetmettau_v5,  "_v5", detail, test);
if (show_steps) { cout << "MC with vertex reweighting as JetMET_2010A..."<<endl; }
//create_unfolding_response(mc_p6_z2, out_p6_z2_jetmet_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmet_v5,  "_v5", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_jetmet_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmet_v5,  "_v5", detail, test);
///create_unfolding_response(mc_hw_eec3, out_hw_eec3_jetmet_unf_allvertex, lumi_hw_eec3, n_files_hw_eec3, "allvertex", vertex_weights_hw_eec3_jetmet_v5,  "_v5", detail, test);
///create_unfolding_response(mc_hw_23, out_hw_23_jetmet_unf_allvertex, lumi_hw_23, n_files_hw_23, "allvertex", vertex_weights_hw_23_jetmet_v5,  "_v5", detail, test);
if (show_steps) { cout << "MC with vertex reweighting as Jet_2010B..."<<endl; }
///create_unfolding_response(mc_p6_z2, out_p6_z2_jet_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jet_v5,  "_v5", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_jet_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jet_v5,  "_v5", detail, test);
if (show_steps) { cout << "The RooUnfold Response Matrixes were been sucessfully created!"<<endl; }
}


//checks the response matrix
if (check_response_matrix)
{
if (show_steps) { cout << "Check the Response Matrixes..."<<endl; }
create_directories(output_dir, check_response_dir, "plots");
if (show_steps) { cout << "Pythia 6 - Tune Z2*..."<<endl; }
//check_response_matrix(out_p6_z2_unf_allvertex, check_response_dir, "pythia6_z2_", detail, test);
//check_response_matrix(out_p6_z2_jetmettau_unf_allvertex, check_response_dir, "pythia6_z2_jetmettau_", detail, test);
//check_response_matrix(out_p6_z2_jetmet_unf_allvertex, check_response_dir, "pythia6_z2_jetmet_", detail, test);
check_response_matrix(out_p6_z2_jet_unf_allvertex, check_response_dir, "pythia6_z2_jet_", detail, test);
if (show_steps) { cout << "Pythia 8 - Tune 4C..."<<endl; }
check_response_matrix(out_p8_4c_unf_allvertex, check_response_dir, "pythia8_4c_", detail, test);
check_response_matrix(out_p8_4c_jetmettau_unf_allvertex, check_response_dir, "pythia8_4c_jetmettau_", detail, test);
check_response_matrix(out_p8_4c_jetmet_unf_allvertex, check_response_dir, "pythia8_4c_jetmet_", detail, test);
//check_response_matrix(out_p8_4c_jet_unf_allvertex, check_response_dir, "pythia8_4c_jet_", detail, test);
if (show_steps) { cout << "Herwig - Tune EEC3..."<<endl; }
check_response_matrix(out_hw_eec3_unf_allvertex, check_response_dir, "herwig_eec3_", detail, test);
check_response_matrix(out_hw_eec3_jetmettau_unf_allvertex, check_response_dir, "herwig_eec3_jetmettau_", detail, test);
check_response_matrix(out_hw_eec3_jetmet_unf_allvertex, check_response_dir, "herwig_eec3_jetmet_", detail, test);
if (show_steps) { cout << "Herwig - Tune 23..."<<endl; }
check_response_matrix(out_hw_23_unf_allvertex, check_response_dir, "herwig_23_", detail, test);
check_response_matrix(out_hw_23_jetmettau_unf_allvertex, check_response_dir, "herwig_23_jetmettau_", detail, test);
check_response_matrix(out_hw_23_jetmet_unf_allvertex, check_response_dir, "herwig_23_jetmet_", detail, test);
if (show_steps) { cout << "The Response Matrixes were been sucessfully checked!"<<endl; }
}


//unfolds the distribution
if (unfold)
{
if (show_steps) { cout << "Unfolding..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, unfolded_dir, "root");
create_directories(output_dir, unfolding_plots, "plots");
double tau[24];
for (int i=0; i < 24; i++)
{ tau[i] = 0.0; }
if (show_steps) { cout << "Unfolding without Vertex Reweight..."<<endl; }

//BinByBin - Data
//unfolding(out_p6_z2_unf_allvertex, data_unc, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_bbb_data, "BinByBin", 0.25, tau[0], unfolding_plots, "data_unfolded_bbb_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, data_unc, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_bbb_data, "BinByBin", 0.25, tau[1], unfolding_plots, "data_unfolded_bbb_pythia8_4c_", false, detail, test);
//unfolding(out_hw_eec3_unf_allvertex, data_unc, out_hw_eec3_gen_allvertex, out_hw_eec3_unfolded_bbb_data, "BinByBin", 0.25, tau[2], unfolding_plots, "data_unfolded_bbb_herwig_eec3_", false, detail, test);
//unfolding(out_hw_23_unf_allvertex, data_unc, out_hw_23_gen_allvertex, out_hw_23_unfolded_bbb_data, "BinByBin", 0.25, tau[3], unfolding_plots, "data_unfolded_bbb_herwig_23_", false, detail, test);

//BinByBin - Closure test
//unfolding(out_p6_z2_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_bbb_p6_z2, "BinByBin", 0.5, tau[4], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_bbb_p8_4c, "BinByBin", 0.5, tau[5], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia8_4c_", false, detail, test);
//unfolding(out_hw_eec3_unf_allvertex, out_hw_eec3_det_allvertex, out_hw_eec3_gen_allvertex, out_p6_z2_unfolded_bbb_hw_eec3, "BinByBin", 0.5, tau[6], unfolding_plots, "herwig_eec3_unfolded_bbb_herwig_eec3_", false, detail, test);
//unfolding(out_hw_23_unf_allvertex, out_hw_23_det_allvertex, out_hw_23_gen_allvertex, out_hw_23_unfolded_bbb_hw_23, "BinByBin", 0.5, tau[7], unfolding_plots, "herwig_23_unfolded_bbb_herwig_23_", false, detail, test);

//BinByBin - Cross-Closure test - Pythia 8 Matrix
//unfolding(out_p8_4c_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p8_4c_unfolded_bbb_p6_z2, "BinByBin", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia8_4c_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_hw_eec3_det_allvertex, out_hw_eec3_gen_allvertex, out_p8_4c_unfolded_bbb_hw_eec3, "BinByBin", 0.5, tau[9], unfolding_plots, "herwig_eec3_unfolded_bbb_pythia8_4c_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_hw_23_det_allvertex, out_hw_23_gen_allvertex, out_p8_4c_unfolded_bbb_hw_23, "BinByBin", 0.5, tau[10], unfolding_plots, "herwig_23_unfolded_bbb_pythia8_4c_", false, detail, test);

//BinByBin - Cross-Closure test - Pythia 6 Matrix
//unfolding(out_p6_z2_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_unfolded_bbb_p8_4c, "BinByBin", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia6_z2_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_hw_eec3_det_allvertex, out_hw_eec3_gen_allvertex, out_p6_z2_unfolded_bbb_hw_eec3, "BinByBin", 0.5, tau[12], unfolding_plots, "herwig_eec3_unfolded_bbb_pythia6_z2_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_hw_23_det_allvertex, out_hw_23_gen_allvertex, out_p6_z2_unfolded_bbb_hw_23, "BinByBin", 0.5, tau[13], unfolding_plots, "herwig_23_unfolded_bbb_pythia6_z2_", false, detail, test);

//BinByBin - Cross-Closure test - Herwig EEC3 Matrix
//unfolding(out_hw_eec3_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_hw_eec3_unfolded_bbb_p8_4c, "BinByBin", 0.5, tau[14], unfolding_plots, "pythia8_4c_unfolded_bbb_herwig_eec3_", false, detail, test);
//unfolding(out_hw_eec3_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_hw_eec3_unfolded_bbb_p6_z2, "BinByBin", 0.5, tau[15], unfolding_plots, "pythia6_z2_unfolded_bbb_herwig_eec3_", false, detail, test);
//unfolding(out_hw_eec3_unf_allvertex, out_hw_23_det_allvertex, out_hw_23_gen_allvertex, out_hw_eec3_unfolded_bbb_hw_23, "BinByBin", 0.5, tau[16], unfolding_plots, "herwig_23_unfolded_bbb_herwig_eec3_", false, detail, test);

//BinByBin - Cross-Closure test - Herwig 23 Matrix
//unfolding(out_hw_23_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_hw_23_unfolded_bbb_p8_4c, "BinByBin", 0.5, tau[17], unfolding_plots, "pythia8_4c_unfolded_bbb_herwig_23_", false, detail, test);
//unfolding(out_hw_23_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_hw_23_unfolded_bbb_p6_z2, "BinByBin", 0.5, tau[18], unfolding_plots, "pythia6_z2_unfolded_bbb_herwig_23_", false, detail, test);
//unfolding(out_hw_23_unf_allvertex, out_hw_eec3_det_allvertex, out_hw_eec3_gen_allvertex, out_hw_23_unfolded_bbb_hw_eec3, "BinByBin", 0.5, tau[19], unfolding_plots, "herwig_23_unfolded_bbb_herwig_23_", false, detail, test);


//SVD
//unfolding(out_p6_z2_unf_allvertex, data_unc, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_svd_data, "SVD", 0.5, tau[0], unfolding_plots, "data_unfolded_svd_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, data_unc, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_svd_data, "SVD", 0.5, tau[1], unfolding_plots, "data_unfolded_svd_pythia8_4c_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_svd_p6_z2, "SVD", 0.5, tau[2], unfolding_plots, "pythia6_z2_unfolded_svd_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_svd_p8_4c, "SVD", 0.5, tau[3], unfolding_plots, "pythia8_4c_unfolded_svd_pythia8_4c_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p8_4c_unfolded_svd_p6_z2, "SVD", 0.5, tau[4], unfolding_plots, "pythia6_z2_unfolded_svd_pythia8_4c_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_unfolded_svd_p8_4c, "SVD", 0.5, tau[5], unfolding_plots, "pythia8_4c_unfolded_svd_pythia6_z2_", false, detail, test);


//Bayes
//unfolding(out_p6_z2_unf_allvertex, data_unc, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_bayes_data, "Bayes", 0.5, tau[0], unfolding_plots, "data_unfolded_bayes_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, data_unc, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_bayes_data, "Bayes", 0.5, tau[1], unfolding_plots, "data_unfolded_bayes_pythia8_4c_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_bayes_p6_z2, "Bayes", 0.5, tau[2], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_bayes_p8_4c, "Bayes", 0.5, tau[3], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia8_4c_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p8_4c_unfolded_bayes_p6_z2, "Bayes", 0.5, tau[4], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia8_4c_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_unfolded_bayes_p8_4c, "Bayes", 0.5, tau[5], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia6_z2_", false, detail, test);


//TUnfold
//unfolding(out_p6_z2_unf_allvertex, data_unc, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_tunfold_data, "TUnfold", 0.5, tau[0], unfolding_plots, "data_unfolded_tunfold_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, data_unc, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_tunfold_data, "TUnfold", 0.5, tau[1], unfolding_plots, "data_unfolded_tunfold_pythia8_4c_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p6_z2_unfolded_tunfold_p6_z2, "TUnfold", 0.5, tau[2], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia6_z2_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p8_4c_unfolded_tunfold_p8_4c, "TUnfold", 0.5, tau[3], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia8_4c_", false, detail, test);
//unfolding(out_p8_4c_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_allvertex, out_p8_4c_unfolded_tunfold_p6_z2, "TUnfold", 0.5, tau[4], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia8_4c_", false, detail, test);
//unfolding(out_p6_z2_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_unfolded_tunfold_p8_4c, "TUnfold", 0.5, tau[5], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia6_z2_", false, detail, test);


if (show_steps) { cout << "Unfolding for JetMETTau_2010A..."<<endl; }

///BinByBin - Data
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p6_z2_jetmettau, out_p6_z2_unfolded_bbb_data_jetmettau, "BinByBin", 0.25, tau[6], unfolding_plots, "data_unfolded_bbb_pythia6_z2_jetmettau_", true, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p8_4c_jetmettau, out_p8_4c_unfolded_bbb_data_jetmettau, "BinByBin", 2.0, tau[7], unfolding_plots, "data_unfolded_bbb_pythia8_4c_jetmettau_", true, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_eec3_jetmettau, out_hw_eec3_unfolded_bbb_data_jetmettau, "BinByBin", 2.0, tau[7], unfolding_plots, "data_unfolded_bbb_herwig_eec3_jetmettau_", true, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_23_jetmettau, out_hw_23_unfolded_bbb_data_jetmettau, "BinByBin", 2.0, tau[7], unfolding_plots, "data_unfolded_bbb_herwig_23_jetmettau_", true, detail, test);


//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p6_z2_unfolded_bbb_p6_z2_jetmettau, "BinByBin", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia6_z2_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_unfolded_bbb_p8_4c_jetmettau, "BinByBin", 2.0, tau[9], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia8_4c_jetmettau_", false, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_eec3_unfolded_bbb_hw_eec3_jetmettau, "BinByBin", 2.0, tau[9], unfolding_plots, "herwig_eec3_unfolded_bbb_herwig_eec3_jetmettau_", false, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_23_unfolded_bbb_hw_23_jetmettau, "BinByBin", 2.0, tau[9], unfolding_plots, "herwig_23_unfolded_bbb_herwig_23_jetmettau_", false, detail, test);


//unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p8_4c_unfolded_bbb_p6_z2_jetmettau, "BinByBin", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia8_4c_jetmettau_", false, detail, test);
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p6_z2_unfolded_bbb_p8_4c_jetmettau, "BinByBin", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia6_z2_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_p8_4c_unfolded_bbb_hw_eec3_jetmettau, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_eec3_unfolded_bbb_pythia8_4c_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_p8_4c_unfolded_bbb_hw_23_jetmettau, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_23_unfolded_bbb_pythia8_4c_jetmettau_", false, detail, test);

///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_eec3_unfolded_bbb_p8_4c_jetmettau, "BinByBin", 2.0, tau[10], unfolding_plots, "pythia8_4c_unfolded_bbb_herwig_eec3_jetmettau_", false, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_eec3_unfolded_bbb_hw_23_jetmettau, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_23_unfolded_bbb_herwig_eec3_jetmettau_", false, detail, test);

///unfolding(out_hw_23_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_23_unfolded_bbb_p8_4c_jetmettau, "BinByBin", 2.0, tau[10], unfolding_plots, "pythia8_4c_unfolded_bbb_herwig_23_jetmettau_", false, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_23_unfolded_bbb_hw_eec3_jetmettau, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_eec3_unfolded_bbb_herwig_23_jetmettau_", false, detail, test);


//SVD
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p6_z2_jetmettau, out_p6_z2_unfolded_svd_data_jetmettau, "SVD", 0.25, tau[6], unfolding_plots, "data_unfolded_svd_pythia6_z2_jetmettau_", true, detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p8_4c_jetmettau, out_p8_4c_unfolded_svd_data_jetmettau, "SVD", 1.25, tau[7], unfolding_plots, "data_unfolded_svd_pythia8_4c_jetmettau_", true, detail, test);
unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_eec3_jetmettau, out_hw_eec3_unfolded_svd_data_jetmettau, "SVD", 1.25, tau[7], unfolding_plots, "data_unfolded_svd_herwig_eec3_jetmettau_", true, detail, test);
unfolding(out_hw_23_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_23_jetmettau, out_hw_23_unfolded_svd_data_jetmettau, "SVD", 1.25, tau[7], unfolding_plots, "data_unfolded_svd_herwig_23_jetmettau_", true, detail, test);


//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p6_z2_unfolded_svd_p6_z2_jetmettau, "SVD", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_svd_pythia6_z2_jetmettau_", false, detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_unfolded_svd_p8_4c_jetmettau, "SVD", 1.25, tau[9], unfolding_plots, "pythia8_4c_unfolded_svd_pythia8_4c_jetmettau_", false, detail, test);
unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_eec3_unfolded_svd_hw_eec3_jetmettau, "SVD", 1.25, tau[9], unfolding_plots, "herwig_eec3_unfolded_svd_herwig_eec3_jetmettau_", false, detail, test);
unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_23_unfolded_svd_hw_23_jetmettau, "SVD", 1.25, tau[9], unfolding_plots, "herwig_23_unfolded_svd_herwig_23_jetmettau_", false, detail, test);


//unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p8_4c_unfolded_svd_p6_z2_jetmettau, "SVD", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_svd_pythia8_4c_jetmettau_", false, detail, test);
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p6_z2_unfolded_svd_p8_4c_jetmettau, "SVD", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_svd_pythia6_z2_jetmettau_", false, detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_p8_4c_unfolded_svd_hw_eec3_jetmettau, "SVD", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_svd_pythia8_4c_jetmettau_", false, detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_p8_4c_unfolded_svd_hw_23_jetmettau, "SVD", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_svd_pythia8_4c_jetmettau_", false, detail, test);

unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_eec3_unfolded_svd_p8_4c_jetmettau, "SVD", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_svd_herwig_eec3_jetmettau_", false, detail, test);
unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_eec3_unfolded_svd_hw_23_jetmettau, "SVD", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_svd_herwig_eec3_jetmettau_", false, detail, test);

unfolding(out_hw_23_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_23_unfolded_svd_p8_4c_jetmettau, "SVD", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_svd_herwig_23_jetmettau_", false, detail, test);
unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_23_unfolded_svd_hw_eec3_jetmettau, "SVD", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_svd_herwig_23_jetmettau_", false, detail, test);


//Bayes
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p6_z2_jetmettau, out_p6_z2_unfolded_bayes_data_jetmettau, "Bayes", 0.25, tau[6], unfolding_plots, "data_unfolded_bayes_pythia6_z2_jetmettau_", true, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p8_4c_jetmettau, out_p8_4c_unfolded_bayes_data_jetmettau, "Bayes", 1.25, tau[7], unfolding_plots, "data_unfolded_bayes_pythia8_4c_jetmettau_", true, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_eec3_jetmettau, out_hw_eec3_unfolded_bayes_data_jetmettau, "Bayes", 1.25, tau[7], unfolding_plots, "data_unfolded_bayes_herwig_eec3_jetmettau_", true, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_23_jetmettau, out_hw_23_unfolded_bayes_data_jetmettau, "Bayes", 1.25, tau[7], unfolding_plots, "data_unfolded_bayes_herwig_23_jetmettau_", true, detail, test);


//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p6_z2_unfolded_bayes_p6_z2_jetmettau, "Bayes", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia6_z2_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_unfolded_bayes_p8_4c_jetmettau, "Bayes", 1.25, tau[9], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia8_4c_jetmettau_", false, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_eec3_unfolded_bayes_hw_eec3_jetmettau, "Bayes", 1.25, tau[9], unfolding_plots, "herwig_eec3_unfolded_bayes_herwig_eec3_jetmettau_", false, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_23_unfolded_bayes_hw_23_jetmettau, "Bayes", 1.25, tau[9], unfolding_plots, "herwig_23_unfolded_bayes_herwig_23_jetmettau_", false, detail, test);


//unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p8_4c_unfolded_bayes_p6_z2_jetmettau, "Bayes", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia8_4c_jetmettau_", false, detail, test);
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p6_z2_unfolded_bayes_p8_4c_jetmettau, "Bayes", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia6_z2_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_p8_4c_unfolded_bayes_hw_eec3_jetmettau, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_bayes_pythia8_4c_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_p8_4c_unfolded_bayes_hw_23_jetmettau, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_bayes_pythia8_4c_jetmettau_", false, detail, test);

///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_eec3_unfolded_bayes_p8_4c_jetmettau, "Bayes", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_bayes_herwig_eec3_jetmettau_", false, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_eec3_unfolded_bayes_hw_23_jetmettau, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_bayes_herwig_eec3_jetmettau_", false, detail, test);

///unfolding(out_hw_23_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_23_unfolded_bayes_p8_4c_jetmettau, "Bayes", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_bayes_herwig_23_jetmettau_", false, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_23_unfolded_bayes_hw_eec3_jetmettau, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_bayes_herwig_23_jetmettau_", false, detail, test);


//TUnfold
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p6_z2_jetmettau, out_p6_z2_unfolded_tunfold_data_jetmettau, "TUnfold", 0.25, tau[6], unfolding_plots, "data_unfolded_tunfold_pythia6_z2_jetmettau_", true, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_p8_4c_jetmettau, out_p8_4c_unfolded_tunfold_data_jetmettau, "TUnfold", 1.25, tau[7], unfolding_plots, "data_unfolded_tunfold_pythia8_4c_jetmettau_", true, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_eec3_jetmettau, out_hw_eec3_unfolded_tunfold_data_jetmettau, "TUnfold", 1.25, tau[7], unfolding_plots, "data_unfolded_tunfold_herwig_eec3_jetmettau_", true, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_jetmettau_allvertex, corrected_hw_23_jetmettau, out_hw_23_unfolded_tunfold_data_jetmettau, "TUnfold", 1.25, tau[7], unfolding_plots, "data_unfolded_tunfold_herwig_23_jetmettau_", true, detail, test);


//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p6_z2_unfolded_tunfold_p6_z2_jetmettau, "TUnfold", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia6_z2_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_unfolded_tunfold_p8_4c_jetmettau, "TUnfold", 1.25, tau[9], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia8_4c_jetmettau_", false, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_eec3_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold", 1.25, tau[9], unfolding_plots, "herwig_eec3_unfolded_tunfold_herwig_eec3_jetmettau_", false, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_23_unfolded_tunfold_hw_23_jetmettau, "TUnfold", 1.25, tau[9], unfolding_plots, "herwig_23_unfolded_tunfold_herwig_23_jetmettau_", false, detail, test);


//unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_allvertex_jetmettau, out_p8_4c_unfolded_tunfold_p6_z2_jetmettau, "TUnfold", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia8_4c_jetmettau_", false, detail, test);
//unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p6_z2_unfolded_tunfold_p8_4c_jetmettau, "TUnfold", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia6_z2_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_p8_4c_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_tunfold_pythia8_4c_jetmettau_", false, detail, test);
///unfolding(out_p8_4c_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_p8_4c_unfolded_tunfold_hw_23_jetmettau, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_tunfold_pythia8_4c_jetmettau_", false, detail, test);

///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_eec3_unfolded_tunfold_p8_4c_jetmettau, "TUnfold", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_tunfold_herwig_eec3_jetmettau_", false, detail, test);
///unfolding(out_hw_eec3_jetmettau_unf_allvertex, out_hw_23_det_allvertex_jetmettau, out_hw_23_gen_allvertex_jetmettau, out_hw_eec3_unfolded_tunfold_hw_23_jetmettau, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_tunfold_herwig_eec3_jetmettau_", false, detail, test);

///unfolding(out_hw_23_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_hw_23_unfolded_tunfold_p8_4c_jetmettau, "TUnfold", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_tunfold_herwig_23_jetmettau_", false, detail, test);
///unfolding(out_hw_23_jetmettau_unf_allvertex, out_hw_eec3_det_allvertex_jetmettau, out_hw_eec3_gen_allvertex_jetmettau, out_hw_23_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_tunfold_herwig_23_jetmettau_", false, detail, test);



if (show_steps) { cout << "Unfolding for JetMET_2010A..."<<endl; }

///BinByBin - Data
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p6_z2_jetmet, out_p6_z2_unfolded_bbb_data_jetmet, "BinByBin", 0.25, tau[6], unfolding_plots, "data_unfolded_bbb_pythia6_z2_jetmet_", true, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p8_4c_jetmet, out_p8_4c_unfolded_bbb_data_jetmet, "BinByBin", 2.0, tau[7], unfolding_plots, "data_unfolded_bbb_pythia8_4c_jetmet_", true, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_eec3_jetmet, out_hw_eec3_unfolded_bbb_data_jetmet, "BinByBin", 2.0, tau[7], unfolding_plots, "data_unfolded_bbb_herwig_eec3_jetmet_", true, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_23_jetmet, out_hw_23_unfolded_bbb_data_jetmet, "BinByBin", 2.0, tau[7], unfolding_plots, "data_unfolded_bbb_herwig_23_jetmet_", true, detail, test);


//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p6_z2_unfolded_bbb_p6_z2_jetmet, "BinByBin", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia6_z2_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p8_4c_unfolded_bbb_p8_4c_jetmet, "BinByBin", 2.0, tau[9], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia8_4c_jetmet_", false, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_eec3_unfolded_bbb_hw_eec3_jetmet, "BinByBin", 2.0, tau[9], unfolding_plots, "herwig_eec3_unfolded_bbb_herwig_eec3_jetmet_", false, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_23_unfolded_bbb_hw_23_jetmet, "BinByBin", 2.0, tau[9], unfolding_plots, "herwig_23_unfolded_bbb_herwig_23_jetmet_", false, detail, test);


//unfolding(out_p8_4c_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p8_4c_unfolded_bbb_p6_z2_jetmet, "BinByBin", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia8_4c_jetmet_", false, detail, test);
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p6_z2_unfolded_bbb_p8_4c_jetmet, "BinByBin", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia6_z2_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_p8_4c_unfolded_bbb_hw_eec3_jetmet, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_eec3_unfolded_bbb_pythia8_4c_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_p8_4c_unfolded_bbb_hw_23_jetmet, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_23_unfolded_bbb_pythia8_4c_jetmet_", false, detail, test);

///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_eec3_unfolded_bbb_p8_4c_jetmet, "BinByBin", 2.0, tau[10], unfolding_plots, "pythia8_4c_unfolded_bbb_herwig_eec3_jetmet_", false, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_eec3_unfolded_bbb_hw_23_jetmet, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_23_unfolded_bbb_herwig_eec3_jetmet_", false, detail, test);

///unfolding(out_hw_23_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_23_unfolded_bbb_p8_4c_jetmet, "BinByBin", 2.0, tau[10], unfolding_plots, "pythia8_4c_unfolded_bbb_herwig_23_jetmet_", false, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_23_unfolded_bbb_hw_eec3_jetmet, "BinByBin", 2.0, tau[10], unfolding_plots, "herwig_eec3_unfolded_bbb_herwig_23_jetmet_", false, detail, test);


//SVD
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p6_z2_jetmet, out_p6_z2_unfolded_svd_data_jetmet, "SVD", 0.25, tau[6], unfolding_plots, "data_unfolded_svd_pythia6_z2_jetmet_", true, detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p8_4c_jetmet, out_p8_4c_unfolded_svd_data_jetmet, "SVD", 1.25, tau[7], unfolding_plots, "data_unfolded_svd_pythia8_4c_jetmet_", true, detail, test);
unfolding(out_hw_eec3_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_eec3_jetmet, out_hw_eec3_unfolded_svd_data_jetmet, "SVD", 1.25, tau[7], unfolding_plots, "data_unfolded_svd_herwig_eec3_jetmet_", true, detail, test);
unfolding(out_hw_23_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_23_jetmet, out_hw_23_unfolded_svd_data_jetmet, "SVD", 1.25, tau[7], unfolding_plots, "data_unfolded_svd_herwig_23_jetmet_", true, detail, test);


//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p6_z2_unfolded_svd_p6_z2_jetmet, "SVD", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_svd_pythia6_z2_jetmet_", false, detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p8_4c_unfolded_svd_p8_4c_jetmet, "SVD", 1.25, tau[9], unfolding_plots, "pythia8_4c_unfolded_svd_pythia8_4c_jetmet_", false, detail, test);
unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_eec3_unfolded_svd_hw_eec3_jetmet, "SVD", 1.25, tau[9], unfolding_plots, "herwig_eec3_unfolded_svd_herwig_eec3_jetmet_", false, detail, test);
unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_23_unfolded_svd_hw_23_jetmet, "SVD", 1.25, tau[9], unfolding_plots, "herwig_23_unfolded_svd_herwig_23_jetmet_", false, detail, test);


//unfolding(out_p8_4c_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p8_4c_unfolded_svd_p6_z2_jetmet, "SVD", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_svd_pythia8_4c_jetmet_", false, detail, test);
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p6_z2_unfolded_svd_p8_4c_jetmet, "SVD", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_svd_pythia6_z2_jetmet_", false, detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_p8_4c_unfolded_svd_hw_eec3_jetmet, "SVD", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_svd_pythia8_4c_jetmet_", false, detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_p8_4c_unfolded_svd_hw_23_jetmet, "SVD", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_svd_pythia8_4c_jetmet_", false, detail, test);

unfolding(out_hw_eec3_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_eec3_unfolded_svd_p8_4c_jetmet, "SVD", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_svd_herwig_eec3_jetmet_", false, detail, test);
unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_eec3_unfolded_svd_hw_23_jetmet, "SVD", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_svd_herwig_eec3_jetmet_", false, detail, test);

unfolding(out_hw_23_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_23_unfolded_svd_p8_4c_jetmet, "SVD", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_svd_herwig_23_jetmet_", false, detail, test);
unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_23_unfolded_svd_hw_eec3_jetmet, "SVD", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_svd_herwig_23_jetmet_", false, detail, test);


//Bayes
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p6_z2_jetmet, out_p6_z2_unfolded_bayes_data_jetmet, "Bayes", 0.25, tau[6], unfolding_plots, "data_unfolded_bayes_pythia6_z2_jetmet_", true, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p8_4c_jetmet, out_p8_4c_unfolded_bayes_data_jetmet, "Bayes", 1.25, tau[7], unfolding_plots, "data_unfolded_bayes_pythia8_4c_jetmet_", true, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_eec3_jetmet, out_hw_eec3_unfolded_bayes_data_jetmet, "Bayes", 1.25, tau[7], unfolding_plots, "data_unfolded_bayes_herwig_eec3_jetmet_", true, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_23_jetmet, out_hw_23_unfolded_bayes_data_jetmet, "Bayes", 1.25, tau[7], unfolding_plots, "data_unfolded_bayes_herwig_23_jetmet_", true, detail, test);


//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p6_z2_unfolded_bayes_p6_z2_jetmet, "Bayes", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia6_z2_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p8_4c_unfolded_bayes_p8_4c_jetmet, "Bayes", 1.25, tau[9], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia8_4c_jetmet_", false, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_eec3_unfolded_bayes_hw_eec3_jetmet, "Bayes", 1.25, tau[9], unfolding_plots, "herwig_eec3_unfolded_bayes_herwig_eec3_jetmet_", false, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_23_unfolded_bayes_hw_23_jetmet, "Bayes", 1.25, tau[9], unfolding_plots, "herwig_23_unfolded_bayes_herwig_23_jetmet_", false, detail, test);


//unfolding(out_p8_4c_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p8_4c_unfolded_bayes_p6_z2_jetmet, "Bayes", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia8_4c_jetmet_", false, detail, test);
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p6_z2_unfolded_bayes_p8_4c_jetmet, "Bayes", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia6_z2_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_p8_4c_unfolded_bayes_hw_eec3_jetmet, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_bayes_pythia8_4c_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_p8_4c_unfolded_bayes_hw_23_jetmet, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_bayes_pythia8_4c_jetmet_", false, detail, test);

///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_eec3_unfolded_bayes_p8_4c_jetmet, "Bayes", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_bayes_herwig_eec3_jetmet_", false, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_eec3_unfolded_bayes_hw_23_jetmet, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_bayes_herwig_eec3_jetmet_", false, detail, test);

///unfolding(out_hw_23_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_23_unfolded_bayes_p8_4c_jetmet, "Bayes", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_bayes_herwig_23_jetmet_", false, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_23_unfolded_bayes_hw_eec3_jetmet, "Bayes", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_bayes_herwig_23_jetmet_", false, detail, test);


//TUnfold
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p6_z2_jetmet, out_p6_z2_unfolded_tunfold_data_jetmet, "TUnfold", 0.25, tau[6], unfolding_plots, "data_unfolded_tunfold_pythia6_z2_jetmet_", true, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_p8_4c_jetmet, out_p8_4c_unfolded_tunfold_data_jetmet, "TUnfold", 1.25, tau[7], unfolding_plots, "data_unfolded_tunfold_pythia8_4c_jetmet_", true, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_eec3_jetmet, out_hw_eec3_unfolded_tunfold_data_jetmet, "TUnfold", 1.25, tau[7], unfolding_plots, "data_unfolded_tunfold_herwig_eec3_jetmet_", true, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_jetmet_allvertex, corrected_hw_23_jetmet, out_hw_23_unfolded_tunfold_data_jetmet, "TUnfold", 1.25, tau[7], unfolding_plots, "data_unfolded_tunfold_herwig_23_jetmet_", true, detail, test);


//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p6_z2_unfolded_tunfold_p6_z2_jetmet, "TUnfold", 0.5, tau[8], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia6_z2_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p8_4c_unfolded_tunfold_p8_4c_jetmet, "TUnfold", 1.25, tau[9], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia8_4c_jetmet_", false, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_eec3_unfolded_tunfold_hw_eec3_jetmet, "TUnfold", 1.25, tau[9], unfolding_plots, "herwig_eec3_unfolded_tunfold_herwig_eec3_jetmet_", false, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_23_unfolded_tunfold_hw_23_jetmet, "TUnfold", 1.25, tau[9], unfolding_plots, "herwig_23_unfolded_tunfold_herwig_23_jetmet_", false, detail, test);


//unfolding(out_p8_4c_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_allvertex_jetmet, out_p8_4c_unfolded_tunfold_p6_z2_jetmet, "TUnfold", 0.5, tau[10], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia8_4c_jetmet_", false, detail, test);
//unfolding(out_p6_z2_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p6_z2_unfolded_tunfold_p8_4c_jetmet, "TUnfold", 0.5, tau[11], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia6_z2_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_p8_4c_unfolded_tunfold_hw_eec3_jetmet, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_tunfold_pythia8_4c_jetmet_", false, detail, test);
///unfolding(out_p8_4c_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_p8_4c_unfolded_tunfold_hw_23_jetmet, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_tunfold_pythia8_4c_jetmet_", false, detail, test);

///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_eec3_unfolded_tunfold_p8_4c_jetmet, "TUnfold", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_tunfold_herwig_eec3_jetmet_", false, detail, test);
///unfolding(out_hw_eec3_jetmet_unf_allvertex, out_hw_23_det_allvertex_jetmet, out_hw_23_gen_allvertex_jetmet, out_hw_eec3_unfolded_tunfold_hw_23_jetmet, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_23_unfolded_tunfold_herwig_eec3_jetmet_", false, detail, test);

///unfolding(out_hw_23_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_hw_23_unfolded_tunfold_p8_4c_jetmet, "TUnfold", 1.25, tau[10], unfolding_plots, "pythia8_4c_unfolded_tunfold_herwig_23_jetmet_", false, detail, test);
///unfolding(out_hw_23_jetmet_unf_allvertex, out_hw_eec3_det_allvertex_jetmet, out_hw_eec3_gen_allvertex_jetmet, out_hw_23_unfolded_tunfold_hw_eec3_jetmet, "TUnfold", 1.25, tau[10], unfolding_plots, "herwig_eec3_unfolded_tunfold_herwig_23_jetmet_", false, detail, test);

if (show_steps) { cout << "Unfolding for Jet_2010B..."<<endl; }

//BinByBin
//unfolding(out_p6_z2_jet_unf_allvertex, out_jet_allvertex, corrected_p6_z2_jet, out_p6_z2_unfolded_bbb_data_jet, "BinByBin", 0.25, tau[18], unfolding_plots, "data_unfolded_bbb_pythia6_z2_jet_", true, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_jet_allvertex, corrected_p8_4c_jet, out_p8_4c_unfolded_bbb_data_jet, "BinByBin", 0.25, tau[19], unfolding_plots, "data_unfolded_bbb_pythia8_4c_jet_", true, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p6_z2_unfolded_bbb_p6_z2_jet, "BinByBin", 0.5, tau[20], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia6_z2_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p8_4c_unfolded_bbb_p8_4c_jet, "BinByBin", 0.5, tau[21], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p8_4c_unfolded_bbb_p6_z2_jet, "BinByBin", 0.5, tau[22], unfolding_plots, "pythia6_z2_unfolded_bbb_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p6_z2_unfolded_bbb_p8_4c_jet, "BinByBin", 0.5, tau[23], unfolding_plots, "pythia8_4c_unfolded_bbb_pythia6_z2_jet_", false, detail, test);


//SVD
//unfolding(out_p6_z2_jet_unf_allvertex, out_jet_allvertex, corrected_p6_z2_jet, out_p6_z2_unfolded_svd_data_jet, "SVD", 0.5, tau[18], unfolding_plots, "data_unfolded_svd_pythia6_z2_jet_", true, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_jet_allvertex, corrected_p8_4c_jet, out_p8_4c_unfolded_svd_data_jet, "SVD", 0.5, tau[19], unfolding_plots, "data_unfolded_svd_pythia8_4c_jet_", true, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p6_z2_unfolded_svd_p6_z2_jet, "SVD", 0.5, tau[20], unfolding_plots, "pythia6_z2_unfolded_svd_pythia6_z2_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p8_4c_unfolded_svd_p8_4c_jet, "SVD", 0.5, tau[21], unfolding_plots, "pythia8_4c_unfolded_svd_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p8_4c_unfolded_svd_p6_z2_jet, "SVD", 0.5, tau[22], unfolding_plots, "pythia6_z2_unfolded_svd_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p6_z2_unfolded_svd_p8_4c_jet, "SVD", 0.5, tau[23], unfolding_plots, "pythia8_4c_unfolded_svd_pythia6_z2_jet_", false, detail, test);


//Bayes
//unfolding(out_p6_z2_jet_unf_allvertex, out_jet_allvertex, corrected_p6_z2_jet, out_p6_z2_unfolded_bayes_data_jet, "Bayes", 0.5, tau[18], unfolding_plots, "data_unfolded_bayes_pythia6_z2_jet_", true, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_jet_allvertex, corrected_p8_4c_jet, out_p8_4c_unfolded_bayes_data_jet, "Bayes", 0.5, tau[19], unfolding_plots, "data_unfolded_bayes_pythia8_4c_jet_", true, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p6_z2_unfolded_bayes_p6_z2_jet, "Bayes", 0.5, tau[20], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia6_z2_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p8_4c_unfolded_bayes_p8_4c_jet, "Bayes", 0.5, tau[21], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p8_4c_unfolded_bayes_p6_z2_jet, "Bayes", 0.5, tau[22], unfolding_plots, "pythia6_z2_unfolded_bayes_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p6_z2_unfolded_bayes_p8_4c_jet, "Bayes", 0.5, tau[23], unfolding_plots, "pythia8_4c_unfolded_bayes_pythia6_z2_jet_", false, detail, test);


//TUnfold
//unfolding(out_p6_z2_jet_unf_allvertex, out_jet_allvertex, corrected_p6_z2_jet, out_p6_z2_unfolded_tunfold_data_jet, "TUnfold", 0.5, tau[18], unfolding_plots, "data_unfolded_tunfold_pythia6_z2_jet_", true, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_jet_allvertex, corrected_p8_4c_jet, out_p8_4c_unfolded_tunfold_data_jet, "TUnfold", 0.5, tau[19], unfolding_plots, "data_unfolded_tunfold_pythia8_4c_jet_", true, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p6_z2_unfolded_tunfold_p6_z2_jet, "TUnfold", 0.5, tau[20], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia6_z2_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p8_4c_unfolded_tunfold_p8_4c_jet, "TUnfold", 0.5, tau[21], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p8_4c_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p8_4c_unfolded_tunfold_p6_z2_jet, "TUnfold", 0.5, tau[22], unfolding_plots, "pythia6_z2_unfolded_tunfold_pythia8_4c_jet_", false, detail, test);
//unfolding(out_p6_z2_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p6_z2_unfolded_tunfold_p8_4c_jet, "TUnfold", 0.5, tau[23], unfolding_plots, "pythia8_4c_unfolded_tunfold_pythia6_z2_jet_", false, detail, test);


if (show_steps) { cout << "Show Tau Values..."<<endl; }
cout << "Tau Values No Reweight:     " << tau[0] << " " << tau[1] << " " << tau[2] << " " << tau[3] << " " << tau[4] << " " << tau[5] << endl;
cout << "Tau Values JetMETTau_2010A: " << tau[6] << " " << tau[7] << " " << tau[8] << " " << tau[9] << " " << tau[10] << " " << tau[11] << endl;
cout << "Tau Values JetMET_2010A:    " << tau[12] << " " << tau[13] << " " << tau[14] << " " << tau[15] << " " << tau[16] << " " << tau[17] << endl;
cout << "Tau Values Jet_2010B:       " << tau[18] << " " << tau[19] << " " << tau[20] << " " << tau[21] << " " << tau[22] << " " << tau[23] << endl;


if (show_steps) { cout << "The Unfolding has sucessful!"<<endl; }
}



//compare the results for the various unfolding methods
if (compare_unfolding_results)
{
if (show_steps) { cout << "Compare the Unfolding Results..."<<endl; }
create_directories(output_dir, compare_unfolding_plots, "plots");

//compare_unfolding_results(out_jetmettau_allvertex, corrected_p6_z2_jetmettau, out_p6_z2_unfolded_bbb_data_jetmettau, out_p6_z2_unfolded_tunfold_data_jetmettau, out_p6_z2_unfolded_svd_data_jetmettau, out_p6_z2_unfolded_bayes_data_jetmettau, compare_unfolding_plots, "data_with_p6_z2_jetmettau_", false, detail);
//compare_unfolding_results(out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_jetmettau_norm, out_p6_z2_unfolded_bbb_p6_z2_jetmettau, out_p6_z2_unfolded_tunfold_p6_z2_jetmettau, out_p6_z2_unfolded_svd_p6_z2_jetmettau, out_p6_z2_unfolded_bayes_p6_z2_jetmettau, compare_unfolding_plots, "p6_z2_with_p6_z2_jetmettau_", true, detail);
//compare_unfolding_results(out_p6_z2_det_allvertex_jetmettau, out_p8_4c_gen_jetmettau_norm, out_p6_z2_unfolded_bbb_p8_4c_jetmettau, out_p6_z2_unfolded_tunfold_p8_4c_jetmettau, out_p6_z2_unfolded_svd_p8_4c_jetmettau, out_p6_z2_unfolded_bayes_p8_4c_jetmettau, compare_unfolding_plots, "p8_4c_with_p6_z2_jetmettau_", true, detail);

//compare_unfolding_results(out_jetmet_allvertex, corrected_p6_z2_jetmet, out_p6_z2_unfolded_bbb_data_jetmet, out_p6_z2_unfolded_tunfold_data_jetmet, out_p6_z2_unfolded_svd_data_jetmet, out_p6_z2_unfolded_bayes_data_jetmet, compare_unfolding_plots, "data_with_p6_z2_jetmet_", false, detail);
//compare_unfolding_results(out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_jetmet_norm, out_p6_z2_unfolded_bbb_p6_z2_jetmet, out_p6_z2_unfolded_tunfold_p6_z2_jetmet, out_p6_z2_unfolded_svd_p6_z2_jetmet, out_p6_z2_unfolded_bayes_p6_z2_jetmet, compare_unfolding_plots, "p6_z2_with_p6_z2_jetmet_", true, detail);
//compare_unfolding_results(out_p6_z2_det_allvertex_jetmet, out_p8_4c_gen_jetmet_norm, out_p6_z2_unfolded_bbb_p8_4c_jetmet, out_p6_z2_unfolded_tunfold_p8_4c_jetmet, out_p6_z2_unfolded_svd_p8_4c_jetmet, out_p6_z2_unfolded_bayes_p8_4c_jetmet, compare_unfolding_plots, "p8_4c_with_p6_z2_jetmet_", true, detail);

//pythia 8 - jetmettau
compare_unfolding_results(out_jetmettau_allvertex, "Uncorrected Data", corrected_p8_4c_jetmettau, "Standard Bin-by-Bin", out_p8_4c_unfolded_bbb_data_jetmettau, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_data_jetmettau, "TUnfold", out_p8_4c_unfolded_svd_data_jetmettau, "SVD", out_p8_4c_unfolded_bayes_data_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-By-Bin", "SVD/Bin-By-Bin", "Bayes/Bin-By-Bin", compare_unfolding_plots, "data_with_p8_4c_jetmettau_", false, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmettau, "Detector Level", out_hw_eec3_gen_allvertex_jetmettau, "Generator Level", out_p8_4c_unfolded_bbb_hw_eec3_jetmettau, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold", out_p8_4c_unfolded_svd_hw_eec3_jetmettau, "SVD", out_p8_4c_unfolded_bayes_hw_eec3_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_eec3_with_p8_4c_jetmettau_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmettau, "Detector Level", out_hw_23_gen_allvertex_jetmettau, "Generator Level", out_p8_4c_unfolded_bbb_hw_23_jetmettau, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_hw_23_jetmettau, "TUnfold", out_p8_4c_unfolded_svd_hw_23_jetmettau, "SVD", out_p8_4c_unfolded_bayes_hw_23_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True",  compare_unfolding_plots, "hw_23_with_p8_4c_jetmettau_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmettau, "Detector Level", out_p8_4c_gen_allvertex_jetmettau, "Generator Level", out_p8_4c_unfolded_bbb_p8_4c_jetmettau, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_p8_4c_jetmettau, "TUnfold", out_p8_4c_unfolded_svd_p8_4c_jetmettau, "SVD", out_p8_4c_unfolded_bayes_p8_4c_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True",  compare_unfolding_plots, "p8_4c_with_p8_4c_jetmettau_", true, detail);

//method test
compare_unfolding_results(out_p8_4c_det_allvertex_jetmettau, "Detector Level", out_p8_4c_gen_allvertex_jetmettau, "Generator Level", out_p8_4c_unfolded_bayes1_p8_4c_jetmettau, "Bayes: iterations = 1", out_p8_4c_unfolded_bayes2_p8_4c_jetmettau, "Bayes: iterations = 2", out_p8_4c_unfolded_bayes3_p8_4c_jetmettau, "Bayes: iterations = 3", out_p8_4c_unfolded_bayes4_p8_4c_jetmettau, "Bayes: iterations = 4", "Bayes: iterations = 1 /True", "Bayes: iterations = 2 /True", "Bayes: iterations = 3 /True", "Bayes: iterations = 4 /True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmettau_bayes_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmettau, "Detector Level", out_p8_4c_gen_allvertex_jetmettau, "Generator Level", out_p8_4c_unfolded_svd1_p8_4c_jetmettau, "SVD: regularization = 1", out_p8_4c_unfolded_svd2_p8_4c_jetmettau, "SVD: regularization = 2", out_p8_4c_unfolded_svd3_p8_4c_jetmettau, "SVD: regularization = 3", out_p8_4c_unfolded_svd4_p8_4c_jetmettau, "SVD: regularization = 4", "SVD: regularization = 1 /True", "SVD: regularization = 2 /True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmettau_svd_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmettau, "Detector Level", out_p8_4c_gen_allvertex_jetmettau, "Generator Level", out_p8_4c_unfolded_tunfold_p8_4c_jetmettau, "TUnfold: Default", out_p8_4c_unfolded_tunfold2_p8_4c_jetmettau, "TUnfold: Size Regularization", out_p8_4c_unfolded_tunfold3_p8_4c_jetmettau, "TUnfold: Derivative Regularization", out_p8_4c_unfolded_tunfold4_p8_4c_jetmettau, "TUnfold: Curvature Regularization", "TUnfold: Default/True", "TUnfold: Size Regularization/True", "TUnfold: Derivative Regularization/True", "TUnfold: Curvature Regularization/True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmettau_tunfold_", true, detail);


//herwig eec3 - jetmettau
compare_unfolding_results(out_jetmettau_allvertex, "Uncorrected Data", corrected_hw_eec3_jetmettau, "Standard Bin-by-Bin", out_hw_eec3_unfolded_bbb_data_jetmettau, "RooUnfold Bin-by-Bin",  out_hw_eec3_unfolded_tunfold_data_jetmettau, "TUnfold", out_hw_eec3_unfolded_svd_data_jetmettau, "SVD",  out_hw_eec3_unfolded_bayes_data_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-By-Bin", "SVD/Bin-By-Bin", "Bayes/Bin-By-Bin", compare_unfolding_plots, "data_with_hw_eec3_jetmettau_", false, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmettau, "Detector Level", out_p8_4c_gen_allvertex_jetmettau, "Generator Level", out_hw_eec3_unfolded_bbb_p8_4c_jetmettau, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_p8_4c_jetmettau, "TUnfold",  out_hw_eec3_unfolded_svd_p8_4c_jetmettau, "SVD", out_hw_eec3_unfolded_bayes_p8_4c_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "p8_4c_with_hw_eec3_jetmettau_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmettau, "Detector Level", out_hw_23_gen_allvertex_jetmettau, "Generator Level", out_hw_eec3_unfolded_bbb_hw_23_jetmettau, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_hw_23_jetmettau, "TUnfold",  out_hw_eec3_unfolded_svd_hw_23_jetmettau, "SVD", out_hw_eec3_unfolded_bayes_hw_23_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_23_with_hw_eec3_jetmettau_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmettau, "Detector Level", out_hw_eec3_gen_allvertex_jetmettau, "Generator Level", out_hw_eec3_unfolded_bbb_hw_eec3_jetmettau, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold",  out_hw_eec3_unfolded_svd_hw_eec3_jetmettau, "SVD", out_hw_eec3_unfolded_bayes_hw_eec3_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True",  compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmettau_", true, detail);

//method test
compare_unfolding_results(out_hw_eec3_det_allvertex_jetmettau, "Detector Level", out_hw_eec3_gen_allvertex_jetmettau, "Generator Level", out_hw_eec3_unfolded_bayes1_hw_eec3_jetmettau, "Bayes: iterations = 1", out_hw_eec3_unfolded_bayes2_hw_eec3_jetmettau, "Bayes: iterations = 2", out_hw_eec3_unfolded_bayes3_hw_eec3_jetmettau, "Bayes: iterations = 3", out_hw_eec3_unfolded_bayes4_hw_eec3_jetmettau, "Bayes: iterations = 4", "Bayes: iterations = 1 /True", "Bayes: iterations = 2 /True", "Bayes: iterations = 3 /True", "Bayes: iterations = 4 /True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmettau_bayes_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmettau, "Detector Level", out_hw_eec3_gen_allvertex_jetmettau, "Generator Level", out_hw_eec3_unfolded_svd1_hw_eec3_jetmettau, "SVD: regularization = 1", out_hw_eec3_unfolded_svd2_hw_eec3_jetmettau, "SVD: regularization = 2", out_hw_eec3_unfolded_svd3_hw_eec3_jetmettau, "SVD: regularization = 3", out_hw_eec3_unfolded_svd4_hw_eec3_jetmettau, "SVD: regularization = 4", "SVD: regularization = 1 /True", "SVD: regularization = 2 /True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmettau_svd_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmettau, "Detector Level", out_hw_eec3_gen_allvertex_jetmettau, "Generator Level", out_hw_eec3_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold: Default", out_hw_eec3_unfolded_tunfold2_hw_eec3_jetmettau, "TUnfold: Size Regularization", out_hw_eec3_unfolded_tunfold3_hw_eec3_jetmettau, "TUnfold: Derivative Regularization", out_hw_eec3_unfolded_tunfold4_hw_eec3_jetmettau, "TUnfold: Curvature Regularization", "TUnfold: Default/True", "TUnfold: No Regularization/True", "TUnfold: Size Regularization/True", "TUnfold: Derivative Regularization/True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmettau_tunfold_", true, detail);


//herwig 23 - jetmettau
compare_unfolding_results(out_jetmettau_allvertex, "Uncorrected Data", corrected_hw_23_jetmettau, "Standard Bin-by-Bin", out_hw_23_unfolded_bbb_data_jetmettau, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_data_jetmettau, "TUnfold", out_hw_23_unfolded_svd_data_jetmettau, "SVD",  out_hw_23_unfolded_bayes_data_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-By-Bin", "SVD/Bin-By-Bin", "Bayes/Bin-By-Bin", compare_unfolding_plots, "data_with_hw_23_jetmettau_", false, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmettau, "Detector Level", out_p8_4c_gen_allvertex_jetmettau, "Generator Level", out_hw_23_unfolded_bbb_p8_4c_jetmettau, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_p8_4c_jetmettau, "TUnfold", out_hw_23_unfolded_svd_p8_4c_jetmettau, "SVD", out_hw_23_unfolded_bayes_p8_4c_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "p8_4c_with_hw_23_jetmettau_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmettau, "Detector Level", out_hw_eec3_gen_allvertex_jetmettau, "Generator Level", out_hw_23_unfolded_bbb_hw_eec3_jetmettau, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_hw_eec3_jetmettau, "TUnfold",  out_hw_23_unfolded_svd_hw_eec3_jetmettau, "SVD", out_hw_23_unfolded_bayes_hw_eec3_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_eec3_with_hw_23_jetmettau_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmettau, "Detector Level",  out_hw_23_gen_allvertex_jetmettau, "Generator Level", out_hw_23_unfolded_bbb_hw_23_jetmettau, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_hw_23_jetmettau, "TUnfold", out_hw_23_unfolded_svd_hw_23_jetmettau, "SVD", out_hw_23_unfolded_bayes_hw_23_jetmettau, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True",  compare_unfolding_plots, "hw_23_with_hw_23_jetmettau_", true, detail);

//method test
compare_unfolding_results(out_hw_23_det_allvertex_jetmettau, "Detector Level", out_hw_23_gen_allvertex_jetmettau, "Generator Level", out_hw_23_unfolded_bayes1_hw_23_jetmettau, "Bayes: iterations = 1", out_hw_23_unfolded_bayes2_hw_23_jetmettau, "Bayes: iterations = 2", out_hw_23_unfolded_bayes3_hw_23_jetmettau, "Bayes: iterations = 3", out_hw_23_unfolded_bayes4_hw_23_jetmettau, "Bayes: iterations = 4", "Bayes: iterations = 1 /True", "Bayes: iterations = 2 /True", "Bayes: iterations = 3 /True", "Bayes: iterations = 4 /True", compare_unfolding_plots, "hw_23_with_hw_23_jetmettau_bayes_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmettau, "Detector Level", out_hw_23_gen_allvertex_jetmettau, "Generator Level", out_hw_23_unfolded_svd1_hw_23_jetmettau, "SVD: regularization = 1", out_hw_23_unfolded_svd2_hw_23_jetmettau, "SVD: regularization = 2", out_hw_23_unfolded_svd3_hw_23_jetmettau, "SVD: regularization = 3", out_hw_23_unfolded_svd4_hw_23_jetmettau, "SVD: regularization = 4", "SVD: regularization = 1 /True", "SVD: regularization = 2 /True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "hw_23_with_hw_23_jetmettau_svd_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmettau, "Detector Level", out_hw_23_gen_allvertex_jetmettau, "Generator Level", out_hw_23_unfolded_tunfold_hw_23_jetmettau, "TUnfold: Default", out_hw_23_unfolded_tunfold2_hw_23_jetmettau, "TUnfold: Size Regularization", out_hw_23_unfolded_tunfold3_hw_23_jetmettau, "TUnfold: Derivative Regularization", out_hw_23_unfolded_tunfold4_hw_23_jetmettau, "TUnfold: Curvature Regularization", "TUnfold: Default/True", "TUnfold: Size Regularization/True", "TUnfold: Derivative Regularization/True", "TUnfold: Curvature Regularization/True", compare_unfolding_plots, "hw_23_with_hw_23_jetmettau_tunfold_", true, detail);


//pythia 8 - jetmet
compare_unfolding_results(out_jetmet_allvertex, "Uncorrected Data", corrected_p8_4c_jetmet, "Standard Bin-by-Bin", out_p8_4c_unfolded_bbb_data_jetmet, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_data_jetmet, "TUnfold", out_p8_4c_unfolded_svd_data_jetmet, "SVD",  out_p8_4c_unfolded_bayes_data_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-by-Bin", "SVD/Bin-by-Bin", "Bayes/Bin-by-Bin", compare_unfolding_plots, "data_with_p8_4c_jetmet_", false, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmet, "Detector Level", out_hw_eec3_gen_allvertex_jetmet, "Generator Level", out_p8_4c_unfolded_bbb_hw_eec3_jetmet, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_hw_eec3_jetmet, "TUnfold", out_p8_4c_unfolded_svd_hw_eec3_jetmet, "SVD", out_p8_4c_unfolded_bayes_hw_eec3_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_eec3_with_p8_4c_jetmet_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_p8_4c_unfolded_bbb_hw_23_jetmet, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_hw_23_jetmet, "TUnfold", out_p8_4c_unfolded_svd_hw_23_jetmet, "SVD",  out_p8_4c_unfolded_bayes_hw_23_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_23_with_p8_4c_jetmet_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmet, "Detector Level", out_p8_4c_gen_allvertex_jetmet, "Generator Level", out_p8_4c_unfolded_bbb_p8_4c_jetmet, "RooUnfold Bin-by-Bin", out_p8_4c_unfolded_tunfold_p8_4c_jetmet, "TUnfold", out_p8_4c_unfolded_svd_p8_4c_jetmet, "SVD",  out_p8_4c_unfolded_bayes_p8_4c_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmet_", true, detail);

//method test
compare_unfolding_results(out_p8_4c_det_allvertex_jetmet, "Detector Level", out_p8_4c_gen_allvertex_jetmet, "Generator Level", out_p8_4c_unfolded_bayes1_p8_4c_jetmet, "Bayes: iterations = 1", out_p8_4c_unfolded_bayes2_p8_4c_jetmet, "Bayes: iterations = 2",  out_p8_4c_unfolded_bayes3_p8_4c_jetmet, "Bayes: iterations = 3", out_p8_4c_unfolded_bayes4_p8_4c_jetmet, "Bayes: iterations = 4", "Bayes: iterations = 1 /True", "Bayes: iterations = 2 /True", "Bayes: iterations = 3 /True", "Bayes: iterations = 4 /True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmet_bayes_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmet, "Detector Level", out_p8_4c_gen_allvertex_jetmet, "Generator Level", out_p8_4c_unfolded_svd1_p8_4c_jetmet, "SVD: regularization = 1", out_p8_4c_unfolded_svd2_p8_4c_jetmet, "SVD: regularization = 2",  out_p8_4c_unfolded_svd3_p8_4c_jetmet, "SVD: regularization = 3", out_p8_4c_unfolded_svd4_p8_4c_jetmet, "SVD: regularization = 4", "SVD: regularization = 1 /True", "SVD: regularization = 2 /True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmet_svd_", true, detail);

compare_unfolding_results(out_p8_4c_det_allvertex_jetmet, "Detector Level", out_p8_4c_gen_allvertex_jetmet, "Generator Level", out_p8_4c_unfolded_tunfold_p8_4c_jetmet, "TUnfold: Default", out_p8_4c_unfolded_tunfold2_p8_4c_jetmet, "TUnfold: Size Regularization",  out_p8_4c_unfolded_tunfold3_p8_4c_jetmet, "TUnfold: Derivative Regularization", out_p8_4c_unfolded_tunfold4_p8_4c_jetmet, "TUnfold: Curvature Regularization", "TUnfold: Default/True", "TUnfold: Size Regularization/True", "TUnfold: Derivative Regularization/True", "TUnfold: Curvature Regularization/True", compare_unfolding_plots, "p8_4c_with_p8_4c_jetmet_tunfold_", true, detail);



//herwig eec3 - jetmet
compare_unfolding_results(out_jetmet_allvertex, "Uncorrected Data", corrected_hw_eec3_jetmet, "Standard Bin-by-Bin", out_hw_eec3_unfolded_bbb_data_jetmet, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_data_jetmet, "TUnfold", out_hw_eec3_unfolded_svd_data_jetmet, "SVD",  out_hw_eec3_unfolded_bayes_data_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-by-Bin", "SVD/Bin-by-Bin", "Bayes/Bin-by-Bin", compare_unfolding_plots, "data_with_hw_eec3_jetmet_", false, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmet, "Detector Level", out_p8_4c_gen_allvertex_jetmet, "Generator Level", out_hw_eec3_unfolded_bbb_p8_4c_jetmet, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_p8_4c_jetmet, "TUnfold", out_hw_eec3_unfolded_svd_p8_4c_jetmet, "SVD", out_hw_eec3_unfolded_bayes_p8_4c_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "p8_4c_with_hw_eec3_jetmet_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_hw_eec3_unfolded_bbb_hw_23_jetmet, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_hw_23_jetmet, "TUnfold", out_hw_eec3_unfolded_svd_hw_23_jetmet, "SVD", out_hw_eec3_unfolded_bayes_hw_23_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_23_with_hw_eec3_jetmet_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmet, "Detector Level", out_hw_eec3_gen_allvertex_jetmet, "Generator Level", out_hw_eec3_unfolded_bbb_hw_eec3_jetmet, "RooUnfold Bin-by-Bin", out_hw_eec3_unfolded_tunfold_hw_eec3_jetmet, "TUnfold", out_hw_eec3_unfolded_svd_hw_eec3_jetmet, "SVD", out_hw_eec3_unfolded_bayes_hw_eec3_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/True", "TUnfold/True", "SVD/True", "Bayes/True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmet_", true, detail);

//Method test
compare_unfolding_results(out_hw_eec3_det_allvertex_jetmet, "Detector Level", out_hw_eec3_gen_allvertex_jetmet, "Generator Level", out_hw_eec3_unfolded_bayes1_hw_eec3_jetmet, "Bayes: iterations = 1", out_hw_eec3_unfolded_bayes2_hw_eec3_jetmet, "Bayes: iterations = 2", out_hw_eec3_unfolded_bayes3_hw_eec3_jetmet, "Bayes: iterations = 3", out_hw_eec3_unfolded_bayes4_hw_eec3_jetmet, "Bayes: iterations = 4", "Bayes: iterations = 1 /True", "Bayes: iterations = 2 /True", "Bayes: iterations = 3 /True", "Bayes: iterations = 4 /True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmet_bayes_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmet, "Detector Level", out_hw_eec3_gen_allvertex_jetmet, "Generator Level", out_hw_eec3_unfolded_svd1_hw_eec3_jetmet, "SVD: regularization = 1", out_hw_eec3_unfolded_svd2_hw_eec3_jetmet, "SVD: regularization = 2", out_hw_eec3_unfolded_svd3_hw_eec3_jetmet, "SVD: regularization = 3", out_hw_eec3_unfolded_svd4_hw_eec3_jetmet, "SVD: regularization = 4", "SVD: regularization = 1 /True", "SVD: regularization = 2 /True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmet_svd_", true, detail);

compare_unfolding_results(out_hw_eec3_det_allvertex_jetmet, "Detector Level", out_hw_eec3_gen_allvertex_jetmet, "Generator Level", out_hw_eec3_unfolded_tunfold_hw_eec3_jetmet, "TUnfold: Default", out_hw_eec3_unfolded_tunfold2_hw_eec3_jetmet, "TUnfold: Size Regularization", out_hw_eec3_unfolded_tunfold3_hw_eec3_jetmet, "TUnfold: Derivative Regularization", out_hw_eec3_unfolded_tunfold4_hw_eec3_jetmet, "TUnfold: Curvature Regularization", "TUnfold: Default/True", "TUnfold: Size Regularization/True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "hw_eec3_with_hw_eec3_jetmet_tunfold_", true, detail);



//herwig 23 - jetmet
compare_unfolding_results(out_jetmet_allvertex, "Uncorrected Data", corrected_hw_23_jetmet, "Standard Bin-by-Bin", out_hw_23_unfolded_bbb_data_jetmet, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_data_jetmet, "TUnfold", out_hw_23_unfolded_svd_data_jetmet, "SVD",  out_hw_23_unfolded_bayes_data_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-by-Bin", "SVD/Bin-by-Bin", "Bayes/Bin-by-Bin", compare_unfolding_plots, "data_with_hw_23_jetmet_", false, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_hw_23_unfolded_bbb_hw_23_jetmet, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_hw_23_jetmet, "TUnfold", out_hw_23_unfolded_svd_hw_23_jetmet, "SVD",  out_hw_23_unfolded_bayes_hw_23_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-by-Bin", "SVD/Bin-by-Bin", "Bayes/Bin-by-Bin", compare_unfolding_plots, "hw_23_with_hw_23_jetmet_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmet, "Detector Level", out_hw_eec3_gen_allvertex_jetmet, "Generator Level", out_hw_23_unfolded_bbb_hw_eec3_jetmet, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_hw_eec3_jetmet, "TUnfold", out_hw_23_unfolded_svd_hw_eec3_jetmet, "SVD", out_hw_23_unfolded_bayes_hw_eec3_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-by-Bin", "SVD/Bin-by-Bin", "Bayes/Bin-by-Bin", compare_unfolding_plots, "hw_eec3_with_hw_23_jetmet_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_hw_23_unfolded_bbb_hw_23_jetmet, "RooUnfold Bin-by-Bin", out_hw_23_unfolded_tunfold_hw_23_jetmet, "TUnfold",  out_hw_23_unfolded_svd_hw_23_jetmet, "SVD", out_hw_23_unfolded_bayes_hw_23_jetmet, "Bayes", "Bin-By-Bin - RooUnfold/Standard", "TUnfold/Bin-by-Bin", "SVD/Bin-by-Bin", "Bayes/Bin-by-Bin", compare_unfolding_plots, "hw_23_with_hw_23_jetmet_", true, detail);

//method test
compare_unfolding_results(out_hw_23_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_hw_23_unfolded_bayes1_hw_23_jetmet, "Bayes: iterations = 1", out_hw_23_unfolded_bayes2_hw_23_jetmet, "Bayes: iterations = 2", out_hw_23_unfolded_bayes3_hw_23_jetmet, "Bayes: iterations = 3", out_hw_23_unfolded_bayes4_hw_23_jetmet, "Bayes: iterations = 4", "Bayes: iterations = 1 /True", "Bayes: iterations = 2 /True", "Bayes: iterations = 3 /True", "Bayes: iterations = 4 /True", compare_unfolding_plots, "hw_23_with_hw_23_jetmet_bayes_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_hw_23_unfolded_svd1_hw_23_jetmet, "SVD: regularization = 1", out_hw_23_unfolded_svd2_hw_23_jetmet, "SVD: regularization = 2", out_hw_23_unfolded_svd3_hw_23_jetmet, "SVD: regularization = 3", out_hw_23_unfolded_svd4_hw_23_jetmet, "SVD: regularization = 4", "SVD: regularization = 1 /True", "SVD: regularization = 2 /True", "SVD: regularization = 3 /True", "SVD: regularization = 4 /True", compare_unfolding_plots, "hw_23_with_hw_23_jetmet_svd_", true, detail);

compare_unfolding_results(out_hw_23_det_allvertex_jetmet, "Detector Level", out_hw_23_gen_allvertex_jetmet, "Generator Level", out_hw_23_unfolded_tunfold_hw_23_jetmet, "TUnfold: Default", out_hw_23_unfolded_tunfold2_hw_23_jetmet, "TUnfold: Size Regularization", out_hw_23_unfolded_tunfold3_hw_23_jetmet, "TUnfold: Derivative Regularization", out_hw_23_unfolded_tunfold4_hw_23_jetmet, "TUnfold: Curvature Regularization", "TUnfold: Default/True", "TUnfold: Size Regularization/True", "TUnfold: Derivative Regularization/True", "TUnfold: Curvature Regularization/True", compare_unfolding_plots, "hw_23_with_hw_23_jetmet_tunfold_", true, detail);


}



//merge the different datasets
if (merge_data)
{
if (show_steps) { cout << "Merging Data..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, merged_data_plots, "plots");
create_directories(output_dir, merged_uncorrected_data_plots, "plots");

//uncorrected data
//merge_data(out_jetmettau_allvertex, out_jetmet_allvertex, out_jet_allvertex, merged_uncorrected_data_plots, merged_uncorrected_data, detail, show_results, test);

//old merging - without unfolding
//merge_data(corrected_jetmettau, corrected_jetmet, corrected_jet, merged_data_plots, merged_data, detail, show_results, test);

//new merging - with unfolding
merge_unfolded_data(out_p6_z2_unfolded_bayes_data_jetmettau, out_p6_z2_unfolded_bayes_data_jetmet, out_p8_4c_unfolded_bayes_data_jetmettau, out_p8_4c_unfolded_bayes_data_jetmet, merged_data_plots, merged_unfolded_data, detail, show_results, test);

if (show_steps) { cout << "The data was sucessfully merged!"<<endl; }
}


//estimate model uncertainty
if (compute_model_uncertainty)
{
if (show_steps) { cout << "Calculating Uncertainties..."<<endl; }
if (show_steps) { cout << "Estimating Model Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, model_unc_dir, "root");
create_directories(output_dir, model_unc_plots, "plots");
//estimate_model_uncertainty(corr_final_p6_z2, corr_final_p8_4c, out_jetmettau_allvertex, model_uncertainty, model_unc_plots, detail, show_results, test); //needed merged - maybe not needed

//old chain - before unfolding
//estimate_model_uncertainty(corr_detector_p6_z2_jetmettau, corr_detector_p8_4c_jetmettau, out_jetmettau_allvertex, model_unc_jetmettau, "JetMETTau_2010A", model_unc_plots, detail, show_results, false, test);
//estimate_model_uncertainty(corr_detector_p6_z2_jetmet, corr_detector_p8_4c_jetmet, out_jetmet_allvertex, model_unc_jetmet, "JetMET_2010A", model_unc_plots, detail, show_results, false, test);
//estimate_model_uncertainty(corr_detector_p6_z2_jet, corr_detector_p8_4c_jet, out_jet_allvertex, model_unc_jet, "Jet_2010B", model_unc_plots, detail, show_results, false, test);


//new chain - after unfolding
estimate_model_uncertainty(out_p6_z2_unfolded_bayes_data_jetmettau, out_p8_4c_unfolded_bayes_data_jetmettau, out_jetmettau_allvertex, model_unc_jetmettau, "JetMETTau_2010A", model_unc_plots, detail, show_results, true, test);
estimate_model_uncertainty(out_p6_z2_unfolded_bayes_data_jetmet, out_p8_4c_unfolded_bayes_data_jetmet, out_jetmet_allvertex, model_unc_jetmet, "JetMET_2010A", model_unc_plots, detail, show_results, true, test);
estimate_model_uncertainty(out_p6_z2_unfolded_bayes_data_jet, out_p8_4c_unfolded_bayes_data_jet, out_jet_allvertex, model_unc_jet, "Jet_2010B", model_unc_plots, detail, show_results, true, test);

//merging
merge_uncertainties(model_unc_jetmettau, model_unc_jetmet, model_unc_jet, model_unc_merged, "model_unc_", "merged_model_unc_", model_unc_plots, detail, show_results, test);
if (show_steps) { cout << "The model uncertainty was sucessfully calculated!"<<endl; }
}

//estimate jes uncertainty
if (compute_jes_uncertainty)
{
if (show_steps) { cout << "Calculating Uncertainties..."<<endl; }
if (show_steps) { cout << "Estimating JES Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, jes_unc_dir, "root");
create_directories(output_dir, jes_unc_plots, "plots");
//estimate_model_uncertainty(corr_final_p6_z2, corr_final_p8_4c, out_jetmettau_allvertex, model_uncertainty, model_unc_plots, detail, show_results, test); //needed merged - maybe not needed

estimate_jes_uncertainty(out_jetmettau_allvertex, out_jetmettau_up, jes_up_unc_jetmettau, "JetMETTau_2010A_up_", "jes_unc_up_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jetmettau_allvertex, out_jetmettau_down, jes_down_unc_jetmettau, "JetMETTau_2010A_down_", "jes_unc_down_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jetmet_allvertex, out_jetmet_up, jes_up_unc_jetmet, "JetMET_2010A_up", "jes_unc_up_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jetmet_allvertex, out_jetmet_down, jes_down_unc_jetmet, "JetMET_2010A_down", "jes_unc_down_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jet_allvertex, out_jet_up, jes_up_unc_jet, "Jet_2010A_up", "jes_unc_up_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jet_allvertex, out_jet_down, jes_down_unc_jet, "Jet_2010A_down", "jes_unc_down_", jes_unc_plots, detail, show_results, test);

//merging
merge_uncertainties(jes_up_unc_jetmettau, jes_up_unc_jetmet, jes_up_unc_jet, jes_unc_up_merged, "jes_unc_up_", "merged_jes_unc_up_", jes_unc_plots, detail,  show_results, test);
merge_uncertainties(jes_down_unc_jetmettau, jes_down_unc_jetmet, jes_down_unc_jet, jes_unc_down_merged, "jes_unc_down_", "merged_jes_unc_down_", jes_unc_plots, detail, show_results, test);
if (show_steps) { cout << "The JES uncertainty was sucessfully calculated!"<<endl; }
}


//estimate jer uncertainty
if (compute_jer_uncertainty)
{
if (show_steps) { cout << "Estimating JER Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, jer_unc_dir, "root");
create_directories(output_dir, jer_unc_plots, "plots");

estimate_jer_uncertainty(out_p6_z2_det_allvertex_jetmettau, out_p6_z2_det_allvertex_jetmettau_up, out_p6_z2_det_allvertex_jetmettau_down, jer_unc_p6_jetmettau, "P6_JetMETTau_2010A_", "jer_unc_", jer_unc_plots, detail, show_results, test);
estimate_jer_uncertainty(out_p8_4c_det_allvertex_jetmettau, out_p8_4c_det_allvertex_jetmettau_up, out_p8_4c_det_allvertex_jetmettau_down, jer_unc_p8_jetmettau, "P8_JetMETTau_2010A_", "jer_unc_", jer_unc_plots, detail, show_results, test);
estimate_jer_uncertainty(out_p6_z2_det_allvertex_jetmet, out_p6_z2_det_allvertex_jetmet_up, out_p6_z2_det_allvertex_jetmet_down, jer_unc_p6_jetmet, "P6_JetMET_2010A_", "jer_unc_", jer_unc_plots, detail, show_results, test);
estimate_jer_uncertainty(out_p8_4c_det_allvertex_jetmet, out_p8_4c_det_allvertex_jetmet_up, out_p8_4c_det_allvertex_jetmet_down, jer_unc_p8_jetmet, "P8_JetMET_2010A_", "jer_unc_", jer_unc_plots, detail, show_results, test);
estimate_jer_uncertainty(out_jet_allvertex, out_p6_z2_det_allvertex_jet_up, out_p6_z2_det_allvertex_jet_down, jer_unc_p6_jet, "P6_Jet_2010B_", "jer_unc_", jer_unc_plots, detail, show_results, test);
estimate_jer_uncertainty(out_jet_allvertex, out_p8_4c_det_allvertex_jet_up, out_p8_4c_det_allvertex_jet_down, jer_unc_p8_jet, "P8_Jet_2010B_", "jer_unc_", jer_unc_plots, detail, show_results, test);

//merging
merge_uncertainties(jer_unc_p6_jetmettau, jer_unc_p6_jetmet, jer_unc_p6_jet, jer_unc_p6_merged, "jer_unc_", "run_merged_jer_unc_", jer_unc_plots, detail,  show_results, test);
merge_uncertainties(jer_unc_p8_jetmettau, jer_unc_p8_jetmet, jer_unc_p8_jet, jer_unc_p8_merged, "jer_unc_", "run_merged_jer_unc_", jer_unc_plots, detail,  show_results, test);

merge_uncertainties(jer_unc_p6_merged, jer_unc_p8_merged, jer_unc_p6_merged, jer_unc_merged, "run_merged_jer_unc_", "merged_jer_unc_", jer_unc_plots, detail,  show_results, test);

if (show_steps) { cout << "The JER uncertainty was sucessfully calculated!"<<endl; }
}

//Compute Correlated Uncertainty
if (compute_corr_uncertainty)
{
if (show_steps) { cout << "Computing Correlated Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, corr_unc_dir, "root");
create_directories(output_dir, corr_unc_plots, "plots");

corr_uncertainty(jes_unc_up_merged, "merged_jes_unc_up_", corr_unc_up_merged, "corr_unc_up_", corr_unc_plots, detail, show_results, test);
corr_uncertainty(jes_unc_down_merged, "merged_jes_unc_down_", corr_unc_down_merged, "corr_unc_down_", corr_unc_plots, detail, show_results, test);
}


//stat_vs_unf_unc
if (run_stat_vs_unf_unc)
{
if (show_steps) { cout << "Computing Statistical VS Unfolding Uncertainty..."<<endl; }
create_directories(output_dir, stat_vs_unf_unc_plots, "plots");

stat_vs_unf_unc(merged_data, merged_unfolded_data, stat_vs_unf_unc_plots, true, true,false);
}


//Compute Uncorrelated Uncertainty
if (compute_uncorr_uncertainty)
{
if (show_steps) { cout << "Computing Uncorrelated Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, uncorr_unc_dir, "root");
create_directories(output_dir, uncorr_unc_plots, "plots");

uncorr_uncertainty(merged_unfolded_data, model_unc_merged, jer_unc_merged, uncorr_unc_merged, "uncorr_unc_", uncorr_unc_plots, detail, show_results, test);
}

//estimate total uncertainty
if (compute_total_uncertainty)
{
if (show_steps) { cout << "Computing Total Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, total_unc_dir, "root");
create_directories(output_dir, total_unc_plots, "plots");

//total_uncertainty(jes_unc_up_merged, model_unc_merged, merged_unfolded_data, total_unc_up, "merged_jes_unc_up_", "total_unc_up_", total_unc_plots, detail, show_results, test);
//total_uncertainty(jes_unc_down_merged, model_unc_merged, merged_unfolded_data, total_unc_down, "merged_jes_unc_down_", "total_unc_down_", total_unc_plots, detail, show_results, test);

total_uncertainty_v2(corr_unc_up_merged, corr_unc_down_merged, uncorr_unc_merged, total_unc, total_unc_plots, detail, show_results, test);

}


//estimate_combination_systematic
if (estimate_combination_systematic)
{
if (show_steps) { cout << "Estimate Combination Systematics..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_syst_dir, "root");
create_directories(output_dir, trigger_syst_plots, "plots");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//estimate_combination_systematic(jetmettaumon, systematics_trigger_jetmettaumon, n_files_jetmettaumon, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//estimate_combination_systematic(jetmettau, systematics_trigger_jetmettau, n_files_jetmettau, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//estimate_combination_systematic(jetmet, systematics_trigger_jetmet, n_files_jetmet, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
estimate_combination_systematic(jet, systematics_trigger_jet, n_files_jet, "allvertex", detail, test);

if (show_steps) { cout << "Merging and ploting"<<endl; }
merge_plot_combination_systematic(systematics_trigger_jetmettau, systematics_trigger_jetmet, systematics_trigger_jet, merged_trig_comb_syst, trigger_syst_plots, detail, test);
if (show_steps) { cout << "The Estimation was sucessfully done!"<<endl; }
}


//apply uncertainty
if (apply_uncertainty)
{
if (show_steps) { cout << "Apply Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, data_unc_plots, "plots");

apply_corr_uncertainty(merged_unfolded_data, corr_unc_up_merged, corr_unc_down_merged, data_unc, data_unc_plots, detail, test);
apply_uncorr_uncertainty(merged_unfolded_data, uncorr_unc_merged, merged_uncorr_unfolded_data, data_unc_plots, detail, test);
if (show_steps) { cout << "The uncertainty was sucessfully applied!"<<endl; }
}


//ratio plots
if (create_ratio)
{
if (show_steps) { cout << "Ratio Plots..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, ratio_dir, "root");
create_directories(output_dir, ratio_plots, "plots");
create_directories(output_dir, correction_ratio_plots, "plots");
create_directories(output_dir, closure_plots, "plots");

create_ratios(merged_uncorr_unfolded_data, merged_uncorr_unfolded_data, "ak5PF_", "ak5PF_", ratio_merged_unfolded_data, ratio_plots, "data_data_", detail, test);
create_ratios(data_unc, merged_uncorr_unfolded_data, "ak5PF_", "ak5PF_", ratio_data_unc, ratio_plots, "dataunc_data_", detail, test);
create_ratios(out_p6_z2_gen_allvertex, merged_uncorr_unfolded_data, "ak5Gen_", "ak5PF_", mc_ratio_py6_z2_off, ratio_plots, "pythia6z2off_data_", detail, test);
create_ratios(out_p8_4c_gen_allvertex, merged_uncorr_unfolded_data, "ak5Gen_", "ak5PF_", mc_ratio_py8_4c_off, ratio_plots, "pythia84coff_data_", detail, test);
create_ratios(mc_pred_herwig6, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_herwig6, ratio_plots, "herwig6_data_", detail, test);
create_ratios(mc_pred_herwigpp, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_herwigpp, ratio_plots, "herwigpp_data_", detail, test);
create_ratios(mc_pred_py8_4c, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py8_4c, ratio_plots, "pythia8_data_", detail, test);
create_ratios(mc_pred_py6_p11, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_p11, ratio_plots, "pythia6p11_data_", detail, test);
create_ratios(mc_pred_py6_z2, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_z2, ratio_plots, "pythia6z2_data_", detail, test);
create_ratios(mc_pred_py6_ambt1, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_ambt1, ratio_plots, "pythia6ambt1_data_", detail, test);
create_ratios(mc_pred_py6_z2_nompi, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_z2_nompi, ratio_plots, "pythia6z2nompi_data_", detail, test);
create_ratios(mc_pred_py6_z2_nompi2, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_z2_nompi2, ratio_plots, "pythia6z2nompi2_data_", detail, test);
create_ratios(mc_pred_py6_z2_nompi_nohad, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_z2_nompi_nohad, ratio_plots, "pythia6z2nompinohad_data_", detail, test);
create_ratios(mc_pred_py6_z2_nompi_nohad_nofsr, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_z2_nompi_nohad_nofsr, ratio_plots, "pythia6z2nompinohadnofsr_data_", detail, test);
create_ratios(mc_pred_py6_z2_nompi_nohad_nofsr_noisr, merged_uncorr_unfolded_data, "", "ak5PF_", mc_ratio_py6_z2_nompi_nohad_nofsr_noisr, ratio_plots, "pythia6z2nompinohadnofsrnoisr_data_", detail, test);


create_ratios(corr_detector_p6_z2_jetmettau, corr_detector_p6_z2, "detector_pythia6_z2_JetMETTau_2010A_", "detector_pythia6_z2_", mc_ratio_py6_z2_jetmettau, ratio_plots, "pythia6z2jetmettau_", detail, test);
create_ratios(corr_detector_p6_z2_jetmet, corr_detector_p6_z2, "detector_pythia6_z2_JetMET_2010A_", "detector_pythia6_z2_", mc_ratio_py6_z2_jetmet, ratio_plots, "pythia6z2jetmet_", detail, test);
create_ratios(corr_detector_p6_z2_jet, corr_detector_p6_z2, "detector_pythia6_z2_Jet_2010B_", "detector_pythia6_z2_", mc_ratio_py6_z2_jet, ratio_plots, "pythia6z2jet_", detail, test);

create_ratios(corr_detector_p8_4c_jetmettau, corr_detector_p8_4c, "detector_pythia8_4c_JetMETTau_2010A_", "detector_pythia8_4c_", mc_ratio_py8_4c_jetmettau, ratio_plots, "pythia84cjetmettau_", detail, test);
create_ratios(corr_detector_p8_4c_jetmet, corr_detector_p8_4c, "detector_pythia8_4c_JetMET_2010A_", "detector_pythia8_4c_", mc_ratio_py8_4c_jetmet, ratio_plots, "pythia84cjetmet_", detail, test);
create_ratios(corr_detector_p8_4c_jet, corr_detector_p8_4c, "detector_pythia8_4c_Jet_2010B_", "detector_pythia8_4c_", mc_ratio_py8_4c_jet, ratio_plots, "pythia84cjet_", detail, test);

create_ratios(corrected_p6_z2_jetmettau, out_p6_z2_gen_allvertex_jetmettau, "ak5PF_", "ak5Gen_", closure_ratio_p6_z2_jetmettau, closure_plots, "closure_p6_z2_jetmettau_", detail, test);

if (show_steps) { cout << "The ratios were sucessfully computed!"<<endl; }
}

//final plots
if (do_final_plots)
{
if (show_steps) { cout << "Final Plots..."<<endl; }
create_directories(output_dir, final_plots, "plots");

final_plots(merged_uncorr_unfolded_data, data_unc, final_mc_list, final_label_list, final_prefix_list, final_plots, "xsec", detail, test);
//final_plots(ratio_merged_unfolded_data, ratio_data_unc, ratio_mc_list, ratio_label_list, ratio_prefix_list, final_plots, "ratio", detail, test);

if (show_steps) { cout << "The final plots were sucessfully generated!"<<endl; }
}

}
