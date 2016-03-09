// Pedro Cipriano, Jun 2013
// DESY, CMS
// Last Update: 17 Jun 2012
//
//


#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "common_methods.h"

void plots_unfolding_ratios(TH1D *hist2, TH1D *hist3, TString ratiolabel1, TH1D *hist4, TString ratiolabel2, TH1D *hist5, TString ratiolabel3, TH1D *hist6, TString ratiolabel4, string plots_path, string name, bool detail)
{

    TH1D *ratio1;
    ratio1 = static_cast<TH1D*>(hist2->Clone());
    ratio1->Divide(hist3,hist2,1.,1.,"B");
    TH1D *ratio2;
    ratio2 = static_cast<TH1D*>(hist2->Clone());
    ratio2->Divide(hist4,hist2,1.,1.,"B");
    TH1D *ratio3;
    ratio3 = static_cast<TH1D*>(hist2->Clone());
    ratio3->Divide(hist5,hist2,1.,1.,"B");
    TH1D *ratio4;
    ratio4 = static_cast<TH1D*>(hist2->Clone());
    ratio4->Divide(hist6,hist2,1.,1.,"B");
    plot_4histograms(ratio1, ratiolabel1, ratio2, ratiolabel2, ratio3, ratiolabel3, ratio4, ratiolabel4, plots_path, name, "top_right", false, detail);
    //plot_3histograms(ratio1, "Bin-By-Bin - RooUnfold/Standard", ratio3, "SVD/Bin-By-Bin", ratio4, "Bayes/Bin-By-Bin", plots_path, prefix+"delta_phi_ratios", "top_right", false, detail);
    plot_histogram(ratio4, plots_path, name + "_simple", ratiolabel4, "top_left", false);

}

void compare_unfolding_results(string infile1, TString label1, string infile2, TString label2, string infile3, TString label3, string infile4, TString label4, string infile5, TString label5, string infile6, TString label6, TString ratiolabel1, TString ratiolabel2, TString ratiolabel3, TString ratiolabel4, string plots_path, string prefix, bool isMC = false, bool detail = false)
{

//output the configuration
   if (detail) { cout<<"Compare Unfolding Results Configuration"<<endl; }
   if (detail) { cout<<"Detector Level File :               "<<infile1<<endl; }
   if (detail) { cout<<"Label 1 :                           "<<label1<<endl; }
   if (detail) { cout<<"Generator Level File :              "<<infile2<<endl; }
   if (detail) { cout<<"Label 2 :                           "<<label2<<endl; }
   if (detail) { cout<<"Method 1 Results File :             "<<infile3<<endl; }
   if (detail) { cout<<"Label 3 :                           "<<label3<<endl; }
   if (detail) { cout<<"Method 2 Results File :             "<<infile4<<endl; }
   if (detail) { cout<<"Label 4 :                           "<<label4<<endl; }
   if (detail) { cout<<"Method 3 Results File :             "<<infile5<<endl; }
   if (detail) { cout<<"Label 5 :                           "<<label5<<endl; }
   if (detail) { cout<<"Method 4 Results File :             "<<infile6<<endl; }
   if (detail) { cout<<"Label 6 :                           "<<label6<<endl; }
   if (detail) { cout<<"Ratio Label 1 :                     "<<ratiolabel1<<endl; }
   if (detail) { cout<<"Ratio Label 2 :                     "<<ratiolabel2<<endl; }
   if (detail) { cout<<"Ratio Label 3 :                     "<<ratiolabel3<<endl; }
   if (detail) { cout<<"Ratio Label 4 :                     "<<ratiolabel4<<endl; }
   if (detail) { cout<<"Plots Path :                        "<<plots_path<<endl; }
   if (detail) { cout<<"Prefix :                            "<<prefix<<endl; }
   if (detail) { cout<<"Detail level :                      "<<detail<<endl; }

    if (detail) { cout << "Opening files..." << endl; }
//opens the MC files
    TFile *file1 = new TFile( infile1.c_str() );
    TFile *file2 = new TFile( infile2.c_str() );
    TFile *file3 = new TFile( infile3.c_str() );
    TFile *file4 = new TFile( infile4.c_str() );
    TFile *file5 = new TFile( infile5.c_str() );
    TFile *file6 = new TFile( infile6.c_str() );

    TString load_prefix = "ak5PF_";
    if (isMC) { load_prefix = "ak5Gen_"; }

//plot delta phi unfolding results
    if (detail) { cout << "Plotting Delta Phi Unfolding Results..." << endl; }

    TH1D *data_delta_phi = 0;
    file1->GetObject("ak5PF_delta_phi",data_delta_phi);
    if (data_delta_phi == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    TH1D *sbbb_delta_phi = 0;
    file2->GetObject(load_prefix+"delta_phi",sbbb_delta_phi);
    if (sbbb_delta_phi == 0) { cout << load_prefix << "delta_phi not found!" << endl; return; }
    TH1D *rbbb_delta_phi = 0;
    file3->GetObject("output_true_delta_phi",rbbb_delta_phi);
    if (rbbb_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }
    TH1D *tunf_delta_phi = 0;
    file4->GetObject("output_true_delta_phi",tunf_delta_phi);
    if (tunf_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }
    TH1D *svd_delta_phi = 0;
    file5->GetObject("output_true_delta_phi",svd_delta_phi);
    if (svd_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }
    TH1D *baye_delta_phi = 0;
    file6->GetObject("output_true_delta_phi",baye_delta_phi);
    if (baye_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }

    data_delta_phi->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi, label1, sbbb_delta_phi, label2, rbbb_delta_phi, label3, tunf_delta_phi, label4, svd_delta_phi, label5, baye_delta_phi, label6, plots_path, prefix, "delta_phi", "top_left", detail);
    plot_3histograms(data_delta_phi, label1, sbbb_delta_phi, label2, baye_delta_phi, label6, plots_path, prefix+"delta_phi_simple", "top_left", detail);  
    data_delta_phi->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi, rbbb_delta_phi, ratiolabel1, tunf_delta_phi, ratiolabel2, svd_delta_phi, ratiolabel3, baye_delta_phi, ratiolabel4, plots_path, prefix+"delta_phi_ratios", detail);
    data_delta_phi->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi norm unfolding results
    if (detail) { cout << "Plotting Delta Phi Norm Unfolding Results..." << endl; }

    TH1D *data_delta_phi_norm = 0;
    file1->GetObject("ak5PF_delta_phi_norm",data_delta_phi_norm);
    if (data_delta_phi_norm == 0) { cout << "ak5PF_delta_phi_norm not found!" << endl; return; }
    TH1D *sbbb_delta_phi_norm = 0;
    file2->GetObject(load_prefix+"delta_phi_norm",sbbb_delta_phi_norm);
    if (sbbb_delta_phi_norm == 0) { cout << load_prefix << "delta_phi_norm not found!" << endl; return; }
    TH1D *rbbb_delta_phi_norm = 0;
    file3->GetObject("output_true_delta_phi_norm",rbbb_delta_phi_norm);
    if (rbbb_delta_phi_norm == 0) { cout << "resp_delta_phi_norm not found!" << endl; return; }
    TH1D *tunf_delta_phi_norm = 0;
    file4->GetObject("output_true_delta_phi_norm",tunf_delta_phi_norm);
    if (tunf_delta_phi_norm == 0) { cout << "resp_delta_phi_norm not found!" << endl; return; }
    TH1D *svd_delta_phi_norm = 0;
    file5->GetObject("output_true_delta_phi_norm",svd_delta_phi_norm);
    if (svd_delta_phi_norm == 0) { cout << "resp_delta_phi_norm not found!" << endl; return; }
    TH1D *baye_delta_phi_norm = 0;
    file6->GetObject("output_true_delta_phi_norm",baye_delta_phi_norm);
    if (baye_delta_phi_norm == 0) { cout << "resp_delta_phi_norm not found!" << endl; return; }

    data_delta_phi_norm->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi} [rad^{-1}]");
    plot_six_dist(data_delta_phi_norm, label1, sbbb_delta_phi_norm, label2, rbbb_delta_phi_norm, label3, tunf_delta_phi_norm, label4, svd_delta_phi_norm, label5, baye_delta_phi_norm, label6, plots_path, prefix, "delta_phi_norm", "top_left", detail);
    plot_3histograms(data_delta_phi_norm, label1, sbbb_delta_phi_norm, "Standard Bin-By-Bin", baye_delta_phi_norm, label6, plots_path, prefix+"delta_phi_norm_simple", "top_left", detail);  
    data_delta_phi_norm->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_norm, rbbb_delta_phi_norm, ratiolabel1, tunf_delta_phi_norm, ratiolabel2, svd_delta_phi_norm, ratiolabel3, baye_delta_phi_norm, ratiolabel4, plots_path, prefix+"delta_phi_norm_ratios", detail);
    data_delta_phi->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi} [rad^{-1}]");


//plot delta phi deta1 unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta1 Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta1 = 0;
    file1->GetObject("ak5PF_delta_phi_deta1",data_delta_phi_deta1);
    if (data_delta_phi_deta1 == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta1 = 0;
    file2->GetObject(load_prefix+"delta_phi_deta1",sbbb_delta_phi_deta1);
    if (sbbb_delta_phi_deta1 == 0) { cout << load_prefix << "delta_phi_deta1 not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta1 = 0;
    file3->GetObject("output_true_delta_phi_deta1",rbbb_delta_phi_deta1);
    if (rbbb_delta_phi_deta1 == 0) { cout << "rbbb_delta_phi_deta1 not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta1 = 0;
    file4->GetObject("output_true_delta_phi_deta1",tunf_delta_phi_deta1);
    if (tunf_delta_phi_deta1 == 0) { cout << "tunf_delta_phi_deta1 not found!" << endl; return; }
    TH1D *svd_delta_phi_deta1 = 0;
    file5->GetObject("output_true_delta_phi_deta1",svd_delta_phi_deta1);
    if (svd_delta_phi_deta1 == 0) { cout << "svd_delta_phi_deta1 not found!" << endl; return; }
    TH1D *baye_delta_phi_deta1 = 0;
    file6->GetObject("output_true_delta_phi_deta1",baye_delta_phi_deta1);
    if (baye_delta_phi_deta1 == 0) { cout << "baye_delta_phi_deta1 not found!" << endl; return; }

    data_delta_phi_deta1->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta1, label1, sbbb_delta_phi_deta1, label2, rbbb_delta_phi_deta1, label3, tunf_delta_phi_deta1, label4, svd_delta_phi_deta1, label5, baye_delta_phi_deta1, label6, plots_path, prefix, "delta_phi_deta1", "top_left", detail);
    plot_3histograms(data_delta_phi_deta1, label1, sbbb_delta_phi_deta1, "Standard Bin-By-Bin", baye_delta_phi_deta1, label6, plots_path, prefix+"delta_phi_deta1_simple", "top_left", detail); 
    data_delta_phi_deta1->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio"); 
    plots_unfolding_ratios(sbbb_delta_phi_deta1, rbbb_delta_phi_deta1, ratiolabel1, tunf_delta_phi_deta1, ratiolabel2, svd_delta_phi_deta1, ratiolabel3, baye_delta_phi_deta1, ratiolabel4, plots_path, prefix+"delta_phi_deta1_ratios", detail);
    data_delta_phi_deta1->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta1 norm unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta1 Norm Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta1_norm = 0;
    file1->GetObject("ak5PF_delta_phi_deta1_norm",data_delta_phi_deta1_norm);
    if (data_delta_phi_deta1_norm == 0) { cout << "ak5PF_delta_phi_deta1_norm not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta1_norm = 0;
    file2->GetObject(load_prefix+"delta_phi_deta1_norm",sbbb_delta_phi_deta1_norm);
    if (sbbb_delta_phi_deta1_norm == 0) { cout << load_prefix << "delta_phi_deta1_norm not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta1_norm = 0;
    file3->GetObject("output_true_delta_phi_deta1_norm",rbbb_delta_phi_deta1_norm);
    if (rbbb_delta_phi_deta1_norm == 0) { cout << "rbbb_delta_phi_deta1_norm not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta1_norm = 0;
    file4->GetObject("output_true_delta_phi_deta1_norm",tunf_delta_phi_deta1_norm);
    if (tunf_delta_phi_deta1_norm == 0) { cout << "tunf_delta_phi_deta1_norm not found!" << endl; return; }
    TH1D *svd_delta_phi_deta1_norm = 0;
    file5->GetObject("output_true_delta_phi_deta1_norm",svd_delta_phi_deta1_norm);
    if (svd_delta_phi_deta1_norm == 0) { cout << "svd_delta_phi_deta1_norm not found!" << endl; return; }
    TH1D *baye_delta_phi_deta1_norm = 0;
    file6->GetObject("output_true_delta_phi_deta1_norm",baye_delta_phi_deta1_norm);
    if (baye_delta_phi_deta1_norm == 0) { cout << "baye_delta_phi_deta1_norm not found!" << endl; return; }

    data_delta_phi_deta1_norm->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^2}{d#Delta#eta d#Delta#phi} [rad^{-1}]");
    plot_six_dist(data_delta_phi_deta1_norm, label1, sbbb_delta_phi_deta1_norm, label2, rbbb_delta_phi_deta1_norm, label3, tunf_delta_phi_deta1_norm, label4, svd_delta_phi_deta1_norm, label5, baye_delta_phi_deta1_norm, label6, plots_path, prefix, "delta_phi_deta1_norm", "top_left", detail);
    plot_3histograms(data_delta_phi_deta1_norm, label1, sbbb_delta_phi_deta1_norm, label2, baye_delta_phi_deta1_norm, label6, plots_path, prefix+"delta_phi_deta1_norm_simple", "top_left", detail); 
    data_delta_phi_deta1_norm->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio"); 
    plots_unfolding_ratios(sbbb_delta_phi_deta1_norm, rbbb_delta_phi_deta1_norm, ratiolabel1, tunf_delta_phi_deta1_norm, ratiolabel2, svd_delta_phi_deta1_norm, ratiolabel3, baye_delta_phi_deta1_norm, ratiolabel4, plots_path, prefix+"delta_phi_deta1_norm_ratios", detail);
    data_delta_phi_deta1->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^2}{d#Delta#eta d#Delta#phi} [rad^{-1}]");


//plot delta phi deta2 unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta2 Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta2 = 0;
    file1->GetObject("ak5PF_delta_phi_deta2",data_delta_phi_deta2);
    if (data_delta_phi_deta2 == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta2 = 0;
    file2->GetObject(load_prefix+"delta_phi_deta2",sbbb_delta_phi_deta2);
    if (sbbb_delta_phi_deta2 == 0) { cout << load_prefix << "delta_phi_deta2 not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta2 = 0;
    file3->GetObject("output_true_delta_phi_deta2",rbbb_delta_phi_deta2);
    if (rbbb_delta_phi_deta2 == 0) { cout << "rbbb_delta_phi_deta2 not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta2 = 0;
    file4->GetObject("output_true_delta_phi_deta2",tunf_delta_phi_deta2);
    if (tunf_delta_phi_deta2 == 0) { cout << "tunf_delta_phi_deta2 not found!" << endl; return; }
    TH1D *svd_delta_phi_deta2 = 0;
    file5->GetObject("output_true_delta_phi_deta2",svd_delta_phi_deta2);
    if (svd_delta_phi_deta2 == 0) { cout << "svd_delta_phi_deta2 not found!" << endl; return; }
    TH1D *baye_delta_phi_deta2 = 0;
    file6->GetObject("output_true_delta_phi_deta2",baye_delta_phi_deta2);
    if (baye_delta_phi_deta2 == 0) { cout << "baye_delta_phi_deta2 not found!" << endl; return; }

    data_delta_phi_deta2->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta2, label1, sbbb_delta_phi_deta2, label2, rbbb_delta_phi_deta2, label3, tunf_delta_phi_deta2, label4, svd_delta_phi_deta2, label5, baye_delta_phi_deta2, label6, plots_path, prefix, "delta_phi_deta2", "top_left", detail);    
    plot_3histograms(data_delta_phi_deta2, label1, sbbb_delta_phi_deta2, label2, baye_delta_phi_deta2, label6, plots_path, prefix+"delta_phi_deta2_simple", "top_left", detail);
    data_delta_phi_deta2->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta2, rbbb_delta_phi_deta2, ratiolabel1, tunf_delta_phi_deta2, ratiolabel2, svd_delta_phi_deta2, ratiolabel3, baye_delta_phi_deta2, ratiolabel4, plots_path, prefix+"delta_phi_deta2_ratios", detail);
    data_delta_phi_deta2->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta3 unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta3 Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta3 = 0;
    file1->GetObject("ak5PF_delta_phi_deta3",data_delta_phi_deta3);
    if (data_delta_phi_deta3 == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta3 = 0;
    file2->GetObject(load_prefix+"delta_phi_deta3",sbbb_delta_phi_deta3);
    if (sbbb_delta_phi_deta3 == 0) { cout << load_prefix << "delta_phi_deta3 not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta3 = 0;
    file3->GetObject("output_true_delta_phi_deta3",rbbb_delta_phi_deta3);
    if (rbbb_delta_phi_deta3 == 0) { cout << "rbbb_delta_phi_deta3 not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta3 = 0;
    file4->GetObject("output_true_delta_phi_deta3",tunf_delta_phi_deta3);
    if (tunf_delta_phi_deta3 == 0) { cout << "tunf_delta_phi_deta3 not found!" << endl; return; }
    TH1D *svd_delta_phi_deta3 = 0;
    file5->GetObject("output_true_delta_phi_deta3",svd_delta_phi_deta3);
    if (svd_delta_phi_deta3 == 0) { cout << "svd_delta_phi_deta3 not found!" << endl; return; }
    TH1D *baye_delta_phi_deta3 = 0;
    file6->GetObject("output_true_delta_phi_deta3",baye_delta_phi_deta3);
    if (baye_delta_phi_deta3 == 0) { cout << "baye_delta_phi_deta3 not found!" << endl; return; }

    data_delta_phi_deta3->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta3, label1, sbbb_delta_phi_deta3, label2, rbbb_delta_phi_deta3, label3, tunf_delta_phi_deta3, label4, svd_delta_phi_deta3, label5, baye_delta_phi_deta3, label6, plots_path, prefix, "delta_phi_deta3", "top_left", detail);
    plot_3histograms(data_delta_phi_deta3, label1, sbbb_delta_phi_deta3, label2, baye_delta_phi_deta3, label6, plots_path, prefix+"delta_phi_deta3_simple", "top_left", detail); 
    data_delta_phi_deta3->SetTitle("#Delta#phi;#Delta#phi [rad];#Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta3, rbbb_delta_phi_deta3, ratiolabel1, tunf_delta_phi_deta3, ratiolabel2, svd_delta_phi_deta3, ratiolabel3, baye_delta_phi_deta3, ratiolabel4, plots_path, prefix+"delta_phi_deta3_ratios", detail);
    data_delta_phi_deta3->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta4 unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta4 Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta4 = 0;
    file1->GetObject("ak5PF_delta_phi_deta4",data_delta_phi_deta4);
    if (data_delta_phi_deta4 == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta4 = 0;
    file2->GetObject(load_prefix+"delta_phi_deta4",sbbb_delta_phi_deta4);
    if (sbbb_delta_phi_deta4 == 0) { cout << load_prefix << "delta_phi_deta4 not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta4 = 0;
    file3->GetObject("output_true_delta_phi_deta4",rbbb_delta_phi_deta4);
    if (rbbb_delta_phi_deta4 == 0) { cout << "rbbb_delta_phi_deta4 not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta4 = 0;
    file4->GetObject("output_true_delta_phi_deta4",tunf_delta_phi_deta4);
    if (tunf_delta_phi_deta4 == 0) { cout << "tunf_delta_phi_deta4 not found!" << endl; return; }
    TH1D *svd_delta_phi_deta4 = 0;
    file5->GetObject("output_true_delta_phi_deta4",svd_delta_phi_deta4);
    if (svd_delta_phi_deta4 == 0) { cout << "svd_delta_phi_deta4 not found!" << endl; return; }
    TH1D *baye_delta_phi_deta4 = 0;
    file6->GetObject("output_true_delta_phi_deta4",baye_delta_phi_deta4);
    if (baye_delta_phi_deta4 == 0) { cout << "baye_delta_phi_deta4 not found!" << endl; return; }

    data_delta_phi_deta4->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta4, label1, sbbb_delta_phi_deta4, label2, rbbb_delta_phi_deta4, label3, tunf_delta_phi_deta4, label4, svd_delta_phi_deta4, label5, baye_delta_phi_deta4, label6, plots_path, prefix, "delta_phi_deta4", "top_left", detail);
    plot_3histograms(data_delta_phi_deta4, label1, sbbb_delta_phi_deta4, label2, baye_delta_phi_deta4, label6, plots_path, prefix+"delta_phi_deta4_simple", "top_left", detail);    
    data_delta_phi_deta4->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta4, rbbb_delta_phi_deta4, ratiolabel1, tunf_delta_phi_deta4, ratiolabel2, svd_delta_phi_deta4, ratiolabel3, baye_delta_phi_deta4, ratiolabel4, plots_path, prefix+"delta_phi_deta4_ratios", detail);
    data_delta_phi_deta4->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi gap unfolding results
    if (detail) { cout << "Plotting Delta Phi Gap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_gap = 0;
    file1->GetObject("ak5PF_delta_phi_gap",data_delta_phi_gap);
    if (data_delta_phi_gap == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_gap = 0;
    file2->GetObject(load_prefix+"delta_phi_gap",sbbb_delta_phi_gap);
    if (sbbb_delta_phi_gap == 0) { cout << load_prefix << "delta_phi_gap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_gap = 0;
    file3->GetObject("output_true_delta_phi_gap",rbbb_delta_phi_gap);
    if (rbbb_delta_phi_gap == 0) { cout << "rbbb_delta_phi_gap not found!" << endl; return; }
    TH1D *tunf_delta_phi_gap = 0;
    file4->GetObject("output_true_delta_phi_gap",tunf_delta_phi_gap);
    if (tunf_delta_phi_gap == 0) { cout << "tunf_delta_phi_gap not found!" << endl; return; }
    TH1D *svd_delta_phi_gap = 0;
    file5->GetObject("output_true_delta_phi_gap",svd_delta_phi_gap);
    if (svd_delta_phi_gap == 0) { cout << "svd_delta_phi_gap not found!" << endl; return; }
    TH1D *baye_delta_phi_gap = 0;
    file6->GetObject("output_true_delta_phi_gap",baye_delta_phi_gap);
    if (baye_delta_phi_gap == 0) { cout << "baye_delta_phi_gap not found!" << endl; return; }

    data_delta_phi_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_gap, label1, sbbb_delta_phi_gap, label2, rbbb_delta_phi_gap, label3, tunf_delta_phi_gap, label4, svd_delta_phi_gap, label5, baye_delta_phi_gap, label6, plots_path, prefix, "delta_phi_gap", "top_left", detail);
    plot_3histograms(data_delta_phi_gap, label1, sbbb_delta_phi_gap, label2, baye_delta_phi_gap, label6, plots_path, prefix+"delta_phi_gap_simple", "top_left", detail);
    data_delta_phi_gap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_gap, rbbb_delta_phi_gap, ratiolabel1, tunf_delta_phi_gap, ratiolabel2, svd_delta_phi_gap, ratiolabel3, baye_delta_phi_gap, ratiolabel4, plots_path, prefix+"delta_phi_gap_ratios", detail);
    data_delta_phi_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta1 gap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta1 Gap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta1_gap = 0;
    file1->GetObject("ak5PF_delta_phi_deta1_gap",data_delta_phi_deta1_gap);
    if (data_delta_phi_deta1_gap == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta1_gap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta1_gap",sbbb_delta_phi_deta1_gap);
    if (sbbb_delta_phi_deta1_gap == 0) { cout << load_prefix << "delta_phi_deta1_gap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta1_gap = 0;
    file3->GetObject("output_true_delta_phi_deta1_gap",rbbb_delta_phi_deta1_gap);
    if (rbbb_delta_phi_deta1_gap == 0) { cout << "rbbb_delta_phi_deta1_gap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta1_gap = 0;
    file4->GetObject("output_true_delta_phi_deta1_gap",tunf_delta_phi_deta1_gap);
    if (tunf_delta_phi_deta1_gap == 0) { cout << "tunf_delta_phi_deta1_gap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta1_gap = 0;
    file5->GetObject("output_true_delta_phi_deta1_gap",svd_delta_phi_deta1_gap);
    if (svd_delta_phi_deta1_gap == 0) { cout << "svd_delta_phi_deta1_gap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta1_gap = 0;
    file6->GetObject("output_true_delta_phi_deta1_gap",baye_delta_phi_deta1_gap);
    if (baye_delta_phi_deta1_gap == 0) { cout << "baye_delta_phi_deta1_gap not found!" << endl; return; }

    data_delta_phi_deta1_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta1_gap, label1, sbbb_delta_phi_deta1_gap, label2, rbbb_delta_phi_deta1_gap, label3, tunf_delta_phi_deta1_gap, label4, svd_delta_phi_deta1_gap, label5, baye_delta_phi_deta1_gap, label6, plots_path, prefix, "delta_phi_deta1_gap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta1_gap, label1, sbbb_delta_phi_deta1_gap, label2, baye_delta_phi_deta1_gap, label6, plots_path, prefix+"delta_phi_deta1_gap_simple", "top_left", detail);
    data_delta_phi_deta1_gap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta1_gap, rbbb_delta_phi_deta1_gap, ratiolabel1, tunf_delta_phi_deta1_gap, ratiolabel2, svd_delta_phi_deta1_gap, ratiolabel3, baye_delta_phi_deta1_gap, ratiolabel4, plots_path, prefix+"delta_phi_deta1_gap_ratios", detail);
    data_delta_phi_deta1_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta2 gap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta2 Gap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta2_gap = 0;
    file1->GetObject("ak5PF_delta_phi_deta2_gap",data_delta_phi_deta2_gap);
    if (data_delta_phi_deta2_gap == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta2_gap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta2_gap",sbbb_delta_phi_deta2_gap);
    if (sbbb_delta_phi_deta2_gap == 0) { cout << load_prefix << "delta_phi_deta2_gap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta2_gap = 0;
    file3->GetObject("output_true_delta_phi_deta2_gap",rbbb_delta_phi_deta2_gap);
    if (rbbb_delta_phi_deta2_gap == 0) { cout << "rbbb_delta_phi_deta2_gap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta2_gap = 0;
    file4->GetObject("output_true_delta_phi_deta2_gap",tunf_delta_phi_deta2_gap);
    if (tunf_delta_phi_deta2_gap == 0) { cout << "tunf_delta_phi_deta2_gap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta2_gap = 0;
    file5->GetObject("output_true_delta_phi_deta2_gap",svd_delta_phi_deta2_gap);
    if (svd_delta_phi_deta2_gap == 0) { cout << "svd_delta_phi_deta2_gap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta2_gap = 0;
    file6->GetObject("output_true_delta_phi_deta2_gap",baye_delta_phi_deta2_gap);
    if (baye_delta_phi_deta2_gap == 0) { cout << "baye_delta_phi_deta2_gap not found!" << endl; return; }

    data_delta_phi_deta2_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta2_gap, label1, sbbb_delta_phi_deta2_gap, label2, rbbb_delta_phi_deta2_gap, label3, tunf_delta_phi_deta2_gap, label4, svd_delta_phi_deta2_gap, label5, baye_delta_phi_deta2_gap, label6, plots_path, prefix, "delta_phi_deta2_gap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta2_gap, label1, sbbb_delta_phi_deta2_gap, label2, baye_delta_phi_deta2_gap, label6, plots_path, prefix+"delta_phi_deta2_gap_simple", "top_left", detail); 
    data_delta_phi_deta2_gap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta2_gap, rbbb_delta_phi_deta2_gap, ratiolabel1, tunf_delta_phi_deta2_gap, ratiolabel2, svd_delta_phi_deta2_gap, ratiolabel3, baye_delta_phi_deta2_gap, ratiolabel4, plots_path, prefix+"delta_phi_deta2_gap_ratios", detail);
    data_delta_phi_deta2_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta3 gap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta3 Gap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta3_gap = 0;
    file1->GetObject("ak5PF_delta_phi_deta3_gap",data_delta_phi_deta3_gap);
    if (data_delta_phi_deta3_gap == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta3_gap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta3_gap",sbbb_delta_phi_deta3_gap);
    if (sbbb_delta_phi_deta3_gap == 0) { cout << load_prefix << "delta_phi_deta3_gap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta3_gap = 0;
    file3->GetObject("output_true_delta_phi_deta3_gap",rbbb_delta_phi_deta3_gap);
    if (rbbb_delta_phi_deta3_gap == 0) { cout << "rbbb_delta_phi_deta3_gap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta3_gap = 0;
    file4->GetObject("output_true_delta_phi_deta3_gap",tunf_delta_phi_deta3_gap);
    if (tunf_delta_phi_deta3_gap == 0) { cout << "tunf_delta_phi_deta3_gap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta3_gap = 0;
    file5->GetObject("output_true_delta_phi_deta3_gap",svd_delta_phi_deta3_gap);
    if (svd_delta_phi_deta3_gap == 0) { cout << "svd_delta_phi_deta3_gap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta3_gap = 0;
    file6->GetObject("output_true_delta_phi_deta3_gap",baye_delta_phi_deta3_gap);
    if (baye_delta_phi_deta3_gap == 0) { cout << "baye_delta_phi_deta3_gap not found!" << endl; return; }

    data_delta_phi_deta3_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta3_gap, label1, sbbb_delta_phi_deta3_gap, label2, rbbb_delta_phi_deta3_gap, label3, tunf_delta_phi_deta3_gap, label4, svd_delta_phi_deta3_gap, label5, baye_delta_phi_deta3_gap, label6, plots_path, prefix, "delta_phi_deta3_gap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta3_gap, label1, sbbb_delta_phi_deta3_gap, label2, baye_delta_phi_deta3_gap, label6, plots_path, prefix+"delta_phi_deta3_gap_simple", "top_left", detail);
    data_delta_phi_deta3_gap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta3_gap, rbbb_delta_phi_deta3_gap, ratiolabel1, tunf_delta_phi_deta3_gap, ratiolabel2, svd_delta_phi_deta3_gap, ratiolabel3, baye_delta_phi_deta3_gap, ratiolabel4, plots_path, prefix+"delta_phi_deta3_gap_ratios", detail);
    data_delta_phi_deta3_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta4 gap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta4 Gap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta4_gap = 0;
    file1->GetObject("ak5PF_delta_phi_deta4_gap",data_delta_phi_deta4_gap);
    if (data_delta_phi_deta4_gap == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta4_gap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta4_gap",sbbb_delta_phi_deta4_gap);
    if (sbbb_delta_phi_deta4_gap == 0) { cout << load_prefix << "delta_phi_deta4_gap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta4_gap = 0;
    file3->GetObject("output_true_delta_phi_deta4_gap",rbbb_delta_phi_deta4_gap);
    if (rbbb_delta_phi_deta4_gap == 0) { cout << "rbbb_delta_phi_deta4_gap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta4_gap = 0;
    file4->GetObject("output_true_delta_phi_deta4_gap",tunf_delta_phi_deta4_gap);
    if (tunf_delta_phi_deta4_gap == 0) { cout << "tunf_delta_phi_deta4_gap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta4_gap = 0;
    file5->GetObject("output_true_delta_phi_deta4_gap",svd_delta_phi_deta4_gap);
    if (svd_delta_phi_deta4_gap == 0) { cout << "svd_delta_phi_deta4_gap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta4_gap = 0;
    file6->GetObject("output_true_delta_phi_deta4_gap",baye_delta_phi_deta4_gap);
    if (baye_delta_phi_deta4_gap == 0) { cout << "baye_delta_phi_deta4_gap not found!" << endl; return; }

    data_delta_phi_deta4_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta4_gap, label1, sbbb_delta_phi_deta4_gap, label2, rbbb_delta_phi_deta4_gap, label3, tunf_delta_phi_deta4_gap, label4, svd_delta_phi_deta4_gap, label5, baye_delta_phi_deta4_gap, label6, plots_path, prefix, "delta_phi_deta4_gap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta4_gap, label1, sbbb_delta_phi_deta4_gap, label2, baye_delta_phi_deta4_gap, label6, plots_path, prefix+"delta_phi_deta4_gap_simple", "top_left", detail);
    data_delta_phi_deta4_gap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta4_gap, rbbb_delta_phi_deta4_gap, ratiolabel1, tunf_delta_phi_deta4_gap, ratiolabel2, svd_delta_phi_deta4_gap, ratiolabel3, baye_delta_phi_deta4_gap, ratiolabel4, plots_path, prefix+"delta_phi_deta4_gap_ratios", detail);
    data_delta_phi_deta4_gap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi nogap unfolding results
    if (detail) { cout << "Plotting Delta Phi Nogap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_nogap = 0;
    file1->GetObject("ak5PF_delta_phi_nogap",data_delta_phi_nogap);
    if (data_delta_phi_nogap == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_nogap = 0;
    file2->GetObject(load_prefix+"delta_phi_nogap",sbbb_delta_phi_nogap);
    if (sbbb_delta_phi_nogap == 0) { cout << load_prefix << "delta_phi_nogap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_nogap = 0;
    file3->GetObject("output_true_delta_phi_nogap",rbbb_delta_phi_nogap);
    if (rbbb_delta_phi_nogap == 0) { cout << "rbbb_delta_phi_nogap not found!" << endl; return; }
    TH1D *tunf_delta_phi_nogap = 0;
    file4->GetObject("output_true_delta_phi_nogap",tunf_delta_phi_nogap);
    if (tunf_delta_phi_nogap == 0) { cout << "tunf_delta_phi_nogap not found!" << endl; return; }
    TH1D *svd_delta_phi_nogap = 0;
    file5->GetObject("output_true_delta_phi_nogap",svd_delta_phi_nogap);
    if (svd_delta_phi_nogap == 0) { cout << "svd_delta_phi_nogap not found!" << endl; return; }
    TH1D *baye_delta_phi_nogap = 0;
    file6->GetObject("output_true_delta_phi_nogap",baye_delta_phi_nogap);
    if (baye_delta_phi_nogap == 0) { cout << "baye_delta_phi_nogap not found!" << endl; return; }

    data_delta_phi_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_nogap, label1, sbbb_delta_phi_nogap, label2, rbbb_delta_phi_nogap, label3, tunf_delta_phi_nogap, label4, svd_delta_phi_nogap, label5, baye_delta_phi_nogap, label6, plots_path, prefix, "delta_phi_nogap", "top_left", detail);
    plot_3histograms(data_delta_phi_nogap, label1, sbbb_delta_phi_nogap, label2, baye_delta_phi_nogap, label6, plots_path, prefix+"delta_phi_nogap_simple", "top_left", detail);
    data_delta_phi_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_nogap, rbbb_delta_phi_nogap, ratiolabel1, tunf_delta_phi_nogap, ratiolabel2, svd_delta_phi_nogap, ratiolabel3, baye_delta_phi_nogap, ratiolabel4, plots_path, prefix+"delta_phi_nogap_ratios", detail);
    data_delta_phi_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta1 nogap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta1 Nogap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta1_nogap = 0;
    file1->GetObject("ak5PF_delta_phi_deta1_nogap",data_delta_phi_deta1_nogap);
    if (data_delta_phi_deta1_nogap == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta1_nogap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta1_nogap",sbbb_delta_phi_deta1_nogap);
    if (sbbb_delta_phi_deta1_nogap == 0) { cout << load_prefix << "delta_phi_deta1_nogap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta1_nogap = 0;
    file3->GetObject("output_true_delta_phi_deta1_nogap",rbbb_delta_phi_deta1_nogap);
    if (rbbb_delta_phi_deta1_nogap == 0) { cout << "rbbb_delta_phi_deta1_nogap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta1_nogap = 0;
    file4->GetObject("output_true_delta_phi_deta1_nogap",tunf_delta_phi_deta1_nogap);
    if (tunf_delta_phi_deta1_nogap == 0) { cout << "tunf_delta_phi_deta1_nogap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta1_nogap = 0;
    file5->GetObject("output_true_delta_phi_deta1_nogap",svd_delta_phi_deta1_nogap);
    if (svd_delta_phi_deta1_nogap == 0) { cout << "svd_delta_phi_deta1_nogap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta1_nogap = 0;
    file6->GetObject("output_true_delta_phi_deta1_nogap",baye_delta_phi_deta1_nogap);
    if (baye_delta_phi_deta1_nogap == 0) { cout << "baye_delta_phi_deta1_nogap not found!" << endl; return; }

    data_delta_phi_deta1_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta1_nogap, label1, sbbb_delta_phi_deta1_nogap, label2, rbbb_delta_phi_deta1_nogap, label3, tunf_delta_phi_deta1_nogap, label4, svd_delta_phi_deta1_nogap, label5, baye_delta_phi_deta1_nogap, label6, plots_path, prefix, "delta_phi_deta1_nogap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta1_nogap, label1, sbbb_delta_phi_deta1_nogap, label2, baye_delta_phi_deta1_nogap, label6, plots_path, prefix+"delta_phi_deta1_nogap_simple", "top_left", detail);
    data_delta_phi_deta1_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta1_nogap, rbbb_delta_phi_deta1_nogap, ratiolabel1, tunf_delta_phi_deta1_nogap, ratiolabel2, svd_delta_phi_deta1_nogap, ratiolabel3, baye_delta_phi_deta1_nogap, ratiolabel4, plots_path, prefix+"delta_phi_deta1_nogap_ratios", detail);
    data_delta_phi_deta1_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta2 nogap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta2 Nogap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta2_nogap = 0;
    file1->GetObject("ak5PF_delta_phi_deta2_nogap",data_delta_phi_deta2_nogap);
    if (data_delta_phi_deta2_nogap == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta2_nogap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta2_nogap",sbbb_delta_phi_deta2_nogap);
    if (sbbb_delta_phi_deta2_nogap == 0) { cout << load_prefix << "delta_phi_deta2_nogap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta2_nogap = 0;
    file3->GetObject("output_true_delta_phi_deta2_nogap",rbbb_delta_phi_deta2_nogap);
    if (rbbb_delta_phi_deta2_nogap == 0) { cout << "rbbb_delta_phi_deta2_nogap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta2_nogap = 0;
    file4->GetObject("output_true_delta_phi_deta2_nogap",tunf_delta_phi_deta2_nogap);
    if (tunf_delta_phi_deta2_nogap == 0) { cout << "tunf_delta_phi_deta2_nogap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta2_nogap = 0;
    file5->GetObject("output_true_delta_phi_deta2_nogap",svd_delta_phi_deta2_nogap);
    if (svd_delta_phi_deta2_nogap == 0) { cout << "svd_delta_phi_deta2_nogap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta2_nogap = 0;
    file6->GetObject("output_true_delta_phi_deta2_nogap",baye_delta_phi_deta2_nogap);
    if (baye_delta_phi_deta2_nogap == 0) { cout << "baye_delta_phi_deta2_nogap not found!" << endl; return; }

    data_delta_phi_deta2_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta2_nogap, label1, sbbb_delta_phi_deta2_nogap, label2, rbbb_delta_phi_deta2_nogap, label3, tunf_delta_phi_deta2_nogap, label4, svd_delta_phi_deta2_nogap, label5, baye_delta_phi_deta2_nogap, label6, plots_path, prefix, "delta_phi_deta2_nogap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta2_nogap, label1, sbbb_delta_phi_deta2_nogap, label2, baye_delta_phi_deta2_nogap, label6, plots_path, prefix+"delta_phi_deta2_nogap_simple", "top_left", detail);
    data_delta_phi_deta2_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta2_nogap, rbbb_delta_phi_deta2_nogap, ratiolabel1, tunf_delta_phi_deta2_nogap, ratiolabel2, svd_delta_phi_deta2_nogap, ratiolabel3, baye_delta_phi_deta2_nogap, ratiolabel4, plots_path, prefix+"delta_phi_deta2_nogap_ratios", detail);
    data_delta_phi_deta2_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta3 nogap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta3 Nogap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta3_nogap = 0;
    file1->GetObject("ak5PF_delta_phi_deta3_nogap",data_delta_phi_deta3_nogap);
    if (data_delta_phi_deta3_nogap == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta3_nogap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta3_nogap",sbbb_delta_phi_deta3_nogap);
    if (sbbb_delta_phi_deta3_nogap == 0) { cout << load_prefix << "delta_phi_deta3_nogap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta3_nogap = 0;
    file3->GetObject("output_true_delta_phi_deta3_nogap",rbbb_delta_phi_deta3_nogap);
    if (rbbb_delta_phi_deta3_nogap == 0) { cout << "rbbb_delta_phi_deta3_nogap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta3_nogap = 0;
    file4->GetObject("output_true_delta_phi_deta3_nogap",tunf_delta_phi_deta3_nogap);
    if (tunf_delta_phi_deta3_nogap == 0) { cout << "tunf_delta_phi_deta3_nogap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta3_nogap = 0;
    file5->GetObject("output_true_delta_phi_deta3_nogap",svd_delta_phi_deta3_nogap);
    if (svd_delta_phi_deta3_nogap == 0) { cout << "svd_delta_phi_deta3_nogap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta3_nogap = 0;
    file6->GetObject("output_true_delta_phi_deta3_nogap",baye_delta_phi_deta3_nogap);
    if (baye_delta_phi_deta3_nogap == 0) { cout << "baye_delta_phi_deta3_nogap not found!" << endl; return; }

    data_delta_phi_deta3_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta3_nogap, label1, sbbb_delta_phi_deta3_nogap, label2, rbbb_delta_phi_deta3_nogap, label3, tunf_delta_phi_deta3_nogap, label4, svd_delta_phi_deta3_nogap, label5, baye_delta_phi_deta3_nogap, label6, plots_path, prefix, "delta_phi_deta3_nogap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta3_nogap, label1, sbbb_delta_phi_deta3_nogap, label2, baye_delta_phi_deta3_nogap, label6, plots_path, prefix+"delta_phi_deta3_nogap_simple", "top_left", detail);
    data_delta_phi_deta3_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta3_nogap, rbbb_delta_phi_deta3_nogap, ratiolabel1, tunf_delta_phi_deta3_nogap, ratiolabel2, svd_delta_phi_deta3_nogap, ratiolabel3, baye_delta_phi_deta3_nogap, ratiolabel4, plots_path, prefix+"delta_phi_deta3_nogap_ratios", detail);
    data_delta_phi_deta3_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot delta phi deta4 nogap unfolding results
    if (detail) { cout << "Plotting Delta Phi Deta4 Nogap Unfolding Results..." << endl; }

    TH1D *data_delta_phi_deta4_nogap = 0;
    file1->GetObject("ak5PF_delta_phi_deta4_nogap",data_delta_phi_deta4_nogap);
    if (data_delta_phi_deta4_nogap == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    TH1D *sbbb_delta_phi_deta4_nogap = 0;
    file2->GetObject(load_prefix+"delta_phi_deta4_nogap",sbbb_delta_phi_deta4_nogap);
    if (sbbb_delta_phi_deta4_nogap == 0) { cout << load_prefix << "delta_phi_deta4_nogap not found!" << endl; return; }
    TH1D *rbbb_delta_phi_deta4_nogap = 0;
    file3->GetObject("output_true_delta_phi_deta4_nogap",rbbb_delta_phi_deta4_nogap);
    if (rbbb_delta_phi_deta4_nogap == 0) { cout << "rbbb_delta_phi_deta4_nogap not found!" << endl; return; }
    TH1D *tunf_delta_phi_deta4_nogap = 0;
    file4->GetObject("output_true_delta_phi_deta4_nogap",tunf_delta_phi_deta4_nogap);
    if (tunf_delta_phi_deta4_nogap == 0) { cout << "tunf_delta_phi_deta4_nogap not found!" << endl; return; }
    TH1D *svd_delta_phi_deta4_nogap = 0;
    file5->GetObject("output_true_delta_phi_deta4_nogap",svd_delta_phi_deta4_nogap);
    if (svd_delta_phi_deta4_nogap == 0) { cout << "svd_delta_phi_deta4_nogap not found!" << endl; return; }
    TH1D *baye_delta_phi_deta4_nogap = 0;
    file6->GetObject("output_true_delta_phi_deta4_nogap",baye_delta_phi_deta4_nogap);
    if (baye_delta_phi_deta4_nogap == 0) { cout << "baye_delta_phi_deta4_nogap not found!" << endl; return; }

    data_delta_phi_deta4_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");
    plot_six_dist(data_delta_phi_deta4_nogap, label1, sbbb_delta_phi_deta4_nogap, label2, rbbb_delta_phi_deta4_nogap, label3, tunf_delta_phi_deta4_nogap, label4, svd_delta_phi_deta4_nogap, label5, baye_delta_phi_deta4_nogap, label6, plots_path, prefix, "delta_phi_deta4_nogap", "top_left", detail);
    plot_3histograms(data_delta_phi_deta4_nogap, label1, sbbb_delta_phi_deta4_nogap, label2, baye_delta_phi_deta4_nogap, label6, plots_path, prefix+"delta_phi_deta4_nogap_simple", "top_left", detail);
    data_delta_phi_deta4_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];Ratio");
    plots_unfolding_ratios(sbbb_delta_phi_deta4_nogap, rbbb_delta_phi_deta4_nogap, ratiolabel1, tunf_delta_phi_deta4_nogap, ratiolabel2, svd_delta_phi_deta4_nogap, ratiolabel3, baye_delta_phi_deta4_nogap, ratiolabel4, plots_path, prefix+"delta_phi_deta4_nogap_ratios", detail);
    data_delta_phi_deta4_nogap->SetTitle("#Delta#phi;#Delta#phi [rad];#frac{d#sigma^2}{d#Delta#eta d#Delta#phi} [#frac{pb}{rad}]");


//plot leading pt inside gap unfolding results
    if (detail) { cout << "Plotting Leading pT Inside Gap Unfolding Results..." << endl; }

    TH1D *data_leading_pt_inside_gap = 0;
    file1->GetObject("ak5PF_leading_pt_inside_gap",data_leading_pt_inside_gap);
    if (data_leading_pt_inside_gap == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    TH1D *sbbb_leading_pt_inside_gap = 0;
    file2->GetObject(load_prefix+"leading_pt_inside_gap",sbbb_leading_pt_inside_gap);
    if (sbbb_leading_pt_inside_gap == 0) { cout << load_prefix << "leading_pt_inside_gap not found!" << endl; return; }
    TH1D *rbbb_leading_pt_inside_gap = 0;
    file3->GetObject("output_true_leading_pt_inside_gap",rbbb_leading_pt_inside_gap);
    if (rbbb_leading_pt_inside_gap == 0) { cout << "rbbb_leading_pt_inside_gap not found!" << endl; return; }
    TH1D *tunf_leading_pt_inside_gap = 0;
    file4->GetObject("output_true_leading_pt_inside_gap",tunf_leading_pt_inside_gap);
    if (tunf_leading_pt_inside_gap == 0) { cout << "tunf_leading_pt_inside_gap not found!" << endl; return; }
    TH1D *svd_leading_pt_inside_gap = 0;
    file5->GetObject("output_true_leading_pt_inside_gap",svd_leading_pt_inside_gap);
    if (svd_leading_pt_inside_gap == 0) { cout << "svd_leading_pt_inside_gap not found!" << endl; return; }
    TH1D *baye_leading_pt_inside_gap = 0;
    file6->GetObject("output_true_leading_pt_inside_gap",baye_leading_pt_inside_gap);
    if (baye_leading_pt_inside_gap == 0) { cout << "baye_leading_pt_inside_gap not found!" << endl; return; }

    data_leading_pt_inside_gap->SetTitle("p_{T}^{inside};p_{T}^{inside} [GeV];#frac{d#sigma}{dp_{T}^{inside}} [#frac{pb}{GeV}]");
    plot_six_dist(data_leading_pt_inside_gap, label1, sbbb_leading_pt_inside_gap, label2, rbbb_leading_pt_inside_gap, label3, tunf_leading_pt_inside_gap, label4, svd_leading_pt_inside_gap, label5, baye_leading_pt_inside_gap, label6, plots_path, prefix, "leading_pt_inside_gap", "top_right", detail); 
    plot_3histograms(data_leading_pt_inside_gap, label1, sbbb_leading_pt_inside_gap, label2, baye_leading_pt_inside_gap, label6, plots_path, prefix+"leading_pt_inside_gap_simple", "top_right", detail);
    data_leading_pt_inside_gap->SetTitle("p_{T}^{inside};p_{T}^{inside} [GeV];Ratio");
    plots_unfolding_ratios(sbbb_leading_pt_inside_gap, rbbb_leading_pt_inside_gap, ratiolabel1,  tunf_leading_pt_inside_gap, ratiolabel2, svd_leading_pt_inside_gap, ratiolabel3, baye_leading_pt_inside_gap, ratiolabel4, plots_path, prefix+"leading_pt_inside_gap_ratios", detail);
    data_leading_pt_inside_gap->SetTitle("p_{T}^{inside};p_{T}^{inside} [GeV];#frac{d#sigma}{dp_{T}^{inside}} [#frac{pb}{GeV}]");

//plot leading eta star inside gap unfolding results
    if (detail) { cout << "Plotting Leading Eta* Inside Gap Unfolding Results..." << endl; }

    TH1D *data_leading_eta_star_inside_gap = 0;
    file1->GetObject("ak5PF_leading_eta_star_inside_gap",data_leading_eta_star_inside_gap);
    if (data_leading_eta_star_inside_gap == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    TH1D *sbbb_leading_eta_star_inside_gap = 0;
    file2->GetObject(load_prefix+"leading_eta_star_inside_gap",sbbb_leading_eta_star_inside_gap);
    if (sbbb_leading_eta_star_inside_gap == 0) { cout << load_prefix << "leading_eta_star_inside_gap not found!" << endl; return; }
    TH1D *rbbb_leading_eta_star_inside_gap = 0;
    file3->GetObject("output_true_leading_eta_star_inside_gap",rbbb_leading_eta_star_inside_gap);
    if (rbbb_leading_eta_star_inside_gap == 0) { cout << "rbbb_leading_eta_star_inside_gap not found!" << endl; return; }
    TH1D *tunf_leading_eta_star_inside_gap = 0;
    file4->GetObject("output_true_leading_eta_star_inside_gap",tunf_leading_eta_star_inside_gap);
    if (tunf_leading_eta_star_inside_gap == 0) { cout << "tunf_leading_eta_star_inside_gap not found!" << endl; return; }
    TH1D *svd_leading_eta_star_inside_gap = 0;
    file5->GetObject("output_true_leading_eta_star_inside_gap",svd_leading_eta_star_inside_gap);
    if (svd_leading_eta_star_inside_gap == 0) { cout << "svd_leading_eta_star_inside_gap not found!" << endl; return; }
    TH1D *baye_leading_eta_star_inside_gap = 0;
    file6->GetObject("output_true_leading_eta_star_inside_gap",baye_leading_eta_star_inside_gap);
    if (baye_leading_eta_star_inside_gap == 0) { cout << "baye_leading_eta_star_inside_gap not found!" << endl; return; }

    data_leading_eta_star_inside_gap->SetTitle("#eta*;#eta*;#frac{d#sigma}{d#eta*} [pb]");
    plot_six_dist(data_leading_eta_star_inside_gap, label1, sbbb_leading_eta_star_inside_gap, label2, rbbb_leading_eta_star_inside_gap, label3, tunf_leading_eta_star_inside_gap, label4, svd_leading_eta_star_inside_gap, label5, baye_leading_eta_star_inside_gap, label6, plots_path, prefix, "leading_eta_star_inside_gap", "bottom_middle", detail);
    plot_3histograms(data_leading_eta_star_inside_gap, label1, sbbb_leading_eta_star_inside_gap, label2, baye_leading_eta_star_inside_gap, label6, plots_path, prefix+"leading_eta_star_inside_gap_simple", "bottom_middle", detail);
    data_leading_eta_star_inside_gap->SetTitle("#eta*;#eta*;Ratio");
    plots_unfolding_ratios(sbbb_leading_eta_star_inside_gap, rbbb_leading_eta_star_inside_gap, ratiolabel1, tunf_leading_eta_star_inside_gap, ratiolabel2, svd_leading_eta_star_inside_gap, ratiolabel3, baye_leading_eta_star_inside_gap, ratiolabel4, plots_path, prefix+"leading_eta_star_inside_gap_ratios", detail);
    data_leading_eta_star_inside_gap->SetTitle("#eta*;#eta*;#frac{d#sigma}{d#eta*} [pb]");

//plot leading pt outside gap unfolding results
    if (detail) { cout << "Plotting Leading pT Outside Gap Unfolding Results..." << endl; }

    TH1D *data_leading_pt_outside_gap = 0;
    file1->GetObject("ak5PF_leading_pt_outside_gap",data_leading_pt_outside_gap);
    if (data_leading_pt_outside_gap == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    TH1D *sbbb_leading_pt_outside_gap = 0;
    file2->GetObject(load_prefix+"leading_pt_outside_gap",sbbb_leading_pt_outside_gap);
    if (sbbb_leading_pt_outside_gap == 0) { cout << load_prefix << "leading_pt_outside_gap not found!" << endl; return; }
    TH1D *rbbb_leading_pt_outside_gap = 0;
    file3->GetObject("output_true_leading_pt_outside_gap",rbbb_leading_pt_outside_gap);
    if (rbbb_leading_pt_outside_gap == 0) { cout << "rbbb_leading_pt_outside_gap not found!" << endl; return; }
    TH1D *tunf_leading_pt_outside_gap = 0;
    file4->GetObject("output_true_leading_pt_outside_gap",tunf_leading_pt_outside_gap);
    if (tunf_leading_pt_outside_gap == 0) { cout << "tunf_leading_pt_outside_gap not found!" << endl; return; }
    TH1D *svd_leading_pt_outside_gap = 0;
    file5->GetObject("output_true_leading_pt_outside_gap",svd_leading_pt_outside_gap);
    if (svd_leading_pt_outside_gap == 0) { cout << "svd_leading_pt_outside_gap not found!" << endl; return; }
    TH1D *baye_leading_pt_outside_gap = 0;
    file6->GetObject("output_true_leading_pt_outside_gap",baye_leading_pt_outside_gap);
    if (baye_leading_pt_outside_gap == 0) { cout << "baye_leading_pt_outside_gap not found!" << endl; return; }

    data_leading_pt_outside_gap->SetTitle("p_{T}^{outside};p_{T}^{outside} [GeV];#frac{d#sigma}{dp_{T}^{outside}} [#frac{pb}{GeV}]");
    plot_six_dist(data_leading_pt_outside_gap, label1, sbbb_leading_pt_outside_gap, label2, rbbb_leading_pt_outside_gap, label3, tunf_leading_pt_outside_gap, label4, svd_leading_pt_outside_gap, label5, baye_leading_pt_outside_gap, label6, plots_path, prefix, "leading_pt_outside_gap", "top_right", detail);    
    plot_3histograms(data_leading_pt_outside_gap, label1, sbbb_leading_pt_outside_gap, label2,  baye_leading_pt_outside_gap, label6, plots_path, prefix+"leading_pt_outside_gap_simple", "top_right", detail);
    data_leading_pt_outside_gap->SetTitle("p_{T}^{inside};p_{T}^{inside} [GeV];Ratio");   
    plots_unfolding_ratios(sbbb_leading_pt_outside_gap, rbbb_leading_pt_outside_gap, ratiolabel1, tunf_leading_pt_outside_gap, ratiolabel2, svd_leading_pt_outside_gap, ratiolabel3, baye_leading_pt_outside_gap, ratiolabel4, plots_path, prefix+"leading_pt_outside_gap_ratios", detail);
    data_leading_pt_outside_gap->SetTitle("p_{T}^{inside};p_{T}^{inside} [GeV];#frac{d#sigma}{dp_{T}^{inside}} [#frac{pb}{GeV}]");


//plot delta eta outside gap unfolding results
    if (detail) { cout << "Plotting Delta Eta Outside Gap Unfolding Results..." << endl; }

    TH1D *data_delta_eta_outside_gap = 0;
    file1->GetObject("ak5PF_delta_eta_outside_gap",data_delta_eta_outside_gap);
    if (data_delta_eta_outside_gap == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    TH1D *sbbb_delta_eta_outside_gap = 0;
    file2->GetObject(load_prefix+"delta_eta_outside_gap",sbbb_delta_eta_outside_gap);
    if (sbbb_delta_eta_outside_gap == 0) { cout << load_prefix << "delta_eta_outside_gap not found!" << endl; return; }
    TH1D *rbbb_delta_eta_outside_gap = 0;
    file3->GetObject("output_true_delta_eta_outside_gap",rbbb_delta_eta_outside_gap);
    if (rbbb_delta_eta_outside_gap == 0) { cout << "rbbb_delta_eta_outside_gap not found!" << endl; return; }
    TH1D *tunf_delta_eta_outside_gap = 0;
    file4->GetObject("output_true_delta_eta_outside_gap",tunf_delta_eta_outside_gap);
    if (tunf_delta_eta_outside_gap == 0) { cout << "tunf_delta_eta_outside_gap not found!" << endl; return; }
    TH1D *svd_delta_eta_outside_gap = 0;
    file5->GetObject("output_true_delta_eta_outside_gap",svd_delta_eta_outside_gap);
    if (svd_delta_eta_outside_gap == 0) { cout << "svd_delta_eta_outside_gap not found!" << endl; return; }
    TH1D *baye_delta_eta_outside_gap = 0;
    file6->GetObject("output_true_delta_eta_outside_gap",baye_delta_eta_outside_gap);
    if (baye_delta_eta_outside_gap == 0) { cout << "baye_delta_eta_outside_gap not found!" << endl; return; }

    data_delta_eta_outside_gap->SetTitle("#Delta#eta^{out};#Delta#eta^{out};#frac{d#sigma}{d#Delta#eta^{out}} [pb]");
    plot_six_dist(data_delta_eta_outside_gap, label1, sbbb_delta_eta_outside_gap, label2, rbbb_delta_eta_outside_gap, label3, tunf_delta_eta_outside_gap, label4, svd_delta_eta_outside_gap, label5, baye_delta_eta_outside_gap, label6, plots_path, prefix, "delta_eta_outside_gap", "bottom_left", detail);
    plot_3histograms(data_delta_eta_outside_gap, label1, sbbb_delta_eta_outside_gap, label2, baye_delta_eta_outside_gap, label6, plots_path, prefix+"delta_eta_outside_gap_simple", "top_right", detail);
    data_delta_eta_outside_gap->SetTitle("#Delta#eta^{out};#Delta#eta^{out};Ratio");
    plots_unfolding_ratios(sbbb_delta_eta_outside_gap, rbbb_delta_eta_outside_gap, ratiolabel1, tunf_delta_eta_outside_gap, ratiolabel2, svd_delta_eta_outside_gap, ratiolabel3, baye_delta_eta_outside_gap, ratiolabel4, plots_path, prefix+"delta_eta_outside_gap_ratios", detail);
    data_delta_eta_outside_gap->SetTitle("#Delta#eta^{out};#Delta#eta^{out};#frac{d#sigma}{d#Delta#eta^{out}} [pb]");

}
