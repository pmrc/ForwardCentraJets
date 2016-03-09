// Pedro Cipriano, Nov 2013
// DESY, CMS
// Last Update: 1 Nov 2013
//
// stat_vs_unf_unc()

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TString.h>
#include <TPaveText.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "common_methods.h"

using namespace std;

void plot_stat_unf(TH1 *hist1 = 0, TH1 *hist2 = 0, TString title = "", TString extra_label = "", string extra_label_position = "", string path= "../output/", string fileout = "test", string legend_position = "top_left", bool detail = false)
{
// plots the two distributions

// check if there are any histograms inputed
if (hist1 == 0) { cout << "Histogram1 is not provided!" << endl; return; } 
if (hist2 == 0) { cout << "Histogram2 is not provided!" << endl; return; } 

//clone the histograms
TH1D *stat = (TH1D*) hist1->Clone();
TH1D *unf = (TH1D*) hist2->Clone();

//set title
stat->SetTitle(title);

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
// declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

//setting the plooting range
    double max = 0.45;
    double min = 0.0;

//set variables
    double stat_val, unf_val, stat_ave, unf_ave;
    double stat_min = 100.0;
    double unf_min = 100.0;
    double stat_max = 0.0;
    double unf_max = 0.0;
    double stat_tot = 0.0;
    double unf_tot = 0.0;
    int stat_nbins = 0.0;
    int unf_nbins = 0.0;

for (int i = 1; i <= hist1->GetNbinsX();i++)
	{
	stat_nbins = stat_nbins + 1;
        stat_val = hist1->GetBinError(i)/hist1->GetBinContent(i);
	stat_tot = stat_tot + stat_val;
	stat->SetBinContent(i,stat_val);
	if (stat_min > stat_val) { stat_min = stat_val; }
	if (stat_max < stat_val) { stat_max = stat_val; }
	unf_nbins = unf_nbins + 1;
	unf_val = hist2->GetBinError(i)/hist2->GetBinContent(i);
	unf_tot = unf_tot + unf_val;
	unf->SetBinContent(i,unf_val);
	if (unf_min > unf_val) { unf_min = unf_val; }
	if (unf_max < unf_val) { unf_max = unf_val; }
	}

//compute stats
	stat_ave = 100.0 * stat_tot / (double) stat_nbins;
	unf_ave = 100.0 * unf_tot / (double) unf_nbins;
	stat_min = stat_min * 100;
	stat_max = stat_max * 100;
	unf_min = unf_min * 100;
	unf_max = unf_max * 100;

//display stats
if (detail)
	{
	cout<<"Statistical Ave = "<<stat_ave<<" Min = "<<stat_min<<" Max = "<<stat_max<<endl;
	cout<<"Unfolding   Ave = "<<unf_ave<<" Min = "<<unf_min<<" Max = "<<unf_max<<endl;
	}

//plooting
    stat->SetMaximum(max);
    stat->SetMinimum(min);
    format_histogram(stat, 1, 1);
    stat->Draw("hist");
    format_histogram(unf, 2, 2);
    unf->Draw("histsame");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 2, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(stat,"Statistical Uncertainty","l");
    leg01->AddEntry(unf,"Unfolding Uncertainty","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    if (extra_label_position == "top_right")    { x1 = 0.79; y1 = 0.77; x2 = 0.98; y2 = 0.92; }
    if (extra_label_position == "bottom_right") { x1 = 0.75; y1 = 0.13; x2 = 0.98; y2 = 0.18; }
    if (extra_label_position == "top_left")     { x1 = 0.13; y1 = 0.77; x2 = 0.49; y2 = 0.92; }
    if (extra_label_position == "middle_left")  { x1 = 0.13; y1 = 0.50; x2 = 0.30; y2 = 0.55; }
    if (extra_label_position == "bottom_left")  { x1 = 0.13; y1 = 0.13; x2 = 0.47; y2 = 0.18; }
    if (extra_label_position == "middle")       { x1 = 0.50; y1 = 0.48; x2 = 0.60; y2 = 0.53; }

    TPaveText *extra = new TPaveText(x1,y1,x2,y2,"NDC"); // NDC sets coords
    extra->SetTextSize(0.05);
    extra->SetBorderSize(0); 
    extra->SetTextAlign(12);
    extra->SetTextFont(42);
    extra->SetLineWidth(1);
    extra->SetLineColor(0);
    extra->SetFillColor(0);
    extra->SetFillStyle(1001);
    extra->AddText(extra_label);
    if (extra_label != "") { extra->Draw(); }

    print_plots(c01, path, fileout);
}


void stat_vs_unf_unc(string path_stat, string path_unf, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Statistical VS Unfolding Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for Statistical: " << path_stat << endl; }
    if (detail) { cout << "Input path for Unfolding:   " << path_unf << endl; }
    if (detail) { cout << "Output Path Plots:          " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:               " << detail << endl; }
    if (detail) { cout << "Display Results:            " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:                  " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_stat = new TFile( path_stat.c_str() );
    TFile *data_unf = new TFile( path_unf.c_str() );


//plot for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_stat = 0;
    TH1D *delta_phi_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi",delta_phi_stat);
    if (delta_phi_stat == 0) { cout << "ak5PF_delta_phi stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi",delta_phi_unf);
    if (delta_phi_unf == 0) { cout << "ak5PF_delta_phi unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_stat, delta_phi_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "", "", output_path_plots, "delta_phi", "top_left", false);


//plot for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi Deta1"<<endl; }

    TH1D *delta_phi_deta1_stat = 0;
    TH1D *delta_phi_deta1_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_stat);
    if (delta_phi_deta1_stat == 0) { cout << "ak5PF_delta_phi_deta1 stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_unf);
    if (delta_phi_deta1_unf == 0) { cout << "ak5PF_delta_phi_deta1 unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta1_stat, delta_phi_deta1_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "0.4 < #Delta#eta < 2.5", "top_right", output_path_plots, "delta_phi_deta1", "top_left", false);


//plot for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi Deta2"<<endl; }

    TH1D *delta_phi_deta2_stat = 0;
    TH1D *delta_phi_deta2_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_stat);
    if (delta_phi_deta2_stat == 0) { cout << "ak5PF_delta_phi_deta2 stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_unf);
    if (delta_phi_deta2_unf == 0) { cout << "ak5PF_delta_phi_deta2 unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta2_stat, delta_phi_deta2_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "2.5 < #Delta#eta < 3.5", "top_right", output_path_plots, "delta_phi_deta2", "top_left", false);


//plot for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi Deta3"<<endl; }

    TH1D *delta_phi_deta3_stat = 0;
    TH1D *delta_phi_deta3_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_stat);
    if (delta_phi_deta3_stat == 0) { cout << "ak5PF_delta_phi_deta3 stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_unf);
    if (delta_phi_deta3_unf == 0) { cout << "ak5PF_delta_phi_deta3 unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta3_stat, delta_phi_deta3_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "3.5 < #Delta#eta < 4.5", "top_right", output_path_plots, "delta_phi_deta3", "top_left", false);


//plot for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi Deta4"<<endl; }

    TH1D *delta_phi_deta4_stat = 0;
    TH1D *delta_phi_deta4_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_stat);
    if (delta_phi_deta4_stat == 0) { cout << "ak5PF_delta_phi_deta4 stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_unf);
    if (delta_phi_deta4_unf == 0) { cout << "ak5PF_delta_phi_deta4 unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta4_stat, delta_phi_deta4_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "3.5 < #Delta#eta < 7.5", "top_right", output_path_plots, "delta_phi_deta4", "top_left", true);


//plot for delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1D *delta_phi_gap_stat = 0;
    TH1D *delta_phi_gap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_stat);
    if (delta_phi_gap_stat == 0) { cout << "ak5PF_delta_phi_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_unf);
    if (delta_phi_gap_unf == 0) { cout << "ak5PF_delta_phi_gap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_gap_stat, delta_phi_gap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "", "", output_path_plots, "delta_phi_gap", "top_left", true);


//plot for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi Deta1 Gap"<<endl; }

    TH1D *delta_phi_deta1_gap_stat = 0;
    TH1D *delta_phi_deta1_gap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_stat);
    if (delta_phi_deta1_gap_stat == 0) { cout << "ak5PF_delta_phi_deta1_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_unf);
    if (delta_phi_deta1_gap_unf == 0) { cout << "ak5PF_delta_phi_deta1_gap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta1_gap_stat, delta_phi_deta1_gap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "0.4 < #Delta#eta < 2.5", "top_right", output_path_plots, "delta_phi_deta1_gap", "top_left", false);


//plot for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi Deta2 Gap"<<endl; }

    TH1D *delta_phi_deta2_gap_stat = 0;
    TH1D *delta_phi_deta2_gap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_stat);
    if (delta_phi_deta2_gap_stat == 0) { cout << "ak5PF_delta_phi_deta2_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_unf);
    if (delta_phi_deta2_gap_unf == 0) { cout << "ak5PF_delta_phi_deta2_gap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta2_gap_stat, delta_phi_deta2_gap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "2.5 < #Delta#eta < 3.5", "top_right", output_path_plots, "delta_phi_deta2_gap", "top_left", false);


//plot for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi Deta3 Gap"<<endl; }

    TH1D *delta_phi_deta3_gap_stat = 0;
    TH1D *delta_phi_deta3_gap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_stat);
    if (delta_phi_deta3_gap_stat == 0) { cout << "ak5PF_delta_phi_deta3_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_unf);
    if (delta_phi_deta3_gap_unf == 0) { cout << "ak5PF_delta_phi_deta3_gap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta3_gap_stat, delta_phi_deta3_gap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "3.5 < #Delta#eta < 4.5", "bottom_left", output_path_plots, "delta_phi_deta3_gap", "top_right", false);


//plot for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi Deta4 Gap"<<endl; }

    TH1D *delta_phi_deta4_gap_stat = 0;
    TH1D *delta_phi_deta4_gap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_stat);
    if (delta_phi_deta4_gap_stat == 0) { cout << "ak5PF_delta_phi_deta4_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_unf);
    if (delta_phi_deta4_gap_unf == 0) { cout << "ak5PF_delta_phi_deta4_gap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta4_gap_stat, delta_phi_deta4_gap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "4.5 < #Delta#eta < 7.5", "bottom_left",  output_path_plots, "delta_phi_deta4_gap", "top_right", false);


//plot for delta phi nogap distribution
    if (detail) { cout<<"Delta phi NoGap"<<endl; }

    TH1D *delta_phi_nogap_stat = 0;
    TH1D *delta_phi_nogap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_stat);
    if (delta_phi_nogap_stat == 0) { cout << "ak5PF_delta_phi_nogap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_unf);
    if (delta_phi_nogap_unf == 0) { cout << "ak5PF_delta_phi_nogap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_nogap_stat, delta_phi_nogap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "", "", output_path_plots, "delta_phi_nogap", "top_right", false);


//plot for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi Deta1 NoGap"<<endl; }

    TH1D *delta_phi_deta1_nogap_stat = 0;
    TH1D *delta_phi_deta1_nogap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_stat);
    if (delta_phi_deta1_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta1_nogap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_unf);
    if (delta_phi_deta1_nogap_unf == 0) { cout << "ak5PF_delta_phi_deta1_nogap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta1_nogap_stat, delta_phi_deta1_nogap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "0.4 < #Delta#eta < 2.5", "top_left",  output_path_plots, "delta_phi_deta1_nogap", "top_right", false);


//plot for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi Deta2 NoGap"<<endl; }

    TH1D *delta_phi_deta2_nogap_stat = 0;
    TH1D *delta_phi_deta2_nogap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_stat);
    if (delta_phi_deta2_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta2_nogap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_unf);
    if (delta_phi_deta2_nogap_unf == 0) { cout << "ak5PF_delta_phi_deta2_nogap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta2_nogap_stat, delta_phi_deta2_nogap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "2.5 < #Delta#eta < 3.5", "top_left",  output_path_plots, "delta_phi_deta2_nogap", "top_right", false);


//plot for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi Deta3 NoGap"<<endl; }

    TH1D *delta_phi_deta3_nogap_stat = 0;
    TH1D *delta_phi_deta3_nogap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_stat);
    if (delta_phi_deta3_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta3_nogap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_unf);
    if (delta_phi_deta3_nogap_unf == 0) { cout << "ak5PF_delta_phi_deta2_nogap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta3_nogap_stat, delta_phi_deta3_nogap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "3.5 < #Delta#eta < 4.5", "top_left",  output_path_plots, "delta_phi_deta3_nogap", "top_right", false);


//plot for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi Deta4 NoGap"<<endl; }

    TH1D *delta_phi_deta4_nogap_stat = 0;
    TH1D *delta_phi_deta4_nogap_unf = 0;

    data_stat->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_stat);
    if (delta_phi_deta4_nogap_stat == 0) { cout << "ak5PF_delta_phi_deta4_nogap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_unf);
    if (delta_phi_deta4_nogap_unf == 0) { cout << "ak5PF_delta_phi_deta4_nogap unf not found!" << endl; return; }

    plot_stat_unf(delta_phi_deta4_nogap_stat, delta_phi_deta4_nogap_unf, "#Delta#phi;#Delta#phi [rad];Relative Uncertainty", "4.5 < #Delta#eta < 7.5", "top_left",  output_path_plots, "delta_phi_deta4_nogap", "top_right", false);


//plot for leading pt inside gap distribution
    if (detail) { cout<<"Leading pT Inside Gap"<<endl; }

    TH1D *leading_pt_inside_gap_stat = 0;
    TH1D *leading_pt_inside_gap_unf = 0;

    data_stat->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_stat);
    if (leading_pt_inside_gap_stat == 0) { cout << "ak5PF_leading_pt_inside_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_unf);
    if (leading_pt_inside_gap_unf == 0) { cout << "ak5PF_leading_pt_inside_gap unf not found!" << endl; return; }

    plot_stat_unf(leading_pt_inside_gap_stat, leading_pt_inside_gap_unf, "p_{T}^{inside};p_{T}^{inside} [GeV];Relative Uncertainty", "", "",  output_path_plots, "leading_pt_inside_gap", "top_right", false);


//plot for leading eta star inside gap distribution
    if (detail) { cout<<"Leading Eta* Inside Gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_stat = 0;
    TH1D *leading_eta_star_inside_gap_unf = 0;

    data_stat->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_stat);
    if (leading_eta_star_inside_gap_stat == 0) { cout << "ak5PF_leading_eta_star_inside_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_unf);
    if (leading_eta_star_inside_gap_unf == 0) { cout << "ak5PF_leading_eta_star_inside_gap unf not found!" << endl; return; }

    plot_stat_unf(leading_eta_star_inside_gap_stat, leading_eta_star_inside_gap_unf, "#eta*;#eta*;Relative Uncertainty", "", "", output_path_plots, "leading_eta_star_inside_gap", "top_right", false);


//plot for leading pt outside gap distribution
    if (detail) { cout<<"Leading pT Outside Gap"<<endl; }

    TH1D *leading_pt_outside_gap_stat = 0;
    TH1D *leading_pt_outside_gap_unf = 0;

    data_stat->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_stat);
    if (leading_pt_outside_gap_stat == 0) { cout << "ak5PF_leading_pt_outside_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_unf);
    if (leading_pt_outside_gap_unf == 0) { cout << "ak5PF_leading_pt_outside_gap unf not found!" << endl; return; }

    plot_stat_unf(leading_pt_outside_gap_stat, leading_pt_outside_gap_unf, "p_{T}^{outside};p_{T}^{outside} [GeV];Relative Uncertainty", "", "", output_path_plots, "leading_pt_outside_gap", "top_right", false);


//plot for delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta Outside Gap"<<endl; }

    TH1D *delta_eta_outside_gap_stat = 0;
    TH1D *delta_eta_outside_gap_unf = 0;

    data_stat->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_stat);
    if (delta_eta_outside_gap_stat == 0) { cout << "ak5PF_delta_eta_outside_gap stat not found!" << endl; return; }
    data_unf->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_unf);
    if (delta_eta_outside_gap_unf == 0) { cout << "ak5PF_delta_eta_outside_gap unf not found!" << endl; return; }

    plot_stat_unf(delta_eta_outside_gap_stat, delta_eta_outside_gap_unf, "#Delta#eta^{outside};#Delta#eta^{outside};Relative Uncertainty", "", "", output_path_plots, "delta_eta_outside_gap", "top_right", false);
}
