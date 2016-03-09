// DESY, CMS
// Last Update: 24 Mar 2013
//
// apply correction()

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TString.h>
#include <TF1.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "common_methods.h"

using namespace std;

void plot_uncertainty(TH1D *data_new, TH1D *data_old, string path, string fileout, string legend_position, TString label)
{
//plots the model uncertainty control plots

//declaring the canvas
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetLogy();
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

//format and ploting the histogram
    data_new->SetFillColor(5);
    data_new->Draw("e2");
    data_old->SetLineStyle(1);
    data_old->SetLineWidth(3);
    data_old->SetLineColor(1);
    data_old->Draw("e1 same");
    
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 2, x1, y1, x2, y2);

//sets and draw the legend
    TLegend *leg00 = new TLegend(x1, y1, x2, y2);
    leg00->AddEntry(data_old,label,"l");
    leg00->AddEntry(data_new,"Total Uncertainty","f");
    leg00->SetFillColor(0);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();
    
 print_plots(c1, path, fileout);
}


void compute_new_bonduaries(TH1D *data, TH1D *up, TH1D *down, TH1D *out, bool detail = false)
{

int nbins = data->GetNbinsX();
double bond_up[nbins], bond_down[nbins];
double val_data = 0.0, central_val = 0.0, error = 0.0;
double val_up = 0.0, val_down = 0.0;

for (int i = 1; i <= nbins; i++)
	{
	bond_up[i-1] = 0.0;
	bond_down[i-1] = 0.0;
	val_data = data->GetBinContent(i);
	val_up = up->GetBinContent(i);
	val_down = down->GetBinContent(i);
	bond_up[i-1] = val_data * (1.0 + val_up);
	bond_down[i-1] = val_data * (1.0- val_down);
	central_val = (bond_up[i-1] + bond_down[i-1]) / 2.0;
	error = bond_up[i-1] - central_val;
	out->SetBinContent(i,central_val);
	out->SetBinError(i,error);
	if (detail) { cout << "Old Value = " << val_data << "+-" << data->GetBinError(i) << " New Value = " << central_val << "+-" << error << endl; }
	//if (detail) { cout << "Old Value = " << val_data << "-" << bond_down[i-1] << "+" << bond_up[i-1] << endl; }
}

out->SetEntries(data->GetEntries());

}


void apply_corr_uncertainty(string path_data, string path_unc_up, string path_unc_down, string path_unc_output, string output_path_plots, bool detail = false, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Apply Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for data:             " << path_data << endl; }
    if (detail) { cout << "Input path for uncertainty up:   " << path_unc_up << endl; }
    if (detail) { cout << "Input path for uncertainty down: " << path_unc_down << endl; }
    if (detail) { cout << "Output path:                     " << path_unc_output << endl; }
    if (detail) { cout << "Output Path Plots:               " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:                    " << detail << endl; }
    if (detail) { cout << "Test Mode:                       " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data = new TFile( path_data.c_str() );
    TFile *up = new TFile( path_unc_up.c_str() );
    TFile *down = new TFile( path_unc_down.c_str() );


//histogram bins
int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};


//apply uncertainty in delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_data = 0;
    TH1D *delta_phi_up = 0;
    TH1D *delta_phi_down = 0;

    data->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi",delta_phi_up);
    if (delta_phi_up == 0) { cout << "corr_unc_up_delta_phi not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi",delta_phi_down);
    if (delta_phi_down == 0) { cout << "corr_unc_down_delta_phi not found!" << endl; return; }

    TH1D *delta_phi;
    delta_phi =  new TH1D("ak5PF_delta_phi","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_data, delta_phi_up, delta_phi_down, delta_phi, detail);
    plot_uncertainty(delta_phi, delta_phi_data, output_path_plots, "data_corr_unc_delta_phi", "top_left", "Data - 2010");


//apply uncertainty in delta phi norm distribution
    if (detail) { cout<<"Delta phi Norm"<<endl; }

    TH1D *delta_phi_norm_data = 0;
    TH1D *delta_phi_norm_up = 0;
    TH1D *delta_phi_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_norm",delta_phi_norm_data);
    if (delta_phi_norm_data == 0) { cout << "ak5PF_delta_phi_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_norm",delta_phi_norm_up);
    if (delta_phi_norm_up == 0) { cout << "corr_unc_up_delta_phi_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_norm",delta_phi_norm_down);
    if (delta_phi_norm_down == 0) { cout << "corr_unc_down_delta_phi_norm not found!" << endl; return; }

    TH1D *delta_phi_norm;
    delta_phi_norm =  new TH1D("ak5PF_delta_phi_norm","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_norm_data, delta_phi_norm_up, delta_phi_norm_down, delta_phi_norm, detail);
    plot_uncertainty(delta_phi_norm, delta_phi_norm_data, output_path_plots, "data_corr_unc_delta_phi_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_data = 0;
    TH1D *delta_phi_deta1_up = 0;
    TH1D *delta_phi_deta1_down = 0;

    data->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta1",delta_phi_deta1_up);
    if (delta_phi_deta1_up == 0) { cout << "corr_unc_up_delta_phi_deta1 not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta1",delta_phi_deta1_down);
    if (delta_phi_deta1_down == 0) { cout << "corr_unc_down_delta_phi_deta1 not found!" << endl; return; }

    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D("ak5PF_delta_phi_deta1","Corrected Data for deta1;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta1_data, delta_phi_deta1_up, delta_phi_deta1_down, delta_phi_deta1, detail);
    plot_uncertainty(delta_phi_deta1, delta_phi_deta1_data, output_path_plots, "data_corr_unc_delta_phi_deta1", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta1 norm distribution
    if (detail) { cout<<"Delta phi deta1 Norm"<<endl; }

    TH1D *delta_phi_deta1_norm_data = 0;
    TH1D *delta_phi_deta1_norm_up = 0;
    TH1D *delta_phi_deta1_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta1_norm",delta_phi_deta1_norm_data);
    if (delta_phi_deta1_norm_data == 0) { cout << "ak5PF_delta_phi_deta1_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta1_norm",delta_phi_deta1_norm_up);
    if (delta_phi_deta1_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta1_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta1_norm",delta_phi_deta1_norm_down);
    if (delta_phi_deta1_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta1_norm not found!" << endl; return; }

    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm =  new TH1D("ak5PF_delta_phi_deta1_norm","Corrected Data for deta1;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta1_norm_data, delta_phi_deta1_norm_up, delta_phi_deta1_norm_down, delta_phi_deta1_norm, detail);
    plot_uncertainty(delta_phi_deta1_norm, delta_phi_deta1_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta1_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_data = 0;
    TH1D *delta_phi_deta2_up = 0;
    TH1D *delta_phi_deta2_down = 0;

    data->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta2",delta_phi_deta2_up);
    if (delta_phi_deta2_up == 0) { cout << "corr_unc_up_delta_phi_deta2 not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta2",delta_phi_deta2_down);
    if (delta_phi_deta2_down == 0) { cout << "corr_unc_down_delta_phi_deta2 not found!" << endl; return; }

    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D("ak5PF_delta_phi_deta2","Corrected Data for deta2;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta2_data, delta_phi_deta2_up, delta_phi_deta2_down, delta_phi_deta2, detail);
    plot_uncertainty(delta_phi_deta2, delta_phi_deta2_data, output_path_plots, "data_corr_unc_delta_phi_deta2", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta2 norm distribution
    if (detail) { cout<<"Delta phi deta2 Norm"<<endl; }

    TH1D *delta_phi_deta2_norm_data = 0;
    TH1D *delta_phi_deta2_norm_up = 0;
    TH1D *delta_phi_deta2_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta2_norm",delta_phi_deta2_norm_data);
    if (delta_phi_deta2_norm_data == 0) { cout << "ak5PF_delta_phi_deta2_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta2_norm",delta_phi_deta2_norm_up);
    if (delta_phi_deta2_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta2_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta2_norm",delta_phi_deta2_norm_down);
    if (delta_phi_deta2_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta2_norm not found!" << endl; return; }

    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm =  new TH1D("ak5PF_delta_phi_deta2_norm","Corrected Data for deta2;#Delta#phi [rad];#frac{d#sigma}{d#Delta#phi} [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta2_norm_data, delta_phi_deta2_norm_up, delta_phi_deta2_norm_down, delta_phi_deta2_norm, detail);
    plot_uncertainty(delta_phi_deta2_norm, delta_phi_deta2_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta2_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_data = 0;
    TH1D *delta_phi_deta3_up = 0;
    TH1D *delta_phi_deta3_down = 0;

    data->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta3",delta_phi_deta3_up);
    if (delta_phi_deta3_up == 0) { cout << "corr_unc_up_delta_phi_deta3 not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta3",delta_phi_deta3_down);
    if (delta_phi_deta3_down == 0) { cout << "corr_unc_down_delta_phi_deta3 not found!" << endl; return; }

    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D("ak5PF_delta_phi_deta3","Corrected Data for deta3;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta3_data, delta_phi_deta3_up, delta_phi_deta3_down, delta_phi_deta3, detail);
    plot_uncertainty(delta_phi_deta3, delta_phi_deta3_data, output_path_plots, "data_corr_unc_delta_phi_deta3", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta3 norm distribution
    if (detail) { cout<<"Delta phi deta3 Norm"<<endl; }

    TH1D *delta_phi_deta3_norm_data = 0;
    TH1D *delta_phi_deta3_norm_up = 0;
    TH1D *delta_phi_deta3_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta3_norm",delta_phi_deta3_norm_data);
    if (delta_phi_deta3_norm_data == 0) { cout << "ak5PF_delta_phi_deta3_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta3_norm",delta_phi_deta3_norm_up);
    if (delta_phi_deta3_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta3_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta3_norm",delta_phi_deta3_norm_down);
    if (delta_phi_deta3_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta3_norm not found!" << endl; return; }

    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm =  new TH1D("ak5PF_delta_phi_deta3_norm","Corrected Data for deta3;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta3_norm_data, delta_phi_deta3_norm_up, delta_phi_deta3_norm_down, delta_phi_deta3_norm, detail);
    plot_uncertainty(delta_phi_deta3_norm, delta_phi_deta3_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta3_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_data = 0;
    TH1D *delta_phi_deta4_up = 0;
    TH1D *delta_phi_deta4_down = 0;

    data->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta4",delta_phi_deta4_up);
    if (delta_phi_deta4_up == 0) { cout << "corr_unc_up_delta_phi_deta4 not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta4",delta_phi_deta4_down);
    if (delta_phi_deta4_down == 0) { cout << "corr_unc_down_delta_phi_deta4 not found!" << endl; return; }

    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D("ak5PF_delta_phi_deta4","Corrected Data for deta4;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta4_data, delta_phi_deta4_up, delta_phi_deta4_down, delta_phi_deta4, detail);
    plot_uncertainty(delta_phi_deta4, delta_phi_deta4_data, output_path_plots, "data_corr_unc_delta_phi_deta4", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta4 norm distribution
    if (detail) { cout<<"Delta phi deta4 Norm"<<endl; }

    TH1D *delta_phi_deta4_norm_data = 0;
    TH1D *delta_phi_deta4_norm_up = 0;
    TH1D *delta_phi_deta4_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta4_norm",delta_phi_deta4_norm_data);
    if (delta_phi_deta4_norm_data == 0) { cout << "ak5PF_delta_phi_deta4_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta4_norm",delta_phi_deta4_norm_up);
    if (delta_phi_deta4_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta4_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta4_norm",delta_phi_deta4_norm_down);
    if (delta_phi_deta4_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta4_norm not found!" << endl; return; }

    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm =  new TH1D("ak5PF_delta_phi_deta4_norm","Corrected Data for deta4;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta4_norm_data, delta_phi_deta4_norm_up, delta_phi_deta4_norm_down, delta_phi_deta4_norm, detail);
    plot_uncertainty(delta_phi_deta4_norm, delta_phi_deta4_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta4_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_data = 0;
    TH1D *delta_phi_gap_up = 0;
    TH1D *delta_phi_gap_down = 0;

    data->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_gap",delta_phi_gap_up);
    if (delta_phi_gap_up == 0) { cout << "corr_unc_up_delta_phi_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_gap",delta_phi_gap_down);
    if (delta_phi_gap_down == 0) { cout << "corr_unc_down_delta_phi_gap not found!" << endl; return; }

    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","Corrected Data gap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_gap_data, delta_phi_gap_up, delta_phi_gap_down, delta_phi_gap, detail);
    plot_uncertainty(delta_phi_gap, delta_phi_gap_data, output_path_plots, "data_corr_unc_delta_phi_gap", "top_left", "Data - 2010");


//apply uncertainty in delta phi gap norm distribution
    if (detail) { cout<<"Delta phi gap Norm"<<endl; }

    TH1D *delta_phi_gap_norm_data = 0;
    TH1D *delta_phi_gap_norm_up = 0;
    TH1D *delta_phi_gap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_gap_norm",delta_phi_gap_norm_data);
    if (delta_phi_gap_norm_data == 0) { cout << "ak5PF_delta_phi_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_gap_norm",delta_phi_gap_norm_up);
    if (delta_phi_gap_norm_up == 0) { cout << "corr_unc_up_delta_phi_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_gap_norm",delta_phi_gap_norm_down);
    if (delta_phi_gap_norm_down == 0) { cout << "corr_unc_down_delta_phi_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D("ak5PF_delta_phi_gap_norm","Corrected Data gap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_gap_norm_data, delta_phi_gap_norm_up, delta_phi_gap_norm_down, delta_phi_gap_norm, detail);
    plot_uncertainty(delta_phi_gap_norm, delta_phi_gap_norm_data, output_path_plots, "data_corr_unc_delta_phi_gap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_data = 0;
    TH1D *delta_phi_deta1_gap_up = 0;
    TH1D *delta_phi_deta1_gap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta1_gap",delta_phi_deta1_gap_up);
    if (delta_phi_deta1_gap_up == 0) { cout << "corr_unc_up_delta_phi_deta1_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta1_gap",delta_phi_deta1_gap_down);
    if (delta_phi_deta1_gap_down == 0) { cout << "corr_unc_down_delta_phi_deta1_gap not found!" << endl; return; }

    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D("ak5PF_delta_phi_deta1_gap","Corrected Data deta1 gap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta1_gap_data, delta_phi_deta1_gap_up, delta_phi_deta1_gap_down, delta_phi_deta1_gap, detail);
    plot_uncertainty(delta_phi_deta1_gap, delta_phi_deta1_gap_data, output_path_plots, "data_corr_unc_delta_phi_deta1_gap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta1 gap norm distribution
    if (detail) { cout<<"Delta phi deta1 gap Norm"<<endl; }

    TH1D *delta_phi_deta1_gap_norm_data = 0;
    TH1D *delta_phi_deta1_gap_norm_up = 0;
    TH1D *delta_phi_deta1_gap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_data);
    if (delta_phi_deta1_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta1_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_up);
    if (delta_phi_deta1_gap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta1_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_down);
    if (delta_phi_deta1_gap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta1_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm =  new TH1D("ak5PF_delta_phi_deta1_gap_norm","Corrected Data deta1 gap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta1_gap_norm_data, delta_phi_deta1_gap_norm_up, delta_phi_deta1_gap_norm_down, delta_phi_deta1_gap_norm, detail);
    plot_uncertainty(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta1_gap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_data = 0;
    TH1D *delta_phi_deta2_gap_up = 0;
    TH1D *delta_phi_deta2_gap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta2_gap",delta_phi_deta2_gap_up);
    if (delta_phi_deta2_gap_up == 0) { cout << "corr_unc_up_delta_phi_deta2_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta2_gap",delta_phi_deta2_gap_down);
    if (delta_phi_deta2_gap_down == 0) { cout << "corr_unc_down_delta_phi_deta2_gap not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D("ak5PF_delta_phi_deta2_gap","Corrected Data deta2 gap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta2_gap_data, delta_phi_deta2_gap_up, delta_phi_deta2_gap_down, delta_phi_deta2_gap, detail);
    plot_uncertainty(delta_phi_deta2_gap, delta_phi_deta2_gap_data, output_path_plots, "data_corr_unc_delta_phi_deta2_gap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta2 gap norm distribution
    if (detail) { cout<<"Delta phi deta2 gap Norm"<<endl; }

    TH1D *delta_phi_deta2_gap_norm_data = 0;
    TH1D *delta_phi_deta2_gap_norm_up = 0;
    TH1D *delta_phi_deta2_gap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_data);
    if (delta_phi_deta2_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta2_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_up);
    if (delta_phi_deta2_gap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta2_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_down);
    if (delta_phi_deta2_gap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta2_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm =  new TH1D("ak5PF_delta_phi_deta2_gap_norm","Corrected Data deta2 gap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta2_gap_norm_data, delta_phi_deta2_gap_norm_up, delta_phi_deta2_gap_norm_down, delta_phi_deta2_gap_norm, detail);
    plot_uncertainty(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta2_gap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_data = 0;
    TH1D *delta_phi_deta3_gap_up = 0;
    TH1D *delta_phi_deta3_gap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta3_gap",delta_phi_deta3_gap_up);
    if (delta_phi_deta3_gap_up == 0) { cout << "corr_unc_up_delta_phi_deta3_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta3_gap",delta_phi_deta3_gap_down);
    if (delta_phi_deta3_gap_down == 0) { cout << "corr_unc_down_delta_phi_deta3_gap not found!" << endl; return; }

    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D("ak5PF_delta_phi_deta3_gap","Corrected Data deta3 gap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta3_gap_data, delta_phi_deta3_gap_up, delta_phi_deta3_gap_down, delta_phi_deta3_gap, detail);
    plot_uncertainty(delta_phi_deta3_gap, delta_phi_deta3_gap_data, output_path_plots, "data_corr_unc_delta_phi_deta3_gap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta3 gap norm distribution
    if (detail) { cout<<"Delta phi deta3 gap Norm"<<endl; }

    TH1D *delta_phi_deta3_gap_norm_data = 0;
    TH1D *delta_phi_deta3_gap_norm_up = 0;
    TH1D *delta_phi_deta3_gap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_data);
    if (delta_phi_deta3_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta3_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_up);
    if (delta_phi_deta3_gap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta3_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_down);
    if (delta_phi_deta3_gap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta3_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm =  new TH1D("ak5PF_delta_phi_deta3_gap_norm","Corrected Data deta3 gap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta3_gap_norm_data, delta_phi_deta3_gap_norm_up, delta_phi_deta3_gap_norm_down, delta_phi_deta3_gap_norm, detail);
    plot_uncertainty(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta3_gap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_data = 0;
    TH1D *delta_phi_deta4_gap_up = 0;
    TH1D *delta_phi_deta4_gap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta4_gap",delta_phi_deta4_gap_up);
    if (delta_phi_deta4_gap_up == 0) { cout << "corr_unc_up_delta_phi_deta4_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta4_gap",delta_phi_deta4_gap_down);
    if (delta_phi_deta4_gap_down == 0) { cout << "corr_unc_down_delta_phi_deta4_gap not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D("ak5PF_delta_phi_deta4_gap","Corrected Data deta4 gap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta4_gap_data, delta_phi_deta4_gap_up, delta_phi_deta4_gap_down, delta_phi_deta4_gap, detail);
    plot_uncertainty(delta_phi_deta4_gap, delta_phi_deta4_gap_data, output_path_plots, "data_corr_unc_delta_phi_deta4_gap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta4 gap Norm distribution
    if (detail) { cout<<"Delta phi deta4 gap Norm"<<endl; }

    TH1D *delta_phi_deta4_gap_norm_data = 0;
    TH1D *delta_phi_deta4_gap_norm_up = 0;
    TH1D *delta_phi_deta4_gap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_data);
    if (delta_phi_deta4_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta4_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_up);
    if (delta_phi_deta4_gap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta4_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_down);
    if (delta_phi_deta4_gap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta4_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm =  new TH1D("ak5PF_delta_phi_deta4_gap_norm","Corrected Data deta4 gap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta4_gap_norm_data, delta_phi_deta4_gap_norm_up, delta_phi_deta4_gap_norm_down, delta_phi_deta4_gap_norm, detail);
    plot_uncertainty(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta4_gap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_data = 0;
    TH1D *delta_phi_nogap_up = 0;
    TH1D *delta_phi_nogap_down = 0;

    data->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_nogap",delta_phi_nogap_up);
    if (delta_phi_nogap_up == 0) { cout << "corr_unc_up_delta_phi_nogap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_nogap",delta_phi_nogap_down);
    if (delta_phi_nogap_down == 0) { cout << "corr_unc_down_delta_phi_nogap not found!" << endl; return; }

    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","Corrected Data nogap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_nogap_data, delta_phi_nogap_up, delta_phi_nogap_down, delta_phi_nogap, detail);
    plot_uncertainty(delta_phi_nogap, delta_phi_nogap_data, output_path_plots, "data_corr_unc_delta_phi_nogap", "top_left", "Data - 2010");


//apply uncertainty in delta phi nogap norm distribution
    if (detail) { cout<<"Delta phi nogap Norm"<<endl; }

    TH1D *delta_phi_nogap_norm_data = 0;
    TH1D *delta_phi_nogap_norm_up = 0;
    TH1D *delta_phi_nogap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_nogap_norm",delta_phi_nogap_norm_data);
    if (delta_phi_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_nogap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_nogap_norm",delta_phi_nogap_norm_up);
    if (delta_phi_nogap_norm_up == 0) { cout << "corr_unc_up_delta_phi_nogap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_nogap_norm",delta_phi_nogap_norm_down);
    if (delta_phi_nogap_norm_down == 0) { cout << "corr_unc_down_delta_phi_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D("ak5PF_delta_phi_nogap_norm","Corrected Data nogap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_nogap_norm_data, delta_phi_nogap_norm_up, delta_phi_nogap_norm_down, delta_phi_nogap_norm, detail);
    plot_uncertainty(delta_phi_nogap_norm, delta_phi_nogap_norm_data, output_path_plots, "data_corr_unc_delta_phi_nogap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_data = 0;
    TH1D *delta_phi_deta1_nogap_up = 0;
    TH1D *delta_phi_deta1_nogap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta1_nogap",delta_phi_deta1_nogap_up);
    if (delta_phi_deta1_nogap_up == 0) { cout << "corr_unc_up_delta_phi_deta1_nogap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta1_nogap",delta_phi_deta1_nogap_down);
    if (delta_phi_deta1_nogap_down == 0) { cout << "corr_unc_down_delta_phi_deta1_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D("ak5PF_delta_phi_deta1_nogap","Corrected Data deta1 nogap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_up, delta_phi_deta1_nogap_down, delta_phi_deta1_nogap, detail);
    plot_uncertainty(delta_phi_deta1_nogap, delta_phi_deta1_nogap_data, output_path_plots, "data_corr_unc_delta_phi_deta1_nogap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta1 nogap norm distribution
    if (detail) { cout<<"Delta phi deta1 nogap Norm"<<endl; }

    TH1D *delta_phi_deta1_nogap_norm_data = 0;
    TH1D *delta_phi_deta1_nogap_norm_up = 0;
    TH1D *delta_phi_deta1_nogap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_data);
    if (delta_phi_deta1_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_up);
    if (delta_phi_deta1_nogap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_down);
    if (delta_phi_deta1_nogap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta1_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm =  new TH1D("ak5PF_delta_phi_deta1_nogap_norm","Corrected Data deta1 nogap norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta1_nogap_norm_data, delta_phi_deta1_nogap_norm_up, delta_phi_deta1_nogap_norm_down, delta_phi_deta1_nogap_norm, detail);
    plot_uncertainty(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta1_nogap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_data = 0;
    TH1D *delta_phi_deta2_nogap_up = 0;
    TH1D *delta_phi_deta2_nogap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta2_nogap",delta_phi_deta2_nogap_up);
    if (delta_phi_deta2_nogap_up == 0) { cout << "corr_unc_up_delta_phi_deta2_nogap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta2_nogap",delta_phi_deta2_nogap_down);
    if (delta_phi_deta2_nogap_down == 0) { cout << "corr_unc_down_delta_phi_deta2_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D("ak5PF_delta_phi_deta2_nogap","Corrected Data deta2 nogap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_up, delta_phi_deta2_nogap_down, delta_phi_deta2_nogap, detail);
    plot_uncertainty(delta_phi_deta2_nogap, delta_phi_deta2_nogap_data, output_path_plots, "data_corr_unc_delta_phi_deta2_nogap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta2 nogap norm distribution
    if (detail) { cout<<"Delta phi deta2 nogap Norm"<<endl; }

    TH1D *delta_phi_deta2_nogap_norm_data = 0;
    TH1D *delta_phi_deta2_nogap_norm_up = 0;
    TH1D *delta_phi_deta2_nogap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_data);
    if (delta_phi_deta2_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_up);
    if (delta_phi_deta2_nogap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_down);
    if (delta_phi_deta2_nogap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta2_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm =  new TH1D("ak5PF_delta_phi_deta2_nogap_norm","Corrected Data deta2 nogap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta2_nogap_norm_data, delta_phi_deta2_nogap_norm_up, delta_phi_deta2_nogap_norm_down, delta_phi_deta2_nogap_norm, detail);
    plot_uncertainty(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta2_nogap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_data = 0;
    TH1D *delta_phi_deta3_nogap_up = 0;
    TH1D *delta_phi_deta3_nogap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta3_nogap",delta_phi_deta3_nogap_up);
    if (delta_phi_deta3_nogap_up == 0) { cout << "corr_unc_up_delta_phi_deta3_nogap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta3_nogap",delta_phi_deta3_nogap_down);
    if (delta_phi_deta3_nogap_down == 0) { cout << "corr_unc_down_delta_phi_deta3_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D("ak5PF_delta_phi_deta3_nogap","Corrected Data deta3 nogap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_up, delta_phi_deta3_nogap_down, delta_phi_deta3_nogap, detail);
    plot_uncertainty(delta_phi_deta3_nogap, delta_phi_deta3_nogap_data, output_path_plots, "data_corr_unc_delta_phi_deta3_nogap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta3 nogap norm distribution
    if (detail) { cout<<"Delta phi deta3 nogap Norm"<<endl; }

    TH1D *delta_phi_deta3_nogap_norm_data = 0;
    TH1D *delta_phi_deta3_nogap_norm_up = 0;
    TH1D *delta_phi_deta3_nogap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_data);
    if (delta_phi_deta3_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_up);
    if (delta_phi_deta3_nogap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_down);
    if (delta_phi_deta3_nogap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta3_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm =  new TH1D("ak5PF_delta_phi_deta3_nogap_norm","Corrected Data deta3 nogap Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta3_nogap_norm_data, delta_phi_deta3_nogap_norm_up, delta_phi_deta3_nogap_norm_down, delta_phi_deta3_nogap_norm, detail);
    plot_uncertainty(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta3_nogap_norm", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_data = 0;
    TH1D *delta_phi_deta4_nogap_up = 0;
    TH1D *delta_phi_deta4_nogap_down = 0;

    data->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta4_nogap",delta_phi_deta4_nogap_up);
    if (delta_phi_deta4_nogap_up == 0) { cout << "corr_unc_up_delta_phi_deta4_nogap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta4_nogap",delta_phi_deta4_nogap_down);
    if (delta_phi_deta4_nogap_down == 0) { cout << "corr_unc_down_delta_phi_deta4_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D("ak5PF_delta_phi_deta4_nogap","Corrected Data deta4 nogap;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_up, delta_phi_deta4_nogap_down, delta_phi_deta4_nogap, detail);
    plot_uncertainty(delta_phi_deta4_nogap, delta_phi_deta4_nogap_data, output_path_plots, "data_corr_unc_delta_phi_deta4_nogap", "top_left", "Data - 2010");


//apply uncertainty in delta phi deta4 nogap norm distribution
    if (detail) { cout<<"Delta phi deta4 nogap Norm"<<endl; }

    TH1D *delta_phi_deta4_nogap_norm_data = 0;
    TH1D *delta_phi_deta4_nogap_norm_up = 0;
    TH1D *delta_phi_deta4_nogap_norm_down = 0;

    data->GetObject("ak5PF_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_data);
    if (delta_phi_deta4_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_up);
    if (delta_phi_deta4_nogap_norm_up == 0) { cout << "corr_unc_up_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_down);
    if (delta_phi_deta4_nogap_norm_down == 0) { cout << "corr_unc_down_delta_phi_deta4_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm =  new TH1D("ak5PF_delta_phi_deta4_nogap_norm","Corrected Data deta4 nogap norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    compute_new_bonduaries(delta_phi_deta4_nogap_norm_data, delta_phi_deta4_nogap_norm_up, delta_phi_deta4_nogap_norm_down, delta_phi_deta4_nogap_norm, detail);
    plot_uncertainty(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_data, output_path_plots, "data_corr_unc_delta_phi_deta4_nogap_norm", "top_left", "Data - 2010");


//apply uncertainty in leading pt inside gap distribution
    if (detail) { cout<<"leading pt inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_data = 0;
    TH1D *leading_pt_inside_gap_up = 0;
    TH1D *leading_pt_inside_gap_down = 0;

    data->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_leading_pt_inside_gap",leading_pt_inside_gap_up);
    if (leading_pt_inside_gap_up == 0) { cout << "corr_unc_up_leading_pt_inside_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_leading_pt_inside_gap",leading_pt_inside_gap_down);
    if (leading_pt_inside_gap_down == 0) { cout << "corr_unc_down_leading_pt_inside_gap not found!" << endl; return; }

    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D("ak5PF_leading_pt_inside_gap","Corrected Data leading_pt_inside_gap;p_{T}^{inside} [GeV];d#sigma/dp_{T}^{inside} [#frac{pb}{GeV}]", in_nbins, in_bins);

    compute_new_bonduaries(leading_pt_inside_gap_data, leading_pt_inside_gap_up, leading_pt_inside_gap_down, leading_pt_inside_gap, detail);
    plot_uncertainty(leading_pt_inside_gap, leading_pt_inside_gap_data, output_path_plots, "data_corr_unc_leading_pt_inside_gap", "top_right", "Data - 2010");


//apply uncertainty in leading pt inside gap norm distribution
    if (detail) { cout<<"leading pt inside gap norm"<<endl; }

    TH1D *leading_pt_inside_gap_norm_data = 0;
    TH1D *leading_pt_inside_gap_norm_up = 0;
    TH1D *leading_pt_inside_gap_norm_down = 0;

    data->GetObject("ak5PF_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_data);
    if (leading_pt_inside_gap_norm_data == 0) { cout << "ak5PF_leading_pt_inside_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_up);
    if (leading_pt_inside_gap_norm_up == 0) { cout << "corr_unc_up_leading_pt_inside_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_down);
    if (leading_pt_inside_gap_norm_down == 0) { cout << "corr_unc_down_leading_pt_inside_gap_norm not found!" << endl; return; }

    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm =  new TH1D("ak5PF_leading_pt_inside_gap_norm","Corrected Data leading_pt_inside_gap Norm;p_{T}^{inside} [GeV];d#sigma/dp_{T}^{inside} [#frac{pb}{GeV}]", in_nbins, in_bins);

    compute_new_bonduaries(leading_pt_inside_gap_norm_data, leading_pt_inside_gap_norm_up, leading_pt_inside_gap_norm_down, leading_pt_inside_gap_norm, detail);
    plot_uncertainty(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_data, output_path_plots, "data_corr_unc_leading_pt_inside_gap_norm", "top_right", "Data - 2010");


//apply uncertainty in leading eta star inside gap distribution
    if (detail) { cout<<"leading eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_data = 0;
    TH1D *leading_eta_star_inside_gap_up = 0;
    TH1D *leading_eta_star_inside_gap_down = 0;

    data->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_leading_eta_star_inside_gap",leading_eta_star_inside_gap_up);
    if (leading_eta_star_inside_gap_up == 0) { cout << "corr_unc_up_leading_eta_star_inside_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_leading_eta_star_inside_gap",leading_eta_star_inside_gap_down);
    if (leading_eta_star_inside_gap_down == 0) { cout << "corr_unc_down_leading_eta_star_inside_gap not found!" << endl; return; }

    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D("ak5PF_leading_eta_star_inside_gap","Corrected Data leading_eta_star_inside_gap;#eta*;d#sigma/d#eta* [pb]", etastar_nbins, etastar_bins);

    compute_new_bonduaries(leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_up, leading_eta_star_inside_gap_down, leading_eta_star_inside_gap, detail);
    plot_uncertainty(leading_eta_star_inside_gap, leading_eta_star_inside_gap_data, output_path_plots, "data_corr_unc_leading_eta_star_inside_gap", "bottom_middle", "Data - 2010");


//apply uncertainty in leading eta star inside gap norm distribution
    if (detail) { cout<<"leading eta star inside gap Norm"<<endl; }

    TH1D *leading_eta_star_inside_gap_norm_data = 0;
    TH1D *leading_eta_star_inside_gap_norm_up = 0;
    TH1D *leading_eta_star_inside_gap_norm_down = 0;

    data->GetObject("ak5PF_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_data);
    if (leading_eta_star_inside_gap_norm_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_up);
    if (leading_eta_star_inside_gap_norm_up == 0) { cout << "corr_unc_up_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_down);
    if (leading_eta_star_inside_gap_norm_down == 0) { cout << "corr_unc_down_leading_eta_star_inside_gap_norm not found!" << endl; return; }

    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm =  new TH1D("ak5PF_leading_eta_star_inside_gap_norm","Corrected Data leading_eta_star_inside_gap norm;#eta*;d#sigma/d#eta* [pb]", etastar_nbins, etastar_bins);

    compute_new_bonduaries(leading_eta_star_inside_gap_norm_data, leading_eta_star_inside_gap_norm_up, leading_eta_star_inside_gap_norm_down, leading_eta_star_inside_gap_norm, detail);
    plot_uncertainty(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_data, output_path_plots, "data_corr_unc_leading_eta_star_inside_gap_norm", "bottom_middle", "Data - 2010");


//apply uncertainty in delta eta outside gap distribution
    if (detail) { cout<<"delta eta outside gap"<<endl; }

    TH1D *delta_eta_outside_gap_data = 0;
    TH1D *delta_eta_outside_gap_up = 0;
    TH1D *delta_eta_outside_gap_down = 0;

    data->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_eta_outside_gap",delta_eta_outside_gap_up);
    if (delta_eta_outside_gap_up == 0) { cout << "corr_unc_up_delta_eta_outside_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_eta_outside_gap",delta_eta_outside_gap_down);
    if (delta_eta_outside_gap_down == 0) { cout << "corr_unc_down_delta_eta_outside_gap not found!" << endl; return; }

    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D("ak5PF_delta_eta_outside_gap","Corrected Data delta_eta_outside_gap;#Delta#eta^{out};d#sigma/d#Delta#eta^{out} [pb]", deta_out_nbins, deta_out_bins);

    compute_new_bonduaries(delta_eta_outside_gap_data, delta_eta_outside_gap_up, delta_eta_outside_gap_down, delta_eta_outside_gap, detail);
    plot_uncertainty(delta_eta_outside_gap, delta_eta_outside_gap_data, output_path_plots, "data_corr_unc_delta_eta_outside_gap", "top_right", "Data - 2010");


//apply uncertainty in delta eta outside gap norm distribution
    if (detail) { cout<<"delta eta outside gap Norm"<<endl; }

    TH1D *delta_eta_outside_gap_norm_data = 0;
    TH1D *delta_eta_outside_gap_norm_up = 0;
    TH1D *delta_eta_outside_gap_norm_down = 0;

    data->GetObject("ak5PF_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_data);
    if (delta_eta_outside_gap_norm_data == 0) { cout << "ak5PF_delta_eta_outside_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_up);
    if (delta_eta_outside_gap_norm_up == 0) { cout << "corr_unc_up_delta_eta_outside_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_down);
    if (delta_eta_outside_gap_norm_down == 0) { cout << "corr_unc_down_delta_eta_outside_gap_norm not found!" << endl; return; }

    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm =  new TH1D("ak5PF_delta_eta_outside_gap_norm","Corrected Data delta_eta_outside_gap_norm;#Delta#eta^{out};d#sigma/d#Delta#eta^{out} [pb]", deta_out_nbins, deta_out_bins);

    compute_new_bonduaries(delta_eta_outside_gap_norm_data, delta_eta_outside_gap_norm_up, delta_eta_outside_gap_norm_down, delta_eta_outside_gap_norm, detail);
    plot_uncertainty(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_data, output_path_plots, "data_corr_unc_delta_eta_outside_gap_norm", "top_right", "Data - 2010");


//apply uncertainty in leading pt outside gap distribution
    if (detail) { cout<<"leading pt outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_data = 0;
    TH1D *leading_pt_outside_gap_up = 0;
    TH1D *leading_pt_outside_gap_down = 0;

    data->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    up->GetObject("corr_unc_up_leading_pt_outside_gap",leading_pt_outside_gap_up);
    if (leading_pt_outside_gap_up == 0) { cout << "corr_unc_up_leading_pt_outside_gap not found!" << endl; return; }
    down->GetObject("corr_unc_down_leading_pt_outside_gap",leading_pt_outside_gap_down);
    if (leading_pt_outside_gap_down == 0) { cout << "corr_unc_down_leading_pt_outside_gap not found!" << endl; return; }

    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D("ak5PF_leading_pt_outside_gap","Corrected Data leading_pt_outside_gap;p_{T}^{outside} [GeV];d#sigma/dp_{T}^{outside} [#frac{pb}{GeV}]", out_nbins, out_bins);

    compute_new_bonduaries(leading_pt_outside_gap_data, leading_pt_outside_gap_up, leading_pt_outside_gap_down, leading_pt_outside_gap, detail);
    plot_uncertainty(leading_pt_outside_gap, leading_pt_outside_gap_data, output_path_plots, "data_corr_unc_leading_pt_outside_gap", "top_right", "Data - 2010");


//apply uncertainty in leading pt outside gap norm distribution
    if (detail) { cout<<"leading pt outside gap Norm"<<endl; }

    TH1D *leading_pt_outside_gap_norm_data = 0;
    TH1D *leading_pt_outside_gap_norm_up = 0;
    TH1D *leading_pt_outside_gap_norm_down = 0;

    data->GetObject("ak5PF_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_data);
    if (leading_pt_outside_gap_norm_data == 0) { cout << "ak5PF_leading_pt_outside_gap_norm not found!" << endl; return; }
    up->GetObject("corr_unc_up_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_up);
    if (leading_pt_outside_gap_norm_up == 0) { cout << "corr_unc_up_leading_pt_outside_gap_norm not found!" << endl; return; }
    down->GetObject("corr_unc_down_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_down);
    if (leading_pt_outside_gap_norm_down == 0) { cout << "corr_unc_down_leading_pt_outside_gap_norm not found!" << endl; return; }

    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm =  new TH1D("ak5PF_leading_pt_outside_gap_norm","Corrected Data leading_pt_outside_gap_norm;p_{T}^{outside} [GeV];d#sigma/dp_{T}^{outside} [#frac{pb}{GeV}]", out_nbins, out_bins);

    compute_new_bonduaries(leading_pt_outside_gap_norm_data, leading_pt_outside_gap_norm_up, leading_pt_outside_gap_norm_down, leading_pt_outside_gap_norm, detail);
    plot_uncertainty(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_data, output_path_plots, "data_corr_unc_leading_pt_outside_gap_norm", "top_right", "Data - 2010");

//Opening the output root file
    if (detail) { cout<<"Creating " << path_unc_output << "..."<<endl; }
    TFile *output = TFile::Open( path_unc_output.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file..."<<endl; }
    delta_phi->Write();
    delta_phi_deta1->Write();
    delta_phi_deta2->Write();
    delta_phi_deta3->Write();
    delta_phi_deta4->Write();
    delta_phi_gap->Write();
    delta_phi_deta1_gap->Write();
    delta_phi_deta2_gap->Write();
    delta_phi_deta3_gap->Write();
    delta_phi_deta4_gap->Write();
    delta_phi_nogap->Write();
    delta_phi_deta1_nogap->Write();
    delta_phi_deta2_nogap->Write();
    delta_phi_deta3_nogap->Write();
    delta_phi_deta4_nogap->Write();
    leading_pt_inside_gap->Write();
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();
    leading_pt_outside_gap->Write();

    delta_phi_norm->Write();
    delta_phi_deta1_norm->Write();
    delta_phi_deta2_norm->Write();
    delta_phi_deta3_norm->Write();
    delta_phi_deta4_norm->Write();
    delta_phi_gap_norm->Write();
    delta_phi_deta1_gap_norm->Write();
    delta_phi_deta2_gap_norm->Write();
    delta_phi_deta3_gap_norm->Write();
    delta_phi_deta4_gap_norm->Write();
    delta_phi_nogap_norm->Write();
    delta_phi_deta1_nogap_norm->Write();
    delta_phi_deta2_nogap_norm->Write();
    delta_phi_deta3_nogap_norm->Write();
    delta_phi_deta4_nogap_norm->Write();
    leading_pt_inside_gap_norm->Write();
    leading_eta_star_inside_gap_norm->Write();
    delta_eta_outside_gap_norm->Write();
    leading_pt_outside_gap_norm->Write();
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data->Close();
    up->Close();
    down->Close();
    output->Close();

//deleting all histograms to avoid memory leak - causes a memory leak!
/*    delete(delta_phi);
    delete(delta_phi_deta1);
    delete(delta_phi_deta2);
    delete(delta_phi_deta3);
    delete(delta_phi_deta4);
    delete(delta_phi_gap);
    delete(delta_phi_deta1_gap);
    delete(delta_phi_deta2_gap);
    delete(delta_phi_deta3_gap);
    delete(delta_phi_deta4_gap);
    delete(delta_phi_nogap);
    delete(delta_phi_deta1_nogap);
    delete(delta_phi_deta2_nogap);
    delete(delta_phi_deta3_nogap);
    delete(delta_phi_deta4_nogap);
    delete(leading_pt_inside_gap);
    delete(leading_pt_outside_gap);
    delete(delta_eta_outside_gap);
    delete(leading_eta_star_inside_gap); */

    if (detail) { cout<<"Done!"<<endl; }

}


void set_uncorr_unc(TH1D *data, TH1D *unc, TH1D *out, bool detail = false)
{

int nbins = data->GetNbinsX();
double val_data = 0.0, error = 0.0, val_unc = 0.0;

for (int i = 1; i <= nbins; i++)
	{
	val_data = data->GetBinContent(i);
	val_unc = unc->GetBinContent(i);
	error = val_unc*val_data;
	out->SetBinContent(i,val_data);
	out->SetBinError(i,error);
	if (detail) { cout << "Old Value = " << data->GetBinError(i) << " New Value = " << error << endl; }
}

out->SetEntries(data->GetEntries());

}

void apply_uncorr_uncertainty(string path_data, string path_unc, string path_out, string output_path_plots, bool detail = false, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Apply Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for data:        " << path_data << endl; }
    if (detail) { cout << "Input path for uncertainty: " << path_unc << endl; }
    if (detail) { cout << "Input path for output:      " << path_out << endl; }
    if (detail) { cout << "Output Path Plots:          " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:               " << detail << endl; }
    if (detail) { cout << "Test Mode:                  " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data = new TFile( path_data.c_str() );
    TFile *unc = new TFile( path_unc.c_str() );



//histogram bins
int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//apply uncertainty in delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_data = 0;
    TH1D *delta_phi_unc = 0;

    data->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi",delta_phi_unc);
    if (delta_phi_unc == 0) { cout << "uncorr_unc_delta_phi not found!" << endl; return; }

    TH1D *delta_phi;
    delta_phi =  new TH1D("ak5PF_delta_phi","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_data, delta_phi_unc, delta_phi, detail);
    plot_2histograms(delta_phi_data, "Only statistical uncertainty", delta_phi, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi", "top_left", true, true);

//apply uncertainty in delta phi norm distribution
    if (detail) { cout<<"Delta phi Norm"<<endl; }

    TH1D *delta_phi_norm_data = 0;
    TH1D *delta_phi_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_norm",delta_phi_norm_data);
    if (delta_phi_norm_data == 0) { cout << "ak5PF_delta_phi_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_norm",delta_phi_norm_unc);
    if (delta_phi_norm_unc == 0) { cout << "uncorr_unc_delta_phi_norm not found!" << endl; return; }

    TH1D *delta_phi_norm;
    delta_phi_norm =  new TH1D("ak5PF_delta_phi_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_norm_data, delta_phi_norm_unc, delta_phi_norm, detail);
    plot_2histograms(delta_phi_norm_data, "Only statistical uncertainty", delta_phi_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_norm", "top_left", true, true);


//apply uncertainty in delta phi deta1 distribution
    if (detail) { cout<<"Delta phi Deta1"<<endl; }

    TH1D *delta_phi_deta1_data = 0;
    TH1D *delta_phi_deta1_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta1",delta_phi_deta1_unc);
    if (delta_phi_deta1_unc == 0) { cout << "uncorr_unc_delta_phi_deta1 not found!" << endl; return; }

    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D("ak5PF_delta_phi_deta1","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta1_data, delta_phi_deta1_unc, delta_phi_deta1, detail);
    plot_2histograms(delta_phi_deta1_data, "Only statistical uncertainty", delta_phi_deta1, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta1", "top_left", true, true);


//apply uncertainty in delta phi deta1 norm distribution
    if (detail) { cout<<"Delta phi Deta1 Norm"<<endl; }

    TH1D *delta_phi_deta1_norm_data = 0;
    TH1D *delta_phi_deta1_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta1_norm",delta_phi_deta1_norm_data);
    if (delta_phi_deta1_norm_data == 0) { cout << "ak5PF_delta_phi_deta1_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta1_norm",delta_phi_deta1_norm_unc);
    if (delta_phi_deta1_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta1_norm not found!" << endl; return; }

    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm =  new TH1D("ak5PF_delta_phi_deta1_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta1_norm_data, delta_phi_deta1_norm_unc, delta_phi_deta1_norm, detail);
    plot_2histograms(delta_phi_deta1_norm_data, "Only statistical uncertainty", delta_phi_deta1_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta1_norm", "top_left", true, true);


//apply uncertainty in delta phi deta2 distribution
    if (detail) { cout<<"Delta phi Deta2"<<endl; }

    TH1D *delta_phi_deta2_data = 0;
    TH1D *delta_phi_deta2_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta2",delta_phi_deta2_unc);
    if (delta_phi_deta2_unc == 0) { cout << "uncorr_unc_delta_phi_deta2 not found!" << endl; return; }

    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D("ak5PF_delta_phi_deta2","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta2_data, delta_phi_deta2_unc, delta_phi_deta2, detail);
    plot_2histograms(delta_phi_deta2_data, "Only statistical uncertainty", delta_phi_deta2, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta2", "top_left", true, true);


//apply uncertainty in delta phi deta2 norm distribution
    if (detail) { cout<<"Delta phi Deta2 Norm"<<endl; }

    TH1D *delta_phi_deta2_norm_data = 0;
    TH1D *delta_phi_deta2_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta2_norm",delta_phi_deta2_norm_data);
    if (delta_phi_deta2_norm_data == 0) { cout << "ak5PF_delta_phi_deta2_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta2_norm",delta_phi_deta2_norm_unc);
    if (delta_phi_deta2_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta2_norm not found!" << endl; return; }

    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm =  new TH1D("ak5PF_delta_phi_deta2_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta2_norm_data, delta_phi_deta2_norm_unc, delta_phi_deta2_norm, detail);
    plot_2histograms(delta_phi_deta2_norm_data, "Only statistical uncertainty", delta_phi_deta2_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta2_norm", "top_left", true, true);


//apply uncertainty in delta phi deta3 distribution
    if (detail) { cout<<"Delta phi Deta3"<<endl; }

    TH1D *delta_phi_deta3_data = 0;
    TH1D *delta_phi_deta3_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta3",delta_phi_deta3_unc);
    if (delta_phi_deta3_unc == 0) { cout << "uncorr_unc_delta_phi_deta3 not found!" << endl; return; }

    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D("ak5PF_delta_phi_deta3","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta3_data, delta_phi_deta3_unc, delta_phi_deta3, detail);
    plot_2histograms(delta_phi_deta3_data, "Only statistical uncertainty", delta_phi_deta3, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta3", "top_left", true, true);


//apply uncertainty in delta phi deta3 norm distribution
    if (detail) { cout<<"Delta phi Deta3 Norm"<<endl; }

    TH1D *delta_phi_deta3_norm_data = 0;
    TH1D *delta_phi_deta3_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta3_norm",delta_phi_deta3_norm_data);
    if (delta_phi_deta3_norm_data == 0) { cout << "ak5PF_delta_phi_deta3_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta3_norm",delta_phi_deta3_norm_unc);
    if (delta_phi_deta3_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta3_norm not found!" << endl; return; }

    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm =  new TH1D("ak5PF_delta_phi_deta3_norm","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta3_norm_data, delta_phi_deta3_norm_unc, delta_phi_deta3_norm, detail);
    plot_2histograms(delta_phi_deta3_norm_data, "Only statistical uncertainty", delta_phi_deta3_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta3_norm", "top_left", true, true);


//apply uncertainty in delta phi deta4 distribution
    if (detail) { cout<<"Delta phi Deta4"<<endl; }

    TH1D *delta_phi_deta4_data = 0;
    TH1D *delta_phi_deta4_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta4",delta_phi_deta4_unc);
    if (delta_phi_deta4_unc == 0) { cout << "uncorr_unc_delta_phi_deta4 not found!" << endl; return; }

    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D("ak5PF_delta_phi_deta4","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta4_data, delta_phi_deta4_unc, delta_phi_deta4, detail);
    plot_2histograms(delta_phi_deta4_data, "Only statistical uncertainty", delta_phi_deta4, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta4", "top_left", true, true);


//apply uncertainty in delta phi deta4 norm distribution
    if (detail) { cout<<"Delta phi Deta4 Norm"<<endl; }

    TH1D *delta_phi_deta4_norm_data = 0;
    TH1D *delta_phi_deta4_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta4_norm",delta_phi_deta4_norm_data);
    if (delta_phi_deta4_norm_data == 0) { cout << "ak5PF_delta_phi_deta4_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta4_norm",delta_phi_deta4_norm_unc);
    if (delta_phi_deta4_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta4_norm not found!" << endl; return; }

    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm =  new TH1D("ak5PF_delta_phi_deta4_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta4_norm_data, delta_phi_deta4_norm_unc, delta_phi_deta4_norm, detail);
    plot_2histograms(delta_phi_deta4_norm_data, "Only statistical uncertainty", delta_phi_deta4_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta4_norm", "top_left", true, true);


//apply uncertainty in delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1D *delta_phi_gap_data = 0;
    TH1D *delta_phi_gap_unc = 0;

    data->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_gap",delta_phi_gap_unc);
    if (delta_phi_gap_unc == 0) { cout << "uncorr_unc_delta_phi_gap not found!" << endl; return; }

    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_gap_data, delta_phi_gap_unc, delta_phi_gap, detail);
    plot_2histograms(delta_phi_gap_data, "Only statistical uncertainty", delta_phi_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_gap", "top_left", true, true);


//apply uncertainty in delta phi gap norm distribution
    if (detail) { cout<<"Delta phi Gap Norm"<<endl; }

    TH1D *delta_phi_gap_norm_data = 0;
    TH1D *delta_phi_gap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_gap_norm",delta_phi_gap_norm_data);
    if (delta_phi_gap_norm_data == 0) { cout << "ak5PF_delta_phi_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_gap_norm",delta_phi_gap_norm_unc);
    if (delta_phi_gap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D("ak5PF_delta_phi_gap_norm","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_gap_norm_data, delta_phi_gap_norm_unc, delta_phi_gap_norm, detail);
    plot_2histograms(delta_phi_gap_norm_data, "Only statistical uncertainty", delta_phi_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_gap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi Deta1 Gap"<<endl; }

    TH1D *delta_phi_deta1_gap_data = 0;
    TH1D *delta_phi_deta1_gap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta1_gap",delta_phi_deta1_gap_unc);
    if (delta_phi_deta1_gap_unc == 0) { cout << "uncorr_unc_delta_phi_deta1_gap not found!" << endl; return; }

    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D("ak5PF_delta_phi_deta1_gap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta1_gap_data, delta_phi_deta1_gap_unc, delta_phi_deta1_gap, detail);
    plot_2histograms(delta_phi_deta1_gap_data, "Only statistical uncertainty", delta_phi_deta1_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta1_gap", "top_left", true, true);


//apply uncertainty in delta phi deta1 gap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Gap Norm"<<endl; }

    TH1D *delta_phi_deta1_gap_norm_data = 0;
    TH1D *delta_phi_deta1_gap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_data);
    if (delta_phi_deta1_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta1_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta1_gap_norm",delta_phi_deta1_gap_norm_unc);
    if (delta_phi_deta1_gap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta1_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm =  new TH1D("ak5PF_delta_phi_deta1_gap_norm","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta1_gap_norm_data, delta_phi_deta1_gap_norm_unc, delta_phi_deta1_gap_norm, detail);
    plot_2histograms(delta_phi_deta1_gap_norm_data, "Only statistical uncertainty", delta_phi_deta1_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta1_gap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi Deta2 Gap"<<endl; }

    TH1D *delta_phi_deta2_gap_data = 0;
    TH1D *delta_phi_deta2_gap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta2_gap",delta_phi_deta2_gap_unc);
    if (delta_phi_deta2_gap_unc == 0) { cout << "uncorr_unc_delta_phi_deta2_gap not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D("ak5PF_delta_phi_deta2_gap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta2_gap_data, delta_phi_deta2_gap_unc, delta_phi_deta2_gap, detail);
    plot_2histograms(delta_phi_deta2_gap_data, "Only statistical uncertainty", delta_phi_deta2_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta2_gap", "top_left", true, true);


//apply uncertainty in delta phi deta2 gap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Gap Norm"<<endl; }

    TH1D *delta_phi_deta2_gap_norm_data = 0;
    TH1D *delta_phi_deta2_gap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_data);
    if (delta_phi_deta2_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta2_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta2_gap_norm",delta_phi_deta2_gap_norm_unc);
    if (delta_phi_deta2_gap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta2_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm =  new TH1D("ak5PF_delta_phi_deta2_gap_norm","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta2_gap_norm_data, delta_phi_deta2_gap_norm_unc, delta_phi_deta2_gap_norm, detail);
    plot_2histograms(delta_phi_deta2_gap_norm_data, "Only statistical uncertainty", delta_phi_deta2_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta2_gap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi Deta3 Gap"<<endl; }

    TH1D *delta_phi_deta3_gap_data = 0;
    TH1D *delta_phi_deta3_gap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta3_gap",delta_phi_deta3_gap_unc);
    if (delta_phi_deta3_gap_unc == 0) { cout << "uncorr_unc_delta_phi_deta3_gap not found!" << endl; return; }

    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D("ak5PF_delta_phi_deta3_gap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta3_gap_data, delta_phi_deta3_gap_unc, delta_phi_deta3_gap, detail);
    plot_2histograms(delta_phi_deta3_gap_data, "Only statistical uncertainty", delta_phi_deta3_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta3_gap", "top_left", true, true);


//apply uncertainty in delta phi deta3 gap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Gap Norm"<<endl; }

    TH1D *delta_phi_deta3_gap_norm_data = 0;
    TH1D *delta_phi_deta3_gap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_data);
    if (delta_phi_deta3_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta3_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta3_gap_norm",delta_phi_deta3_gap_norm_unc);
    if (delta_phi_deta3_gap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta3_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm =  new TH1D("ak5PF_delta_phi_deta3_gap_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta3_gap_norm_data, delta_phi_deta3_gap_norm_unc, delta_phi_deta3_gap_norm, detail);
    plot_2histograms(delta_phi_deta3_gap_norm_data, "Only statistical uncertainty", delta_phi_deta3_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta3_gap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi Deta4 Gap"<<endl; }

    TH1D *delta_phi_deta4_gap_data = 0;
    TH1D *delta_phi_deta4_gap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta4_gap",delta_phi_deta4_gap_unc);
    if (delta_phi_deta4_gap_unc == 0) { cout << "uncorr_unc_delta_phi_deta4_gap not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D("ak5PF_delta_phi_deta4_gap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta4_gap_data, delta_phi_deta4_gap_unc, delta_phi_deta4_gap, detail);
    plot_2histograms(delta_phi_deta4_gap_data, "Only statistical uncertainty", delta_phi_deta4_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta4_gap", "top_left", true, true);


//apply uncertainty in delta phi deta4 gap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Gap Norm"<<endl; }

    TH1D *delta_phi_deta4_gap_norm_data = 0;
    TH1D *delta_phi_deta4_gap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_data);
    if (delta_phi_deta4_gap_norm_data == 0) { cout << "ak5PF_delta_phi_deta4_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta4_gap_norm",delta_phi_deta4_gap_norm_unc);
    if (delta_phi_deta4_gap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta4_gap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm =  new TH1D("ak5PF_delta_phi_deta4_gap_norm","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta4_gap_norm_data, delta_phi_deta4_gap_norm_unc, delta_phi_deta4_gap_norm, detail);
    plot_2histograms(delta_phi_deta4_gap_norm_data, "Only statistical uncertainty", delta_phi_deta4_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta4_gap_norm", "top_left", true, true);


//apply uncertainty in delta phi nogap distribution
    if (detail) { cout<<"Delta phi Nogap"<<endl; }

    TH1D *delta_phi_nogap_data = 0;
    TH1D *delta_phi_nogap_unc = 0;

    data->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_nogap",delta_phi_nogap_unc);
    if (delta_phi_nogap_unc == 0) { cout << "uncorr_unc_delta_phi_nogap not found!" << endl; return; }

    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_nogap_data, delta_phi_nogap_unc, delta_phi_nogap, detail);
    plot_2histograms(delta_phi_nogap_data, "Only statistical uncertainty", delta_phi_nogap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_nogap", "top_left", true, true);


//apply uncertainty in delta phi nogap norm distribution
    if (detail) { cout<<"Delta phi Nogap Norm"<<endl; }

    TH1D *delta_phi_nogap_norm_data = 0;
    TH1D *delta_phi_nogap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_nogap_norm",delta_phi_nogap_norm_data);
    if (delta_phi_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_nogap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_nogap_norm",delta_phi_nogap_norm_unc);
    if (delta_phi_nogap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D("ak5PF_delta_phi_nogap_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_nogap_norm_data, delta_phi_nogap_norm_unc, delta_phi_nogap_norm, detail);
    plot_2histograms(delta_phi_nogap_norm_data, "Only statistical uncertainty", delta_phi_nogap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_nogap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_data = 0;
    TH1D *delta_phi_deta1_nogap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta1_nogap",delta_phi_deta1_nogap_unc);
    if (delta_phi_deta1_nogap_unc == 0) { cout << "uncorr_unc_delta_phi_deta1_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D("ak5PF_delta_phi_deta1_nogap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_unc, delta_phi_deta1_nogap, detail);
    plot_2histograms(delta_phi_deta1_nogap_data, "Only statistical uncertainty", delta_phi_deta1_nogap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta1_nogap", "top_left", true, true);


//apply uncertainty in delta phi deta1 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta1 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta1_nogap_norm_data = 0;
    TH1D *delta_phi_deta1_nogap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_data);
    if (delta_phi_deta1_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta1_nogap_norm",delta_phi_deta1_nogap_norm_unc);
    if (delta_phi_deta1_nogap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta1_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm =  new TH1D("ak5PF_delta_phi_deta1_nogap_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta1_nogap_norm_data, delta_phi_deta1_nogap_norm_unc, delta_phi_deta1_nogap_norm, detail);
    plot_2histograms(delta_phi_deta1_nogap_norm_data, "Only statistical uncertainty", delta_phi_deta1_nogap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta1_nogap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_data = 0;
    TH1D *delta_phi_deta2_nogap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta2_nogap",delta_phi_deta2_nogap_unc);
    if (delta_phi_deta2_nogap_unc == 0) { cout << "uncorr_unc_delta_phi_deta2_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D("ak5PF_delta_phi_deta2_nogap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_unc, delta_phi_deta2_nogap, detail);
    plot_2histograms(delta_phi_deta2_nogap_data, "Only statistical uncertainty", delta_phi_deta2_nogap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta2_nogap", "top_left", true, true);


//apply uncertainty in delta phi deta2 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta2 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta2_nogap_norm_data = 0;
    TH1D *delta_phi_deta2_nogap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_data);
    if (delta_phi_deta2_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta2_nogap_norm",delta_phi_deta2_nogap_norm_unc);
    if (delta_phi_deta2_nogap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta2_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm =  new TH1D("ak5PF_delta_phi_deta2_nogap_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta2_nogap_norm_data, delta_phi_deta2_nogap_norm_unc, delta_phi_deta2_nogap_norm, detail);
    plot_2histograms(delta_phi_deta2_nogap_norm_data, "Only statistical uncertainty", delta_phi_deta2_nogap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta2_nogap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_data = 0;
    TH1D *delta_phi_deta3_nogap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta3_nogap",delta_phi_deta3_nogap_unc);
    if (delta_phi_deta3_nogap_unc == 0) { cout << "uncorr_unc_delta_phi_deta3_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D("ak5PF_delta_phi_deta3_nogap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_unc, delta_phi_deta3_nogap, detail);
    plot_2histograms(delta_phi_deta3_nogap_data, "Only statistical uncertainty", delta_phi_deta3_nogap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta3_nogap", "top_left", true, true);


//apply uncertainty in delta phi deta3 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta3 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta3_nogap_norm_data = 0;
    TH1D *delta_phi_deta3_nogap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_data);
    if (delta_phi_deta3_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta3_nogap_norm",delta_phi_deta3_nogap_norm_unc);
    if (delta_phi_deta3_nogap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta3_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm =  new TH1D("ak5PF_delta_phi_deta3_nogap_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta3_nogap_norm_data, delta_phi_deta3_nogap_norm_unc, delta_phi_deta3_nogap_norm, detail);
    plot_2histograms(delta_phi_deta3_nogap_norm_data, "Only statistical uncertainty", delta_phi_deta3_nogap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta3_nogap_norm", "top_left", true, true);


//apply uncertainty in delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_data = 0;
    TH1D *delta_phi_deta4_nogap_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta4_nogap",delta_phi_deta4_nogap_unc);
    if (delta_phi_deta4_nogap_unc == 0) { cout << "uncorr_unc_delta_phi_deta4_nogap not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D("ak5PF_delta_phi_deta4_nogap","Corrected Data;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_unc, delta_phi_deta4_nogap, detail);
    plot_2histograms(delta_phi_deta4_nogap_data, "Only statistical uncertainty", delta_phi_deta4_nogap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta4_nogap", "top_left", true, true);


//apply uncertainty in delta phi deta4 nogap norm distribution
    if (detail) { cout<<"Delta phi Deta4 Nogap Norm"<<endl; }

    TH1D *delta_phi_deta4_nogap_norm_data = 0;
    TH1D *delta_phi_deta4_nogap_norm_unc = 0;

    data->GetObject("ak5PF_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_data);
    if (delta_phi_deta4_nogap_norm_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_phi_deta4_nogap_norm",delta_phi_deta4_nogap_norm_unc);
    if (delta_phi_deta4_nogap_norm_unc == 0) { cout << "uncorr_unc_delta_phi_deta4_nogap_norm not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm =  new TH1D("ak5PF_delta_phi_deta4_nogap_norm","Corrected Data Norm;#Delta#phi [rad];d#sigma/d#Delta#phi [pb/rad]", dphi_nbins, dphi_bins);

    set_uncorr_unc(delta_phi_deta4_nogap_norm_data, delta_phi_deta4_nogap_norm_unc, delta_phi_deta4_nogap_norm, detail);
    plot_2histograms(delta_phi_deta4_nogap_norm_data, "Only statistical uncertainty", delta_phi_deta4_nogap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_phi_deta4_nogap_norm", "top_left", true, true);


//apply uncertainty in leading pt inside gap distribution
    if (detail) { cout<<"Leading pT Inside Gap"<<endl; }

    TH1D *leading_pt_inside_gap_data = 0;
    TH1D *leading_pt_inside_gap_unc = 0;

    data->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_leading_pt_inside_gap",leading_pt_inside_gap_unc);
    if (leading_pt_inside_gap_unc == 0) { cout << "uncorr_unc_leading_pt_inside_gap not found!" << endl; return; }

    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D("ak5PF_leading_pt_inside_gap","Corrected Data;p_{T}^{inside} [rad];d#sigma/dp_{T}^{inside} [pb/GeV]", in_nbins, in_bins);

    set_uncorr_unc(leading_pt_inside_gap_data, leading_pt_inside_gap_unc, leading_pt_inside_gap, detail);
    plot_2histograms(leading_pt_inside_gap_data, "Only statistical uncertainty", leading_pt_inside_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_leading_pt_inside_gap", "top_right", true, true);


//apply uncertainty in leading pt inside gap norm distribution
    if (detail) { cout<<"Leading pT Inside Gap Norm"<<endl; }

    TH1D *leading_pt_inside_gap_norm_data = 0;
    TH1D *leading_pt_inside_gap_norm_unc = 0;

    data->GetObject("ak5PF_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_data);
    if (leading_pt_inside_gap_norm_data == 0) { cout << "ak5PF_leading_pt_inside_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_leading_pt_inside_gap_norm",leading_pt_inside_gap_norm_unc);
    if (leading_pt_inside_gap_norm_unc == 0) { cout << "uncorr_unc_leading_pt_inside_gap_norm not found!" << endl; return; }

    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm =  new TH1D("ak5PF_leading_pt_inside_gap_norm","Corrected Data Norm;p_{T}^{inside} [rad];d#sigma/dp_{T}^{inside} [pb/GeV]", in_nbins, in_bins);

    set_uncorr_unc(leading_pt_inside_gap_norm_data, leading_pt_inside_gap_norm_unc, leading_pt_inside_gap_norm, detail);
    plot_2histograms(leading_pt_inside_gap_norm_data, "Only statistical uncertainty", leading_pt_inside_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_leading_pt_inside_gap_norm", "top_right", true, true);


//apply uncertainty in leading eta star inside gap distribution
    if (detail) { cout<<"Leading Eta Star Inside Gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_data = 0;
    TH1D *leading_eta_star_inside_gap_unc = 0;

    data->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_leading_eta_star_inside_gap",leading_eta_star_inside_gap_unc);
    if (leading_eta_star_inside_gap_unc == 0) { cout << "uncorr_unc_leading_eta_star_inside_gap not found!" << endl; return; }

    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D("ak5PF_leading_eta_star_inside_gap","Corrected Data;#Eta*;d#sigma/d#Eta* [pb]", etastar_nbins, etastar_bins);

    set_uncorr_unc(leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_unc, leading_eta_star_inside_gap, detail);
    plot_2histograms(leading_eta_star_inside_gap_data, "Only statistical uncertainty", leading_eta_star_inside_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_leading_eta_star_inside_gap", "top_right", true, true);


//apply uncertainty in leading eta star inside gap norm distribution
    if (detail) { cout<<"Leading Eta Star Inside Gap Norm"<<endl; }

    TH1D *leading_eta_star_inside_gap_norm_data = 0;
    TH1D *leading_eta_star_inside_gap_norm_unc = 0;

    data->GetObject("ak5PF_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_data);
    if (leading_eta_star_inside_gap_norm_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_leading_eta_star_inside_gap_norm",leading_eta_star_inside_gap_norm_unc);
    if (leading_eta_star_inside_gap_norm_unc == 0) { cout << "uncorr_unc_leading_eta_star_inside_gap_norm not found!" << endl; return; }

    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm =  new TH1D("ak5PF_leading_eta_star_inside_gap_norm","Corrected Data Norm;#Eta*;d#sigma/d#Eta* [pb]", etastar_nbins, etastar_bins);

    set_uncorr_unc(leading_eta_star_inside_gap_norm_data, leading_eta_star_inside_gap_norm_unc, leading_eta_star_inside_gap_norm, detail);
    plot_2histograms(leading_eta_star_inside_gap_norm_data, "Only statistical uncertainty", leading_eta_star_inside_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_leading_eta_star_inside_gap_norm", "top_right", true, true);


//apply uncertainty in delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta Outside Gap"<<endl; }

    TH1D *delta_eta_outside_gap_data = 0;
    TH1D *delta_eta_outside_gap_unc = 0;

    data->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_eta_outside_gap",delta_eta_outside_gap_unc);
    if (delta_eta_outside_gap_unc == 0) { cout << "uncorr_unc_delta_eta_outside_gap not found!" << endl; return; }

    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D("ak5PF_delta_eta_outside_gap","Corrected Data;#Delta#eta^{out};d#sigma/d#Delta#eta^{out} [pb]", deta_out_nbins, deta_out_bins);

    set_uncorr_unc(delta_eta_outside_gap_data, delta_eta_outside_gap_unc, delta_eta_outside_gap, detail);
    plot_2histograms(delta_eta_outside_gap_data, "Only statistical uncertainty", delta_eta_outside_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_eta_outside_gap", "top_right", true, true);


//apply uncertainty in delta eta outside gap norm distribution
    if (detail) { cout<<"Delta Eta Outside Gap Norm"<<endl; }

    TH1D *delta_eta_outside_gap_norm_data = 0;
    TH1D *delta_eta_outside_gap_norm_unc = 0;

    data->GetObject("ak5PF_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_data);
    if (delta_eta_outside_gap_norm_data == 0) { cout << "ak5PF_delta_eta_outside_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_delta_eta_outside_gap_norm",delta_eta_outside_gap_norm_unc);
    if (delta_eta_outside_gap_norm_unc == 0) { cout << "uncorr_unc_delta_eta_outside_gap_norm not found!" << endl; return; }

    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm =  new TH1D("ak5PF_delta_eta_outside_gap_norm","Corrected Data Norm;#Delta#eta^{out};d#sigma/d#Delta#eta^{out} [pb]", deta_out_nbins, deta_out_bins);

    set_uncorr_unc(delta_eta_outside_gap_norm_data, delta_eta_outside_gap_norm_unc, delta_eta_outside_gap_norm, detail);
    plot_2histograms(delta_eta_outside_gap_norm_data, "Only statistical uncertainty", delta_eta_outside_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_delta_eta_outside_gap_norm", "top_right", true, true);


//apply uncertainty in leading pt outside gap distribution
    if (detail) { cout<<"Leading pT Outside Gap"<<endl; }

    TH1D *leading_pt_outside_gap_data = 0;
    TH1D *leading_pt_outside_gap_unc = 0;

    data->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    unc->GetObject("uncorr_unc_leading_pt_outside_gap",leading_pt_outside_gap_unc);
    if (leading_pt_outside_gap_unc == 0) { cout << "uncorr_unc_leading_pt_outside_gap not found!" << endl; return; }

    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D("ak5PF_leading_pt_outside_gap","Corrected Data;p_{T}^{outside} [GeV];d#sigma/dp_{T}^{outside} [pb/GeV]", out_nbins, out_bins);

    set_uncorr_unc(leading_pt_outside_gap_data, leading_pt_outside_gap_unc, leading_pt_outside_gap, detail);
    plot_2histograms(leading_pt_outside_gap_data, "Only statistical uncertainty", leading_pt_outside_gap, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_leading_pt_outside_gap", "top_right", true, true);


//apply uncertainty in leading pt outside gap norm distribution
    if (detail) { cout<<"Leading pT Outside Gap Norm"<<endl; }

    TH1D *leading_pt_outside_gap_norm_data = 0;
    TH1D *leading_pt_outside_gap_norm_unc = 0;

    data->GetObject("ak5PF_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_data);
    if (leading_pt_outside_gap_norm_data == 0) { cout << "ak5PF_leading_pt_outside_gap_norm not found!" << endl; return; }
    unc->GetObject("uncorr_unc_leading_pt_outside_gap_norm",leading_pt_outside_gap_norm_unc);
    if (leading_pt_outside_gap_norm_unc == 0) { cout << "uncorr_unc_leading_pt_outside_gap_norm not found!" << endl; return; }

    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm =  new TH1D("ak5PF_leading_pt_outside_gap_norm","Corrected Data Norm;p_{T}^{outside} [GeV];d#sigma/dp_{T}^{outside} [pb/GeV]", out_nbins, out_bins);

    set_uncorr_unc(leading_pt_outside_gap_norm_data, leading_pt_outside_gap_norm_unc, leading_pt_outside_gap_norm, detail);
    plot_2histograms(leading_pt_outside_gap_norm_data, "Only statistical uncertainty", leading_pt_outside_gap_norm, "Uncorrelated uncertainty", output_path_plots, "data_uncorr_unc_leading_pt_outside_gap_norm", "top_right", true, true);


//Opening the output root file
    if (detail) { cout<<"Creating " << path_out << "..."<<endl; }
    TFile *output = TFile::Open( path_out.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file..."<<endl; }
    delta_phi->Write();
    delta_phi_deta1->Write();
    delta_phi_deta2->Write();
    delta_phi_deta3->Write();
    delta_phi_deta4->Write();
    delta_phi_gap->Write();
    delta_phi_deta1_gap->Write();
    delta_phi_deta2_gap->Write();
    delta_phi_deta3_gap->Write();
    delta_phi_deta4_gap->Write();
    delta_phi_nogap->Write();
    delta_phi_deta1_nogap->Write();
    delta_phi_deta2_nogap->Write();
    delta_phi_deta3_nogap->Write();
    delta_phi_deta4_nogap->Write();
    leading_pt_inside_gap->Write();
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();
    leading_pt_outside_gap->Write();

    delta_phi_norm->Write();
    delta_phi_deta1_norm->Write();
    delta_phi_deta2_norm->Write();
    delta_phi_deta3_norm->Write();
    delta_phi_deta4_norm->Write();
    delta_phi_gap_norm->Write();
    delta_phi_deta1_gap_norm->Write();
    delta_phi_deta2_gap_norm->Write();
    delta_phi_deta3_gap_norm->Write();
    delta_phi_deta4_gap_norm->Write();
    delta_phi_nogap_norm->Write();
    delta_phi_deta1_nogap_norm->Write();
    delta_phi_deta2_nogap_norm->Write();
    delta_phi_deta3_nogap_norm->Write();
    delta_phi_deta4_nogap_norm->Write();
    leading_pt_inside_gap_norm->Write();
    leading_eta_star_inside_gap_norm->Write();
    delta_eta_outside_gap_norm->Write();
    leading_pt_outside_gap_norm->Write();
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data->Close();
    unc->Close();
    output->Close();

//deleting all histograms to avoid memory leak - causes a memory leak!
/*    delete(delta_phi);
    delete(delta_phi_deta1);
    delete(delta_phi_deta2);
    delete(delta_phi_deta3);
    delete(delta_phi_deta4);
    delete(delta_phi_gap);
    delete(delta_phi_deta1_gap);
    delete(delta_phi_deta2_gap);
    delete(delta_phi_deta3_gap);
    delete(delta_phi_deta4_gap);
    delete(delta_phi_nogap);
    delete(delta_phi_deta1_nogap);
    delete(delta_phi_deta2_nogap);
    delete(delta_phi_deta3_nogap);
    delete(delta_phi_deta4_nogap);
    delete(leading_pt_inside_gap);
    delete(leading_pt_outside_gap);
    delete(delta_eta_outside_gap);
    delete(leading_eta_star_inside_gap); */

    if (detail) { cout<<"Done!"<<endl; }

}
