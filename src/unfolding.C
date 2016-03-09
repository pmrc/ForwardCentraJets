// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 17 Jun 2013
//
// unfolding(string response_file, string path_data, string mc_file, string root_out, bool detail = false)
// unfolds the delta_phi distribution present in the path_data with TUnfold method using the response matrix suplied by response_file and saves the result in root_out

#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TPad.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TSpline.h>
#include "TLegend.h"
#include "TDecompSVD.h"

#include "../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldTUnfold.h" 

#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include "common_methods.h"

// Global definitions
const int n_loop = 10000;
const int seed = 5;
//TUnfold::ERegMode regmode=TUnfold::kRegModeNone;
TUnfold::ERegMode regmode=TUnfold::kRegModeSize;
//TUnfold::ERegMode regmode=TUnfold::kRegModeDerivative;
//TUnfold::ERegMode regmode=TUnfold::kRegModeCurvature;

void plot_covariance(TMatrixD cov, string path, string prefix, string name, bool detail = false)
{
//declaring the canvas
    if (detail) { cout << "Ploting " << name << endl; }
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogz();

    cov.Draw("colz");
    cov.Draw("text same");


//setting the output files
   string fileout = prefix + name;
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c1->Print( out_png.c_str() );
    c1->Print( out_c.c_str() );
    c1->Print( out_eps.c_str() );
    c1->Close();

}


void plot_lcurve(TGraph *lcurve, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{
//declaring the canvas
    if (detail) { cout << "Ploting " << name << endl; }
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    //gPad->SetLogy();

    lcurve->SetTitle("L-curve;log(#chi^2);log(reg_cond)");
    lcurve->Draw("AC*");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1, y1, x2, y2);
    leg00->AddEntry(lcurve,"L-curve","l");
    leg00->SetFillColor(0);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();


//setting the output files
   string fileout = prefix + name;
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c1->Print( out_png.c_str() );
    c1->Print( out_c.c_str() );
    c1->Print( out_eps.c_str() );
    c1->Close();

}


void plot_taucurves(TSpline *logTauX, TSpline *logTauY, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{
    if (detail) { cout << "Ploting " << name << endl; }
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    //gPad->SetLogy();
    c1->Divide(1,2);

    c1->cd(1);
    logTauX->SetTitle("#tau as function of #chi^2;log(#tau);log(#chi^2)");
    logTauX->Draw("");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1, y1, x2, y2);
    leg00->AddEntry(logTauX,"tau-curves","l");
    leg00->SetFillColor(0);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();


    c1->cd(2);
    logTauY->SetTitle("#tau as function of regularization condition;log(#tau);log(#reg_cond)");
    logTauY->Draw("");


//setting the output files
   string fileout = prefix + name;
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c1->Print( out_png.c_str() );
    c1->Print( out_c.c_str() );
    c1->Print( out_eps.c_str() );
    c1->Close();
}

void plot_unfolded_result(TH1D *unfolded, TH1D *true_dist, TH1D *measured_dist, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{

//declaring the canvas
    if (detail) { cout << "Ploting " << name << endl; }
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogy();

for (int j = 1; j <= true_dist->GetNbinsX(); j++)
	{
	//cout << j << " : " << true_dist->GetBinContent(j) << " +- " << true_dist->GetBinError(j) << endl;
	}

for (int j = 1; j <= measured_dist->GetNbinsX(); j++)
	{
	//cout << j << " : " << measured_dist->GetBinContent(j) << " +- " << measured_dist->GetBinError(j) << endl;
	}

//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 0.0;
    double max = unfolded->GetMaximum();
    if (unfolded->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(unfolded,detail);
    }
    else
    {
    min = unfolded->GetMinimum();
    } 

    set_histogram_min_max(true_dist, min, max, detail);
    set_histogram_min_max(measured_dist, min, max, detail);
//    set_histogram_min_max(measured_plot, min, max, detail);

    max = 1.3 * max;
    min = 0.7 * min;

//format and ploting the histogram
    if (detail) { cout << "Drawning on the canvas..." << endl; }
    unfolded->SetMaximum(max);
    unfolded->SetMinimum(min);
    unfolded->SetLineColor(2);
    unfolded->SetLineStyle(1);
    unfolded->SetLineWidth(3);
    unfolded->Draw("e1");
    true_dist->SetLineColor(3);
    true_dist->SetLineStyle(2);
    true_dist->SetLineWidth(3);
    true_dist->Draw("e1 same");
    measured_dist->SetLineColor(4);
    measured_dist->SetLineStyle(3);
    measured_dist->SetLineWidth(3);
    measured_dist->Draw("e1 same");
//    measured_plot->SetLineColor(6);
//    measured_plot->SetLineStyle(4);
//    measured_plot->SetLineWidth(3);
//    measured_plot->Draw("e1 same");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1,y1,x2,y2);
    leg00->AddEntry(unfolded,"Unfolded result","l");
    leg00->AddEntry(true_dist,"True distribution","l");
    leg00->AddEntry(measured_dist,"Measured distribution","l");
//    leg00->AddEntry(measured_plot,"Measured distribution","l");
    leg00->SetFillColor(0);
    leg00->SetLineStyle(1);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();

//setting the output files
   print_plots(c1, path, prefix + name);

}


void set_input(TH1 * detected, TH1 *input_detected)
{

input_detected->SetBinContent(1,0.0);
input_detected->SetBinError(1,0.0);

for (int a = 1; a <= detected->GetNbinsX();a++)
	{
	input_detected->SetBinContent(a+1,detected->GetBinContent(a));
	input_detected->SetBinError(a+1,detected->GetBinError(a));
	cout << "Bin : " << a << " -> " << detected->GetBinContent(a) << " +- " << detected->GetBinError(a) << endl;
	}

}

void unfold_distribution(TH1D *detected, TH1D *input_detected, TH1D *result, TH1D *output_true, TH1D *generated, RooUnfoldResponse *response, RooUnfoldResponse *response_all, string plots_path, string prefix, string name, string method, double scale, double& tau, bool detail = false)
{

//setup
// 0 - empty dist / 1 - too high at low bins / 2 - too low at lower bins
// 3 - too low at lower bins / 4 - too low at lower bins / 5 - too low at lower bins
// 6 - too low at lower bins / 7 - too low at lower bins / 8 - too low at lower bins 
// 9 - too low at lower bins / 10 - too low at lower bins / 11 - too low at lower bins 
int bayes_iterations = 4; //small int bigger than 0
// 1 - too low at high bins / 2 - too low at high bins / 3 - 
int svd_regularization = 1; //between 1 and number of bins

//get entries
int num = detected->GetEntries();

//initialize the output objects
TMatrixD cov(7,7);
TMatrixD invertida(7,7);
TMatrixD chk(7,7);
TGraph *lcurve = 0;
TSpline *logTauX = 0, *logTauY = 0;
int fac = 0;

//  Creating the RooUnfoldTUnfold object
if (method == "TUnfold")
	{
	RooUnfoldTUnfold unfold(response, detected, regmode);
	result = (TH1D*) unfold.Hreco();
	tau = unfold.GetTau();
	cout<<"Tau : "<<tau<<endl;
	lcurve = static_cast<TGraph*>(unfold.GetLCurve()->Clone());
	logTauX = static_cast<TSpline*>(unfold.GetLogTauX()->Clone());
	logTauY = static_cast<TSpline*>(unfold.GetLogTauY()->Clone());
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (method == "BinByBin")
	{
	RooUnfoldBinByBin unfold(response_all, input_detected);
	result = (TH1D*) unfold.Hreco();
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (method == "SVD")
	{
	RooUnfoldSvd unfold(response_all, input_detected, svd_regularization);
	result = (TH1D*) unfold.Hreco();
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (method == "Bayes")
	{
	RooUnfoldBayes unfold(response, detected, bayes_iterations);
	result = (TH1D*) unfold.Hreco();
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (result == 0) { cout << "Invalid unfolding method! Routine termination!!!" << endl; return; }

// unfolding


//because we input the double of the bins the result is the double of normal, we need to scale the histogram
result->Scale(scale);

//matrix investion and cross-check multiplication
invertida = cov;
invertida.Invert();
chk = cov * invertida;

for (int a = 1; a <= output_true->GetNbinsX();a++)
	{
	if (method == "TUnfold" or method == "Bayes") { fac = a;} else { fac = a+1; }
	output_true->SetBinContent(a,result->GetBinContent(fac));
	output_true->SetBinError(a,result->GetBinError(fac));
	if (detail)
	  { cout << "Bin : " << a << " -> " << output_true->GetBinContent(a) << " +- " << output_true->GetBinError(a) << endl; }
	}

//set the correct number of entries for merging porpuses
result->SetEntries(num);
output_true->SetEntries(num);

//ploting the result
if (method == "TUnfold" or method == "Bayes")
	{
	plot_unfolded_result(result, generated, detected, plots_path, prefix, name, "top_left", detail);

	}
else
	{
	plot_unfolded_result(output_true, generated, detected, plots_path, prefix, name, "top_left", detail);
	}

if (method == "TUnfold")
	{
	plot_lcurve(lcurve, plots_path, prefix, name+"lcurve", "top_left", detail);
	plot_taucurves(logTauX, logTauY, plots_path, prefix, name+"taucurves", "top_left", detail);
	}

plot_covariance(cov, plots_path, prefix, name+"covariance", detail);
plot_covariance(invertida, plots_path, prefix, name+"inv_covariance", detail);
plot_covariance(chk, plots_path, prefix, name+"inv_covariance_mult_convariance", detail);


}

void unfolding(string response_file, string data_file, string mc_file, string root_out, string method, double scale, double& tau, string plots_path, string prefix, bool true_is_corr_data = false, bool detail = false, bool test = false)
{

if (detail)
{
cout << "Unfolding Configuration" << endl;
cout << "Response File :  " << response_file << endl;
cout << "Data File :      " << data_file << endl;
cout << "MC File :        " << mc_file << endl;
cout << "Root output :    " << root_out << endl;
cout << "Method :         " << method << endl;
cout << "Scaling Factor : " << scale << endl;
cout << "Plots Path :     " << plots_path << endl;
cout << "Prefix :         " << prefix << endl;
cout << "Detail :         " << detail << endl;
cout << "Test Mode :      " << test << endl;
}

//opening the files
TFile *respfile = new TFile( response_file.c_str() );
TFile *datafile = new TFile( data_file.c_str() );
TFile *mcfile = new TFile( mc_file.c_str() );

//binning
int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};
int dphi_nbins_meas = 15;
double dphi_bins_meas[16] = {-1.0, 0.0, 0.225, 0.45, 0.675, 0.9, 1.125, 1.35, 1.575, 1.8, 2.025, 2.25, 2.475, 2.7, 2.925, 3.15};

int in_nbins = 9;
double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
int in_nbins_meas = 19;
double in_bins_meas[20] = {0, 20, 24, 27, 31, 35, 40, 45, 51, 57, 65, 72, 80, 90, 105, 120, 130, 150, 170, 200};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};
int etastar_nbins_meas = 25;
double etastar_bins_meas[26] = {-10,-3.6, -3.0, -2.5, -2.3, -2.0, -1.7, -1.5, -1.3, -1.0, -0.8, -0.5, -0.3, 0.0, 0.3, 0.5, 0.8, 1.0, 1.3, 1.5, 1.8, 2.0, 2.3, 2.5, 3.0, 3.6};

int out_nbins = 9;
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
int out_nbins_meas = 19;
double out_bins_meas[20] = {0, 20, 24, 27, 31, 35, 40, 45, 51, 57, 65, 72, 80, 90, 105, 120, 130, 150, 170, 200};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};
int deta_out_nbins_meas = 13;
double deta_out_bins_meas[14] = {-1.0, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7.5};

//Unfolding Delta Phi
TString delta_phi_hname = "resp_delta_phi";
TString delta_phi_hname_all = "resp_delta_phi_all";
TString delta_phi_measured = "ak5PF_delta_phi_fine";
TString delta_phi_measured_plot = "ak5PF_delta_phi";
TString delta_phi_true;
if (true_is_corr_data)
	{
	delta_phi_true = "ak5PF_delta_phi";
	}
else
	{
	delta_phi_true = "ak5Gen_delta_phi";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_response = 0;
respfile->GetObject(delta_phi_hname,delta_phi_response);
if (delta_phi_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_response_all = 0;
respfile->GetObject(delta_phi_hname_all,delta_phi_response_all);
if (delta_phi_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_detected = 0;
datafile->GetObject(delta_phi_measured,delta_phi_detected);
//respfile->GetObject("ak5PF_delta_phi_measured",detected);
if (delta_phi_detected == 0) { cout<<"Measured data " << delta_phi_measured << " not found!"<<endl; return; }

TH1D *delta_phi_detected_plot = 0;
datafile->GetObject(delta_phi_measured_plot,delta_phi_detected_plot);
if (delta_phi_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_generated = 0;
mcfile->GetObject(delta_phi_true,delta_phi_generated);
//respfile->GetObject("ak5Gen_delta_phi_truth",generated);
if (delta_phi_generated == 0) { cout<< delta_phi_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_input_detected;
  delta_phi_input_detected =  new TH1D("input_detected_delta_phi","#Delta#phi;|#Delta#phi|;Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_detected, delta_phi_input_detected);

  TH1D *delta_phi_output_true;
  delta_phi_output_true =  new TH1D("output_true_delta_phi","#Delta#phi;|#Delta#phi| [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_result = 0;

  unfold_distribution(delta_phi_detected, delta_phi_input_detected, delta_phi_result, delta_phi_output_true, delta_phi_generated, delta_phi_response, delta_phi_response_all, plots_path, prefix, "delta_phi", method, scale, tau, detail);


//Unfolding Delta Phi Norm
TString delta_phi_norm_measured = "ak5PF_delta_phi_norm_fine";
TString delta_phi_norm_measured_plot = "ak5PF_delta_phi_norm";
TString delta_phi_norm_true;
if (true_is_corr_data)
	{
	delta_phi_norm_true = "ak5PF_delta_phi_norm";
	}
else
	{
	delta_phi_norm_true = "ak5Gen_delta_phi_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_norm_detected = 0;
datafile->GetObject(delta_phi_norm_measured,delta_phi_norm_detected);
if (delta_phi_norm_detected == 0) { cout<<"Measured data " << delta_phi_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_norm_detected_plot = 0;
datafile->GetObject(delta_phi_norm_measured_plot,delta_phi_norm_detected_plot);
if (delta_phi_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_norm_generated = 0;
mcfile->GetObject(delta_phi_norm_true,delta_phi_norm_generated);
if (delta_phi_norm_generated == 0) { cout<< delta_phi_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_norm_input_detected;
  delta_phi_norm_input_detected =  new TH1D("input_detected_delta_phi_norm","#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi} [rad^{-1}]", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_norm_detected, delta_phi_norm_input_detected);

  TH1D *delta_phi_norm_output_true;
  delta_phi_norm_output_true =  new TH1D("output_true_delta_phi_norm","#Delta#phi;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi} [rad^{-1}]", dphi_nbins, dphi_bins);

  TH1D *delta_phi_norm_result = 0;

  unfold_distribution(delta_phi_norm_detected, delta_phi_norm_input_detected, delta_phi_norm_result, delta_phi_norm_output_true, delta_phi_norm_generated, delta_phi_response, delta_phi_response_all, plots_path, prefix, "delta_phi_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta1
TString delta_phi_deta1_hname = "resp_delta_phi_deta1";
TString delta_phi_deta1_hname_all = "resp_delta_phi_deta1_all";
TString delta_phi_deta1_measured = "ak5PF_delta_phi_deta1_fine";
TString delta_phi_deta1_measured_plot = "ak5PF_delta_phi_deta1";
TString delta_phi_deta1_true;
if (true_is_corr_data)
	{
	delta_phi_deta1_true = "ak5PF_delta_phi_deta1";
	}
else
	{
	delta_phi_deta1_true = "ak5Gen_delta_phi_deta1";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta1_response = 0;
respfile->GetObject(delta_phi_deta1_hname,delta_phi_deta1_response);
if (delta_phi_deta1_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta1_response_all = 0;
respfile->GetObject(delta_phi_deta1_hname_all,delta_phi_deta1_response_all);
if (delta_phi_deta1_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta1_detected = 0;
datafile->GetObject(delta_phi_deta1_measured,delta_phi_deta1_detected);
//respfile->GetObject("ak5PF_delta_phi_measured",detected);
if (delta_phi_deta1_detected == 0) { cout<<"Measured data " << delta_phi_deta1_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta1_detected_plot = 0;
datafile->GetObject(delta_phi_deta1_measured_plot,delta_phi_deta1_detected_plot);
if (delta_phi_deta1_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta1_generated = 0;
mcfile->GetObject(delta_phi_deta1_true,delta_phi_deta1_generated);
if (delta_phi_deta1_generated == 0) { cout<< delta_phi_deta1_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta1_input_detected;
  delta_phi_deta1_input_detected =  new TH1D("input_detected_delta_phi_deta1","#Delta#phi Deta1;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta1_detected, delta_phi_deta1_input_detected);

  TH1D *delta_phi_deta1_output_true;
  delta_phi_deta1_output_true =  new TH1D("output_true_delta_phi_deta1","#Delta#phi Deta1;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta1_result = 0;

  unfold_distribution(delta_phi_deta1_detected, delta_phi_deta1_input_detected, delta_phi_deta1_result, delta_phi_deta1_output_true, delta_phi_deta1_generated, delta_phi_deta1_response, delta_phi_deta1_response_all, plots_path, prefix, "delta_phi_deta1", method, scale, tau, detail);


//Unfolding Delta Phi Deta1 Norm
TString delta_phi_deta1_norm_measured = "ak5PF_delta_phi_deta1_norm_fine";
TString delta_phi_deta1_norm_measured_plot = "ak5PF_delta_phi_deta1_norm";
TString delta_phi_deta1_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta1_norm_true = "ak5PF_delta_phi_deta1_norm";
	}
else
	{
	delta_phi_deta1_norm_true = "ak5Gen_delta_phi_deta1_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta1_norm_detected = 0;
datafile->GetObject(delta_phi_deta1_norm_measured,delta_phi_deta1_norm_detected);
if (delta_phi_deta1_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta1_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta1_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta1_norm_measured_plot,delta_phi_deta1_norm_detected_plot);
if (delta_phi_deta1_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta1_norm_generated = 0;
mcfile->GetObject(delta_phi_deta1_norm_true,delta_phi_deta1_norm_generated);
if (delta_phi_deta1_norm_generated == 0) { cout<< delta_phi_deta1_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta1_norm_input_detected;
  delta_phi_deta1_norm_input_detected =  new TH1D("input_detected_delta_phi_deta1_norm","#Delta#phi Deta1;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta1_norm_detected, delta_phi_deta1_norm_input_detected);

  TH1D *delta_phi_deta1_norm_output_true;
  delta_phi_deta1_norm_output_true =  new TH1D("output_true_delta_phi_deta1_norm","#Delta#phi Deta1;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta1_norm_result = 0;

  unfold_distribution(delta_phi_deta1_norm_detected, delta_phi_deta1_norm_input_detected, delta_phi_deta1_norm_result, delta_phi_deta1_norm_output_true, delta_phi_deta1_norm_generated, delta_phi_deta1_response, delta_phi_deta1_response_all, plots_path, prefix, "delta_phi_deta1_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta2
TString delta_phi_deta2_hname = "resp_delta_phi_deta2";
TString delta_phi_deta2_hname_all = "resp_delta_phi_deta2_all";
TString delta_phi_deta2_measured = "ak5PF_delta_phi_deta2_fine";
TString delta_phi_deta2_measured_plot = "ak5PF_delta_phi_deta2";
TString delta_phi_deta2_true;
if (true_is_corr_data)
	{
	delta_phi_deta2_true = "ak5PF_delta_phi_deta2";
	}
else
	{
	delta_phi_deta2_true = "ak5Gen_delta_phi_deta2";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta2_response = 0;
respfile->GetObject(delta_phi_deta2_hname,delta_phi_deta2_response);
if (delta_phi_deta2_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta2_response_all = 0;
respfile->GetObject(delta_phi_deta2_hname_all,delta_phi_deta2_response_all);
if (delta_phi_deta2_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta2_detected = 0;
datafile->GetObject(delta_phi_deta2_measured,delta_phi_deta2_detected);
if (delta_phi_deta2_detected == 0) { cout<<"Measured data " << delta_phi_deta2_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta2_detected_plot = 0;
datafile->GetObject(delta_phi_deta2_measured_plot,delta_phi_deta2_detected_plot);
if (delta_phi_deta2_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta2_generated = 0;
mcfile->GetObject(delta_phi_deta2_true,delta_phi_deta2_generated);
if (delta_phi_deta2_generated == 0) { cout<< delta_phi_deta2_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta2_input_detected;
  delta_phi_deta2_input_detected =  new TH1D("input_detected_delta_phi_deta2","#Delta#phi Deta2;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta2_detected, delta_phi_deta2_input_detected);

  TH1D *delta_phi_deta2_output_true;
  delta_phi_deta2_output_true =  new TH1D("output_true_delta_phi_deta2","#Delta#phi Deta2;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta2_result = 0;

  unfold_distribution(delta_phi_deta2_detected, delta_phi_deta2_input_detected, delta_phi_deta2_result, delta_phi_deta2_output_true, delta_phi_deta2_generated, delta_phi_deta2_response, delta_phi_deta2_response_all, plots_path, prefix, "delta_phi_deta2", method, scale, tau, detail);


//Unfolding Delta Phi Deta2 Norm
TString delta_phi_deta2_norm_measured = "ak5PF_delta_phi_deta2_norm_fine";
TString delta_phi_deta2_norm_measured_plot = "ak5PF_delta_phi_deta2_norm";
TString delta_phi_deta2_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta2_norm_true = "ak5PF_delta_phi_deta2_norm";
	}
else
	{
	delta_phi_deta2_norm_true = "ak5Gen_delta_phi_deta2_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta2_norm_detected = 0;
datafile->GetObject(delta_phi_deta2_norm_measured,delta_phi_deta2_norm_detected);
if (delta_phi_deta2_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta2_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta2_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta2_norm_measured_plot,delta_phi_deta2_norm_detected_plot);
if (delta_phi_deta2_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta2_norm_generated = 0;
mcfile->GetObject(delta_phi_deta2_norm_true,delta_phi_deta2_norm_generated);
if (delta_phi_deta2_norm_generated == 0) { cout<< delta_phi_deta2_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta2_norm_input_detected;
  delta_phi_deta2_norm_input_detected =  new TH1D("input_detected_delta_phi_deta2_norm","#Delta#phi Deta2;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta2_norm_detected, delta_phi_deta2_norm_input_detected);

  TH1D *delta_phi_deta2_norm_output_true;
  delta_phi_deta2_norm_output_true =  new TH1D("output_true_delta_phi_deta2_norm","#Delta#phi Deta2;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta2_norm_result = 0;

  unfold_distribution(delta_phi_deta2_norm_detected, delta_phi_deta2_norm_input_detected, delta_phi_deta2_norm_result, delta_phi_deta2_norm_output_true, delta_phi_deta2_norm_generated, delta_phi_deta2_response, delta_phi_deta2_response_all, plots_path, prefix, "delta_phi_deta2_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta3
TString delta_phi_deta3_hname = "resp_delta_phi_deta3";
TString delta_phi_deta3_hname_all = "resp_delta_phi_deta3_all";
TString delta_phi_deta3_measured = "ak5PF_delta_phi_deta3_fine";
TString delta_phi_deta3_measured_plot = "ak5PF_delta_phi_deta3";
TString delta_phi_deta3_true;
if (true_is_corr_data)
	{
	delta_phi_deta3_true = "ak5PF_delta_phi_deta3";
	}
else
	{
	delta_phi_deta3_true = "ak5Gen_delta_phi_deta3";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta3_response = 0;
respfile->GetObject(delta_phi_deta3_hname,delta_phi_deta3_response);
if (delta_phi_deta3_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta3_response_all = 0;
respfile->GetObject(delta_phi_deta3_hname_all,delta_phi_deta3_response_all);
if (delta_phi_deta3_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta3_detected = 0;
datafile->GetObject(delta_phi_deta3_measured,delta_phi_deta3_detected);
if (delta_phi_deta3_detected == 0) { cout<<"Measured data " << delta_phi_deta3_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta3_detected_plot = 0;
datafile->GetObject(delta_phi_deta3_measured_plot,delta_phi_deta3_detected_plot);
if (delta_phi_deta3_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta3_generated = 0;
mcfile->GetObject(delta_phi_deta3_true,delta_phi_deta3_generated);
if (delta_phi_deta3_generated == 0) { cout<< delta_phi_deta3_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta3_input_detected;
  delta_phi_deta3_input_detected =  new TH1D("input_detected_delta_phi_deta3","#Delta#phi Deta3;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta3_detected, delta_phi_deta3_input_detected);

  TH1D *delta_phi_deta3_output_true;
  delta_phi_deta3_output_true =  new TH1D("output_true_delta_phi_deta3","#Delta#phi Deta3;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta3_result = 0;

  unfold_distribution(delta_phi_deta3_detected, delta_phi_deta3_input_detected, delta_phi_deta3_result, delta_phi_deta3_output_true, delta_phi_deta3_generated, delta_phi_deta3_response, delta_phi_deta3_response_all, plots_path, prefix, "delta_phi_deta3", method, scale, tau, detail);


//Unfolding Delta Phi Deta3 Norm
TString delta_phi_deta3_norm_measured = "ak5PF_delta_phi_deta3_norm_fine";
TString delta_phi_deta3_norm_measured_plot = "ak5PF_delta_phi_deta3_norm";
TString delta_phi_deta3_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta3_norm_true = "ak5PF_delta_phi_deta3_norm";
	}
else
	{
	delta_phi_deta3_norm_true = "ak5Gen_delta_phi_deta3_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta3_norm_detected = 0;
datafile->GetObject(delta_phi_deta3_norm_measured,delta_phi_deta3_norm_detected);
if (delta_phi_deta3_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta3_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta3_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta3_norm_measured_plot,delta_phi_deta3_norm_detected_plot);
if (delta_phi_deta3_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta3_norm_generated = 0;
mcfile->GetObject(delta_phi_deta3_norm_true,delta_phi_deta3_norm_generated);
if (delta_phi_deta3_norm_generated == 0) { cout<< delta_phi_deta3_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta3_norm_input_detected;
  delta_phi_deta3_norm_input_detected =  new TH1D("input_detected_delta_phi_deta3_norm","#Delta#phi Deta3;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta3_norm_detected, delta_phi_deta3_norm_input_detected);

  TH1D *delta_phi_deta3_norm_output_true;
  delta_phi_deta3_norm_output_true =  new TH1D("output_true_delta_phi_deta3_norm","#Delta#phi Deta3;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta3_norm_result = 0;

  unfold_distribution(delta_phi_deta3_norm_detected, delta_phi_deta3_norm_input_detected, delta_phi_deta3_norm_result, delta_phi_deta3_norm_output_true, delta_phi_deta3_norm_generated, delta_phi_deta3_response, delta_phi_deta3_response_all, plots_path, prefix, "delta_phi_deta3_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta4
TString delta_phi_deta4_hname = "resp_delta_phi_deta4";
TString delta_phi_deta4_hname_all = "resp_delta_phi_deta4_all";
TString delta_phi_deta4_measured = "ak5PF_delta_phi_deta4_fine";
TString delta_phi_deta4_measured_plot = "ak5PF_delta_phi_deta4";
TString delta_phi_deta4_true;
if (true_is_corr_data)
	{
	delta_phi_deta4_true = "ak5PF_delta_phi_deta4";
	}
else
	{
	delta_phi_deta4_true = "ak5Gen_delta_phi_deta4";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta4_response = 0;
respfile->GetObject(delta_phi_deta4_hname,delta_phi_deta4_response);
if (delta_phi_deta4_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta4_response_all = 0;
respfile->GetObject(delta_phi_deta4_hname_all,delta_phi_deta4_response_all);
if (delta_phi_deta4_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta4_detected = 0;
datafile->GetObject(delta_phi_deta4_measured,delta_phi_deta4_detected);
if (delta_phi_deta4_detected == 0) { cout<<"Measured data " << delta_phi_deta4_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta4_detected_plot = 0;
datafile->GetObject(delta_phi_deta4_measured_plot,delta_phi_deta4_detected_plot);
if (delta_phi_deta4_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta4_generated = 0;
mcfile->GetObject(delta_phi_deta4_true,delta_phi_deta4_generated);
if (delta_phi_deta4_generated == 0) { cout<< delta_phi_deta4_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta4_input_detected;
  delta_phi_deta4_input_detected =  new TH1D("input_detected_delta_phi_deta4","#Delta#phi Deta4;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta4_detected, delta_phi_deta4_input_detected);

  TH1D *delta_phi_deta4_output_true;
  delta_phi_deta4_output_true =  new TH1D("output_true_delta_phi_deta4","#Delta#phi Deta4;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta4_result = 0;

  unfold_distribution(delta_phi_deta4_detected, delta_phi_deta4_input_detected, delta_phi_deta4_result, delta_phi_deta4_output_true, delta_phi_deta4_generated, delta_phi_deta4_response, delta_phi_deta4_response_all, plots_path, prefix, "delta_phi_deta4", method, scale, tau, detail);


//Unfolding Delta Phi Deta4 Norm
TString delta_phi_deta4_norm_measured = "ak5PF_delta_phi_deta4_norm_fine";
TString delta_phi_deta4_norm_measured_plot = "ak5PF_delta_phi_deta4_norm";
TString delta_phi_deta4_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta4_norm_true = "ak5PF_delta_phi_deta4_norm";
	}
else
	{
	delta_phi_deta4_norm_true = "ak5Gen_delta_phi_deta4_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta4_norm_detected = 0;
datafile->GetObject(delta_phi_deta4_norm_measured,delta_phi_deta4_norm_detected);
if (delta_phi_deta4_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta4_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta4_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta4_norm_measured_plot,delta_phi_deta4_norm_detected_plot);
if (delta_phi_deta4_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta4_norm_generated = 0;
mcfile->GetObject(delta_phi_deta4_norm_true,delta_phi_deta4_norm_generated);
if (delta_phi_deta4_norm_generated == 0) { cout<< delta_phi_deta4_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta4_norm_input_detected;
  delta_phi_deta4_norm_input_detected =  new TH1D("input_detected_delta_phi_deta4_norm","#Delta#phi Deta4;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta4_norm_detected, delta_phi_deta4_norm_input_detected);

  TH1D *delta_phi_deta4_norm_output_true;
  delta_phi_deta4_norm_output_true =  new TH1D("output_true_delta_phi_deta4_norm","#Delta#phi Deta4;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta4_norm_result = 0;

  unfold_distribution(delta_phi_deta4_norm_detected, delta_phi_deta4_norm_input_detected, delta_phi_deta4_norm_result, delta_phi_deta4_norm_output_true, delta_phi_deta4_norm_generated, delta_phi_deta4_response, delta_phi_deta4_response_all, plots_path, prefix, "delta_phi_deta4_norm", method, scale, tau, detail);


//Unfolding Delta Phi Gap
TString delta_phi_gap_hname = "resp_delta_phi_gap";
TString delta_phi_gap_hname_all = "resp_delta_phi_gap_all";
TString delta_phi_gap_measured = "ak5PF_delta_phi_gap_fine";
TString delta_phi_gap_measured_plot = "ak5PF_delta_phi_gap";
TString delta_phi_gap_true;
if (true_is_corr_data)
	{
	delta_phi_gap_true = "ak5PF_delta_phi_gap";
	}
else
	{
	delta_phi_gap_true = "ak5Gen_delta_phi_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_gap_response = 0;
respfile->GetObject(delta_phi_gap_hname,delta_phi_gap_response);
if (delta_phi_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_gap_response_all = 0;
respfile->GetObject(delta_phi_gap_hname_all,delta_phi_gap_response_all);
if (delta_phi_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_gap_detected = 0;
datafile->GetObject(delta_phi_gap_measured,delta_phi_gap_detected);
if (delta_phi_gap_detected == 0) { cout<<"Measured data " << delta_phi_gap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_gap_detected_plot = 0;
datafile->GetObject(delta_phi_gap_measured_plot,delta_phi_gap_detected_plot);
if (delta_phi_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_gap_generated = 0;
mcfile->GetObject(delta_phi_gap_true,delta_phi_gap_generated);
if (delta_phi_gap_generated == 0) { cout<< delta_phi_gap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_gap_input_detected;
  delta_phi_gap_input_detected =  new TH1D("input_detected_delta_phi_gap","#Delta#phi Gap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_gap_detected, delta_phi_gap_input_detected);

  TH1D *delta_phi_gap_output_true;
  delta_phi_gap_output_true =  new TH1D("output_true_delta_phi_gap","#Delta#phi Gap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_gap_result = 0;

  unfold_distribution(delta_phi_gap_detected, delta_phi_gap_input_detected, delta_phi_gap_result, delta_phi_gap_output_true, delta_phi_gap_generated, delta_phi_gap_response, delta_phi_gap_response_all, plots_path, prefix, "delta_phi_gap", method, scale, tau, detail);


//Unfolding Delta Phi Gap Norm
TString delta_phi_gap_norm_measured = "ak5PF_delta_phi_gap_norm_fine";
TString delta_phi_gap_norm_measured_plot = "ak5PF_delta_phi_gap_norm";
TString delta_phi_gap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_gap_norm_true = "ak5PF_delta_phi_gap_norm";
	}
else
	{
	delta_phi_gap_norm_true = "ak5Gen_delta_phi_gap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_gap_norm_detected = 0;
datafile->GetObject(delta_phi_gap_norm_measured,delta_phi_gap_norm_detected);
if (delta_phi_gap_norm_detected == 0) { cout<<"Measured data " << delta_phi_gap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_gap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_gap_norm_measured_plot,delta_phi_gap_norm_detected_plot);
if (delta_phi_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_gap_norm_generated = 0;
mcfile->GetObject(delta_phi_gap_norm_true,delta_phi_gap_norm_generated);
if (delta_phi_gap_norm_generated == 0) { cout<< delta_phi_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_gap_norm_input_detected;
  delta_phi_gap_norm_input_detected =  new TH1D("input_detected_delta_phi_gap_norm","#Delta#phi Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_gap_norm_detected, delta_phi_gap_norm_input_detected);

  TH1D *delta_phi_gap_norm_output_true;
  delta_phi_gap_norm_output_true =  new TH1D("output_true_delta_phi_gap_norm","#Delta#phi Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_gap_norm_result = 0;

  unfold_distribution(delta_phi_gap_norm_detected, delta_phi_gap_norm_input_detected, delta_phi_gap_norm_result, delta_phi_gap_norm_output_true, delta_phi_gap_norm_generated, delta_phi_gap_response, delta_phi_gap_response_all, plots_path, prefix, "delta_phi_gap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta1 Gap
TString delta_phi_deta1_gap_hname = "resp_delta_phi_deta1_gap";
TString delta_phi_deta1_gap_hname_all = "resp_delta_phi_deta1_gap_all";
TString delta_phi_deta1_gap_measured = "ak5PF_delta_phi_deta1_gap_fine";
TString delta_phi_deta1_gap_measured_plot = "ak5PF_delta_phi_deta1_gap";
TString delta_phi_deta1_gap_true;
if (true_is_corr_data)
	{
	delta_phi_deta1_gap_true = "ak5PF_delta_phi_deta1_gap";
	}
else
	{
	delta_phi_deta1_gap_true = "ak5Gen_delta_phi_deta1_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta1_gap_response = 0;
respfile->GetObject(delta_phi_deta1_gap_hname,delta_phi_deta1_gap_response);
if (delta_phi_deta1_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta1_gap_response_all = 0;
respfile->GetObject(delta_phi_deta1_gap_hname_all,delta_phi_deta1_gap_response_all);
if (delta_phi_deta1_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta1_gap_detected = 0;
datafile->GetObject(delta_phi_deta1_gap_measured,delta_phi_deta1_gap_detected);
if (delta_phi_deta1_gap_detected == 0) { cout<<"Measured data " << delta_phi_deta1_gap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta1_gap_detected_plot = 0;
datafile->GetObject(delta_phi_deta1_gap_measured_plot,delta_phi_deta1_gap_detected_plot);
if (delta_phi_deta1_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta1_gap_generated = 0;
mcfile->GetObject(delta_phi_deta1_gap_true,delta_phi_deta1_gap_generated);
if (delta_phi_deta1_gap_generated == 0) { cout<< delta_phi_deta1_gap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta1_gap_input_detected;
  delta_phi_deta1_gap_input_detected =  new TH1D("input_detected_delta_phi_deta1_gap","#Delta#phi Deta1 Gap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta1_gap_detected, delta_phi_deta1_gap_input_detected);

  TH1D *delta_phi_deta1_gap_output_true;
  delta_phi_deta1_gap_output_true =  new TH1D("output_true_delta_phi_deta1_gap","#Delta#phi Deta1 Gap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta1_gap_result = 0;

  unfold_distribution(delta_phi_deta1_gap_detected, delta_phi_deta1_gap_input_detected, delta_phi_deta1_gap_result, delta_phi_deta1_gap_output_true, delta_phi_deta1_gap_generated, delta_phi_deta1_gap_response, delta_phi_deta1_gap_response_all, plots_path, prefix, "delta_phi_deta1_gap", method, scale, tau, detail);


//Unfolding Delta Phi Deta1 Gap Norm
TString delta_phi_deta1_gap_norm_measured = "ak5PF_delta_phi_deta1_gap_norm_fine";
TString delta_phi_deta1_gap_norm_measured_plot = "ak5PF_delta_phi_deta1_gap_norm";
TString delta_phi_deta1_gap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta1_gap_norm_true = "ak5PF_delta_phi_deta1_gap_norm";
	}
else
	{
	delta_phi_deta1_gap_norm_true = "ak5Gen_delta_phi_deta1_gap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta1_gap_norm_detected = 0;
datafile->GetObject(delta_phi_deta1_gap_norm_measured,delta_phi_deta1_gap_norm_detected);
if (delta_phi_deta1_gap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta1_gap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta1_gap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta1_gap_norm_measured_plot,delta_phi_deta1_gap_norm_detected_plot);
if (delta_phi_deta1_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta1_gap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta1_gap_norm_true,delta_phi_deta1_gap_norm_generated);
if (delta_phi_deta1_gap_norm_generated == 0) { cout<< delta_phi_deta1_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta1_gap_norm_input_detected;
  delta_phi_deta1_gap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta1_gap_norm","#Delta#phi Deta1 Gap;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta1_gap_norm_detected, delta_phi_deta1_gap_norm_input_detected);

  TH1D *delta_phi_deta1_gap_norm_output_true;
  delta_phi_deta1_gap_norm_output_true =  new TH1D("output_true_delta_phi_deta1_gap_norm","#Delta#phi Deta1 Gap;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta1_gap_norm_result = 0;

  unfold_distribution(delta_phi_deta1_gap_norm_detected, delta_phi_deta1_gap_norm_input_detected, delta_phi_deta1_gap_norm_result, delta_phi_deta1_gap_norm_output_true, delta_phi_deta1_gap_norm_generated, delta_phi_deta1_gap_response, delta_phi_deta1_gap_response_all, plots_path, prefix, "delta_phi_deta1_gap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta2 Gap
TString delta_phi_deta2_gap_hname = "resp_delta_phi_deta2_gap";
TString delta_phi_deta2_gap_hname_all = "resp_delta_phi_deta2_gap_all";
TString delta_phi_deta2_gap_measured = "ak5PF_delta_phi_deta2_gap_fine";
TString delta_phi_deta2_gap_measured_plot = "ak5PF_delta_phi_deta2_gap";
TString delta_phi_deta2_gap_true;
if (true_is_corr_data)
	{
	delta_phi_deta2_gap_true = "ak5PF_delta_phi_deta2_gap";
	}
else
	{
	delta_phi_deta2_gap_true = "ak5Gen_delta_phi_deta2_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta2_gap_response = 0;
respfile->GetObject(delta_phi_deta2_gap_hname,delta_phi_deta2_gap_response);
if (delta_phi_deta2_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta2_gap_response_all = 0;
respfile->GetObject(delta_phi_deta2_gap_hname_all,delta_phi_deta2_gap_response_all);
if (delta_phi_deta2_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta2_gap_detected = 0;
datafile->GetObject(delta_phi_deta2_gap_measured,delta_phi_deta2_gap_detected);
if (delta_phi_deta2_gap_detected == 0) { cout<<"Measured data " << delta_phi_deta2_gap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta2_gap_detected_plot = 0;
datafile->GetObject(delta_phi_deta2_gap_measured_plot,delta_phi_deta2_gap_detected_plot);
if (delta_phi_deta2_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta2_gap_generated = 0;
mcfile->GetObject(delta_phi_deta2_gap_true,delta_phi_deta2_gap_generated);
if (delta_phi_deta2_gap_generated == 0) { cout<< delta_phi_deta2_gap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta2_gap_input_detected;
  delta_phi_deta2_gap_input_detected =  new TH1D("input_detected_delta_phi_deta2_gap","#Delta#phi Deta2 Gap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta2_gap_detected, delta_phi_deta2_gap_input_detected);

  TH1D *delta_phi_deta2_gap_output_true;
  delta_phi_deta2_gap_output_true =  new TH1D("output_true_delta_phi_deta2_gap","#Delta#phi Deta2 Gap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta2_gap_result = 0;

  unfold_distribution(delta_phi_deta2_gap_detected, delta_phi_deta2_gap_input_detected, delta_phi_deta2_gap_result, delta_phi_deta2_gap_output_true, delta_phi_deta2_gap_generated, delta_phi_deta2_gap_response, delta_phi_deta2_gap_response_all, plots_path, prefix, "delta_phi_deta2_gap", method, scale, tau, detail);


//Unfolding Delta Phi Deta2 Gap Norm
TString delta_phi_deta2_gap_norm_measured = "ak5PF_delta_phi_deta2_gap_norm_fine";
TString delta_phi_deta2_gap_norm_measured_plot = "ak5PF_delta_phi_deta2_gap_norm";
TString delta_phi_deta2_gap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta2_gap_norm_true = "ak5PF_delta_phi_deta2_gap_norm";
	}
else
	{
	delta_phi_deta2_gap_norm_true = "ak5Gen_delta_phi_deta2_gap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta2_gap_norm_detected = 0;
datafile->GetObject(delta_phi_deta2_gap_norm_measured,delta_phi_deta2_gap_norm_detected);
if (delta_phi_deta2_gap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta2_gap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta2_gap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta2_gap_norm_measured_plot,delta_phi_deta2_gap_norm_detected_plot);
if (delta_phi_deta2_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta2_gap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta2_gap_norm_true,delta_phi_deta2_gap_norm_generated);
if (delta_phi_deta2_gap_norm_generated == 0) { cout<< delta_phi_deta2_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta2_gap_norm_input_detected;
  delta_phi_deta2_gap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta2_gap_norm","#Delta#phi Deta2 Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta2_gap_norm_detected, delta_phi_deta2_gap_norm_input_detected);

  TH1D *delta_phi_deta2_gap_norm_output_true;
  delta_phi_deta2_gap_norm_output_true =  new TH1D("output_true_delta_phi_deta2_gap_norm","#Delta#phi Deta2 Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta2_gap_norm_result = 0;

  unfold_distribution(delta_phi_deta2_gap_norm_detected, delta_phi_deta2_gap_norm_input_detected, delta_phi_deta2_gap_norm_result, delta_phi_deta2_gap_norm_output_true, delta_phi_deta2_gap_norm_generated, delta_phi_deta2_gap_response, delta_phi_deta2_gap_response_all, plots_path, prefix, "delta_phi_deta2_gap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta3 Gap
TString delta_phi_deta3_gap_hname = "resp_delta_phi_deta3_gap";
TString delta_phi_deta3_gap_hname_all = "resp_delta_phi_deta3_gap_all";
TString delta_phi_deta3_gap_measured = "ak5PF_delta_phi_deta3_gap_fine";
TString delta_phi_deta3_gap_measured_plot = "ak5PF_delta_phi_deta3_gap";
TString delta_phi_deta3_gap_true;
if (true_is_corr_data)
	{
	delta_phi_deta3_gap_true = "ak5PF_delta_phi_deta3_gap";
	}
else
	{
	delta_phi_deta3_gap_true = "ak5Gen_delta_phi_deta3_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta3_gap_response = 0;
respfile->GetObject(delta_phi_deta3_gap_hname,delta_phi_deta3_gap_response);
if (delta_phi_deta3_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta3_gap_response_all = 0;
respfile->GetObject(delta_phi_deta3_gap_hname_all,delta_phi_deta3_gap_response_all);
if (delta_phi_deta3_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta3_gap_detected = 0;
datafile->GetObject(delta_phi_deta3_gap_measured,delta_phi_deta3_gap_detected);
if (delta_phi_deta3_gap_detected == 0) { cout<<"Measured data " << delta_phi_deta3_gap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta3_gap_detected_plot = 0;
datafile->GetObject(delta_phi_deta3_gap_measured_plot,delta_phi_deta3_gap_detected_plot);
if (delta_phi_deta3_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta3_gap_generated = 0;
mcfile->GetObject(delta_phi_deta3_gap_true,delta_phi_deta3_gap_generated);
if (delta_phi_deta3_gap_generated == 0) { cout<< delta_phi_deta3_gap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta3_gap_input_detected;
  delta_phi_deta3_gap_input_detected =  new TH1D("input_detected_delta_phi_deta3_gap","#Delta#phi Deta3 Gap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta3_gap_detected, delta_phi_deta3_gap_input_detected);

  TH1D *delta_phi_deta3_gap_output_true;
  delta_phi_deta3_gap_output_true =  new TH1D("output_true_delta_phi_deta3_gap","#Delta#phi Deta3 Gap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta3_gap_result = 0;

  unfold_distribution(delta_phi_deta3_gap_detected, delta_phi_deta3_gap_input_detected, delta_phi_deta3_gap_result, delta_phi_deta3_gap_output_true, delta_phi_deta3_gap_generated, delta_phi_deta3_gap_response, delta_phi_deta3_gap_response_all, plots_path, prefix, "delta_phi_deta3_gap", method, scale, tau, detail);


//Unfolding Delta Phi Deta3 Gap Norm
TString delta_phi_deta3_gap_norm_measured = "ak5PF_delta_phi_deta3_gap_norm_fine";
TString delta_phi_deta3_gap_norm_measured_plot = "ak5PF_delta_phi_deta3_gap_norm";
TString delta_phi_deta3_gap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta3_gap_norm_true = "ak5PF_delta_phi_deta3_gap_norm";
	}
else
	{
	delta_phi_deta3_gap_norm_true = "ak5Gen_delta_phi_deta3_gap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta3_gap_norm_detected = 0;
datafile->GetObject(delta_phi_deta3_gap_norm_measured,delta_phi_deta3_gap_norm_detected);
if (delta_phi_deta3_gap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta3_gap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta3_gap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta3_gap_norm_measured_plot,delta_phi_deta3_gap_norm_detected_plot);
if (delta_phi_deta3_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta3_gap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta3_gap_norm_true,delta_phi_deta3_gap_norm_generated);
if (delta_phi_deta3_gap_norm_generated == 0) { cout<< delta_phi_deta3_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta3_gap_norm_input_detected;
  delta_phi_deta3_gap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta3_gap_norm","#Delta#phi Deta3 Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta3_gap_norm_detected, delta_phi_deta3_gap_norm_input_detected);

  TH1D *delta_phi_deta3_gap_norm_output_true;
  delta_phi_deta3_gap_norm_output_true =  new TH1D("output_true_delta_phi_deta3_gap_norm","#Delta#phi Deta3 Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta3_gap_norm_result = 0;

  unfold_distribution(delta_phi_deta3_gap_norm_detected, delta_phi_deta3_gap_norm_input_detected, delta_phi_deta3_gap_norm_result, delta_phi_deta3_gap_norm_output_true, delta_phi_deta3_gap_norm_generated, delta_phi_deta3_gap_response, delta_phi_deta3_gap_response_all, plots_path, prefix, "delta_phi_deta3_gap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta4 Gap
TString delta_phi_deta4_gap_hname = "resp_delta_phi_deta4_gap";
TString delta_phi_deta4_gap_hname_all = "resp_delta_phi_deta4_gap_all";
TString delta_phi_deta4_gap_measured = "ak5PF_delta_phi_deta4_gap_fine";
TString delta_phi_deta4_gap_measured_plot = "ak5PF_delta_phi_deta4_gap";
TString delta_phi_deta4_gap_true;
if (true_is_corr_data)
	{
	delta_phi_deta4_gap_true = "ak5PF_delta_phi_deta4_gap";
	}
else
	{
	delta_phi_deta4_gap_true = "ak5Gen_delta_phi_deta4_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta4_gap_response = 0;
respfile->GetObject(delta_phi_deta4_gap_hname,delta_phi_deta4_gap_response);
if (delta_phi_deta4_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta4_gap_response_all = 0;
respfile->GetObject(delta_phi_deta4_gap_hname_all,delta_phi_deta4_gap_response_all);
if (delta_phi_deta4_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta4_gap_detected = 0;
datafile->GetObject(delta_phi_deta4_gap_measured,delta_phi_deta4_gap_detected);
if (delta_phi_deta4_gap_detected == 0) { cout<<"Measured data " << delta_phi_deta4_gap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta4_gap_detected_plot = 0;
datafile->GetObject(delta_phi_deta4_gap_measured_plot,delta_phi_deta4_gap_detected_plot);
if (delta_phi_deta4_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta4_gap_generated = 0;
mcfile->GetObject(delta_phi_deta4_gap_true,delta_phi_deta4_gap_generated);
if (delta_phi_deta4_gap_generated == 0) { cout<< delta_phi_deta4_gap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta4_gap_input_detected;
  delta_phi_deta4_gap_input_detected =  new TH1D("input_detected_delta_phi_deta4_gap","#Delta#phi Deta4 Gap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta4_gap_detected, delta_phi_deta4_gap_input_detected);

  TH1D *delta_phi_deta4_gap_output_true;
  delta_phi_deta4_gap_output_true =  new TH1D("output_true_delta_phi_deta4_gap","#Delta#phi Deta4 Gap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta4_gap_result = 0;

  unfold_distribution(delta_phi_deta4_gap_detected, delta_phi_deta4_gap_input_detected, delta_phi_deta4_gap_result, delta_phi_deta4_gap_output_true, delta_phi_deta4_gap_generated, delta_phi_deta4_gap_response, delta_phi_deta4_gap_response_all, plots_path, prefix, "delta_phi_deta4_gap", method, scale, tau, detail);


//Unfolding Delta Phi Deta4 Gap Norm
TString delta_phi_deta4_gap_norm_measured = "ak5PF_delta_phi_deta4_gap_norm_fine";
TString delta_phi_deta4_gap_norm_measured_plot = "ak5PF_delta_phi_deta4_gap_norm";
TString delta_phi_deta4_gap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta4_gap_norm_true = "ak5PF_delta_phi_deta4_gap_norm";
	}
else
	{
	delta_phi_deta4_gap_norm_true = "ak5Gen_delta_phi_deta4_gap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta4_gap_norm_detected = 0;
datafile->GetObject(delta_phi_deta4_gap_norm_measured,delta_phi_deta4_gap_norm_detected);
if (delta_phi_deta4_gap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta4_gap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta4_gap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta4_gap_norm_measured_plot,delta_phi_deta4_gap_norm_detected_plot);
if (delta_phi_deta4_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta4_gap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta4_gap_norm_true,delta_phi_deta4_gap_norm_generated);
if (delta_phi_deta4_gap_norm_generated == 0) { cout<< delta_phi_deta4_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta4_gap_norm_input_detected;
  delta_phi_deta4_gap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta4_gap_norm","#Delta#phi Deta4 Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta4_gap_norm_detected, delta_phi_deta4_gap_norm_input_detected);

  TH1D *delta_phi_deta4_gap_norm_output_true;
  delta_phi_deta4_gap_norm_output_true =  new TH1D("output_true_delta_phi_deta4_gap_norm","#Delta#phi Deta4 Gap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta4_gap_norm_result = 0;

  unfold_distribution(delta_phi_deta4_gap_norm_detected, delta_phi_deta4_gap_norm_input_detected, delta_phi_deta4_gap_norm_result, delta_phi_deta4_gap_norm_output_true, delta_phi_deta4_gap_norm_generated, delta_phi_deta4_gap_response, delta_phi_deta4_gap_response_all, plots_path, prefix, "delta_phi_deta4_gap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Nogap
TString delta_phi_nogap_hname = "resp_delta_phi_nogap";
TString delta_phi_nogap_hname_all = "resp_delta_phi_nogap_all";
TString delta_phi_nogap_measured = "ak5PF_delta_phi_nogap_fine";
TString delta_phi_nogap_measured_plot = "ak5PF_delta_phi_nogap";
TString delta_phi_nogap_true;
if (true_is_corr_data)
	{
	delta_phi_nogap_true = "ak5PF_delta_phi_nogap";
	}
else
	{
	delta_phi_nogap_true = "ak5Gen_delta_phi_nogap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_nogap_response = 0;
respfile->GetObject(delta_phi_nogap_hname,delta_phi_nogap_response);
if (delta_phi_nogap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_nogap_response_all = 0;
respfile->GetObject(delta_phi_nogap_hname_all,delta_phi_nogap_response_all);
if (delta_phi_nogap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_nogap_detected = 0;
datafile->GetObject(delta_phi_nogap_measured,delta_phi_nogap_detected);
if (delta_phi_nogap_detected == 0) { cout<<"Measured data " << delta_phi_nogap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_nogap_detected_plot = 0;
datafile->GetObject(delta_phi_nogap_measured_plot,delta_phi_nogap_detected_plot);
if (delta_phi_nogap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_nogap_generated = 0;
mcfile->GetObject(delta_phi_nogap_true,delta_phi_nogap_generated);
if (delta_phi_nogap_generated == 0) { cout<< delta_phi_nogap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_nogap_input_detected;
  delta_phi_nogap_input_detected =  new TH1D("input_detected_delta_phi_nogap","#Delta#phi Nogap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_nogap_detected, delta_phi_nogap_input_detected);

  TH1D *delta_phi_nogap_output_true;
  delta_phi_nogap_output_true =  new TH1D("output_true_delta_phi_nogap","#Delta#phi Nogap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_nogap_result = 0;

  unfold_distribution(delta_phi_nogap_detected, delta_phi_nogap_input_detected, delta_phi_nogap_result, delta_phi_nogap_output_true, delta_phi_nogap_generated, delta_phi_nogap_response, delta_phi_nogap_response_all, plots_path, prefix, "delta_phi_nogap", method, scale, tau, detail);


//Unfolding Delta Phi Nogap Norm
TString delta_phi_nogap_norm_measured = "ak5PF_delta_phi_nogap_norm_fine";
TString delta_phi_nogap_norm_measured_plot = "ak5PF_delta_phi_nogap_norm";
TString delta_phi_nogap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_nogap_norm_true = "ak5PF_delta_phi_nogap_norm";
	}
else
	{
	delta_phi_nogap_norm_true = "ak5Gen_delta_phi_nogap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_nogap_norm_detected = 0;
datafile->GetObject(delta_phi_nogap_norm_measured,delta_phi_nogap_norm_detected);
if (delta_phi_nogap_norm_detected == 0) { cout<<"Measured data " << delta_phi_nogap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_nogap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_nogap_norm_measured_plot,delta_phi_nogap_norm_detected_plot);
if (delta_phi_nogap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_nogap_norm_generated = 0;
mcfile->GetObject(delta_phi_nogap_norm_true,delta_phi_nogap_norm_generated);
if (delta_phi_nogap_norm_generated == 0) { cout<< delta_phi_nogap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_nogap_norm_input_detected;
  delta_phi_nogap_norm_input_detected =  new TH1D("input_detected_delta_phi_nogap_norm","#Delta#phi Nogap;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_nogap_norm_detected, delta_phi_nogap_norm_input_detected);

  TH1D *delta_phi_nogap_norm_output_true;
  delta_phi_nogap_norm_output_true =  new TH1D("output_true_delta_phi_nogap_norm","#Delta#phi Nogap;#Delta#phi [rad];#Delta#phi [rad];#frac{#sigma^{-1} d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_nogap_norm_result = 0;

  unfold_distribution(delta_phi_nogap_norm_detected, delta_phi_nogap_norm_input_detected, delta_phi_nogap_norm_result, delta_phi_nogap_norm_output_true, delta_phi_nogap_norm_generated, delta_phi_nogap_response, delta_phi_nogap_response_all, plots_path, prefix, "delta_phi_nogap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta1 Nogap
TString delta_phi_deta1_nogap_hname = "resp_delta_phi_deta1_nogap";
TString delta_phi_deta1_nogap_hname_all = "resp_delta_phi_deta1_nogap_all";
TString delta_phi_deta1_nogap_measured = "ak5PF_delta_phi_deta1_nogap_fine";
TString delta_phi_deta1_nogap_measured_plot = "ak5PF_delta_phi_deta1_nogap";
TString delta_phi_deta1_nogap_true;
if (true_is_corr_data)
	{
	delta_phi_deta1_nogap_true = "ak5PF_delta_phi_deta1_nogap";
	}
else
	{
	delta_phi_deta1_nogap_true = "ak5Gen_delta_phi_deta1_nogap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta1_nogap_response = 0;
respfile->GetObject(delta_phi_deta1_nogap_hname,delta_phi_deta1_nogap_response);
if (delta_phi_deta1_nogap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta1_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta1_nogap_hname_all,delta_phi_deta1_nogap_response_all);
if (delta_phi_deta1_nogap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta1_nogap_detected = 0;
datafile->GetObject(delta_phi_deta1_nogap_measured,delta_phi_deta1_nogap_detected);
if (delta_phi_deta1_nogap_detected == 0) { cout<<"Measured data " << delta_phi_deta1_nogap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta1_nogap_detected_plot = 0;
datafile->GetObject(delta_phi_deta1_nogap_measured_plot,delta_phi_deta1_nogap_detected_plot);
if (delta_phi_deta1_nogap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta1_nogap_generated = 0;
mcfile->GetObject(delta_phi_deta1_nogap_true,delta_phi_deta1_nogap_generated);
if (delta_phi_deta1_nogap_generated == 0) { cout<< delta_phi_deta1_nogap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta1_nogap_input_detected;
  delta_phi_deta1_nogap_input_detected =  new TH1D("input_detected_delta_phi_deta1_nogap","#Delta#phi Deta1 Nogap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta1_nogap_detected, delta_phi_deta1_nogap_input_detected);

  TH1D *delta_phi_deta1_nogap_output_true;
  delta_phi_deta1_nogap_output_true =  new TH1D("output_true_delta_phi_deta1_nogap","#Delta#phi Deta1 Nogap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta1_nogap_result = 0;

  unfold_distribution(delta_phi_deta1_nogap_detected, delta_phi_deta1_nogap_input_detected, delta_phi_deta1_nogap_result, delta_phi_deta1_nogap_output_true, delta_phi_deta1_nogap_generated, delta_phi_deta1_nogap_response, delta_phi_deta1_nogap_response_all, plots_path, prefix, "delta_phi_deta1_nogap", method, scale, tau, detail);


//Unfolding Delta Phi Deta1 Nogap Norm
TString delta_phi_deta1_nogap_norm_measured = "ak5PF_delta_phi_deta1_nogap_norm_fine";
TString delta_phi_deta1_nogap_norm_measured_plot = "ak5PF_delta_phi_deta1_nogap_norm";
TString delta_phi_deta1_nogap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta1_nogap_norm_true = "ak5PF_delta_phi_deta1_nogap_norm";
	}
else
	{
	delta_phi_deta1_nogap_norm_true = "ak5Gen_delta_phi_deta1_nogap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta1_nogap_norm_detected = 0;
datafile->GetObject(delta_phi_deta1_nogap_norm_measured,delta_phi_deta1_nogap_norm_detected);
if (delta_phi_deta1_nogap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta1_nogap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta1_nogap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta1_nogap_norm_measured_plot,delta_phi_deta1_nogap_norm_detected_plot);
if (delta_phi_deta1_nogap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta1_nogap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta1_nogap_norm_true,delta_phi_deta1_nogap_norm_generated);
if (delta_phi_deta1_nogap_norm_generated == 0) { cout<< delta_phi_deta1_nogap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta1_nogap_norm_input_detected;
  delta_phi_deta1_nogap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta1_nogap_norm","#Delta#phi Deta1 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta1_nogap_norm_detected, delta_phi_deta1_nogap_norm_input_detected);

  TH1D *delta_phi_deta1_nogap_norm_output_true;
  delta_phi_deta1_nogap_norm_output_true =  new TH1D("output_true_delta_phi_deta1_nogap_norm","#Delta#phi Deta1 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta1_nogap_norm_result = 0;

  unfold_distribution(delta_phi_deta1_nogap_norm_detected, delta_phi_deta1_nogap_norm_input_detected, delta_phi_deta1_nogap_norm_result, delta_phi_deta1_nogap_norm_output_true, delta_phi_deta1_nogap_norm_generated, delta_phi_deta1_nogap_response, delta_phi_deta1_nogap_response_all, plots_path, prefix, "delta_phi_deta1_nogap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta2 Nogap
TString delta_phi_deta2_nogap_hname = "resp_delta_phi_deta2_nogap";
TString delta_phi_deta2_nogap_hname_all = "resp_delta_phi_deta2_nogap_all";
TString delta_phi_deta2_nogap_measured = "ak5PF_delta_phi_deta2_nogap_fine";
TString delta_phi_deta2_nogap_measured_plot = "ak5PF_delta_phi_deta2_nogap";
TString delta_phi_deta2_nogap_true;
if (true_is_corr_data)
	{
	delta_phi_deta2_nogap_true = "ak5PF_delta_phi_deta2_nogap";
	}
else
	{
	delta_phi_deta2_nogap_true = "ak5Gen_delta_phi_deta2_nogap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta2_nogap_response = 0;
respfile->GetObject(delta_phi_deta2_nogap_hname,delta_phi_deta2_nogap_response);
if (delta_phi_deta2_nogap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta2_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta2_nogap_hname_all,delta_phi_deta2_nogap_response_all);
if (delta_phi_deta2_nogap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta2_nogap_detected = 0;
datafile->GetObject(delta_phi_deta2_nogap_measured,delta_phi_deta2_nogap_detected);
if (delta_phi_deta2_nogap_detected == 0) { cout<<"Measured data " << delta_phi_deta2_nogap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta2_nogap_detected_plot = 0;
datafile->GetObject(delta_phi_deta2_nogap_measured_plot,delta_phi_deta2_nogap_detected_plot);
if (delta_phi_deta2_nogap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta2_nogap_generated = 0;
mcfile->GetObject(delta_phi_deta2_nogap_true,delta_phi_deta2_nogap_generated);
if (delta_phi_deta2_nogap_generated == 0) { cout<< delta_phi_deta2_nogap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta2_nogap_input_detected;
  delta_phi_deta2_nogap_input_detected =  new TH1D("input_detected_delta_phi_deta2_nogap","#Delta#phi Deta2 Nogap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta2_nogap_detected, delta_phi_deta2_nogap_input_detected);

  TH1D *delta_phi_deta2_nogap_output_true;
  delta_phi_deta2_nogap_output_true =  new TH1D("output_true_delta_phi_deta2_nogap","#Delta#phi Deta2 Nogap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta2_nogap_result = 0;

  unfold_distribution(delta_phi_deta2_nogap_detected, delta_phi_deta2_nogap_input_detected, delta_phi_deta2_nogap_result, delta_phi_deta2_nogap_output_true, delta_phi_deta2_nogap_generated, delta_phi_deta2_nogap_response, delta_phi_deta2_nogap_response_all, plots_path, prefix, "delta_phi_deta2_nogap", method, scale, tau, detail);


//Unfolding Delta Phi Deta2 Nogap Norm
TString delta_phi_deta2_nogap_norm_measured = "ak5PF_delta_phi_deta2_nogap_norm_fine";
TString delta_phi_deta2_nogap_norm_measured_plot = "ak5PF_delta_phi_deta2_nogap_norm";
TString delta_phi_deta2_nogap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta2_nogap_norm_true = "ak5PF_delta_phi_deta2_nogap_norm";
	}
else
	{
	delta_phi_deta2_nogap_norm_true = "ak5Gen_delta_phi_deta2_nogap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta2_nogap_norm_detected = 0;
datafile->GetObject(delta_phi_deta2_nogap_norm_measured,delta_phi_deta2_nogap_norm_detected);
if (delta_phi_deta2_nogap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta2_nogap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta2_nogap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta2_nogap_norm_measured_plot,delta_phi_deta2_nogap_norm_detected_plot);
if (delta_phi_deta2_nogap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta2_nogap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta2_nogap_norm_true,delta_phi_deta2_nogap_norm_generated);
if (delta_phi_deta2_nogap_norm_generated == 0) { cout<< delta_phi_deta2_nogap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta2_nogap_norm_input_detected;
  delta_phi_deta2_nogap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta2_nogap_norm","#Delta#phi Deta2 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta2_nogap_norm_detected, delta_phi_deta2_nogap_norm_input_detected);

  TH1D *delta_phi_deta2_nogap_norm_output_true;
  delta_phi_deta2_nogap_norm_output_true =  new TH1D("output_true_delta_phi_deta2_nogap_norm","#Delta#phi Deta2 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta2_nogap_norm_result = 0;

  unfold_distribution(delta_phi_deta2_nogap_norm_detected, delta_phi_deta2_nogap_norm_input_detected, delta_phi_deta2_nogap_norm_result, delta_phi_deta2_nogap_norm_output_true, delta_phi_deta2_nogap_norm_generated, delta_phi_deta2_nogap_response, delta_phi_deta2_nogap_response_all, plots_path, prefix, "delta_phi_deta2_nogap_norm", method, scale, tau, detail);


//Unfolding Delta Phi Deta3 Nogap
TString delta_phi_deta3_nogap_hname = "resp_delta_phi_deta3_nogap";
TString delta_phi_deta3_nogap_hname_all = "resp_delta_phi_deta3_nogap_all";
TString delta_phi_deta3_nogap_measured = "ak5PF_delta_phi_deta3_nogap_fine";
TString delta_phi_deta3_nogap_measured_plot = "ak5PF_delta_phi_deta3_nogap";
TString delta_phi_deta3_nogap_true;
if (true_is_corr_data)
	{
	delta_phi_deta3_nogap_true = "ak5PF_delta_phi_deta3_nogap";
	}
else
	{
	delta_phi_deta3_nogap_true = "ak5Gen_delta_phi_deta3_nogap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta3_nogap_response = 0;
respfile->GetObject(delta_phi_deta3_nogap_hname,delta_phi_deta3_nogap_response);
if (delta_phi_deta3_nogap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta3_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta3_nogap_hname_all,delta_phi_deta3_nogap_response_all);
if (delta_phi_deta3_nogap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta3_nogap_detected = 0;
datafile->GetObject(delta_phi_deta3_nogap_measured,delta_phi_deta3_nogap_detected);
if (delta_phi_deta3_nogap_detected == 0) { cout<<"Measured data " << delta_phi_deta3_nogap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta3_nogap_detected_plot = 0;
datafile->GetObject(delta_phi_deta3_nogap_measured_plot,delta_phi_deta3_nogap_detected_plot);
if (delta_phi_deta3_nogap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta3_nogap_generated = 0;
mcfile->GetObject(delta_phi_deta3_nogap_true,delta_phi_deta3_nogap_generated);
if (delta_phi_deta3_nogap_generated == 0) { cout<< delta_phi_deta3_nogap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta3_nogap_input_detected;
  delta_phi_deta3_nogap_input_detected =  new TH1D("input_detected_delta_phi_deta3_nogap","#Delta#phi Deta3 Nogap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta3_nogap_detected, delta_phi_deta3_nogap_input_detected);

  TH1D *delta_phi_deta3_nogap_output_true;
  delta_phi_deta3_nogap_output_true =  new TH1D("output_true_delta_phi_deta3_nogap","#Delta#phi Deta3 Nogap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta3_nogap_result = 0;

  unfold_distribution(delta_phi_deta3_nogap_detected, delta_phi_deta3_nogap_input_detected, delta_phi_deta3_nogap_result, delta_phi_deta3_nogap_output_true, delta_phi_deta3_nogap_generated, delta_phi_deta3_nogap_response, delta_phi_deta3_nogap_response_all, plots_path, prefix, "delta_phi_deta3_nogap", method, scale, tau, detail);


//Unfolding Delta Phi Deta3 Nogap Norm
TString delta_phi_deta3_nogap_norm_measured = "ak5PF_delta_phi_deta3_nogap_norm_fine";
TString delta_phi_deta3_nogap_norm_measured_plot = "ak5PF_delta_phi_deta3_nogap_norm";
TString delta_phi_deta3_nogap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta3_nogap_norm_true = "ak5PF_delta_phi_deta3_nogap_norm";
	}
else
	{
	delta_phi_deta3_nogap_norm_true = "ak5Gen_delta_phi_deta3_nogap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta3_nogap_norm_detected = 0;
datafile->GetObject(delta_phi_deta3_nogap_norm_measured,delta_phi_deta3_nogap_norm_detected);
if (delta_phi_deta3_nogap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta3_nogap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta3_nogap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta3_nogap_norm_measured_plot,delta_phi_deta3_nogap_norm_detected_plot);
if (delta_phi_deta3_nogap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta3_nogap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta3_nogap_norm_true,delta_phi_deta3_nogap_norm_generated);
if (delta_phi_deta3_nogap_norm_generated == 0) { cout<< delta_phi_deta3_nogap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta3_nogap_norm_input_detected;
  delta_phi_deta3_nogap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta3_nogap_norm","#Delta#phi Deta3 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta3_nogap_norm_detected, delta_phi_deta3_nogap_norm_input_detected);

  TH1D *delta_phi_deta3_nogap_norm_output_true;
  delta_phi_deta3_nogap_norm_output_true =  new TH1D("output_true_delta_phi_deta3_nogap_norm","#Delta#phi Deta3 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta3_nogap_norm_result = 0;

  unfold_distribution(delta_phi_deta3_nogap_norm_detected, delta_phi_deta3_nogap_norm_input_detected, delta_phi_deta3_nogap_norm_result, delta_phi_deta3_nogap_norm_output_true, delta_phi_deta3_nogap_norm_generated, delta_phi_deta3_nogap_response, delta_phi_deta3_nogap_response_all, plots_path, prefix, "delta_phi_deta3_nogap_norm", method, scale, tau, detail);



//Unfolding Delta Phi Deta4 Nogap
TString delta_phi_deta4_nogap_hname = "resp_delta_phi_deta4_nogap";
TString delta_phi_deta4_nogap_hname_all = "resp_delta_phi_deta4_nogap_all";
TString delta_phi_deta4_nogap_measured = "ak5PF_delta_phi_deta4_nogap_fine";
TString delta_phi_deta4_nogap_measured_plot = "ak5PF_delta_phi_deta4_nogap";
TString delta_phi_deta4_nogap_true;
if (true_is_corr_data)
	{
	delta_phi_deta4_nogap_true = "ak5PF_delta_phi_deta4_nogap";
	}
else
	{
	delta_phi_deta4_nogap_true = "ak5Gen_delta_phi_deta4_nogap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_phi_deta4_nogap_response = 0;
respfile->GetObject(delta_phi_deta4_nogap_hname,delta_phi_deta4_nogap_response);
if (delta_phi_deta4_nogap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_deta4_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta4_nogap_hname_all,delta_phi_deta4_nogap_response_all);
if (delta_phi_deta4_nogap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_phi_deta4_nogap_detected = 0;
datafile->GetObject(delta_phi_deta4_nogap_measured,delta_phi_deta4_nogap_detected);
if (delta_phi_deta4_nogap_detected == 0) { cout<<"Measured data " << delta_phi_deta4_nogap_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta4_nogap_detected_plot = 0;
datafile->GetObject(delta_phi_deta4_nogap_measured_plot,delta_phi_deta4_nogap_detected_plot);
if (delta_phi_deta4_nogap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta4_nogap_generated = 0;
mcfile->GetObject(delta_phi_deta4_nogap_true,delta_phi_deta4_nogap_generated);
if (delta_phi_deta4_nogap_generated == 0) { cout<< delta_phi_deta4_nogap_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta4_nogap_input_detected;
  delta_phi_deta4_nogap_input_detected =  new TH1D("input_detected_delta_phi_deta4_nogap","#Delta#phi Deta4 Nogap;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta4_nogap_detected, delta_phi_deta4_nogap_input_detected);

  TH1D *delta_phi_deta4_nogap_output_true;
  delta_phi_deta4_nogap_output_true =  new TH1D("output_true_delta_phi_deta4_nogap","#Delta#phi Deta4 Nogap;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta4_nogap_result = 0;

  unfold_distribution(delta_phi_deta4_nogap_detected, delta_phi_deta4_nogap_input_detected, delta_phi_deta4_nogap_result, delta_phi_deta4_nogap_output_true, delta_phi_deta4_nogap_generated, delta_phi_deta4_nogap_response, delta_phi_deta4_nogap_response_all, plots_path, prefix, "delta_phi_deta4_nogap", method, scale, tau, detail);


//Unfolding Delta Phi Deta4 Nogap Norm
TString delta_phi_deta4_nogap_norm_measured = "ak5PF_delta_phi_deta4_nogap_norm_fine";
TString delta_phi_deta4_nogap_norm_measured_plot = "ak5PF_delta_phi_deta4_nogap_norm";
TString delta_phi_deta4_nogap_norm_true;
if (true_is_corr_data)
	{
	delta_phi_deta4_nogap_norm_true = "ak5PF_delta_phi_deta4_nogap_norm";
	}
else
	{
	delta_phi_deta4_nogap_norm_true = "ak5Gen_delta_phi_deta4_nogap_norm";
	}

//get the objects on a safe way
TH1D *delta_phi_deta4_nogap_norm_detected = 0;
datafile->GetObject(delta_phi_deta4_nogap_norm_measured,delta_phi_deta4_nogap_norm_detected);
if (delta_phi_deta4_nogap_norm_detected == 0) { cout<<"Measured data " << delta_phi_deta4_nogap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_phi_deta4_nogap_norm_detected_plot = 0;
datafile->GetObject(delta_phi_deta4_nogap_norm_measured_plot,delta_phi_deta4_nogap_norm_detected_plot);
if (delta_phi_deta4_nogap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_phi_deta4_nogap_norm_generated = 0;
mcfile->GetObject(delta_phi_deta4_nogap_norm_true,delta_phi_deta4_nogap_norm_generated);
if (delta_phi_deta4_nogap_norm_generated == 0) { cout<< delta_phi_deta4_nogap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_phi_deta4_nogap_norm_input_detected;
  delta_phi_deta4_nogap_norm_input_detected =  new TH1D("input_detected_delta_phi_deta4_nogap_norm","#Delta#phi Deta4 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins_meas, dphi_bins_meas);
  set_input(delta_phi_deta4_nogap_norm_detected, delta_phi_deta4_nogap_norm_input_detected);

  TH1D *delta_phi_deta4_nogap_norm_output_true;
  delta_phi_deta4_nogap_norm_output_true =  new TH1D("output_true_delta_phi_deta4_nogap_norm","#Delta#phi Deta4 Nogap;#Delta#phi [rad];#frac{#sigma^{-1} d#sigma^{2}}{d#Delta#eta d#Delta#phi}", dphi_nbins, dphi_bins);

  TH1D *delta_phi_deta4_nogap_norm_result = 0;

  unfold_distribution(delta_phi_deta4_nogap_norm_detected, delta_phi_deta4_nogap_norm_input_detected, delta_phi_deta4_nogap_norm_result, delta_phi_deta4_nogap_norm_output_true, delta_phi_deta4_nogap_norm_generated, delta_phi_deta4_nogap_response, delta_phi_deta4_nogap_response_all, plots_path, prefix, "delta_phi_deta4_nogap_norm", method, scale, tau, detail);



//Unfolding Leading pT Inside Gap
TString leading_pt_inside_gap_hname = "resp_leading_pt_inside_gap";
TString leading_pt_inside_gap_hname_all = "resp_leading_pt_inside_gap_all";
TString leading_pt_inside_gap_measured = "ak5PF_leading_pt_inside_gap_fine";
TString leading_pt_inside_gap_measured_plot = "ak5PF_leading_pt_inside_gap";
TString leading_pt_inside_gap_true;
if (true_is_corr_data)
	{
	leading_pt_inside_gap_true = "ak5PF_leading_pt_inside_gap";
	}
else
	{
	leading_pt_inside_gap_true = "ak5Gen_leading_pt_inside_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *leading_pt_inside_gap_response = 0;
respfile->GetObject(leading_pt_inside_gap_hname,leading_pt_inside_gap_response);
if (leading_pt_inside_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *leading_pt_inside_gap_response_all = 0;
respfile->GetObject(leading_pt_inside_gap_hname_all,leading_pt_inside_gap_response_all);
if (leading_pt_inside_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *leading_pt_inside_gap_detected = 0;
datafile->GetObject(leading_pt_inside_gap_measured,leading_pt_inside_gap_detected);
if (leading_pt_inside_gap_detected == 0) { cout<<"Measured data " << leading_pt_inside_gap_measured << " not found!"<<endl; return; }

TH1D *leading_pt_inside_gap_detected_plot = 0;
datafile->GetObject(leading_pt_inside_gap_measured_plot,leading_pt_inside_gap_detected_plot);
if (leading_pt_inside_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *leading_pt_inside_gap_generated = 0;
mcfile->GetObject(leading_pt_inside_gap_true,leading_pt_inside_gap_generated);
if (leading_pt_inside_gap_generated == 0) { cout<< leading_pt_inside_gap_true << " true data not found!"<<endl; return; }

  TH1D *leading_pt_inside_gap_input_detected;
  leading_pt_inside_gap_input_detected =  new TH1D("input_detected_leading_pt_inside_gap","Leading pT Inside Gap;p_{T}^{inside};Events", in_nbins_meas, in_bins_meas);
  set_input(leading_pt_inside_gap_detected, leading_pt_inside_gap_input_detected);

  TH1D *leading_pt_inside_gap_output_true;
  leading_pt_inside_gap_output_true =  new TH1D("output_true_leading_pt_inside_gap","Leading pT Inside Gap;p_{T}^{inside} [GeV];Events", in_nbins, in_bins);

  TH1D *leading_pt_inside_gap_result = 0;

  unfold_distribution(leading_pt_inside_gap_detected, leading_pt_inside_gap_input_detected, leading_pt_inside_gap_result, leading_pt_inside_gap_output_true, leading_pt_inside_gap_generated, leading_pt_inside_gap_response, leading_pt_inside_gap_response_all, plots_path, prefix, "leading_pt_inside_gap", method, scale, tau, detail);


//Unfolding Leading pT Inside Gap Norm
TString leading_pt_inside_gap_norm_measured = "ak5PF_leading_pt_inside_gap_norm_fine";
TString leading_pt_inside_gap_norm_measured_plot = "ak5PF_leading_pt_inside_gap_norm";
TString leading_pt_inside_gap_norm_true;
if (true_is_corr_data)
	{
	leading_pt_inside_gap_norm_true = "ak5PF_leading_pt_inside_gap_norm";
	}
else
	{
	leading_pt_inside_gap_norm_true = "ak5Gen_leading_pt_inside_gap_norm";
	}

//get the objects on a safe way
TH1D *leading_pt_inside_gap_norm_detected = 0;
datafile->GetObject(leading_pt_inside_gap_norm_measured,leading_pt_inside_gap_norm_detected);
if (leading_pt_inside_gap_norm_detected == 0) { cout<<"Measured data " << leading_pt_inside_gap_norm_measured << " not found!"<<endl; return; }

TH1D *leading_pt_inside_gap_norm_detected_plot = 0;
datafile->GetObject(leading_pt_inside_gap_norm_measured_plot,leading_pt_inside_gap_norm_detected_plot);
if (leading_pt_inside_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *leading_pt_inside_gap_norm_generated = 0;
mcfile->GetObject(leading_pt_inside_gap_norm_true,leading_pt_inside_gap_norm_generated);
if (leading_pt_inside_gap_norm_generated == 0) { cout<< leading_pt_inside_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *leading_pt_inside_gap_norm_input_detected;
  leading_pt_inside_gap_norm_input_detected =  new TH1D("input_detected_leading_pt_inside_gap_norm","Leading pT Inside Gap;p_{T}^{inside};#frac{#sigma^{-1} d#sigma}{dp_{T}^{inside}}", in_nbins_meas, in_bins_meas);
  set_input(leading_pt_inside_gap_norm_detected, leading_pt_inside_gap_norm_input_detected);

  TH1D *leading_pt_inside_gap_norm_output_true;
  leading_pt_inside_gap_norm_output_true =  new TH1D("output_true_leading_pt_inside_gap_norm","Leading pT Inside Gap;p_{T}^{inside} [GeV];#frac{#sigma^{-1} d#sigma}{dp_{T}^{inside}}", in_nbins, in_bins);

  TH1D *leading_pt_inside_gap_norm_result = 0;

  unfold_distribution(leading_pt_inside_gap_norm_detected, leading_pt_inside_gap_norm_input_detected, leading_pt_inside_gap_norm_result, leading_pt_inside_gap_norm_output_true, leading_pt_inside_gap_norm_generated, leading_pt_inside_gap_response, leading_pt_inside_gap_response_all, plots_path, prefix, "leading_pt_inside_gap_norm", method, scale, tau, detail);


//Unfolding Leading Eta* Inside Gap
TString leading_eta_star_inside_gap_hname = "resp_leading_eta_star_inside_gap";
TString leading_eta_star_inside_gap_hname_all = "resp_leading_eta_star_inside_gap_all";
TString leading_eta_star_inside_gap_measured = "ak5PF_leading_eta_star_inside_gap_fine";
TString leading_eta_star_inside_gap_measured_plot = "ak5PF_leading_eta_star_inside_gap";
TString leading_eta_star_inside_gap_true;
if (true_is_corr_data)
	{
	leading_eta_star_inside_gap_true = "ak5PF_leading_eta_star_inside_gap";
	}
else
	{
	leading_eta_star_inside_gap_true = "ak5Gen_leading_eta_star_inside_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *leading_eta_star_inside_gap_response = 0;
respfile->GetObject(leading_eta_star_inside_gap_hname,leading_eta_star_inside_gap_response);
if (leading_eta_star_inside_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *leading_eta_star_inside_gap_response_all = 0;
respfile->GetObject(leading_eta_star_inside_gap_hname_all,leading_eta_star_inside_gap_response_all);
if (leading_eta_star_inside_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *leading_eta_star_inside_gap_detected = 0;
datafile->GetObject(leading_eta_star_inside_gap_measured,leading_eta_star_inside_gap_detected);
if (leading_eta_star_inside_gap_detected == 0) { cout<<"Measured data " << leading_eta_star_inside_gap_measured << " not found!"<<endl; return; }

TH1D *leading_eta_star_inside_gap_detected_plot = 0;
datafile->GetObject(leading_eta_star_inside_gap_measured_plot,leading_eta_star_inside_gap_detected_plot);
if (leading_eta_star_inside_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *leading_eta_star_inside_gap_generated = 0;
mcfile->GetObject(leading_eta_star_inside_gap_true,leading_eta_star_inside_gap_generated);
if (leading_eta_star_inside_gap_generated == 0) { cout<< leading_eta_star_inside_gap_true << " true data not found!"<<endl; return; }

  TH1D *leading_eta_star_inside_gap_input_detected;
  leading_eta_star_inside_gap_input_detected =  new TH1D("input_detected_leading_eta_star_inside_gap","Leading Eta* Inside Gap;#eta*^{inside};Events", etastar_nbins_meas, etastar_bins_meas);
  set_input(leading_eta_star_inside_gap_detected, leading_eta_star_inside_gap_input_detected);

  TH1D *leading_eta_star_inside_gap_output_true;
  leading_eta_star_inside_gap_output_true =  new TH1D("output_true_leading_eta_star_inside_gap","Leading Eta* Inside Gap;#eta*^{inside};Events", etastar_nbins, etastar_bins);

  TH1D *leading_eta_star_inside_gap_result = 0;

  unfold_distribution(leading_eta_star_inside_gap_detected, leading_eta_star_inside_gap_input_detected, leading_eta_star_inside_gap_result, leading_eta_star_inside_gap_output_true, leading_eta_star_inside_gap_generated, leading_eta_star_inside_gap_response, leading_eta_star_inside_gap_response_all, plots_path, prefix, "leading_eta_star_inside_gap", method, scale, tau, detail);


//Unfolding Leading Eta* Inside Gap Norm
TString leading_eta_star_inside_gap_norm_measured = "ak5PF_leading_eta_star_inside_gap_norm_fine";
TString leading_eta_star_inside_gap_norm_measured_plot = "ak5PF_leading_eta_star_inside_gap_norm";
TString leading_eta_star_inside_gap_norm_true;
if (true_is_corr_data)
	{
	leading_eta_star_inside_gap_norm_true = "ak5PF_leading_eta_star_inside_gap_norm";
	}
else
	{
	leading_eta_star_inside_gap_norm_true = "ak5Gen_leading_eta_star_inside_gap_norm";
	}

//get the objects on a safe way
TH1D *leading_eta_star_inside_gap_norm_detected = 0;
datafile->GetObject(leading_eta_star_inside_gap_norm_measured,leading_eta_star_inside_gap_norm_detected);
if (leading_eta_star_inside_gap_norm_detected == 0) { cout<<"Measured data " << leading_eta_star_inside_gap_norm_measured << " not found!"<<endl; return; }

TH1D *leading_eta_star_inside_gap_norm_detected_plot = 0;
datafile->GetObject(leading_eta_star_inside_gap_norm_measured_plot,leading_eta_star_inside_gap_norm_detected_plot);
if (leading_eta_star_inside_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *leading_eta_star_inside_gap_norm_generated = 0;
mcfile->GetObject(leading_eta_star_inside_gap_norm_true,leading_eta_star_inside_gap_norm_generated);
if (leading_eta_star_inside_gap_norm_generated == 0) { cout<< leading_eta_star_inside_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *leading_eta_star_inside_gap_norm_input_detected;
  leading_eta_star_inside_gap_norm_input_detected =  new TH1D("input_detected_leading_eta_star_inside_gap_norm","Leading Eta* Inside Gap;#eta*^{inside};#frac{#sigma^{-1} d#sigma}{d#eta*}", etastar_nbins_meas, etastar_bins_meas);
  set_input(leading_eta_star_inside_gap_norm_detected, leading_eta_star_inside_gap_norm_input_detected);

  TH1D *leading_eta_star_inside_gap_norm_output_true;
  leading_eta_star_inside_gap_norm_output_true =  new TH1D("output_true_leading_eta_star_inside_gap_norm","Leading Eta* Inside Gap;#eta*^{inside};#frac{#sigma^{-1} d#sigma}{d#eta*}", etastar_nbins, etastar_bins);

  TH1D *leading_eta_star_inside_gap_norm_result = 0;

  unfold_distribution(leading_eta_star_inside_gap_norm_detected, leading_eta_star_inside_gap_norm_input_detected, leading_eta_star_inside_gap_norm_result, leading_eta_star_inside_gap_norm_output_true, leading_eta_star_inside_gap_norm_generated, leading_eta_star_inside_gap_response, leading_eta_star_inside_gap_response_all, plots_path, prefix, "leading_eta_star_inside_gap_norm", method, scale, tau, detail);

//Unfolding Leading pT Outside Gap
TString leading_pt_outside_gap_hname = "resp_leading_pt_outside_gap";
TString leading_pt_outside_gap_hname_all = "resp_leading_pt_outside_gap_all";
TString leading_pt_outside_gap_measured = "ak5PF_leading_pt_outside_gap_fine";
TString leading_pt_outside_gap_measured_plot = "ak5PF_leading_pt_outside_gap";
TString leading_pt_outside_gap_true;
if (true_is_corr_data)
	{
	leading_pt_outside_gap_true = "ak5PF_leading_pt_outside_gap";
	}
else
	{
	leading_pt_outside_gap_true = "ak5Gen_leading_pt_outside_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *leading_pt_outside_gap_response = 0;
respfile->GetObject(leading_pt_outside_gap_hname,leading_pt_outside_gap_response);
if (leading_pt_outside_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *leading_pt_outside_gap_response_all = 0;
respfile->GetObject(leading_pt_outside_gap_hname_all,leading_pt_outside_gap_response_all);
if (leading_pt_outside_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *leading_pt_outside_gap_detected = 0;
datafile->GetObject(leading_pt_outside_gap_measured,leading_pt_outside_gap_detected);
if (leading_pt_outside_gap_detected == 0) { cout<<"Measured data " << leading_pt_outside_gap_measured << " not found!"<<endl; return; }

TH1D *leading_pt_outside_gap_detected_plot = 0;
datafile->GetObject(leading_pt_outside_gap_measured_plot,leading_pt_outside_gap_detected_plot);
if (leading_pt_outside_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *leading_pt_outside_gap_generated = 0;
mcfile->GetObject(leading_pt_outside_gap_true,leading_pt_outside_gap_generated);
if (leading_pt_outside_gap_generated == 0) { cout<< leading_pt_outside_gap_true << " true data not found!"<<endl; return; }

  TH1D *leading_pt_outside_gap_input_detected;
  leading_pt_outside_gap_input_detected =  new TH1D("input_detected_leading_pt_outside_gap","Leading pT Outside Gap;p_{T}^{outside};Events", out_nbins_meas, out_bins_meas);
  set_input(leading_pt_outside_gap_detected, leading_pt_outside_gap_input_detected);

  TH1D *leading_pt_outside_gap_output_true;
  leading_pt_outside_gap_output_true =  new TH1D("output_true_leading_pt_outside_gap","Leading pT Outside Gap;p_{T}^{outside} [GeV];Events", out_nbins, out_bins);

  TH1D *leading_pt_outside_gap_result = 0;

  unfold_distribution(leading_pt_outside_gap_detected, leading_pt_outside_gap_input_detected, leading_pt_outside_gap_result, leading_pt_outside_gap_output_true, leading_pt_outside_gap_generated, leading_pt_outside_gap_response, leading_pt_outside_gap_response_all, plots_path, prefix, "leading_pt_outside_gap", method, scale, tau, detail);


//Unfolding Leading pT Outside Gap Norm
TString leading_pt_outside_gap_norm_measured = "ak5PF_leading_pt_outside_gap_norm_fine";
TString leading_pt_outside_gap_norm_measured_plot = "ak5PF_leading_pt_outside_gap_norm";
TString leading_pt_outside_gap_norm_true;
if (true_is_corr_data)
	{
	leading_pt_outside_gap_norm_true = "ak5PF_leading_pt_outside_gap_norm";
	}
else
	{
	leading_pt_outside_gap_norm_true = "ak5Gen_leading_pt_outside_gap_norm";
	}

//get the objects on a safe way
TH1D *leading_pt_outside_gap_norm_detected = 0;
datafile->GetObject(leading_pt_outside_gap_norm_measured,leading_pt_outside_gap_norm_detected);
if (leading_pt_outside_gap_norm_detected == 0) { cout<<"Measured data " << leading_pt_outside_gap_norm_measured << " not found!"<<endl; return; }

TH1D *leading_pt_outside_gap_norm_detected_plot = 0;
datafile->GetObject(leading_pt_outside_gap_norm_measured_plot,leading_pt_outside_gap_norm_detected_plot);
if (leading_pt_outside_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *leading_pt_outside_gap_norm_generated = 0;
mcfile->GetObject(leading_pt_outside_gap_norm_true,leading_pt_outside_gap_norm_generated);
if (leading_pt_outside_gap_norm_generated == 0) { cout<< leading_pt_outside_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *leading_pt_outside_gap_norm_input_detected;
  leading_pt_outside_gap_norm_input_detected =  new TH1D("input_detected_leading_pt_outside_gap_norm","Leading pT Outside Gap;p_{T}^{outside};#frac{#sigma^{-1} d#sigma}{dp_{T}^{ouside}}", out_nbins_meas, out_bins_meas);
  set_input(leading_pt_outside_gap_norm_detected, leading_pt_outside_gap_norm_input_detected);

  TH1D *leading_pt_outside_gap_norm_output_true;
  leading_pt_outside_gap_norm_output_true =  new TH1D("output_true_leading_pt_outside_gap_norm","Leading pT Outside Gap;p_{T}^{outside} [GeV];#frac{#sigma^{-1} d#sigma}{dp_{T}^{ouside}}", out_nbins, out_bins);

  TH1D *leading_pt_outside_gap_norm_result = 0;

  unfold_distribution(leading_pt_outside_gap_norm_detected, leading_pt_outside_gap_norm_input_detected, leading_pt_outside_gap_norm_result, leading_pt_outside_gap_norm_output_true, leading_pt_outside_gap_norm_generated, leading_pt_outside_gap_response, leading_pt_outside_gap_response_all, plots_path, prefix, "leading_pt_outside_gap_norm", method, scale, tau, detail);


//Unfolding Deta Eta Outside Gap
TString delta_eta_outside_gap_hname = "resp_delta_eta_outside_gap";
TString delta_eta_outside_gap_hname_all = "resp_delta_eta_outside_gap_all";
TString delta_eta_outside_gap_measured = "ak5PF_delta_eta_outside_gap_fine";
TString delta_eta_outside_gap_measured_plot = "ak5PF_delta_eta_outside_gap";
TString delta_eta_outside_gap_true;
if (true_is_corr_data)
	{
	delta_eta_outside_gap_true = "ak5PF_delta_eta_outside_gap";
	}
else
	{
	delta_eta_outside_gap_true = "ak5Gen_delta_eta_outside_gap";
	}

//get the objects on a safe way
RooUnfoldResponse *delta_eta_outside_gap_response = 0;
respfile->GetObject(delta_eta_outside_gap_hname,delta_eta_outside_gap_response);
if (delta_eta_outside_gap_response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *delta_eta_outside_gap_response_all = 0;
respfile->GetObject(delta_eta_outside_gap_hname_all,delta_eta_outside_gap_response_all);
if (delta_eta_outside_gap_response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *delta_eta_outside_gap_detected = 0;
datafile->GetObject(delta_eta_outside_gap_measured,delta_eta_outside_gap_detected);
if (delta_eta_outside_gap_detected == 0) { cout<<"Measured data " << delta_eta_outside_gap_measured << " not found!"<<endl; return; }

TH1D *delta_eta_outside_gap_detected_plot = 0;
datafile->GetObject(delta_eta_outside_gap_measured_plot,delta_eta_outside_gap_detected_plot);
if (delta_eta_outside_gap_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_eta_outside_gap_generated = 0;
mcfile->GetObject(delta_eta_outside_gap_true,delta_eta_outside_gap_generated);
if (delta_eta_outside_gap_generated == 0) { cout<< delta_eta_outside_gap_true << " true data not found!"<<endl; return; }

  TH1D *delta_eta_outside_gap_input_detected;
  delta_eta_outside_gap_input_detected =  new TH1D("input_detected_delta_eta_outside_gap","#Delta#eta Outside Gap;#Delta#eta^{outside};Events", deta_out_nbins_meas, deta_out_bins_meas);
  set_input(delta_eta_outside_gap_detected, delta_eta_outside_gap_input_detected);

  TH1D *delta_eta_outside_gap_output_true;
  delta_eta_outside_gap_output_true =  new TH1D("output_true_delta_eta_outside_gap","#Delta#Eta Outside Gap;#Delta#eta^{outside};Events", deta_out_nbins, deta_out_bins);

  TH1D *delta_eta_outside_gap_result = 0;

  unfold_distribution(delta_eta_outside_gap_detected, delta_eta_outside_gap_input_detected, delta_eta_outside_gap_result, delta_eta_outside_gap_output_true, delta_eta_outside_gap_generated, delta_eta_outside_gap_response, delta_eta_outside_gap_response_all, plots_path, prefix, "delta_eta_outside_gap", method, scale, tau, detail);


//Unfolding Deta Eta Outside Gap Norm
TString delta_eta_outside_gap_norm_measured = "ak5PF_delta_eta_outside_gap_norm_fine";
TString delta_eta_outside_gap_norm_measured_plot = "ak5PF_delta_eta_outside_gap_norm";
TString delta_eta_outside_gap_norm_true;
if (true_is_corr_data)
	{
	delta_eta_outside_gap_norm_true = "ak5PF_delta_eta_outside_gap_norm";
	}
else
	{
	delta_eta_outside_gap_norm_true = "ak5Gen_delta_eta_outside_gap_norm";
	}

//get the objects on a safe way
TH1D *delta_eta_outside_gap_norm_detected = 0;
datafile->GetObject(delta_eta_outside_gap_norm_measured,delta_eta_outside_gap_norm_detected);
if (delta_eta_outside_gap_norm_detected == 0) { cout<<"Measured data " << delta_eta_outside_gap_norm_measured << " not found!"<<endl; return; }

TH1D *delta_eta_outside_gap_norm_detected_plot = 0;
datafile->GetObject(delta_eta_outside_gap_norm_measured_plot,delta_eta_outside_gap_norm_detected_plot);
if (delta_eta_outside_gap_norm_detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *delta_eta_outside_gap_norm_generated = 0;
mcfile->GetObject(delta_eta_outside_gap_norm_true,delta_eta_outside_gap_norm_generated);
if (delta_eta_outside_gap_norm_generated == 0) { cout<< delta_eta_outside_gap_norm_true << " true data not found!"<<endl; return; }

  TH1D *delta_eta_outside_gap_norm_input_detected;
  delta_eta_outside_gap_norm_input_detected =  new TH1D("input_detected_delta_eta_outside_gap_norm","#Delta#eta Outside Gap;#Delta#eta^{outside};#frac{#sigma^{/1} d#sigma}{d#Delta#eta}", deta_out_nbins_meas, deta_out_bins_meas);
  set_input(delta_eta_outside_gap_norm_detected, delta_eta_outside_gap_norm_input_detected);

  TH1D *delta_eta_outside_gap_norm_output_true;
  delta_eta_outside_gap_norm_output_true =  new TH1D("output_true_delta_eta_outside_gap_norm","#Delta#Eta Outside Gap;#Delta#eta^{outside};#frac{#sigma^{/1} d#sigma}{d#Delta#eta}", deta_out_nbins, deta_out_bins);

  TH1D *delta_eta_outside_gap_norm_result = 0;

  unfold_distribution(delta_eta_outside_gap_norm_detected, delta_eta_outside_gap_norm_input_detected, delta_eta_outside_gap_norm_result, delta_eta_outside_gap_norm_output_true, delta_eta_outside_gap_norm_generated, delta_eta_outside_gap_response, delta_eta_outside_gap_response_all, plots_path, prefix, "delta_eta_outside_gap_norm", method, scale, tau, detail);

//normalized
normalize_histogram(delta_phi_norm_output_true, "Delta Phi Norm", true, true);
normalize_histogram(delta_phi_gap_norm_output_true, "Delta Phi Gap Norm", true, true);
normalize_histogram(delta_phi_nogap_norm_output_true, "Delta Phi noGap Norm", true, true);
normalize_histogram(delta_phi_deta1_norm_output_true, "Delta Phi Deta1 Norm", true, true);
normalize_histogram(delta_phi_deta2_norm_output_true, "Delta Phi Deta2 Norm", true, true);
normalize_histogram(delta_phi_deta3_norm_output_true, "Delta Phi Deta3 Norm", true, true);
normalize_histogram(delta_phi_deta4_norm_output_true, "Delta Phi Deta4 Norm", true, true);
normalize_histogram(delta_phi_deta1_gap_norm_output_true, "Delta Phi Deta1 Gap Norm", true, true);
normalize_histogram(delta_phi_deta2_gap_norm_output_true, "Delta Phi Deta2 Gap Norm", true, true);
normalize_histogram(delta_phi_deta3_gap_norm_output_true, "Delta Phi Deta3 Gap Norm", true, true);
normalize_histogram(delta_phi_deta4_gap_norm_output_true, "Delta Phi Deta4 Gap Norm", true, true);
normalize_histogram(delta_phi_deta1_nogap_norm_output_true, "Delta Phi Deta1 noGap Norm", true, true);
normalize_histogram(delta_phi_deta2_nogap_norm_output_true, "Delta Phi Deta2 noGap Norm", true, true);
normalize_histogram(delta_phi_deta3_nogap_norm_output_true, "Delta Phi Deta3 noGap Norm", true, true);
normalize_histogram(delta_phi_deta4_nogap_norm_output_true, "Delta Phi Deta4 noGap Norm", true, true);
normalize_histogram(leading_pt_inside_gap_norm_output_true, "Leading pT Inside Gap Norm", true, true);
normalize_histogram(leading_eta_star_inside_gap_norm_output_true, "Leading Eta* Outside Gap Norm", true, true);
normalize_histogram(delta_eta_outside_gap_norm_output_true, "Delta Eta Outside Gap Norm", true, true);
normalize_histogram(leading_pt_outside_gap_norm_output_true, "Leading pT Outside Gap Norm", true, true);


// recreating the output file
TFile *f = TFile::Open(root_out.c_str(), "RECREATE");

// save the unfolded result in the output file
delta_phi_output_true->Write();
delta_phi_deta1_output_true->Write();
delta_phi_deta2_output_true->Write();
delta_phi_deta3_output_true->Write();
delta_phi_deta4_output_true->Write();
delta_phi_gap_output_true->Write();
delta_phi_deta1_gap_output_true->Write();
delta_phi_deta2_gap_output_true->Write();
delta_phi_deta3_gap_output_true->Write();
delta_phi_deta4_gap_output_true->Write();
delta_phi_nogap_output_true->Write();
delta_phi_deta1_nogap_output_true->Write();
delta_phi_deta2_nogap_output_true->Write();
delta_phi_deta3_nogap_output_true->Write();
delta_phi_deta4_nogap_output_true->Write();
leading_pt_inside_gap_output_true->Write();
leading_eta_star_inside_gap_output_true->Write();
leading_pt_outside_gap_output_true->Write();
delta_eta_outside_gap_output_true->Write();
delta_phi_norm_output_true->Write();
delta_phi_deta1_norm_output_true->Write();
delta_phi_deta2_norm_output_true->Write();
delta_phi_deta3_norm_output_true->Write();
delta_phi_deta4_norm_output_true->Write();
delta_phi_gap_norm_output_true->Write();
delta_phi_deta1_gap_norm_output_true->Write();
delta_phi_deta2_gap_norm_output_true->Write();
delta_phi_deta3_gap_norm_output_true->Write();
delta_phi_deta4_gap_norm_output_true->Write();
delta_phi_nogap_norm_output_true->Write();
delta_phi_deta1_nogap_norm_output_true->Write();
delta_phi_deta2_nogap_norm_output_true->Write();
delta_phi_deta3_nogap_norm_output_true->Write();
delta_phi_deta4_nogap_norm_output_true->Write();
leading_pt_inside_gap_norm_output_true->Write();
leading_eta_star_inside_gap_norm_output_true->Write();
leading_pt_outside_gap_norm_output_true->Write();
delta_eta_outside_gap_norm_output_true->Write();
cout<<"Unfolding result saved to : "<<root_out<<endl;

//close output file
f->Close();

//delete the variables to avoid memory leak
delete f;
delete delta_phi_output_true;
delete delta_phi_result;
delete delta_phi_deta1_output_true;
delete delta_phi_deta1_result;
delete delta_phi_deta2_output_true;
delete delta_phi_deta2_result;
delete delta_phi_deta3_output_true;
delete delta_phi_deta3_result;
delete delta_phi_deta4_output_true;
delete delta_phi_deta4_result;
delete delta_phi_gap_output_true;
delete delta_phi_gap_result;
delete delta_phi_deta1_gap_output_true;
delete delta_phi_deta1_gap_result;
delete delta_phi_deta2_gap_output_true;
delete delta_phi_deta2_gap_result;
delete delta_phi_deta3_gap_output_true;
delete delta_phi_deta3_gap_result;
delete delta_phi_deta4_gap_output_true;
delete delta_phi_deta4_gap_result;
delete delta_phi_nogap_output_true;
delete delta_phi_nogap_result;
delete delta_phi_deta1_nogap_output_true;
delete delta_phi_deta1_nogap_result;
delete delta_phi_deta2_nogap_output_true;
delete delta_phi_deta2_nogap_result;
delete delta_phi_deta3_nogap_output_true;
delete delta_phi_deta3_nogap_result;
delete delta_phi_deta4_nogap_output_true;
delete delta_phi_deta4_nogap_result;
delete leading_pt_inside_gap_output_true;
delete leading_pt_inside_gap_result;
delete leading_eta_star_inside_gap_output_true;
delete leading_eta_star_inside_gap_result;
delete leading_pt_outside_gap_output_true;
delete leading_pt_outside_gap_result;
delete delta_eta_outside_gap_output_true;
delete delta_eta_outside_gap_result;
delete delta_phi_norm_output_true;
delete delta_phi_norm_result;
delete delta_phi_deta1_norm_output_true;
delete delta_phi_deta1_norm_result;
delete delta_phi_deta2_norm_output_true;
delete delta_phi_deta2_norm_result;
delete delta_phi_deta3_norm_output_true;
delete delta_phi_deta3_norm_result;
delete delta_phi_deta4_norm_output_true;
delete delta_phi_deta4_norm_result;
delete delta_phi_gap_norm_output_true;
delete delta_phi_gap_norm_result;
delete delta_phi_deta1_gap_norm_output_true;
delete delta_phi_deta1_gap_norm_result;
delete delta_phi_deta2_gap_norm_output_true;
delete delta_phi_deta2_gap_norm_result;
delete delta_phi_deta3_gap_norm_output_true;
delete delta_phi_deta3_gap_norm_result;
delete delta_phi_deta4_gap_norm_output_true;
delete delta_phi_deta4_gap_norm_result;
delete delta_phi_nogap_norm_output_true;
delete delta_phi_nogap_norm_result;
delete delta_phi_deta1_nogap_norm_output_true;
delete delta_phi_deta1_nogap_norm_result;
delete delta_phi_deta2_nogap_norm_output_true;
delete delta_phi_deta2_nogap_norm_result;
delete delta_phi_deta3_nogap_norm_output_true;
delete delta_phi_deta3_nogap_norm_result;
delete delta_phi_deta4_nogap_norm_output_true;
delete delta_phi_deta4_nogap_norm_result;
delete leading_pt_inside_gap_norm_output_true;
delete leading_pt_inside_gap_norm_result;
delete leading_eta_star_inside_gap_norm_output_true;
delete leading_eta_star_inside_gap_norm_result;
delete leading_pt_outside_gap_norm_output_true;
delete leading_pt_outside_gap_norm_result;
delete delta_eta_outside_gap_norm_output_true;
delete delta_eta_outside_gap_norm_result;
}
