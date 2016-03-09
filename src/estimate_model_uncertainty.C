// Pedro Cipriano, Nov 2011
// DESY, CMS
// Last Update: 22 Mar 2013
//
// estimate_model_uncertainty()
// calculate the model uncertainty for pythia 8 - tune 1, pyhtia 6 - tune z2 and herwig 6

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

void calc_model_unc(TH1D *histogram, TH1D *hist_data, TH1D *hist_p8, TH1D *hist_p6_z2, double *final_model_unc, int index, bool detail)
{
// calculates the model uncertainty for a given distribution taking the difference between p6_z2 and h6

//declaring variables
double min = 0.0;
double max = 0.0;
double tot = 0.0;
double ave = 0.0;

double sist = 0.0;
double sist_error = 0.0;

double cont = 0.0;

double corr_p8 = 0.0;
double er_p8 = 0.0;
double sist_p8 = 0.0;

double corr_p6_z2 = 0.0;
double er_p6_z2 = 0.0;
double sist_p6_z2 = 0.0;

//loops over the bins, gets and saves the deviations
    for(Int_t i=1;i<=hist_data->GetNbinsX();i++)
    {
//loads the relevant histograms contents for a given bin
    cont = hist_data->GetBinContent(i);

    corr_p8 = hist_p8->GetBinContent(i);
    er_p8 = hist_p8->GetBinError(i);
    sist_p8 = cont*corr_p8;

    corr_p6_z2 = hist_p6_z2->GetBinContent(i);
    er_p6_z2 = hist_p6_z2->GetBinError(i);
    sist_p6_z2 = cont*corr_p6_z2;

//computes the model uncertainty using only p6 and p8
	cont = cont*.5*(corr_p6_z2+corr_p8);
	sist = abs(sist_p6_z2 - sist_p8)/2.0;
	sist_error = (er_p6_z2+er_p8)/2.0;

//makes the sistematic relative
    sist = sist/cont;

//control output of the calculation
//    if (detail) { cout<<"cont before correction = "<<hist_data->GetBinContent(i)<<", after correction = "<<cont<<", sistematic = "<<sist<<", sistematic error = "<<sist_error<<endl; }

//checks for minimum, maximum and sums the uncertainties
    if (sist > max) {max = sist;}
    tot = tot + sist;
    if (sist < min || i == 1) { min = sist;}

//save the results in a histogram
    histogram->SetBinError(i,sist_error);
    histogram->SetBinContent(i,sist);
    }

//calculates the average model uncertainty
    ave = tot/hist_data->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    final_model_unc[index*3+0] = ave*100;
    final_model_unc[index*3+1] = min*100;
    final_model_unc[index*3+2] = max*100;

}


void calc_model_unc2(TH1D *histogram, TH1D *hist_p8, TH1D *hist_p6_z2, double *final_model_unc, int index, bool detail)
{
// calculates the model uncertainty for a given distribution taking the difference between p6_z2 and h6

//declaring variables
double min = 0.0;
double max = 0.0;
double tot = 0.0;
double ave = 0.0;

double sist = 0.0;
double sist_error = 0.0;

double cont = 0.0;

double cont_p8 = 0.0;
double er_p8 = 0.0;

double cont_p6_z2 = 0.0;
double er_p6_z2 = 0.0;

//loops over the bins, gets and saves the deviations
    for(Int_t i=1;i<=hist_p8->GetNbinsX();i++)
    {
//loads the relevant histograms contents for a given bin
    cont_p8 = hist_p8->GetBinContent(i);
    er_p8 = hist_p8->GetBinError(i);

    cont_p6_z2 = hist_p6_z2->GetBinContent(i);
    er_p6_z2 = hist_p6_z2->GetBinError(i);

//computes the model uncertainty using only p6 and p8
	cont = .5*(cont_p6_z2+cont_p8);
	sist = abs(cont_p6_z2 - cont_p8)/2.0;
	sist_error = (er_p6_z2+er_p8)/(2.0*cont);

//makes the sistematic relative
    sist = sist/cont;

//control output of the calculation
//    if (detail) { cout<<"cont before correction = "<<hist_data->GetBinContent(i)<<", after correction = "<<cont<<", sistematic = "<<sist<<", sistematic error = "<<sist_error<<endl; }

//checks for minimum, maximum and sums the uncertainties
    if (sist > max) {max = sist;}
    tot = tot + sist;
    if (sist < min || i == 1) { min = sist;}

//save the results in a histogram
    histogram->SetBinError(i,sist_error);
    histogram->SetBinContent(i,sist);
    }

//calculates the average model uncertainty
    ave = tot/hist_p8->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    final_model_unc[index*3+0] = ave*100;
    final_model_unc[index*3+1] = min*100;
    final_model_unc[index*3+2] = max*100;

}

void show_model_uncertainties(double *final_model_unc)
{
//shows the computed model uncertainties
    cout<<" "<<endl;
    cout<<"Model Uncertainty"<<endl;
    cout<<"Observable                  	    Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   	    "<<final_model_unc[0]<<"  "<<final_model_unc[1]<<" "<<final_model_unc[2]<<endl;
    cout<<"Delta phi deta1             	    "<<final_model_unc[3]<<"  "<<final_model_unc[4]<<"  "<<final_model_unc[5]<<endl;
    cout<<"Delta phi deta2             	    "<<final_model_unc[6]<<"   "<<final_model_unc[7]<<" "<<final_model_unc[8]<<endl;
    cout<<"Delta phi deta3             	    "<<final_model_unc[9]<<"  "<<final_model_unc[10]<<" "<<final_model_unc[11]<<endl;
    cout<<"Delta phi deta4             	    "<<final_model_unc[12]<<"   "<<final_model_unc[13]<<" "<<final_model_unc[14]<<endl;
    cout<<"Delta phi gap               	    "<<final_model_unc[15]<<"  "<<final_model_unc[16]<<" "<<final_model_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         	    "<<final_model_unc[18]<<"  "<<final_model_unc[19]<<"  "<<final_model_unc[20]<<endl;
    cout<<"Delta phi deta2 gap        	    "<<final_model_unc[21]<<"  "<<final_model_unc[22]<<" "<<final_model_unc[23]<<endl;
    cout<<"Delta phi deta3 gap        	    "<<final_model_unc[24]<<"  "<<final_model_unc[25]<<" "<<final_model_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         	    "<<final_model_unc[27]<<"  "<<final_model_unc[28]<<" "<<final_model_unc[29]<<endl;
    cout<<"Delta phi nogap                  "<<final_model_unc[30]<<"  "<<final_model_unc[31]<<" "<<final_model_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap            "<<final_model_unc[33]<<"  "<<final_model_unc[34]<<"  "<<final_model_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap            "<<final_model_unc[36]<<"   "<<final_model_unc[37]<<"  "<<final_model_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap            "<<final_model_unc[39]<<"  "<<final_model_unc[40]<<"  "<<final_model_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap            "<<final_model_unc[42]<<"  "<<final_model_unc[43]<<"  "<<final_model_unc[44]<<endl;
    cout<<"Leading pT inside gap            "<<final_model_unc[45]<<"  "<<final_model_unc[46]<<" "<<final_model_unc[47]<<endl;
    cout<<"Leading eta star inside gap      "<<final_model_unc[48]<<"  "<<final_model_unc[49]<<"  "<<final_model_unc[50]<<endl;
    cout<<"Delta eta outside gap            "<<final_model_unc[51]<<"  "<<final_model_unc[52]<<" "<<final_model_unc[53]<<endl;
    cout<<"Leading pT outside gap           "<<final_model_unc[54]<<"  "<<final_model_unc[55]<<" "<<final_model_unc[56]<<endl;
    cout<<"Delta phi Norm                   "<<final_model_unc[57]<<"  "<<final_model_unc[58]<<" "<<final_model_unc[59]<<endl;
    cout<<"Delta phi deta1 Norm             "<<final_model_unc[60]<<"  "<<final_model_unc[61]<<"  "<<final_model_unc[62]<<endl;
    cout<<"Delta phi deta2 Norm             "<<final_model_unc[63]<<"   "<<final_model_unc[64]<<" "<<final_model_unc[65]<<endl;
    cout<<"Delta phi deta3 Norm             "<<final_model_unc[66]<<"  "<<final_model_unc[67]<<" "<<final_model_unc[68]<<endl;
    cout<<"Delta phi deta4 Norm             "<<final_model_unc[69]<<"   "<<final_model_unc[70]<<" "<<final_model_unc[71]<<endl;
    cout<<"Delta phi gap Norm               "<<final_model_unc[72]<<"  "<<final_model_unc[73]<<" "<<final_model_unc[74]<<endl;
    cout<<"Delta phi deta1 gap Norm         "<<final_model_unc[75]<<"  "<<final_model_unc[76]<<"  "<<final_model_unc[77]<<endl;
    cout<<"Delta phi deta2 gap Norm         "<<final_model_unc[78]<<"  "<<final_model_unc[79]<<" "<<final_model_unc[80]<<endl;
    cout<<"Delta phi deta3 gap Norm         "<<final_model_unc[81]<<"  "<<final_model_unc[82]<<" "<<final_model_unc[83]<<endl;
    cout<<"Delta phi deta4 gap Norm         "<<final_model_unc[84]<<"  "<<final_model_unc[85]<<" "<<final_model_unc[86]<<endl;
    cout<<"Delta phi nogap Norm             "<<final_model_unc[87]<<"  "<<final_model_unc[88]<<" "<<final_model_unc[89]<<endl;
    cout<<"Delta phi deta1 nogap Norm       "<<final_model_unc[90]<<"  "<<final_model_unc[91]<<"  "<<final_model_unc[92]<<endl;
    cout<<"Delta phi deta2 nogap Norm       "<<final_model_unc[93]<<"   "<<final_model_unc[94]<<"  "<<final_model_unc[95]<<endl;
    cout<<"Delta phi deta3 nogap Norm       "<<final_model_unc[96]<<"  "<<final_model_unc[97]<<"  "<<final_model_unc[98]<<endl;
    cout<<"Delta phi deta4 nogap Norm       "<<final_model_unc[99]<<"  "<<final_model_unc[100]<<"  "<<final_model_unc[101]<<endl;
    cout<<"Leading pT inside gap Norm       "<<final_model_unc[102]<<"  "<<final_model_unc[103]<<" "<<final_model_unc[104]<<endl;
    cout<<"Leading eta star inside gap Norm "<<final_model_unc[105]<<"  "<<final_model_unc[106]<<"  "<<final_model_unc[107]<<endl;
    cout<<"Delta eta outside gap Norm       "<<final_model_unc[108]<<"  "<<final_model_unc[109]<<" "<<final_model_unc[110]<<endl;
    cout<<"Leading pT outside gap Norm      "<<final_model_unc[111]<<"  "<<final_model_unc[112]<<" "<<final_model_unc[113]<<endl;
}

void estimate_model_uncertainty(string path_pythia6, string path_pythia8, string path_data, string out_model_uncertainty, string label, string output_path_plots = "../output/model_uncertainty/", bool detail = false, bool disp_uncertainty = true, bool with_unfold = false, bool test = false)
{
//main estimate model uncertainty routine
//output detail: true-output all the major steps; false-run quietly
//disp_uncertainty: true-display the model uncartainty estimation values; false-dont show anything

//outputs the configuration
    if (detail) { cout << "Estimate Model Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for Pythia 6 - Tune Z2*: " << path_pythia6 << endl; }
    if (detail) { cout << "Input path for Pythia 8 - Tune 4C:  " << path_pythia8 << endl; }
    if (detail) { cout << "Input path for Data:                " << path_data << endl; }
    if (detail) { cout << "Histogram label:                    " << label << endl; }
    if (detail) { cout << "Output path:                        " << out_model_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:                  " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:                       " << detail << endl; }
    if (detail) { cout << "Display Results:                    " << disp_uncertainty << endl; }
    if (detail) { cout << "With Unfolded Data:                 " << with_unfold << endl; }
    if (detail) { cout << "Test Mode:                          " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the correction files" << endl; }
    TFile *pythia6_file = new TFile( path_pythia6.c_str() );
    TFile *pythia8_file = new TFile( path_pythia8.c_str() );
    if (detail) { cout << "Opening the data file" << endl; }
    TFile *data_file = new TFile( path_data.c_str() );


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

//declaring and initializing the uncertainties array
    double final_model_unc[38*3];

    for (int i=0; i<= 38*3-1;i++)
    {
    final_model_unc[i] = 0.0;
    }

//starts the estimation of the model uncertainty

//estimate the model uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_data = 0;
    TH1D *delta_phi_p6_z2 = 0;
    TH1D *delta_phi_p8_4c = 0;
    TString delta_phi_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi";
    TString delta_phi_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi";
    TString delta_phi_generic_name = "ak5PF_delta_phi";
    TString delta_phi_unfold_name = "output_true_delta_phi";

    data_file->GetObject(delta_phi_generic_name,delta_phi_data);
    if (delta_phi_data == 0) { cout << delta_phi_generic_name << " not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_unfold_name,delta_phi_p6_z2);
	if (delta_phi_p6_z2 == 0) { cout << delta_phi_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_p6_z2_name,delta_phi_p6_z2);
	if (delta_phi_p6_z2 == 0) { cout << delta_phi_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_unfold_name,delta_phi_p8_4c);
	if (delta_phi_p8_4c == 0) { cout << delta_phi_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_p8_4c_name,delta_phi_p8_4c);
	if (delta_phi_p8_4c == 0) { cout << delta_phi_p8_4c_name << " not found!" << endl; return; }
	}

    
    TH1D *delta_phi;
    delta_phi =  new TH1D("model_unc_delta_phi","Model Uncertainty;#Delta#phi;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi, delta_phi_p6_z2, delta_phi_p8_4c, final_model_unc, 0, detail);
	}
    else
	{
	calc_model_unc(delta_phi, delta_phi_data, delta_phi_p6_z2, delta_phi_p8_4c, final_model_unc, 0, detail);
	}
    delta_phi->SetEntries(delta_phi_data->GetEntries());
    plot_histogram(delta_phi, output_path_plots, "model_uncertainty_" + label + "_delta_phi", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi norm distribution
    if (detail) { cout<<"Delta phi Norm"<<endl; }

    TH1D *delta_phi_norm_data = 0;
    TH1D *delta_phi_norm_p6_z2 = 0;
    TH1D *delta_phi_norm_p8_4c = 0;
    TString delta_phi_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_norm";
    TString delta_phi_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_norm";
    TString delta_phi_norm_generic_name = "ak5PF_delta_phi_norm";
    TString delta_phi_norm_unfold_name = "output_true_delta_phi_norm";

    data_file->GetObject(delta_phi_norm_generic_name,delta_phi_norm_data);
    if (delta_phi_norm_data == 0) { cout << delta_phi_norm_generic_name << " not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_norm_unfold_name,delta_phi_norm_p6_z2);
	if (delta_phi_norm_p6_z2 == 0) { cout << delta_phi_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_norm_p6_z2_name,delta_phi_norm_p6_z2);
	if (delta_phi_norm_p6_z2 == 0) { cout << delta_phi_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_norm_unfold_name,delta_phi_norm_p8_4c);
	if (delta_phi_norm_p8_4c == 0) { cout << delta_phi_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_norm_p8_4c_name,delta_phi_norm_p8_4c);
	if (delta_phi_norm_p8_4c == 0) { cout << delta_phi_norm_p8_4c_name << " not found!" << endl; return; }
	}

    
    TH1D *delta_phi_norm;
    delta_phi_norm =  new TH1D("model_unc_delta_phi_norm","Model Uncertainty;#Delta#phi [rad];Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_norm, delta_phi_norm_p6_z2, delta_phi_norm_p8_4c, final_model_unc, 19, detail);
	}
    else
	{
	calc_model_unc(delta_phi_norm, delta_phi_norm_data, delta_phi_norm_p6_z2, delta_phi_norm_p8_4c, final_model_unc, 19, detail);
	}
    delta_phi_norm->SetEntries(delta_phi_norm_data->GetEntries());
    plot_histogram(delta_phi_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_data = 0;
    TH1D *delta_phi_deta1_p6_z2 = 0;
    TH1D *delta_phi_deta1_p8_4c = 0;
    TString delta_phi_deta1_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta1";
    TString delta_phi_deta1_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta1";
    TString delta_phi_deta1_generic_name = "ak5PF_delta_phi_deta1";
    TString delta_phi_deta1_unfold_name = "output_true_delta_phi_deta1";

    data_file->GetObject(delta_phi_deta1_generic_name,delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << delta_phi_deta1_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta1_unfold_name,delta_phi_deta1_p6_z2);
	if (delta_phi_deta1_p6_z2 == 0) { cout << delta_phi_deta1_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta1_p6_z2_name,delta_phi_deta1_p6_z2);
	if (delta_phi_deta1_p6_z2 == 0) { cout << delta_phi_deta1_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta1_unfold_name,delta_phi_deta1_p8_4c);
	if (delta_phi_deta1_p8_4c == 0) { cout << delta_phi_deta1_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta1_p8_4c_name,delta_phi_deta1_p8_4c);
	if (delta_phi_deta1_p8_4c == 0) { cout << delta_phi_deta1_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D("model_unc_delta_phi_deta1","Model Uncertainty for 0.4 > #Delta#eta > 1.5;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta1, delta_phi_deta1_p6_z2, delta_phi_deta1_p8_4c, final_model_unc, 1, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta1, delta_phi_deta1_data, delta_phi_deta1_p6_z2, delta_phi_deta1_p8_4c, final_model_unc, 1, detail);
	}
    delta_phi_deta1->SetEntries(delta_phi_deta1_data->GetEntries());
    plot_histogram(delta_phi_deta1, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta1", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta1 norm distribution
    if (detail) { cout<<"Delta phi deta1 Norm"<<endl; }

    TH1D *delta_phi_deta1_norm_data = 0;
    TH1D *delta_phi_deta1_norm_p6_z2 = 0;
    TH1D *delta_phi_deta1_norm_p8_4c = 0;
    TString delta_phi_deta1_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta1_norm";
    TString delta_phi_deta1_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta1_norm";
    TString delta_phi_deta1_norm_generic_name = "ak5PF_delta_phi_deta1_norm";
    TString delta_phi_deta1_norm_unfold_name = "output_true_delta_phi_deta1_norm";

    data_file->GetObject(delta_phi_deta1_norm_generic_name,delta_phi_deta1_norm_data);
    if (delta_phi_deta1_norm_data == 0) { cout << delta_phi_deta1_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta1_norm_unfold_name,delta_phi_deta1_norm_p6_z2);
	if (delta_phi_deta1_norm_p6_z2 == 0) { cout << delta_phi_deta1_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta1_norm_p6_z2_name,delta_phi_deta1_norm_p6_z2);
	if (delta_phi_deta1_norm_p6_z2 == 0) { cout << delta_phi_deta1_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta1_norm_unfold_name,delta_phi_deta1_norm_p8_4c);
	if (delta_phi_deta1_p8_4c == 0) { cout << delta_phi_deta1_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta1_norm_p8_4c_name,delta_phi_deta1_norm_p8_4c);
	if (delta_phi_deta1_norm_p8_4c == 0) { cout << delta_phi_deta1_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta1_norm;
    delta_phi_deta1_norm =  new TH1D("model_unc_delta_phi_deta1_norm","Model Uncertainty for 0.4 > #Delta#eta > 1.5;#Delta#phi [rad] for 0.4 > |#Delta#eta| > 2.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta1_norm, delta_phi_deta1_norm_p6_z2, delta_phi_deta1_norm_p8_4c, final_model_unc, 20, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta1_norm, delta_phi_deta1_norm_data, delta_phi_deta1_norm_p6_z2, delta_phi_deta1_norm_p8_4c, final_model_unc, 20, detail);
	}

    delta_phi_deta1_norm->SetEntries(delta_phi_deta1_norm_data->GetEntries());
    plot_histogram(delta_phi_deta1_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta1_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_data = 0;
    TH1D *delta_phi_deta2_p6_z2 = 0;
    TH1D *delta_phi_deta2_p8_4c = 0;
    TString delta_phi_deta2_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta2";
    TString delta_phi_deta2_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta2";
    TString delta_phi_deta2_generic_name = "ak5PF_delta_phi_deta2";
    TString delta_phi_deta2_unfold_name = "output_true_delta_phi_deta2";

    data_file->GetObject(delta_phi_deta2_generic_name,delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << delta_phi_deta2_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta2_unfold_name,delta_phi_deta2_p6_z2);
	if (delta_phi_deta2_p6_z2 == 0) { cout << delta_phi_deta2_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta2_p6_z2_name,delta_phi_deta2_p6_z2);
	if (delta_phi_deta2_p6_z2 == 0) { cout << delta_phi_deta2_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta2_unfold_name,delta_phi_deta2_p8_4c);
	if (delta_phi_deta2_p8_4c == 0) { cout << delta_phi_deta2_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta2_p8_4c_name,delta_phi_deta2_p8_4c);
	if (delta_phi_deta2_p8_4c == 0) { cout << delta_phi_deta2_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D("model_unc_delta_phi_deta2","Model Uncertainty for 2.5 > #Delta#eta > 3.5;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta2, delta_phi_deta2_p6_z2, delta_phi_deta2_p8_4c, final_model_unc, 2, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta2, delta_phi_deta2_data, delta_phi_deta2_p6_z2, delta_phi_deta2_p8_4c, final_model_unc, 2, detail);
	}
    delta_phi_deta2->SetEntries(delta_phi_deta2_data->GetEntries());
    plot_histogram(delta_phi_deta2, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta2", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta2 norm distribution
    if (detail) { cout<<"Delta phi deta2 Norm"<<endl; }

    TH1D *delta_phi_deta2_norm_data = 0;
    TH1D *delta_phi_deta2_norm_p6_z2 = 0;
    TH1D *delta_phi_deta2_norm_p8_4c = 0;
    TString delta_phi_deta2_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta2_norm";
    TString delta_phi_deta2_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta2_norm";
    TString delta_phi_deta2_norm_generic_name = "ak5PF_delta_phi_deta2_norm";
    TString delta_phi_deta2_norm_unfold_name = "output_true_delta_phi_deta2_norm";

    data_file->GetObject(delta_phi_deta2_norm_generic_name,delta_phi_deta2_norm_data);
    if (delta_phi_deta2_norm_data == 0) { cout << delta_phi_deta2_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta2_norm_unfold_name,delta_phi_deta2_norm_p6_z2);
	if (delta_phi_deta2_norm_p6_z2 == 0) { cout << delta_phi_deta2_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta2_norm_p6_z2_name,delta_phi_deta2_norm_p6_z2);
	if (delta_phi_deta2_norm_p6_z2 == 0) { cout << delta_phi_deta2_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta2_norm_unfold_name,delta_phi_deta2_norm_p8_4c);
	if (delta_phi_deta2_norm_p8_4c == 0) { cout << delta_phi_deta2_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta2_norm_p8_4c_name,delta_phi_deta2_norm_p8_4c);
	if (delta_phi_deta2_norm_p8_4c == 0) { cout << delta_phi_deta2_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta2_norm;
    delta_phi_deta2_norm =  new TH1D("model_unc_delta_phi_deta2_norm","Model Uncertainty for 2.5 > #Delta#eta > 3.5;#Delta#phi [rad] for 2.5 > #Delta#eta > 3.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta2_norm, delta_phi_deta2_norm_p6_z2, delta_phi_deta2_norm_p8_4c, final_model_unc, 21, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta2_norm, delta_phi_deta2_norm_data, delta_phi_deta2_norm_p6_z2, delta_phi_deta2_norm_p8_4c, final_model_unc, 21, detail);
	}
    delta_phi_deta2_norm->SetEntries(delta_phi_deta2_norm_data->GetEntries());
    plot_histogram(delta_phi_deta2_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta2_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_data = 0;
    TH1D *delta_phi_deta3_p6_z2 = 0;
    TH1D *delta_phi_deta3_p8_4c = 0;
    TString delta_phi_deta3_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta3";
    TString delta_phi_deta3_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta3";
    TString delta_phi_deta3_generic_name = "ak5PF_delta_phi_deta3";
    TString delta_phi_deta3_unfold_name = "output_true_delta_phi_deta3";

    data_file->GetObject(delta_phi_deta3_generic_name,delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << delta_phi_deta3_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta3_unfold_name,delta_phi_deta3_p6_z2);
	if (delta_phi_deta3_p6_z2 == 0) { cout << delta_phi_deta3_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta3_p6_z2_name,delta_phi_deta3_p6_z2);
	if (delta_phi_deta3_p6_z2 == 0) { cout << delta_phi_deta3_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta3_unfold_name,delta_phi_deta3_p8_4c);
	if (delta_phi_deta3_p8_4c == 0) { cout << delta_phi_deta3_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta3_p8_4c_name,delta_phi_deta3_p8_4c);
	if (delta_phi_deta3_p8_4c == 0) { cout << delta_phi_deta3_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D("model_unc_delta_phi_deta3","Model Uncertainty for 3.5 > #Delta#eta > 4.5;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta3, delta_phi_deta3_p6_z2, delta_phi_deta3_p8_4c, final_model_unc, 3, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta3, delta_phi_deta3_data, delta_phi_deta3_p6_z2, delta_phi_deta3_p8_4c, final_model_unc, 3, detail);
	}
    delta_phi_deta3->SetEntries(delta_phi_deta3_data->GetEntries());
    plot_histogram(delta_phi_deta3, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta3", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta3 norm distribution
    if (detail) { cout<<"Delta phi deta3 Norm"<<endl; }

    TH1D *delta_phi_deta3_norm_data = 0;
    TH1D *delta_phi_deta3_norm_p6_z2 = 0;
    TH1D *delta_phi_deta3_norm_p8_4c = 0;
    TString delta_phi_deta3_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta3_norm";
    TString delta_phi_deta3_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta3_norm";
    TString delta_phi_deta3_norm_generic_name = "ak5PF_delta_phi_deta3_norm";
    TString delta_phi_deta3_norm_unfold_name = "output_true_delta_phi_deta3_norm";

    data_file->GetObject(delta_phi_deta3_norm_generic_name,delta_phi_deta3_norm_data);
    if (delta_phi_deta3_norm_data == 0) { cout << delta_phi_deta3_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta3_norm_unfold_name,delta_phi_deta3_norm_p6_z2);
	if (delta_phi_deta3_norm_p6_z2 == 0) { cout << delta_phi_deta3_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta3_norm_p6_z2_name,delta_phi_deta3_norm_p6_z2);
	if (delta_phi_deta3_norm_p6_z2 == 0) { cout << delta_phi_deta3_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta3_norm_unfold_name,delta_phi_deta3_norm_p8_4c);
	if (delta_phi_deta3_norm_p8_4c == 0) { cout << delta_phi_deta3_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta3_norm_p8_4c_name,delta_phi_deta3_norm_p8_4c);
	if (delta_phi_deta3_norm_p8_4c == 0) { cout << delta_phi_deta3_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta3_norm;
    delta_phi_deta3_norm =  new TH1D("model_unc_delta_phi_deta3_norm","Model Uncertainty for 3.5 > #Delta#eta > 4.5;#Delta#phi [rad] for 3.5 > #Delta#eta > 4.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta3_norm, delta_phi_deta3_norm_p6_z2, delta_phi_deta3_norm_p8_4c, final_model_unc, 22, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta3_norm, delta_phi_deta3_norm_data, delta_phi_deta3_norm_p6_z2, delta_phi_deta3_norm_p8_4c, final_model_unc, 22, detail);
	}
    delta_phi_deta3_norm->SetEntries(delta_phi_deta3_norm_data->GetEntries());
    plot_histogram(delta_phi_deta3_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta3_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_data = 0;
    TH1D *delta_phi_deta4_p6_z2 = 0;
    TH1D *delta_phi_deta4_p8_4c = 0;
    TString delta_phi_deta4_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta4";
    TString delta_phi_deta4_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta4";
    TString delta_phi_deta4_generic_name = "ak5PF_delta_phi_deta4";
    TString delta_phi_deta4_unfold_name = "output_true_delta_phi_deta4";

    data_file->GetObject(delta_phi_deta4_generic_name,delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << delta_phi_deta4_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta4_unfold_name,delta_phi_deta4_p6_z2);
	if (delta_phi_deta4_p6_z2 == 0) { cout << delta_phi_deta4_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta4_p6_z2_name,delta_phi_deta4_p6_z2);
	if (delta_phi_deta4_p6_z2 == 0) { cout << delta_phi_deta4_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta4_unfold_name,delta_phi_deta4_p8_4c);
	if (delta_phi_deta4_p8_4c == 0) { cout << delta_phi_deta4_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta4_p8_4c_name,delta_phi_deta4_p8_4c);
	if (delta_phi_deta4_p8_4c == 0) { cout << delta_phi_deta4_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D("model_unc_delta_phi_deta4","Model Uncertainty for 4.5 > #Delta#eta > 7.5;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta4, delta_phi_deta4_p6_z2, delta_phi_deta4_p8_4c, final_model_unc, 4, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta4, delta_phi_deta4_data, delta_phi_deta4_p6_z2, delta_phi_deta4_p8_4c, final_model_unc, 4, detail);
	}
    delta_phi_deta4->SetEntries(delta_phi_deta4_data->GetEntries());
    plot_histogram(delta_phi_deta4, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta4", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta4 norm distribution
    if (detail) { cout<<"Delta phi deta4 Norm"<<endl; }

    TH1D *delta_phi_deta4_norm_data = 0;
    TH1D *delta_phi_deta4_norm_p6_z2 = 0;
    TH1D *delta_phi_deta4_norm_p8_4c = 0;
    TString delta_phi_deta4_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta4_norm";
    TString delta_phi_deta4_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta4_norm";
    TString delta_phi_deta4_norm_generic_name = "ak5PF_delta_phi_deta4_norm";
    TString delta_phi_deta4_norm_unfold_name = "output_true_delta_phi_deta4_norm";

    data_file->GetObject(delta_phi_deta4_norm_generic_name,delta_phi_deta4_norm_data);
    if (delta_phi_deta4_norm_data == 0) { cout << delta_phi_deta4_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta4_norm_unfold_name,delta_phi_deta4_norm_p6_z2);
	if (delta_phi_deta4_norm_p6_z2 == 0) { cout << delta_phi_deta4_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta4_norm_p6_z2_name,delta_phi_deta4_norm_p6_z2);
	if (delta_phi_deta4_norm_p6_z2 == 0) { cout << delta_phi_deta4_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta4_norm_unfold_name,delta_phi_deta4_norm_p8_4c);
	if (delta_phi_deta4_norm_p8_4c == 0) { cout << delta_phi_deta4_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta4_norm_p8_4c_name,delta_phi_deta4_norm_p8_4c);
	if (delta_phi_deta4_norm_p8_4c == 0) { cout << delta_phi_deta4_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta4_norm;
    delta_phi_deta4_norm =  new TH1D("model_unc_delta_phi_deta4_norm","Model Uncertainty for 4.5 > #Delta#eta > 7.5;#Delta#phi [rad] for 4.5 > #Delta#eta > 7.5;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta4_norm, delta_phi_deta4_norm_p6_z2, delta_phi_deta4_norm_p8_4c, final_model_unc, 23, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta4, delta_phi_deta4_norm_data, delta_phi_deta4_p6_z2, delta_phi_deta4_p8_4c, final_model_unc, 23, detail);
	}
    delta_phi_deta4_norm->SetEntries(delta_phi_deta4_data->GetEntries());
    plot_histogram(delta_phi_deta4_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta4_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_data = 0;
    TH1D *delta_phi_gap_p6_z2 = 0;
    TH1D *delta_phi_gap_p8_4c = 0;
    TString delta_phi_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_gap";
    TString delta_phi_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_gap";
    TString delta_phi_gap_generic_name = "ak5PF_delta_phi_gap";
    TString delta_phi_gap_unfold_name = "output_true_delta_phi_gap";

    data_file->GetObject(delta_phi_gap_generic_name,delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << delta_phi_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_gap_unfold_name,delta_phi_gap_p6_z2);
	if (delta_phi_gap_p6_z2 == 0) { cout << delta_phi_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_gap_p6_z2_name,delta_phi_gap_p6_z2);
	if (delta_phi_gap_p6_z2 == 0) { cout << delta_phi_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_gap_unfold_name,delta_phi_gap_p8_4c);
	if (delta_phi_gap_p8_4c == 0) { cout << delta_phi_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_gap_p8_4c_name,delta_phi_gap_p8_4c);
	if (delta_phi_gap_p8_4c == 0) { cout << delta_phi_gap_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("model_unc_delta_phi_gap","Model Uncartainty;|#Delta#phi| [rad];Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_gap, delta_phi_gap_p6_z2, delta_phi_gap_p8_4c, final_model_unc, 5, detail);
	}
    else
	{
	calc_model_unc(delta_phi_gap, delta_phi_gap_data, delta_phi_gap_p6_z2, delta_phi_gap_p8_4c, final_model_unc, 5, detail);
	}
    delta_phi_gap->SetEntries(delta_phi_gap_data->GetEntries());
    plot_histogram(delta_phi_gap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi gap norm distribution
    if (detail) { cout<<"Delta phi gap Norm"<<endl; }

    TH1D *delta_phi_gap_norm_data = 0;
    TH1D *delta_phi_gap_norm_p6_z2 = 0;
    TH1D *delta_phi_gap_norm_p8_4c = 0;
    TString delta_phi_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_gap_norm";
    TString delta_phi_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_gap_norm";
    TString delta_phi_gap_norm_generic_name = "ak5PF_delta_phi_gap_norm";
    TString delta_phi_gap_norm_unfold_name = "output_true_delta_phi_gap_norm";

    data_file->GetObject(delta_phi_gap_norm_generic_name,delta_phi_gap_norm_data);
    if (delta_phi_gap_norm_data == 0) { cout << delta_phi_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_gap_norm_unfold_name,delta_phi_gap_norm_p6_z2);
	if (delta_phi_gap_norm_p6_z2 == 0) { cout << delta_phi_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_gap_norm_p6_z2_name,delta_phi_gap_norm_p6_z2);
	if (delta_phi_gap_norm_p6_z2 == 0) { cout << delta_phi_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_gap_norm_unfold_name,delta_phi_gap_norm_p8_4c);
	if (delta_phi_gap_norm_p8_4c == 0) { cout << delta_phi_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_gap_norm_p8_4c_name,delta_phi_gap_norm_p8_4c);
	if (delta_phi_gap_norm_p8_4c == 0) { cout << delta_phi_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_gap_norm;
    delta_phi_gap_norm =  new TH1D("model_unc_delta_phi_gap_norm","Model Uncartainty;#Delta#phi [rad];Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_gap_norm, delta_phi_gap_norm_p6_z2, delta_phi_gap_norm_p8_4c, final_model_unc, 24, detail);
	}
    else
	{
	calc_model_unc(delta_phi_gap_norm, delta_phi_gap_norm_data, delta_phi_gap_norm_p6_z2, delta_phi_gap_norm_p8_4c, final_model_unc, 24, detail);
	}
    delta_phi_gap_norm->SetEntries(delta_phi_gap_norm_data->GetEntries());
    plot_histogram(delta_phi_gap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }
    
    TH1D *delta_phi_deta1_gap_data = 0;
    TH1D *delta_phi_deta1_gap_p6_z2 = 0;
    TH1D *delta_phi_deta1_gap_p8_4c = 0;
    TString delta_phi_deta1_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_generic_name = "ak5PF_delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_unfold_name = "output_true_delta_phi_deta1_gap";

    data_file->GetObject(delta_phi_deta1_gap_generic_name,delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << delta_phi_deta1_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta1_gap_unfold_name,delta_phi_deta1_gap_p6_z2);
	if (delta_phi_deta1_gap_p6_z2 == 0) { cout << delta_phi_deta1_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta1_gap_p6_z2_name,delta_phi_deta1_gap_p6_z2);
	if (delta_phi_deta1_gap_p6_z2 == 0) { cout << delta_phi_deta1_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta1_gap_unfold_name,delta_phi_deta1_gap_p8_4c);
	if (delta_phi_deta1_gap_p8_4c == 0) { cout << delta_phi_deta1_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta1_gap_p8_4c_name,delta_phi_deta1_gap_p8_4c);
	if (delta_phi_deta1_gap_p8_4c == 0) { cout << delta_phi_deta1_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D("model_unc_delta_phi_deta1_gap","Model Uncertanty for 0.4 > #Delta#eta > 2.5 when requiring a gap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta1_gap, delta_phi_deta1_gap_p6_z2, delta_phi_deta1_gap_p8_4c, final_model_unc, 6, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta1_gap, delta_phi_deta1_gap_data, delta_phi_deta1_gap_p6_z2, delta_phi_deta1_gap_p8_4c, final_model_unc, 6, detail);
	}
    delta_phi_deta1_gap->SetEntries(delta_phi_deta1_gap_data->GetEntries());
    plot_histogram(delta_phi_deta1_gap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta1_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta1 gap norm distribution
    if (detail) { cout<<"Delta phi deta1 gap Norm"<<endl; }
    
    TH1D *delta_phi_deta1_gap_norm_data = 0;
    TH1D *delta_phi_deta1_gap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta1_gap_norm_p8_4c = 0;
    TString delta_phi_deta1_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta1_gap_norm";
    TString delta_phi_deta1_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta1_gap_norm";
    TString delta_phi_deta1_gap_norm_generic_name = "ak5PF_delta_phi_deta1_gap_norm";
    TString delta_phi_deta1_gap_norm_unfold_name = "output_true_delta_phi_deta1_gap_norm";

    data_file->GetObject(delta_phi_deta1_gap_norm_generic_name,delta_phi_deta1_gap_norm_data);
    if (delta_phi_deta1_gap_norm_data == 0) { cout << delta_phi_deta1_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta1_gap_norm_unfold_name,delta_phi_deta1_gap_norm_p6_z2);
	if (delta_phi_deta1_gap_norm_p6_z2 == 0) { cout << delta_phi_deta1_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta1_gap_norm_p6_z2_name,delta_phi_deta1_gap_norm_p6_z2);
	if (delta_phi_deta1_gap_norm_p6_z2 == 0) { cout << delta_phi_deta1_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta1_gap_norm_unfold_name,delta_phi_deta1_gap_norm_p8_4c);
	if (delta_phi_deta1_gap_norm_p8_4c == 0) { cout << delta_phi_deta1_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta1_gap_norm_p8_4c_name,delta_phi_deta1_gap_norm_p8_4c);
	if (delta_phi_deta1_gap_norm_p8_4c == 0) { cout << delta_phi_deta1_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta1_gap_norm;
    delta_phi_deta1_gap_norm =  new TH1D("model_unc_delta_phi_deta1_gap_norm","Model Uncertanty for 0.4 > #Delta#eta > 2.5 when requiring a gap;#Delta#phi [rad] for 0.4 > #Delta#eta > 2.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_p6_z2, delta_phi_deta1_gap_norm_p8_4c, final_model_unc, 25, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta1_gap_norm, delta_phi_deta1_gap_norm_data, delta_phi_deta1_gap_norm_p6_z2, delta_phi_deta1_gap_norm_p8_4c, final_model_unc, 25, detail);
	}
    delta_phi_deta1_gap_norm->SetEntries(delta_phi_deta1_gap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta1_gap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta1_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);

 
//estimate the model uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }
    
    TH1D *delta_phi_deta2_gap_data = 0;
    TH1D *delta_phi_deta2_gap_p6_z2 = 0;
    TH1D *delta_phi_deta2_gap_p8_4c = 0;
    TString delta_phi_deta2_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_generic_name = "ak5PF_delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_unfold_name = "output_true_delta_phi_deta2_gap";

    data_file->GetObject(delta_phi_deta2_gap_generic_name,delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << delta_phi_deta2_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta2_gap_unfold_name,delta_phi_deta2_gap_p6_z2);
	if (delta_phi_deta2_gap_p6_z2 == 0) { cout << delta_phi_deta2_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta2_gap_p6_z2_name,delta_phi_deta2_gap_p6_z2);
	if (delta_phi_deta2_gap_p6_z2 == 0) { cout << delta_phi_deta2_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta2_gap_unfold_name,delta_phi_deta2_gap_p8_4c);
	if (delta_phi_deta2_gap_p8_4c == 0) { cout << delta_phi_deta2_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta2_gap_p8_4c_name,delta_phi_deta2_gap_p8_4c);
	if (delta_phi_deta2_gap_p8_4c == 0) { cout << delta_phi_deta2_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D("model_unc_delta_phi_deta2_gap","Model Uncertanty for 2.5 > #Delta#eta > 3.5 when requiring a gap;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta2_gap, delta_phi_deta2_gap_p6_z2, delta_phi_deta2_gap_p8_4c, final_model_unc, 7, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta2_gap, delta_phi_deta2_gap_data, delta_phi_deta2_gap_p6_z2, delta_phi_deta2_gap_p8_4c, final_model_unc, 7, detail);
	}
    delta_phi_deta2_gap->SetEntries(delta_phi_deta2_gap_data->GetEntries());
    plot_histogram(delta_phi_deta2_gap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta2_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta2 gap norm distribution
    if (detail) { cout<<"Delta phi deta2 gap Norm"<<endl; }
    
    TH1D *delta_phi_deta2_gap_norm_data = 0;
    TH1D *delta_phi_deta2_gap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta2_gap_norm_p8_4c = 0;
    TString delta_phi_deta2_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta2_gap_norm";
    TString delta_phi_deta2_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta2_gap_norm";
    TString delta_phi_deta2_gap_norm_generic_name = "ak5PF_delta_phi_deta2_gap_norm";
    TString delta_phi_deta2_gap_norm_unfold_name = "output_true_delta_phi_deta2_gap_norm";

    data_file->GetObject(delta_phi_deta2_gap_norm_generic_name,delta_phi_deta2_gap_norm_data);
    if (delta_phi_deta2_gap_norm_data == 0) { cout << delta_phi_deta2_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta2_gap_norm_unfold_name,delta_phi_deta2_gap_norm_p6_z2);
	if (delta_phi_deta2_gap_norm_p6_z2 == 0) { cout << delta_phi_deta2_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta2_gap_norm_p6_z2_name,delta_phi_deta2_gap_norm_p6_z2);
	if (delta_phi_deta2_gap_norm_p6_z2 == 0) { cout << delta_phi_deta2_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta2_gap_norm_unfold_name,delta_phi_deta2_gap_norm_p8_4c);
	if (delta_phi_deta2_gap_norm_p8_4c == 0) { cout << delta_phi_deta2_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta2_gap_norm_p8_4c_name,delta_phi_deta2_gap_norm_p8_4c);
	if (delta_phi_deta2_gap_norm_p8_4c == 0) { cout << delta_phi_deta2_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta2_gap_norm;
    delta_phi_deta2_gap_norm =  new TH1D("model_unc_delta_phi_deta2_gap_norm","Model Uncertanty for 2.5 > #Delta#eta > 3.5 when requiring a gap;#Delta#phi [rad] for 2.5 > #Delta#eta > 3.5;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_p6_z2, delta_phi_deta2_gap_norm_p8_4c, final_model_unc, 26, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta2_gap_norm, delta_phi_deta2_gap_norm_data, delta_phi_deta2_gap_norm_p6_z2, delta_phi_deta2_gap_norm_p8_4c, final_model_unc, 26, detail);
	}
    delta_phi_deta2_gap_norm->SetEntries(delta_phi_deta2_gap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta2_gap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta2_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }
    
    TH1D *delta_phi_deta3_gap_data = 0;
    TH1D *delta_phi_deta3_gap_p6_z2 = 0;
    TH1D *delta_phi_deta3_gap_p8_4c = 0;
    TString delta_phi_deta3_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_generic_name = "ak5PF_delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_unfold_name = "output_true_delta_phi_deta3_gap";

    data_file->GetObject(delta_phi_deta3_gap_generic_name,delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << delta_phi_deta3_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta3_gap_unfold_name,delta_phi_deta3_gap_p6_z2);
	if (delta_phi_deta3_gap_p6_z2 == 0) { cout << delta_phi_deta3_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta3_gap_p6_z2_name,delta_phi_deta3_gap_p6_z2);
	if (delta_phi_deta3_gap_p6_z2 == 0) { cout << delta_phi_deta3_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta3_gap_unfold_name,delta_phi_deta3_gap_p8_4c);
	if (delta_phi_deta3_gap_p8_4c == 0) { cout << delta_phi_deta3_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta3_gap_p8_4c_name,delta_phi_deta3_gap_p8_4c);
	if (delta_phi_deta3_gap_p8_4c == 0) { cout << delta_phi_deta3_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D("model_unc_delta_phi_deta3_gap","Model Uncertainty for 3.5 > #Delta#eta > 4.5 when requiring a gap;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta3_gap, delta_phi_deta3_gap_p6_z2, delta_phi_deta3_gap_p8_4c, final_model_unc, 8, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta3_gap, delta_phi_deta3_gap_data, delta_phi_deta3_gap_p6_z2, delta_phi_deta3_gap_p8_4c, final_model_unc, 8, detail);
	}
    delta_phi_deta3_gap->SetEntries(delta_phi_deta3_gap_data->GetEntries());
    plot_histogram(delta_phi_deta3_gap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta3_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta3 gap norm distribution
    if (detail) { cout<<"Delta phi deta3 gap Norm"<<endl; }
    
    TH1D *delta_phi_deta3_gap_norm_data = 0;
    TH1D *delta_phi_deta3_gap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta3_gap_norm_p8_4c = 0;
    TString delta_phi_deta3_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta3_gap_norm";
    TString delta_phi_deta3_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta3_gap_norm";
    TString delta_phi_deta3_gap_norm_generic_name = "ak5PF_delta_phi_deta3_gap_norm";
    TString delta_phi_deta3_gap_norm_unfold_name = "output_true_delta_phi_deta3_gap_norm";

    data_file->GetObject(delta_phi_deta3_gap_norm_generic_name,delta_phi_deta3_gap_norm_data);
    if (delta_phi_deta3_gap_norm_data == 0) { cout << delta_phi_deta3_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta3_gap_norm_unfold_name,delta_phi_deta3_gap_norm_p6_z2);
	if (delta_phi_deta3_gap_norm_p6_z2 == 0) { cout << delta_phi_deta3_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta3_gap_p6_z2_name,delta_phi_deta3_gap_p6_z2);
	if (delta_phi_deta3_gap_p6_z2 == 0) { cout << delta_phi_deta3_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta3_gap_norm_unfold_name,delta_phi_deta3_gap_norm_p8_4c);
	if (delta_phi_deta3_gap_norm_p8_4c == 0) { cout << delta_phi_deta3_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta3_gap_norm_p8_4c_name,delta_phi_deta3_gap_norm_p8_4c);
	if (delta_phi_deta3_gap_norm_p8_4c == 0) { cout << delta_phi_deta3_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta3_gap_norm;
    delta_phi_deta3_gap_norm =  new TH1D("model_unc_delta_phi_deta3_gap_norm","Model Uncertainty for 3.5 > #Delta#eta > 4.5 when requiring a gap;#Delta#phi [rad] for 3.5 > #Delta#eta > 4.5;Model Uncertainty", dphi_nbins, dphi_bins);
    
    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_p6_z2, delta_phi_deta3_gap_norm_p8_4c, final_model_unc, 27, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta3_gap_norm, delta_phi_deta3_gap_norm_data, delta_phi_deta3_gap_norm_p6_z2, delta_phi_deta3_gap_norm_p8_4c, final_model_unc, 27, detail);
	}
    delta_phi_deta3_gap_norm->SetEntries(delta_phi_deta3_gap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta3_gap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta3_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }
    
    TH1D *delta_phi_deta4_gap_data = 0;
    TH1D *delta_phi_deta4_gap_p6_z2 = 0;
    TH1D *delta_phi_deta4_gap_p8_4c = 0;
    TString delta_phi_deta4_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_generic_name = "ak5PF_delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_unfold_name = "output_true_delta_phi_deta4_gap";

    data_file->GetObject(delta_phi_deta4_gap_generic_name,delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << delta_phi_deta4_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta4_gap_unfold_name,delta_phi_deta4_gap_p6_z2);
	if (delta_phi_deta4_gap_p6_z2 == 0) { cout << delta_phi_deta4_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta4_gap_p6_z2_name,delta_phi_deta4_gap_p6_z2);
	if (delta_phi_deta4_gap_p6_z2 == 0) { cout << delta_phi_deta4_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta4_gap_unfold_name,delta_phi_deta4_gap_p8_4c);
	if (delta_phi_deta4_gap_p8_4c == 0) { cout << delta_phi_deta4_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta4_gap_p8_4c_name,delta_phi_deta4_gap_p8_4c);
	if (delta_phi_deta4_gap_p8_4c == 0) { cout << delta_phi_deta4_gap_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D("model_unc_delta_phi_deta4_gap","Model Uncertainty for 4.5 > #Delta#eta > 7.5 when requiring a gap;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta4_gap, delta_phi_deta4_gap_p6_z2, delta_phi_deta4_gap_p8_4c, final_model_unc, 9, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta4_gap, delta_phi_deta4_gap_data, delta_phi_deta4_gap_p6_z2, delta_phi_deta4_gap_p8_4c, final_model_unc, 9, detail);
	}
    delta_phi_deta4_gap->SetEntries(delta_phi_deta4_gap_data->GetEntries());
    plot_histogram(delta_phi_deta4_gap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta4_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta4 gap Norm distribution
    if (detail) { cout<<"Delta phi deta4 gap Norm"<<endl; }
    
    TH1D *delta_phi_deta4_gap_norm_data = 0;
    TH1D *delta_phi_deta4_gap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta4_gap_norm_p8_4c = 0;
    TString delta_phi_deta4_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta4_gap_norm";
    TString delta_phi_deta4_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta4_gap_norm";
    TString delta_phi_deta4_gap_norm_generic_name = "ak5PF_delta_phi_deta4_gap_norm";
    TString delta_phi_deta4_gap_norm_unfold_name = "output_true_delta_phi_deta4_gap_norm";

    data_file->GetObject(delta_phi_deta4_gap_norm_generic_name,delta_phi_deta4_gap_norm_data);
    if (delta_phi_deta4_gap_norm_data == 0) { cout << delta_phi_deta4_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta4_gap_norm_unfold_name,delta_phi_deta4_gap_norm_p6_z2);
	if (delta_phi_deta4_gap_norm_p6_z2 == 0) { cout << delta_phi_deta4_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta4_gap_norm_p6_z2_name,delta_phi_deta4_gap_norm_p6_z2);
	if (delta_phi_deta4_gap_norm_p6_z2 == 0) { cout << delta_phi_deta4_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta4_gap_norm_unfold_name,delta_phi_deta4_gap_norm_p8_4c);
	if (delta_phi_deta4_gap_norm_p8_4c == 0) { cout << delta_phi_deta4_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta4_gap_norm_p8_4c_name,delta_phi_deta4_gap_norm_p8_4c);
	if (delta_phi_deta4_gap_norm_p8_4c == 0) { cout << delta_phi_deta4_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta4_gap_norm;
    delta_phi_deta4_gap_norm =  new TH1D("model_unc_delta_phi_deta4_gap_norm","Model Uncertainty for 4.5 > #Delta#eta > 7.5 when requiring a gap;#Delta#phi [rad] for 4.5 > #Delta#eta > 7.5;Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_p6_z2, delta_phi_deta4_gap_norm_p8_4c, final_model_unc, 28, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta4_gap_norm, delta_phi_deta4_gap_norm_data, delta_phi_deta4_gap_norm_p6_z2, delta_phi_deta4_gap_norm_p8_4c, final_model_unc, 28, detail);
	}
    delta_phi_deta4_gap_norm->SetEntries(delta_phi_deta4_gap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta4_gap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta4_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }
    
    TH1D *delta_phi_nogap_data = 0;
    TH1D *delta_phi_nogap_p6_z2 = 0;
    TH1D *delta_phi_nogap_p8_4c = 0;
    TString delta_phi_nogap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_nogap";
    TString delta_phi_nogap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_nogap";
    TString delta_phi_nogap_generic_name = "ak5PF_delta_phi_nogap";
    TString delta_phi_nogap_unfold_name = "output_true_delta_phi_nogap";

    data_file->GetObject(delta_phi_nogap_generic_name,delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << delta_phi_nogap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_nogap_unfold_name,delta_phi_nogap_p6_z2);
	if (delta_phi_nogap_p6_z2 == 0) { cout << delta_phi_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_nogap_p6_z2_name,delta_phi_nogap_p6_z2);
	if (delta_phi_nogap_p6_z2 == 0) { cout << delta_phi_nogap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_nogap_unfold_name,delta_phi_nogap_p8_4c);
	if (delta_phi_nogap_p8_4c == 0) { cout << delta_phi_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_nogap_p8_4c_name,delta_phi_nogap_p8_4c);
	if (delta_phi_nogap_p8_4c == 0) { cout << delta_phi_nogap_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("model_unc_delta_phi_nogap","Model Uncertainty;|#Delta#phi| [rad];Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_nogap, delta_phi_nogap_p6_z2, delta_phi_nogap_p8_4c, final_model_unc, 10, detail);
	}
    else
	{
	calc_model_unc(delta_phi_nogap, delta_phi_nogap_data, delta_phi_nogap_p6_z2, delta_phi_nogap_p8_4c, final_model_unc, 10, detail);
	}
    delta_phi_nogap->SetEntries(delta_phi_nogap_data->GetEntries());
    plot_histogram(delta_phi_nogap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_nogap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);



//estimate the model uncertainty for delta phi nogap norm distribution
    if (detail) { cout<<"Delta phi nogap Norm"<<endl; }
    
    TH1D *delta_phi_nogap_norm_data = 0;
    TH1D *delta_phi_nogap_norm_p6_z2 = 0;
    TH1D *delta_phi_nogap_norm_p8_4c = 0;
    TString delta_phi_nogap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_nogap_norm";
    TString delta_phi_nogap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_nogap_norm";
    TString delta_phi_nogap_norm_generic_name = "ak5PF_delta_phi_nogap_norm";
    TString delta_phi_nogap_norm_unfold_name = "output_true_delta_phi_nogap_norm";

    data_file->GetObject(delta_phi_nogap_norm_generic_name,delta_phi_nogap_norm_data);
    if (delta_phi_nogap_norm_data == 0) { cout << delta_phi_nogap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_nogap_norm_unfold_name,delta_phi_nogap_norm_p6_z2);
	if (delta_phi_nogap_norm_p6_z2 == 0) { cout << delta_phi_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_nogap_norm_p6_z2_name,delta_phi_nogap_norm_p6_z2);
	if (delta_phi_nogap_norm_p6_z2 == 0) { cout << delta_phi_nogap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_nogap_norm_unfold_name,delta_phi_nogap_norm_p8_4c);
	if (delta_phi_nogap_norm_p8_4c == 0) { cout << delta_phi_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_nogap_norm_p8_4c_name,delta_phi_nogap_norm_p8_4c);
	if (delta_phi_nogap_norm_p8_4c == 0) { cout << delta_phi_nogap_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_nogap_norm;
    delta_phi_nogap_norm =  new TH1D("model_unc_delta_phi_nogap_norm","Model Uncertainty;#Delta#phi [rad];Model Uncertainty", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_nogap_norm, delta_phi_nogap_norm_p6_z2, delta_phi_nogap_norm_p8_4c, final_model_unc, 29, detail);
	}
    else
	{
	calc_model_unc(delta_phi_nogap_norm, delta_phi_nogap_norm_data, delta_phi_nogap_norm_p6_z2, delta_phi_nogap_norm_p8_4c, final_model_unc, 29, detail);
	}
    delta_phi_nogap_norm->SetEntries(delta_phi_nogap_norm_data->GetEntries());
    plot_histogram(delta_phi_nogap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_nogap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_data = 0;
    TH1D *delta_phi_deta1_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta1_nogap_p8_4c = 0;
    TString delta_phi_deta1_nogap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_generic_name = "ak5PF_delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_unfold_name = "output_true_delta_phi_deta1_nogap";

    data_file->GetObject(delta_phi_deta1_nogap_generic_name,delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << delta_phi_deta1_nogap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta1_nogap_unfold_name,delta_phi_deta1_nogap_p6_z2);
	if (delta_phi_deta1_nogap_p6_z2 == 0) { cout << delta_phi_deta1_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta1_nogap_p6_z2_name,delta_phi_deta1_nogap_p6_z2);
	if (delta_phi_deta1_nogap_p6_z2 == 0) { cout << delta_phi_deta1_nogap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta1_nogap_unfold_name,delta_phi_deta1_nogap_p8_4c);
	if (delta_phi_deta1_nogap_p8_4c == 0) { cout << delta_phi_deta1_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta1_nogap_p8_4c_name,delta_phi_deta1_nogap_p8_4c);
	if (delta_phi_deta1_nogap_p8_4c == 0) { cout << delta_phi_deta1_nogap_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D("model_unc_delta_phi_deta1_nogap","Model Uncertainty for 0.4 > #Delta#eta > 2.5 when requiring a gap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta1_nogap, delta_phi_deta1_nogap_p6_z2, delta_phi_deta1_nogap_p8_4c, final_model_unc, 11, detail);
	}
	else
	{
	calc_model_unc(delta_phi_deta1_nogap, delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_p6_z2, delta_phi_deta1_nogap_p8_4c, final_model_unc, 11, detail);
	}
    delta_phi_deta1_nogap->SetEntries(delta_phi_deta1_nogap_data->GetEntries());
    plot_histogram(delta_phi_deta1_nogap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta1_nogap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta1 nogap norm distribution
    if (detail) { cout<<"Delta phi deta1 nogap Norm"<<endl; }

    TH1D *delta_phi_deta1_nogap_norm_data = 0;
    TH1D *delta_phi_deta1_nogap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta1_nogap_norm_p8_4c = 0;
    TString delta_phi_deta1_nogap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta1_nogap_norm";
    TString delta_phi_deta1_nogap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta1_nogap_norm";
    TString delta_phi_deta1_nogap_norm_generic_name = "ak5PF_delta_phi_deta1_nogap_norm";
    TString delta_phi_deta1_nogap_norm_unfold_name = "output_true_delta_phi_deta1_nogap_norm";

    data_file->GetObject(delta_phi_deta1_nogap_norm_generic_name,delta_phi_deta1_nogap_norm_data);
    if (delta_phi_deta1_nogap_norm_data == 0) { cout << delta_phi_deta1_nogap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta1_nogap_norm_unfold_name,delta_phi_deta1_nogap_norm_p6_z2);
	if (delta_phi_deta1_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta1_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta1_nogap_norm_p6_z2_name,delta_phi_deta1_nogap_norm_p6_z2);
	if (delta_phi_deta1_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta1_nogap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta1_nogap_norm_unfold_name,delta_phi_deta1_nogap_norm_p8_4c);
	if (delta_phi_deta1_nogap_norm_p8_4c == 0) { cout << delta_phi_deta1_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta1_nogap_norm_p8_4c_name,delta_phi_deta1_nogap_norm_p8_4c);
	if (delta_phi_deta1_nogap_norm_p8_4c == 0) { cout << delta_phi_deta1_nogap_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta1_nogap_norm;
    delta_phi_deta1_nogap_norm =  new TH1D("model_unc_delta_phi_deta1_nogap_norm","Model Uncertainty for 0.4 > #Delta#eta > 2.5 when requiring a gap;#Delta#phi [rad] for 0.4 > #Delta#eta > 2.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_p6_z2, delta_phi_deta1_nogap_norm_p8_4c, final_model_unc, 30, detail);
	}
	else
	{
	calc_model_unc(delta_phi_deta1_nogap_norm, delta_phi_deta1_nogap_norm_data, delta_phi_deta1_nogap_norm_p6_z2, delta_phi_deta1_nogap_norm_p8_4c, final_model_unc, 30, detail);
	}
    delta_phi_deta1_nogap_norm->SetEntries(delta_phi_deta1_nogap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta1_nogap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta1_nogap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_data = 0;
    TH1D *delta_phi_deta2_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta2_nogap_p8_4c = 0;
    TString delta_phi_deta2_nogap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_generic_name = "ak5PF_delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_unfold_name = "output_true_delta_phi_deta2_nogap";

    data_file->GetObject(delta_phi_deta2_nogap_generic_name,delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << delta_phi_deta2_nogap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta2_nogap_unfold_name,delta_phi_deta2_nogap_p6_z2);
	if (delta_phi_deta2_nogap_p6_z2 == 0) { cout << delta_phi_deta2_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta2_nogap_p6_z2_name,delta_phi_deta2_nogap_p6_z2);
	if (delta_phi_deta2_nogap_p6_z2 == 0) { cout << delta_phi_deta2_nogap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta2_nogap_unfold_name,delta_phi_deta2_nogap_p8_4c);
	if (delta_phi_deta2_nogap_p8_4c == 0) { cout << delta_phi_deta2_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta2_nogap_p8_4c_name,delta_phi_deta2_nogap_p8_4c);
	if (delta_phi_deta2_nogap_p8_4c == 0) { cout << delta_phi_deta2_nogap_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D("model_unc_delta_phi_deta2_nogap","Model Uncertainty for 2.5 > #Delta#eta > 3.5 when requiring a gap;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta2_nogap, delta_phi_deta2_nogap_p6_z2, delta_phi_deta2_nogap_p8_4c, final_model_unc, 12, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta2_nogap, delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_p6_z2, delta_phi_deta2_nogap_p8_4c, final_model_unc, 12, detail);
	}
    delta_phi_deta2_nogap->SetEntries(delta_phi_deta2_nogap_data->GetEntries());
    plot_histogram(delta_phi_deta2_nogap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta2_nogap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);   


//estimate the model uncertainty for delta phi deta2 nogap norm distribution
    if (detail) { cout<<"Delta phi deta2 nogap Norm"<<endl; }

    TH1D *delta_phi_deta2_nogap_norm_data = 0;
    TH1D *delta_phi_deta2_nogap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta2_nogap_norm_p8_4c = 0;
    TString delta_phi_deta2_nogap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta2_nogap_norm";
    TString delta_phi_deta2_nogap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta2_nogap_norm";
    TString delta_phi_deta2_nogap_norm_generic_name = "ak5PF_delta_phi_deta2_nogap_norm";
    TString delta_phi_deta2_nogap_norm_unfold_name = "output_true_delta_phi_deta2_nogap_norm";

    data_file->GetObject(delta_phi_deta2_nogap_norm_generic_name,delta_phi_deta2_nogap_norm_data);
    if (delta_phi_deta2_nogap_norm_data == 0) { cout << delta_phi_deta2_nogap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta2_nogap_norm_unfold_name,delta_phi_deta2_nogap_norm_p6_z2);
	if (delta_phi_deta2_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta2_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta2_nogap_norm_p6_z2_name,delta_phi_deta2_nogap_norm_p6_z2);
	if (delta_phi_deta2_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta2_nogap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta2_nogap_norm_unfold_name,delta_phi_deta2_nogap_norm_p8_4c);
	if (delta_phi_deta2_nogap_norm_p8_4c == 0) { cout << delta_phi_deta2_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta2_nogap_norm_p8_4c_name,delta_phi_deta2_nogap_norm_p8_4c);
	if (delta_phi_deta2_nogap_norm_p8_4c == 0) { cout << delta_phi_deta2_nogap_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta2_nogap_norm;
    delta_phi_deta2_nogap_norm =  new TH1D("model_unc_delta_phi_deta2_nogap_norm","Model Uncertainty for 2.5 > #Delta#eta > 3.5 when requiring a gap;#Delta#phi [rad] for 2.5 > #Delta#eta > 3.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_p6_z2, delta_phi_deta2_nogap_norm_p8_4c, final_model_unc, 31, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta2_nogap_norm, delta_phi_deta2_nogap_norm_data, delta_phi_deta2_nogap_norm_p6_z2, delta_phi_deta2_nogap_norm_p8_4c, final_model_unc, 31, detail);
	}
    delta_phi_deta2_nogap_norm->SetEntries(delta_phi_deta2_nogap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta2_nogap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta2_nogap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_data = 0;
    TH1D *delta_phi_deta3_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta3_nogap_p8_4c = 0;
    TString delta_phi_deta3_nogap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_generic_name = "ak5PF_delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_unfold_name = "output_true_delta_phi_deta3_nogap";

    data_file->GetObject(delta_phi_deta3_nogap_generic_name,delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << delta_phi_deta3_nogap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta3_nogap_unfold_name,delta_phi_deta3_nogap_p6_z2);
	if (delta_phi_deta3_nogap_p6_z2 == 0) { cout << delta_phi_deta3_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta3_nogap_p6_z2_name,delta_phi_deta3_nogap_p6_z2);
	if (delta_phi_deta3_nogap_p6_z2 == 0) { cout << delta_phi_deta3_nogap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta3_nogap_unfold_name,delta_phi_deta3_nogap_p8_4c);
	if (delta_phi_deta3_nogap_p8_4c == 0) { cout << delta_phi_deta3_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta3_nogap_p8_4c_name,delta_phi_deta3_nogap_p8_4c);
	if (delta_phi_deta3_nogap_p8_4c == 0) { cout << delta_phi_deta3_nogap_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D("model_unc_delta_phi_deta3_nogap","Model Uncertainty for 3.5 > #Delta#eta > 4.5 when requiring a gap;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta3_nogap, delta_phi_deta3_nogap_p6_z2, delta_phi_deta3_nogap_p8_4c, final_model_unc, 13, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta3_nogap, delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_p6_z2, delta_phi_deta3_nogap_p8_4c, final_model_unc, 13, detail);
	}
    delta_phi_deta3_nogap->SetEntries(delta_phi_deta3_nogap_data->GetEntries());
    plot_histogram(delta_phi_deta3_nogap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta3_nogap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta3 nogap norm distribution
    if (detail) { cout<<"Delta phi deta3 nogap norm"<<endl; }

    TH1D *delta_phi_deta3_nogap_norm_data = 0;
    TH1D *delta_phi_deta3_nogap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta3_nogap_norm_p8_4c = 0;
    TString delta_phi_deta3_nogap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta3_nogap_norm";
    TString delta_phi_deta3_nogap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta3_nogap_norm";
    TString delta_phi_deta3_nogap_norm_generic_name = "ak5PF_delta_phi_deta3_nogap_norm";
    TString delta_phi_deta3_nogap_norm_unfold_name = "output_true_delta_phi_deta3_nogap_norm";

    data_file->GetObject(delta_phi_deta3_nogap_norm_generic_name,delta_phi_deta3_nogap_norm_data);
    if (delta_phi_deta3_nogap_norm_data == 0) { cout << delta_phi_deta3_nogap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta3_nogap_norm_unfold_name,delta_phi_deta3_nogap_norm_p6_z2);
	if (delta_phi_deta3_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta3_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta3_nogap_norm_p6_z2_name,delta_phi_deta3_nogap_norm_p6_z2);
	if (delta_phi_deta3_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta3_nogap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta3_nogap_norm_unfold_name,delta_phi_deta3_nogap_norm_p8_4c);
	if (delta_phi_deta3_nogap_norm_p8_4c == 0) { cout << delta_phi_deta3_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta3_nogap_norm_p8_4c_name,delta_phi_deta3_nogap_norm_p8_4c);
	if (delta_phi_deta3_nogap_norm_p8_4c == 0) { cout << delta_phi_deta3_nogap_norm_p8_4c_name << " not found!" << endl; return; }
	}

    TH1D *delta_phi_deta3_nogap_norm;
    delta_phi_deta3_nogap_norm =  new TH1D("model_unc_delta_phi_deta3_nogap_norm","Model Uncertainty for 3.5 > #Delta#eta > 4.5 when requiring a gap;#Delta#phi [rad] for 3.5 > #Delta#eta > 4.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_p6_z2, delta_phi_deta3_nogap_norm_p8_4c, final_model_unc, 32, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta3_nogap_norm, delta_phi_deta3_nogap_norm_data, delta_phi_deta3_nogap_norm_p6_z2, delta_phi_deta3_nogap_norm_p8_4c, final_model_unc, 32, detail);
	}
    delta_phi_deta3_nogap_norm->SetEntries(delta_phi_deta3_nogap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta3_nogap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta3_nogap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_data = 0;
    TH1D *delta_phi_deta4_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta4_nogap_p8_4c = 0;
    TString delta_phi_deta4_nogap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_generic_name = "ak5PF_delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_unfold_name = "output_true_delta_phi_deta4_nogap";

    data_file->GetObject(delta_phi_deta4_nogap_generic_name,delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << delta_phi_deta4_nogap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta4_nogap_unfold_name,delta_phi_deta4_nogap_p6_z2);
	if (delta_phi_deta4_nogap_p6_z2 == 0) { cout << delta_phi_deta4_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta4_nogap_p6_z2_name,delta_phi_deta4_nogap_p6_z2);
	if (delta_phi_deta4_nogap_p6_z2 == 0) { cout << delta_phi_deta4_nogap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta4_nogap_unfold_name,delta_phi_deta4_nogap_p8_4c);
	if (delta_phi_deta4_nogap_p8_4c == 0) { cout << delta_phi_deta4_nogap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta4_nogap_p8_4c_name,delta_phi_deta4_nogap_p8_4c);
	if (delta_phi_deta4_nogap_p8_4c == 0) { cout << delta_phi_deta4_nogap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D("model_unc_delta_phi_deta4_nogap","Model Uncertainty for 4.5 > #Delta#eta > 7.5 when requiring a gap;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta4_nogap, delta_phi_deta4_nogap_p6_z2, delta_phi_deta4_nogap_p8_4c, final_model_unc, 14, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta4_nogap, delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_p6_z2, delta_phi_deta4_nogap_p8_4c, final_model_unc, 14, detail);
	}
    delta_phi_deta4_nogap->SetEntries(delta_phi_deta4_nogap_data->GetEntries());
    plot_histogram(delta_phi_deta4_nogap, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta4_nogap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta phi deta4 nogap norm distribution
    if (detail) { cout<<"Delta phi deta4 nogap Norm"<<endl; }

    TH1D *delta_phi_deta4_nogap_norm_data = 0;
    TH1D *delta_phi_deta4_nogap_norm_p6_z2 = 0;
    TH1D *delta_phi_deta4_nogap_norm_p8_4c = 0;
    TString delta_phi_deta4_nogap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_phi_deta4_nogap_norm";
    TString delta_phi_deta4_nogap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_phi_deta4_nogap_norm";
    TString delta_phi_deta4_nogap_norm_generic_name = "ak5PF_delta_phi_deta4_nogap_norm";
    TString delta_phi_deta4_nogap_norm_unfold_name = "output_true_delta_phi_deta4_nogap_norm";

    data_file->GetObject(delta_phi_deta4_nogap_norm_generic_name,delta_phi_deta4_nogap_norm_data);
    if (delta_phi_deta4_nogap_norm_data == 0) { cout << delta_phi_deta4_nogap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_phi_deta4_nogap_norm_unfold_name,delta_phi_deta4_nogap_norm_p6_z2);
	if (delta_phi_deta4_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta4_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_phi_deta4_nogap_norm_p6_z2_name,delta_phi_deta4_nogap_norm_p6_z2);
	if (delta_phi_deta4_nogap_norm_p6_z2 == 0) { cout << delta_phi_deta4_nogap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_phi_deta4_nogap_norm_unfold_name,delta_phi_deta4_nogap_norm_p8_4c);
	if (delta_phi_deta4_nogap_norm_p8_4c == 0) { cout << delta_phi_deta4_nogap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_phi_deta4_nogap_norm_p8_4c_name,delta_phi_deta4_nogap_norm_p8_4c);
	if (delta_phi_deta4_nogap_norm_p8_4c == 0) { cout << delta_phi_deta4_nogap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_phi_deta4_nogap_norm;
    delta_phi_deta4_nogap_norm =  new TH1D("model_unc_delta_phi_deta4_nogap_norm","Model Uncertainty for 4.5 > #Delta#eta > 7.5 when requiring a gap;#Delta#phi [rad] for 4.5 > #Delta#eta > 7.5;Model Dependency", dphi_nbins, dphi_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_p6_z2, delta_phi_deta4_nogap_norm_p8_4c, final_model_unc, 33, detail);
	}
    else
	{
	calc_model_unc(delta_phi_deta4_nogap_norm, delta_phi_deta4_nogap_norm_data, delta_phi_deta4_nogap_norm_p6_z2, delta_phi_deta4_nogap_norm_p8_4c, final_model_unc, 33, detail);
	}
    delta_phi_deta4_nogap_norm->SetEntries(delta_phi_deta4_nogap_norm_data->GetEntries());
    plot_histogram(delta_phi_deta4_nogap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_phi_deta4_nogap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for leading pT inside gap distribution
    if (detail) { cout<<"Leading pT inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_data = 0;
    TH1D *leading_pt_inside_gap_p6_z2 = 0;
    TH1D *leading_pt_inside_gap_p8_4c = 0;
    TString leading_pt_inside_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_leading_pt_inside_gap";
    TString leading_pt_inside_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_leading_pt_inside_gap";
    TString leading_pt_inside_gap_generic_name = "ak5PF_leading_pt_inside_gap";
    TString leading_pt_inside_gap_unfold_name = "output_true_leading_pt_inside_gap";

    data_file->GetObject(leading_pt_inside_gap_generic_name,leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << leading_pt_inside_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(leading_pt_inside_gap_unfold_name,leading_pt_inside_gap_p6_z2);
	if (leading_pt_inside_gap_p6_z2 == 0) { cout << leading_pt_inside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(leading_pt_inside_gap_p6_z2_name,leading_pt_inside_gap_p6_z2);
	if (leading_pt_inside_gap_p6_z2 == 0) { cout << leading_pt_inside_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(leading_pt_inside_gap_unfold_name,leading_pt_inside_gap_p8_4c);
	if (leading_pt_inside_gap_p8_4c == 0) { cout << leading_pt_inside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(leading_pt_inside_gap_p8_4c_name,leading_pt_inside_gap_p8_4c);
	if (leading_pt_inside_gap_p8_4c == 0) { cout << leading_pt_inside_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D("model_unc_leading_pt_inside_gap","Model Uncertainty;p_{T} [GeV];Model Uncertainty", in_nbins, in_bins);

    if (with_unfold)
	{
	calc_model_unc2(leading_pt_inside_gap, leading_pt_inside_gap_p6_z2, leading_pt_inside_gap_p8_4c, final_model_unc, 15, detail);
	}
    else
	{
	calc_model_unc(leading_pt_inside_gap, leading_pt_inside_gap_data, leading_pt_inside_gap_p6_z2, leading_pt_inside_gap_p8_4c, final_model_unc, 15, detail);
	}
    leading_pt_inside_gap->SetEntries(leading_pt_inside_gap_data->GetEntries());
    plot_histogram(leading_pt_inside_gap, output_path_plots, "model_uncertainty_" + label + "_leading_pt_inside_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for leading pT inside gap norm distribution
    if (detail) { cout<<"Leading pT inside gap Norm"<<endl; }

    TH1D *leading_pt_inside_gap_norm_data = 0;
    TH1D *leading_pt_inside_gap_norm_p6_z2 = 0;
    TH1D *leading_pt_inside_gap_norm_p8_4c = 0;
    TString leading_pt_inside_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_leading_pt_inside_gap_norm";
    TString leading_pt_inside_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_leading_pt_inside_gap_norm";
    TString leading_pt_inside_gap_norm_generic_name = "ak5PF_leading_pt_inside_gap_norm";
    TString leading_pt_inside_gap_norm_unfold_name = "output_true_leading_pt_inside_gap_norm";

    data_file->GetObject(leading_pt_inside_gap_norm_generic_name,leading_pt_inside_gap_norm_data);
    if (leading_pt_inside_gap_norm_data == 0) { cout << leading_pt_inside_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(leading_pt_inside_gap_norm_unfold_name,leading_pt_inside_gap_norm_p6_z2);
	if (leading_pt_inside_gap_norm_p6_z2 == 0) { cout << leading_pt_inside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(leading_pt_inside_gap_norm_p6_z2_name,leading_pt_inside_gap_norm_p6_z2);
	if (leading_pt_inside_gap_norm_p6_z2 == 0) { cout << leading_pt_inside_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(leading_pt_inside_gap_norm_unfold_name,leading_pt_inside_gap_norm_p8_4c);
	if (leading_pt_inside_gap_norm_p8_4c == 0) { cout << leading_pt_inside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(leading_pt_inside_gap_norm_p8_4c_name,leading_pt_inside_gap_norm_p8_4c);
	if (leading_pt_inside_gap_norm_p8_4c == 0) { cout << leading_pt_inside_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *leading_pt_inside_gap_norm;
    leading_pt_inside_gap_norm =  new TH1D("model_unc_leading_pt_inside_gap_norm","Model Uncertainty;p_{T}^{inside} [GeV];Model Uncertainty", in_nbins, in_bins);

    if (with_unfold)
	{
	calc_model_unc2(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_p6_z2, leading_pt_inside_gap_norm_p8_4c, final_model_unc, 34, detail);
	}
    else
	{
	calc_model_unc(leading_pt_inside_gap_norm, leading_pt_inside_gap_norm_data, leading_pt_inside_gap_norm_p6_z2, leading_pt_inside_gap_norm_p8_4c, final_model_unc, 34, detail);
	}
    leading_pt_inside_gap_norm->SetEntries(leading_pt_inside_gap_norm_data->GetEntries());
    plot_histogram(leading_pt_inside_gap_norm, output_path_plots, "model_uncertainty_" + label + "_leading_pt_inside_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading Eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_data = 0;
    TH1D *leading_eta_star_inside_gap_p6_z2 = 0;
    TH1D *leading_eta_star_inside_gap_p8_4c = 0;
    TString leading_eta_star_inside_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_generic_name = "ak5PF_leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_unfold_name = "output_true_leading_eta_star_inside_gap";

    data_file->GetObject(leading_eta_star_inside_gap_generic_name,leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << leading_eta_star_inside_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(leading_eta_star_inside_gap_unfold_name,leading_eta_star_inside_gap_p6_z2);
	if (leading_eta_star_inside_gap_p6_z2 == 0) { cout << leading_eta_star_inside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(leading_eta_star_inside_gap_p6_z2_name,leading_eta_star_inside_gap_p6_z2);
	if (leading_eta_star_inside_gap_p6_z2 == 0) { cout << leading_eta_star_inside_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(leading_eta_star_inside_gap_unfold_name,leading_eta_star_inside_gap_p8_4c);
	if (leading_eta_star_inside_gap_p8_4c == 0) { cout << leading_eta_star_inside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(leading_eta_star_inside_gap_p8_4c_name,leading_eta_star_inside_gap_p8_4c);
	if (leading_eta_star_inside_gap_p8_4c == 0) { cout << leading_eta_star_inside_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D("model_unc_leading_eta_star_inside_gap","Model Uncertainty;#eta*;Model Uncertainty", etastar_nbins, etastar_bins);

    if (with_unfold)
	{
	calc_model_unc2(leading_eta_star_inside_gap, leading_eta_star_inside_gap_p6_z2, leading_eta_star_inside_gap_p8_4c, final_model_unc, 16, detail);
	}
    else
	{
	calc_model_unc(leading_eta_star_inside_gap, leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_p6_z2, leading_eta_star_inside_gap_p8_4c, final_model_unc, 16, detail);
	}
    leading_eta_star_inside_gap->SetEntries(leading_eta_star_inside_gap_data->GetEntries());
    plot_histogram(leading_eta_star_inside_gap, output_path_plots, "model_uncertainty_" + label + "_leading_eta_star_inside_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);



//estimate the model uncertainty for leading eta star inside gap norm distribution
    if (detail) { cout<<"Leading Eta star inside gap Norm"<<endl; }

    TH1D *leading_eta_star_inside_gap_norm_data = 0;
    TH1D *leading_eta_star_inside_gap_norm_p6_z2 = 0;
    TH1D *leading_eta_star_inside_gap_norm_p8_4c = 0;
    TString leading_eta_star_inside_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_leading_eta_star_inside_gap_norm";
    TString leading_eta_star_inside_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_leading_eta_star_inside_gap_norm";
    TString leading_eta_star_inside_gap_norm_generic_name = "ak5PF_leading_eta_star_inside_gap_norm";
    TString leading_eta_star_inside_gap_norm_unfold_name = "output_true_leading_eta_star_inside_gap_norm";

    data_file->GetObject(leading_eta_star_inside_gap_norm_generic_name,leading_eta_star_inside_gap_norm_data);
    if (leading_eta_star_inside_gap_norm_data == 0) { cout << leading_eta_star_inside_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(leading_eta_star_inside_gap_norm_unfold_name,leading_eta_star_inside_gap_norm_p6_z2);
	if (leading_eta_star_inside_gap_norm_p6_z2 == 0) { cout << leading_eta_star_inside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(leading_eta_star_inside_gap_norm_p6_z2_name,leading_eta_star_inside_gap_norm_p6_z2);
	if (leading_eta_star_inside_gap_norm_p6_z2 == 0) { cout << leading_eta_star_inside_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(leading_eta_star_inside_gap_norm_unfold_name,leading_eta_star_inside_gap_norm_p8_4c);
	if (leading_eta_star_inside_gap_norm_p8_4c == 0) { cout << leading_eta_star_inside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(leading_eta_star_inside_gap_norm_p8_4c_name,leading_eta_star_inside_gap_norm_p8_4c);
	if (leading_eta_star_inside_gap_norm_p8_4c == 0) { cout << leading_eta_star_inside_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *leading_eta_star_inside_gap_norm;
    leading_eta_star_inside_gap_norm =  new TH1D("model_unc_leading_eta_star_inside_gap_norm","Model Uncertainty;#eta*;Model Uncertainty", etastar_nbins, etastar_bins);

    if (with_unfold)
	{
	calc_model_unc2(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_p6_z2, leading_eta_star_inside_gap_norm_p8_4c, final_model_unc, 35, detail);
	}
    else
	{
	calc_model_unc(leading_eta_star_inside_gap_norm, leading_eta_star_inside_gap_norm_data, leading_eta_star_inside_gap_norm_p6_z2, leading_eta_star_inside_gap_norm_p8_4c, final_model_unc, 35, detail);
	}
    leading_eta_star_inside_gap_norm->SetEntries(leading_eta_star_inside_gap_norm_data->GetEntries());
    plot_histogram(leading_eta_star_inside_gap_norm, output_path_plots, "model_uncertainty_" + label + "_leading_eta_star_inside_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Leading Eta star inside gap"<<endl; }

    TH1D *delta_eta_outside_gap_data = 0;
    TH1D *delta_eta_outside_gap_p6_z2 = 0;
    TH1D *delta_eta_outside_gap_p8_4c = 0;
    TString delta_eta_outside_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_eta_outside_gap";
    TString delta_eta_outside_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_eta_outside_gap";
    TString delta_eta_outside_gap_generic_name = "ak5PF_delta_eta_outside_gap";
    TString delta_eta_outside_gap_unfold_name = "output_true_delta_eta_outside_gap";

    data_file->GetObject(delta_eta_outside_gap_generic_name,delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << delta_eta_outside_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_eta_outside_gap_unfold_name,delta_eta_outside_gap_p6_z2);
	if (delta_eta_outside_gap_p6_z2 == 0) { cout << delta_eta_outside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_eta_outside_gap_p6_z2_name,delta_eta_outside_gap_p6_z2);
	if (delta_eta_outside_gap_p6_z2 == 0) { cout << delta_eta_outside_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_eta_outside_gap_unfold_name,delta_eta_outside_gap_p8_4c);
	if (delta_eta_outside_gap_p8_4c == 0) { cout << delta_eta_outside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_eta_outside_gap_p8_4c_name,delta_eta_outside_gap_p8_4c);
	if (delta_eta_outside_gap_p8_4c == 0) { cout << delta_eta_outside_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D("model_unc_delta_eta_outside_gap","Model Uncertainty;#Delta#eta;Model Uncertainty", deta_out_nbins, deta_out_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_eta_outside_gap, delta_eta_outside_gap_p6_z2, delta_eta_outside_gap_p8_4c, final_model_unc, 17, detail);
	}
    else
	{
	calc_model_unc(delta_eta_outside_gap, delta_eta_outside_gap_data, delta_eta_outside_gap_p6_z2, delta_eta_outside_gap_p8_4c, final_model_unc, 17, detail);
	}
    delta_eta_outside_gap->SetEntries(delta_eta_outside_gap_data->GetEntries());
    plot_histogram(delta_eta_outside_gap, output_path_plots, "model_uncertainty_" + label + "_delta_eta_outside_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for delta eta outside gap norm distribution
    if (detail) { cout<<"Leading Eta star inside gap Norm"<<endl; }

    TH1D *delta_eta_outside_gap_norm_data = 0;
    TH1D *delta_eta_outside_gap_norm_p6_z2 = 0;
    TH1D *delta_eta_outside_gap_norm_p8_4c = 0;
    TString delta_eta_outside_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_delta_eta_outside_gap_norm";
    TString delta_eta_outside_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_delta_eta_outside_gap_norm";
    TString delta_eta_outside_gap_norm_generic_name = "ak5PF_delta_eta_outside_gap_norm";
    TString delta_eta_outside_gap_norm_unfold_name = "output_true_delta_eta_outside_gap_norm";

    data_file->GetObject(delta_eta_outside_gap_norm_generic_name,delta_eta_outside_gap_norm_data);
    if (delta_eta_outside_gap_norm_data == 0) { cout << delta_eta_outside_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(delta_eta_outside_gap_norm_unfold_name,delta_eta_outside_gap_norm_p6_z2);
	if (delta_eta_outside_gap_norm_p6_z2 == 0) { cout << delta_eta_outside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(delta_eta_outside_gap_norm_p6_z2_name,delta_eta_outside_gap_norm_p6_z2);
	if (delta_eta_outside_gap_norm_p6_z2 == 0) { cout << delta_eta_outside_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(delta_eta_outside_gap_norm_unfold_name,delta_eta_outside_gap_norm_p8_4c);
	if (delta_eta_outside_gap_norm_p8_4c == 0) { cout << delta_eta_outside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(delta_eta_outside_gap_norm_p8_4c_name,delta_eta_outside_gap_norm_p8_4c);
	if (delta_eta_outside_gap_norm_p8_4c == 0) { cout << delta_eta_outside_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *delta_eta_outside_gap_norm;
    delta_eta_outside_gap_norm =  new TH1D("model_unc_delta_eta_outside_gap_norm","Model Uncertainty;#Delta#eta;Model Uncertainty", deta_out_nbins, deta_out_bins);

    if (with_unfold)
	{
	calc_model_unc2(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_p6_z2, delta_eta_outside_gap_norm_p8_4c, final_model_unc, 36, detail);
	}
    else
	{
	calc_model_unc(delta_eta_outside_gap_norm, delta_eta_outside_gap_norm_data, delta_eta_outside_gap_norm_p6_z2, delta_eta_outside_gap_norm_p8_4c, final_model_unc, 36, detail);
	}
    delta_eta_outside_gap_norm->SetEntries(delta_eta_outside_gap_norm_data->GetEntries());
    plot_histogram(delta_eta_outside_gap_norm, output_path_plots, "model_uncertainty_" + label + "_delta_eta_outside_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for deta eta outside gap distribution
    if (detail) { cout<<"Leading pT outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_data = 0;
    TH1D *leading_pt_outside_gap_p6_z2 = 0;
    TH1D *leading_pt_outside_gap_p8_4c = 0;
    TString leading_pt_outside_gap_p6_z2_name = "detector_pythia6_z2_" + label + "_leading_pt_outside_gap";
    TString leading_pt_outside_gap_p8_4c_name = "detector_pythia8_4c_" + label + "_leading_pt_outside_gap";
    TString leading_pt_outside_gap_generic_name = "ak5PF_leading_pt_outside_gap";
    TString leading_pt_outside_gap_unfold_name = "output_true_leading_pt_outside_gap";

    data_file->GetObject(leading_pt_outside_gap_generic_name,leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << leading_pt_outside_gap_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(leading_pt_outside_gap_unfold_name,leading_pt_outside_gap_p6_z2);
	if (leading_pt_outside_gap_p6_z2 == 0) { cout << leading_pt_outside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(leading_pt_outside_gap_p6_z2_name,leading_pt_outside_gap_p6_z2);
	if (leading_pt_outside_gap_p6_z2 == 0) { cout << leading_pt_outside_gap_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(leading_pt_outside_gap_unfold_name,leading_pt_outside_gap_p8_4c);
	if (leading_pt_outside_gap_p8_4c == 0) { cout << leading_pt_outside_gap_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(leading_pt_outside_gap_p8_4c_name,leading_pt_outside_gap_p8_4c);
	if (leading_pt_outside_gap_p8_4c == 0) { cout << leading_pt_outside_gap_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D("model_unc_leading_pt_outside_gap","Model Uncertainty;p_{T} [GeV];Model Uncertainty", out_nbins, out_bins);

    if (with_unfold)
	{
	calc_model_unc2(leading_pt_outside_gap, leading_pt_outside_gap_p6_z2, leading_pt_outside_gap_p8_4c, final_model_unc, 18, detail);
	}
    else
	{
	calc_model_unc(leading_pt_outside_gap, leading_pt_outside_gap_data, leading_pt_outside_gap_p6_z2, leading_pt_outside_gap_p8_4c, final_model_unc, 18, detail);
	}
    leading_pt_outside_gap->SetEntries(leading_pt_outside_gap_data->GetEntries());
    plot_histogram(leading_pt_outside_gap, output_path_plots, "model_uncertainty_" + label + "_leading_pt_outside_gap", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//estimate the model uncertainty for deta eta outside gap norm distribution
    if (detail) { cout<<"Leading pT outside gap Norm"<<endl; }

    TH1D *leading_pt_outside_gap_norm_data = 0;
    TH1D *leading_pt_outside_gap_norm_p6_z2 = 0;
    TH1D *leading_pt_outside_gap_norm_p8_4c = 0;
    TString leading_pt_outside_gap_norm_p6_z2_name = "detector_pythia6_z2_" + label + "_leading_pt_outside_gap_norm";
    TString leading_pt_outside_gap_norm_p8_4c_name = "detector_pythia8_4c_" + label + "_leading_pt_outside_gap_norm";
    TString leading_pt_outside_gap_norm_generic_name = "ak5PF_leading_pt_outside_gap_norm";
    TString leading_pt_outside_gap_norm_unfold_name = "output_true_leading_pt_outside_gap_norm";

    data_file->GetObject(leading_pt_outside_gap_norm_generic_name,leading_pt_outside_gap_norm_data);
    if (leading_pt_outside_gap_norm_data == 0) { cout << leading_pt_outside_gap_norm_generic_name << "not found!" << endl; return; }
    if (with_unfold)
	{
	pythia6_file->GetObject(leading_pt_outside_gap_norm_unfold_name,leading_pt_outside_gap_norm_p6_z2);
	if (leading_pt_outside_gap_norm_p6_z2 == 0) { cout << leading_pt_outside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia6_file->GetObject(leading_pt_outside_gap_norm_p6_z2_name,leading_pt_outside_gap_norm_p6_z2);
	if (leading_pt_outside_gap_norm_p6_z2 == 0) { cout << leading_pt_outside_gap_norm_p6_z2_name << " not found!" << endl; return; }
	}

    if (with_unfold)
	{
	pythia8_file->GetObject(leading_pt_outside_gap_norm_unfold_name,leading_pt_outside_gap_norm_p8_4c);
	if (leading_pt_outside_gap_norm_p8_4c == 0) { cout << leading_pt_outside_gap_norm_unfold_name << " not found!" << endl; return; }
	}
    else
	{
	pythia8_file->GetObject(leading_pt_outside_gap_norm_p8_4c_name,leading_pt_outside_gap_norm_p8_4c);
	if (leading_pt_outside_gap_norm_p8_4c == 0) { cout << leading_pt_outside_gap_norm_p8_4c_name << " not found!" << endl; return; }
	}
    
    TH1D *leading_pt_outside_gap_norm;
    leading_pt_outside_gap_norm =  new TH1D("model_unc_leading_pt_outside_gap_norm","Model Uncertainty;p_{T}^{outside} [GeV];Model Uncertainty", out_nbins, out_bins);

    if (with_unfold)
	{
	calc_model_unc2(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_p6_z2, leading_pt_outside_gap_norm_p8_4c, final_model_unc, 37, detail);
	}
    else
	{
	calc_model_unc(leading_pt_outside_gap_norm, leading_pt_outside_gap_norm_data, leading_pt_outside_gap_norm_p6_z2, leading_pt_outside_gap_norm_p8_4c, final_model_unc, 37, detail);
	}
    leading_pt_outside_gap_norm->SetEntries(leading_pt_outside_gap_norm_data->GetEntries());
    plot_histogram(leading_pt_outside_gap_norm, output_path_plots, "model_uncertainty_" + label + "_leading_pt_outside_gap_norm", "Model Dependecy from Pythia8-4C and Pythia6-Z2", "top_left", true);


//output the error variation
    if (detail) { cout<<"Display the final model uncertainties..."<<endl; }
    if (disp_uncertainty) { show_model_uncertainties(final_model_unc); }

//Opening the output root file
    if (detail) { cout<<"Creating " << out_model_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( out_model_uncertainty.c_str() , "RECREATE");

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
    leading_pt_outside_gap->Write();
    delta_eta_outside_gap->Write();
    leading_eta_star_inside_gap->Write();

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
    leading_pt_outside_gap_norm->Write();
    delta_eta_outside_gap_norm->Write();
    leading_eta_star_inside_gap_norm->Write();
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data_file->Close();
    pythia6_file->Close();
    pythia8_file->Close();
    data_output->Close();

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
