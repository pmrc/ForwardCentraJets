// Pedro Cipriano, Jul 2012
// DESY, CMS
// Last Update: 30 Jul 2012
//
//

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TF1.h>

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "common_methods.h"

void fit_resolution(TH1D *hist = 0, double range_min = 0.0, double range_max = 0.0, double *result = 0, string path_plots = "../output/", string name = "test", TString label = "test", bool detail = false)
{
    //check if histogram is valid
    if (hist == 0) { cout << "The input histogram is invalid!"<<endl; return; }

    if (detail) { cout << "Fiting and Ploting " << name << endl; }
    //declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    //gPad->SetLogy();
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    //drawning and fitting the function
    format_histogram(hist, 2, 2);
    hist->Draw("e1");
    hist->GetXaxis()->SetRangeUser(range_min,range_max);
    hist->Fit("gaus", "", "", range_min, range_max);
    TF1 *tfit = 0;
    tfit = hist->GetFunction("gaus");

    //saving the results
    result[0] = tfit->GetParameter(1);
    result[1] = tfit->GetParError(1);
    result[2] = tfit->GetParameter(2);
    result[3] = tfit->GetParError(2);

    //assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position("bottom_middle", 2, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hist,label,"l");
    leg01->AddEntry(tfit,"fit","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    //print output to file
    print_plots(c01, path_plots, name);
}

void fit_and_plot(TH1D *p6_jetmettau_hist = 0, TH1D *p6_jetmet_hist = 0, TH1D *p8_jetmettau_hist = 0, TH1D *p8_jetmet_hist = 0, double range_min = 0.0, double range_max = 0.0, string path_plots = "../output/", string name = "test", TString label = "test", TString variable = "var", double *result = 0, int index = 0, bool detail = false)
{
//checking if histograms are valid
if (p6_jetmettau_hist == 0) { cout << "The input p6_jetmettau histogram is invalid!" << endl; return; }
if (p6_jetmet_hist == 0) { cout << "The input p6_jetmet histogram is invalid!" << endl; return; }
if (p8_jetmettau_hist == 0) { cout << "The input p8_jetmettau histogram is invalid!" << endl; return; }
if (p8_jetmet_hist == 0) { cout << "The input p8_jetmet histogram is invalid!" << endl; return; }

//ploting the four histograms
p6_jetmettau_hist->SetMinimum(p6_jetmettau_hist->GetMinimum()*0.01);
plot_4histograms(p6_jetmettau_hist, "Herwig pp - EEC3 with JetMETTau_2010A", p6_jetmet_hist, "Herwig pp - EEC3 with JetMET_2010A", p8_jetmettau_hist, "Pythia 8 - 4C with JetMETTau_2010A", p8_jetmet_hist, "Pythia 8 - 4C with JetMET_2010A", path_plots, "comparison_" + name, "bottom_right", true, detail);


//declaring variables
double result1[4] = {0.0, 0.0, 0.0, 0.0};
double result2[4] = {0.0, 0.0, 0.0, 0.0};
double result3[4] = {0.0, 0.0, 0.0, 0.0};
double result4[4] = {0.0, 0.0, 0.0, 0.0};

fit_resolution(p6_jetmettau_hist, range_min, range_max, result1, path_plots, "fit_hw_eec3_jetmettau_" + name, label, detail);
fit_resolution(p6_jetmet_hist, range_min, range_max, result2, path_plots, "fit_hw_eec3_jetmet_" + name, label, detail);
fit_resolution(p8_jetmettau_hist, range_min, range_max, result3, path_plots, "fit_p8_4c_jetmettau_" + name, label, detail);
fit_resolution(p8_jetmet_hist, range_min, range_max, result4, path_plots, "fit_p8_4c_jetmet_" + name, label, detail);


    //declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.20);
    gPad->SetFrameBorderMode(0);

    TH1D *bias;
    bias =  new TH1D("deviation","Bias "+variable+";;Bias "+variable, 4,-0,4);
    bias->Fill("Herwig pp - EEC3 with JetMETTau_2010A",0);
    bias->Fill("Herwig pp - EEC3 with JetMET_2010A",0);
    bias->Fill("Pythia 8 - 4C with JetMETTau_2010A",0);
    bias->Fill("Pythia 8 - 4C with JetMET_2010A",0);
    
    bias->SetBinContent(1,result1[0]);
    bias->SetBinError(1,result1[1]);
    bias->SetBinContent(2,result2[0]);
    bias->SetBinError(2,result2[1]);
    bias->SetBinContent(3,result3[0]);
    bias->SetBinError(3,result3[1]);
    bias->SetBinContent(4,result4[0]);
    bias->SetBinError(4,result4[1]);

    format_histogram(bias,2,1);
    bias->Draw("e1");

    //print output to file
    print_plots(c01, path_plots, "bias_"+name);

    //declare and configure the canvas
    TCanvas *c02 = new TCanvas("c02","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.20);
    gPad->SetFrameBorderMode(0);

    TH1D *spread;
    spread =  new TH1D("spread","Spread "+variable+";;Spread "+variable, 4,-0,4);
    spread->Fill("Herwig pp - EEC3 with JetMETTau_2010A",0);
    spread->Fill("Herwig pp - EEC3 with JetMET_2010A",0);
    spread->Fill("Pythia 8 - 4C with JetMETTau_2010A",0);
    spread->Fill("Pythia 8 - 4C with JetMET_2010A",0);

    spread->SetBinContent(1,result1[2]);
    spread->SetBinError(1,result1[3]);
    spread->SetBinContent(2,result2[2]);
    spread->SetBinError(2,result2[3]);
    spread->SetBinContent(3,result3[2]);
    spread->SetBinError(3,result3[3]);
    spread->SetBinContent(4,result4[2]);
    spread->SetBinError(4,result4[3]);

    format_histogram(spread,4,2);
    spread->Draw("e1");

    //print output to file
    print_plots(c02, path_plots, "spread_"+name);

    //declare and configure the canvas
    TCanvas *c03 = new TCanvas("c03","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.20);
    gPad->SetFrameBorderMode(0);

    TH1D *resolution;
    resolution =  new TH1D("resolution","Resolution "+variable+";;Resolution "+variable, 4,-0,4);
    resolution->Fill("Herwig pp - EEC3 with JetMETTau_2010A",0);
    resolution->Fill("Herwig pp - EEC3 with JetMET_2010A",0);
    resolution->Fill("Pythia 8 - 4C with JetMETTau_2010A",0);
    resolution->Fill("Pythia 8 - 4C with JetMET_2010A",0);

    resolution->SetBinContent(1,result1[0]);
    resolution->SetBinError(1,result1[2]);
    resolution->SetBinContent(2,result2[0]);
    resolution->SetBinError(2,result2[2]);
    resolution->SetBinContent(3,result3[0]);
    resolution->SetBinError(3,result3[2]);
    resolution->SetBinContent(4,result4[0]);
    resolution->SetBinError(4,result4[2]);

    format_histogram(resolution,9,3);
    resolution->Draw("e1");

    //print output to file
    print_plots(c03, path_plots, "resolution_"+name);

    //average results
    for (int i = 0; i < 4; i++)
	{
	result[i+4*index] = (result1[i] + result2[i] + result3[i] + result4[i])/4.0;
	}

    //delete histograms to avoid memory leak
    delete(bias);
    delete(spread);
    delete(resolution);
}


void show_results(double *results)
{

    cout<<" "<<endl;
    cout<<setiosflags(ios::fixed)<<setprecision(3)<<setfill('0'); //setprecision(3);
    cout<<"Merging data errors"<<endl;
    cout<<"Observable                Bias                Spread"<<endl;
    cout<<"Central pT                "<<setw(6)<<results[0]/10.0<<"+/-"<<setw(6)<<results[1]/10.0<<"     "<<setw(6)<<results[2]/10.0<<"+/-"<<setw(6)<<results[3]/10.0<<endl;
    cout<<"Central Eta               "<<setw(6)<<results[4]<<"+/-"<<setw(6)<<results[5]<<"     "<<setw(6)<<results[6]<<"+/-"<<setw(6)<<results[7]<<endl;
    cout<<"Central Phi               "<<setw(6)<<results[8]<<"+/-"<<setw(6)<<results[9]<<"     "<<setw(6)<<results[10]<<"+/-"<<setw(6)<<results[11]<<endl;
    cout<<"Forward pT                "<<setw(6)<<results[12]/10.0<<"+/-"<<setw(6)<<results[13]/10.0<<"     "<<setw(6)<<results[14]/10.0<<"+/-"<<setw(6)<<results[15]/10.0<<endl;
    cout<<"Forward Eta               "<<setw(6)<<results[16]<<"+/-"<<setw(6)<<results[17]<<"     "<<setw(6)<<results[18]<<"+/-"<<setw(6)<<results[19]<<endl;
    cout<<"Forward Phi               "<<setw(6)<<results[20]<<"+/-"<<setw(6)<<results[21]<<"     "<<setw(6)<<results[22]<<"+/-"<<setw(6)<<results[23]<<endl;
    cout<<"Delta Phi                 "<<setw(6)<<results[24]<<"+/-"<<setw(6)<<results[25]<<"     "<<setw(6)<<results[26]<<"+/-"<<setw(6)<<results[27]<<endl;
    cout<<"Delta Eta                 "<<setw(6)<<results[28]<<"+/-"<<setw(6)<<results[29]<<"     "<<setw(6)<<results[30]<<"+/-"<<setw(6)<<results[31]<<endl;
    cout<<"Delta Phi Deta1           "<<setw(6)<<results[32]<<"+/-"<<setw(6)<<results[33]<<"     "<<setw(6)<<results[34]<<"+/-"<<setw(6)<<results[35]<<endl;
    cout<<"Delta Phi Deta2           "<<setw(6)<<results[36]<<"+/-"<<setw(6)<<results[37]<<"     "<<setw(6)<<results[38]<<"+/-"<<setw(6)<<results[39]<<endl;
    cout<<"Delta Phi Deta3           "<<setw(6)<<results[40]<<"+/-"<<setw(6)<<results[41]<<"     "<<setw(6)<<results[42]<<"+/-"<<setw(6)<<results[43]<<endl;
    cout<<"Delta Phi Deta4           "<<setw(6)<<results[44]<<"+/-"<<setw(6)<<results[45]<<"     "<<setw(6)<<results[46]<<"+/-"<<setw(6)<<results[47]<<endl;
    cout<<"Delta Phi Gap             "<<setw(6)<<results[48]<<"+/-"<<setw(6)<<results[49]<<"     "<<setw(6)<<results[50]<<"+/-"<<setw(6)<<results[51]<<endl;
    cout<<"Delta Eta Gap             "<<setw(6)<<results[52]<<"+/-"<<setw(6)<<results[53]<<"     "<<setw(6)<<results[54]<<"+/-"<<setw(6)<<results[55]<<endl;
    cout<<"Delta Phi Deta1 Gap       "<<setw(6)<<results[56]<<"+/-"<<setw(6)<<results[57]<<"     "<<setw(6)<<results[58]<<"+/-"<<setw(6)<<results[59]<<endl;
    cout<<"Delta Phi Deta2 Gap       "<<setw(6)<<results[60]<<"+/-"<<setw(6)<<results[61]<<"     "<<setw(6)<<results[62]<<"+/-"<<setw(6)<<results[63]<<endl;
    cout<<"Delta Phi Deta3 Gap       "<<setw(6)<<results[64]<<"+/-"<<setw(6)<<results[65]<<"     "<<setw(6)<<results[66]<<"+/-"<<setw(6)<<results[67]<<endl;
    cout<<"Delta Phi Deta4 Gap       "<<setw(6)<<results[68]<<"+/-"<<setw(6)<<results[69]<<"     "<<setw(6)<<results[70]<<"+/-"<<setw(6)<<results[71]<<endl;
    cout<<"Delta Phi Nogap           "<<setw(6)<<results[72]<<"+/-"<<setw(6)<<results[73]<<"     "<<setw(6)<<results[74]<<"+/-"<<setw(6)<<results[75]<<endl;
    cout<<"Delta Eta Nogap           "<<setw(6)<<results[76]<<"+/-"<<setw(6)<<results[77]<<"     "<<setw(6)<<results[78]<<"+/-"<<setw(6)<<results[79]<<endl;
    cout<<"Delta Phi Deta1 Nogap     "<<setw(6)<<results[80]<<"+/-"<<setw(6)<<results[81]<<"     "<<setw(6)<<results[82]<<"+/-"<<setw(6)<<results[83]<<endl;
    cout<<"Delta Phi Deta2 Nogap     "<<setw(6)<<results[84]<<"+/-"<<setw(6)<<results[85]<<"     "<<setw(6)<<results[86]<<"+/-"<<setw(6)<<results[87]<<endl;
    cout<<"Delta Phi Deta3 Nogap     "<<setw(6)<<results[88]<<"+/-"<<setw(6)<<results[89]<<"     "<<setw(6)<<results[90]<<"+/-"<<setw(6)<<results[91]<<endl;
    cout<<"Delta Phi Deta4 Nogap     "<<setw(6)<<results[92]<<"+/-"<<setw(6)<<results[93]<<"     "<<setw(6)<<results[94]<<"+/-"<<setw(6)<<results[95]<<endl;
    cout<<"Leading pT Inside         "<<setw(6)<<results[96]/10.0<<"+/-"<<setw(6)<<results[97]/10.0<<"     "<<setw(6)<<results[98]/10.0<<"+/-"<<setw(6)<<results[99]/10.0<<endl;
    cout<<"Leading Eta Inside        "<<setw(6)<<results[100]<<"+/-"<<setw(6)<<results[101]<<"     "<<setw(6)<<results[102]<<"+/-"<<setw(6)<<results[103]<<endl;
    cout<<"Leading Phi Inside        "<<setw(6)<<results[104]<<"+/-"<<setw(6)<<results[105]<<"     "<<setw(6)<<results[106]<<"+/-"<<setw(6)<<results[107]<<endl;
    cout<<"Leading Eta Star Inside   "<<setw(6)<<results[108]<<"+/-"<<setw(6)<<results[109]<<"     "<<setw(6)<<results[110]<<"+/-"<<setw(6)<<results[111]<<endl;
    cout<<"Leading pT Outside        "<<setw(6)<<results[112]/10.0<<"+/-"<<setw(6)<<results[113]/10.0<<"     "<<setw(6)<<results[114]/10.0<<"+/-"<<setw(6)<<results[115]/10.0<<endl;
    cout<<"Leading Eta Outside       "<<setw(6)<<results[116]<<"+/-"<<setw(6)<<results[117]<<"     "<<setw(6)<<results[118]<<"+/-"<<setw(6)<<results[119]<<endl;
    cout<<"Leading Phi Outside       "<<setw(6)<<results[120]<<"+/-"<<setw(6)<<results[121]<<"     "<<setw(6)<<results[122]<<"+/-"<<setw(6)<<results[123]<<endl;
    cout<<"Delta Eta Outside         "<<setw(6)<<results[124]<<"+/-"<<setw(6)<<results[125]<<"     "<<setw(6)<<results[126]<<"+/-"<<setw(6)<<results[127]<<endl;
}


void plot_resolution(string p6_jetmettau, string p6_jetmet, string p8_jetmettau, string p8_jetmet, string path_plots = "../output/", string plot_prefix = "test_", bool show_result = false, bool detail = false)
{

//output the configuration
   if (detail) { cout<<"Plot Resolution Configuration"<<endl; }
   if (detail) { cout<<"Input Herwig pp - EEC3 with JetMETTau_2010A weights : "<<p6_jetmettau<<endl; }
   if (detail) { cout<<"Input Herwig pp - EEC3 with JetMET_2010A weights :    "<<p6_jetmet<<endl; }
   if (detail) { cout<<"Input Pythia 8 - 4C with JetMETTau_2010A weights :  "<<p8_jetmettau<<endl; }
   if (detail) { cout<<"Input Pythia 8 - 4C with JetMET_2010A weights :     "<<p8_jetmet<<endl; }
   if (detail) { cout<<"Path Plots :                                        "<<path_plots<<endl; }
   if (detail) { cout<<"Plot Prefix :                                       "<<plot_prefix<<endl; }
   if (detail) { cout<<"Show Result :                                       "<<show_result<<endl; }
   if (detail) { cout<<"Detail level :                                      "<<detail<<endl; }

   //opening the files
   if (detail) { cout << "Opening files..." << endl; }
   TFile *p6_z2_jetmettau = new TFile( p6_jetmettau.c_str() );
   TFile *p6_z2_jetmet = new TFile( p6_jetmet.c_str() );
   TFile *p8_4c_jetmettau = new TFile( p8_jetmettau.c_str() );
   TFile *p8_4c_jetmet = new TFile( p8_jetmet.c_str() );
   if (detail) { cout << "All files opened sucessfully!" << endl; }

   //inicializar as variaveis
   int dist = 32;
   double resultado[4*dist];

   for (int i = 0; i < 4*dist; i++)
	{
	resultado[i] = 0.0;
	}


   //central jet pt
   if (detail) { cout << "Central Jet pT..." << endl; }
   TH1D *leading_central_pt_p6_jetmettau = 0;
   TH1D *leading_central_pt_p6_jetmet = 0;
   TH1D *leading_central_pt_p8_jetmettau = 0;
   TH1D *leading_central_pt_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_central_pt",leading_central_pt_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_central_pt",leading_central_pt_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_central_pt",leading_central_pt_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_central_pt",leading_central_pt_p8_jetmet);


   if (leading_central_pt_p6_jetmettau == 0) { cout << "leading_central_pt_p6_jetmettau not found!" << endl; return; }
   if (leading_central_pt_p6_jetmet == 0) { cout << "leading_central_pt_p6_jetmet not found!" << endl; return; }
   if (leading_central_pt_p8_jetmettau == 0) { cout << "leading_central_pt_p8_jetmettau not found!" << endl; return; }
   if (leading_central_pt_p8_jetmet == 0) { cout << "leading_central_pt_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_central_pt_p6_jetmettau, leading_central_pt_p6_jetmet, leading_central_pt_p8_jetmettau, leading_central_pt_p8_jetmet, -0.1, 0.1, path_plots, plot_prefix+"leading_central_pt", "hadron-detector/hadron", "p_{T}^{central}", resultado, 0, detail);


   //central jet eta
   if (detail) { cout << "Central Jet Eta..." << endl; }
   TH1D *leading_central_eta_p6_jetmettau = 0;
   TH1D *leading_central_eta_p6_jetmet = 0;
   TH1D *leading_central_eta_p8_jetmettau = 0;
   TH1D *leading_central_eta_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_central_eta",leading_central_eta_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_central_eta",leading_central_eta_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_central_eta",leading_central_eta_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_central_eta",leading_central_eta_p8_jetmet);


   if (leading_central_eta_p6_jetmettau == 0) { cout << "leading_central_eta_p6_jetmettau not found!" << endl; return; }
   if (leading_central_eta_p6_jetmet == 0) { cout << "leading_central_eta_p6_jetmet not found!" << endl; return; }
   if (leading_central_eta_p8_jetmettau == 0) { cout << "leading_central_eta_p8_jetmettau not found!" << endl; return; }
   if (leading_central_eta_p8_jetmet == 0) { cout << "leading_central_eta_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_central_eta_p6_jetmettau, leading_central_eta_p6_jetmet, leading_central_eta_p8_jetmettau, leading_central_eta_p8_jetmet, -0.018, 0.018, path_plots, plot_prefix+"leading_central_eta", "hadron-detector", "#eta^{central}", resultado, 1, detail);


   //central jet phi
   if (detail) { cout << "Central Jet Phi..." << endl; }
   TH1D *leading_central_phi_p6_jetmettau = 0;
   TH1D *leading_central_phi_p6_jetmet = 0;
   TH1D *leading_central_phi_p8_jetmettau = 0;
   TH1D *leading_central_phi_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_central_phi",leading_central_phi_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_central_phi",leading_central_phi_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_central_phi",leading_central_phi_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_central_phi",leading_central_phi_p8_jetmet);


   if (leading_central_phi_p6_jetmettau == 0) { cout << "leading_central_phi_p6_jetmettau not found!" << endl; return; }
   if (leading_central_phi_p6_jetmet == 0) { cout << "leading_central_phi_p6_jetmet not found!" << endl; return; }
   if (leading_central_phi_p8_jetmettau == 0) { cout << "leading_central_phi_p8_jetmettau not found!" << endl; return; }
   if (leading_central_phi_p8_jetmet == 0) { cout << "leading_central_phi_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_central_phi_p6_jetmettau, leading_central_phi_p6_jetmet, leading_central_phi_p8_jetmettau, leading_central_phi_p8_jetmet, -0.020, 0.020, path_plots, plot_prefix+"leading_central_phi", "hadron-detector", "#phi^{central}", resultado, 2, detail);


   //forward jet pt
   if (detail) { cout << "Forward Jet pT..." << endl; }
   TH1D *leading_forward_pt_p6_jetmettau = 0;
   TH1D *leading_forward_pt_p6_jetmet = 0;
   TH1D *leading_forward_pt_p8_jetmettau = 0;
   TH1D *leading_forward_pt_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_forward_pt",leading_forward_pt_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_forward_pt",leading_forward_pt_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_forward_pt",leading_forward_pt_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_forward_pt",leading_forward_pt_p8_jetmet);


   if (leading_forward_pt_p6_jetmettau == 0) { cout << "leading_forward_pt_p6_jetmettau not found!" << endl; return; }
   if (leading_forward_pt_p6_jetmet == 0) { cout << "leading_forward_pt_p6_jetmet not found!" << endl; return; }
   if (leading_forward_pt_p8_jetmettau == 0) { cout << "leading_forward_pt_p8_jetmettau not found!" << endl; return; }
   if (leading_forward_pt_p8_jetmet == 0) { cout << "leading_forward_pt_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_forward_pt_p6_jetmettau, leading_forward_pt_p6_jetmet, leading_forward_pt_p8_jetmettau, leading_forward_pt_p8_jetmet, -0.08, 0.08, path_plots, plot_prefix+"leading_forward_pt", "hadron-detector/hadron", "p_{T}^{forward}", resultado, 3, detail);


   //forward jet eta
   if (detail) { cout << "Forward Jet Eta..." << endl; }
   TH1D *leading_forward_eta_p6_jetmettau = 0;
   TH1D *leading_forward_eta_p6_jetmet = 0;
   TH1D *leading_forward_eta_p8_jetmettau = 0;
   TH1D *leading_forward_eta_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_forward_eta",leading_forward_eta_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_forward_eta",leading_forward_eta_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_forward_eta",leading_forward_eta_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_forward_eta",leading_forward_eta_p8_jetmet);


   if (leading_forward_eta_p6_jetmettau == 0) { cout << "leading_forward_eta_p6_jetmettau not found!" << endl; return; }
   if (leading_forward_eta_p6_jetmet == 0) { cout << "leading_forward_eta_p6_jetmet not found!" << endl; return; }
   if (leading_forward_eta_p8_jetmettau == 0) { cout << "leading_forward_eta_p8_jetmettau not found!" << endl; return; }
   if (leading_forward_eta_p8_jetmet == 0) { cout << "leading_forward_eta_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_forward_eta_p6_jetmettau, leading_forward_eta_p6_jetmet, leading_forward_eta_p8_jetmettau, leading_forward_eta_p8_jetmet, -0.05, 0.05, path_plots, plot_prefix+"leading_forward_eta", "hadron-detector", "#eta^{forward}", resultado, 4, detail);


   //forward jet phi
   if (detail) { cout << "Forward Jet Phi..." << endl; }
   TH1D *leading_forward_phi_p6_jetmettau = 0;
   TH1D *leading_forward_phi_p6_jetmet = 0;
   TH1D *leading_forward_phi_p8_jetmettau = 0;
   TH1D *leading_forward_phi_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_forward_phi",leading_forward_phi_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_forward_phi",leading_forward_phi_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_forward_phi",leading_forward_phi_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_forward_phi",leading_forward_phi_p8_jetmet);


   if (leading_forward_phi_p6_jetmettau == 0) { cout << "leading_forward_phi_p6_jetmettau not found!" << endl; return; }
   if (leading_forward_phi_p6_jetmet == 0) { cout << "leading_forward_phi_p6_jetmet not found!" << endl; return; }
   if (leading_forward_phi_p8_jetmettau == 0) { cout << "leading_forward_phi_p8_jetmettau not found!" << endl; return; }
   if (leading_forward_phi_p8_jetmet == 0) { cout << "leading_forward_phi_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_forward_phi_p6_jetmettau, leading_forward_phi_p6_jetmet, leading_forward_phi_p8_jetmettau, leading_forward_phi_p8_jetmet, -0.03, 0.03, path_plots, plot_prefix+"leading_forward_phi", "hadron-detector", "#phi^{forward}", resultado, 5, detail);


   //delta phi
   if (detail) { cout << "Delta Phi..." << endl; }
   TH1D *delta_phi_p6_jetmettau = 0;
   TH1D *delta_phi_p6_jetmet = 0;
   TH1D *delta_phi_p8_jetmettau = 0;
   TH1D *delta_phi_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi",delta_phi_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi",delta_phi_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi",delta_phi_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi",delta_phi_p8_jetmet);


   if (delta_phi_p6_jetmettau == 0) { cout << "delta_phi_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_p6_jetmet == 0) { cout << "delta_phi_p6_jetmet not found!" << endl; return; }
   if (delta_phi_p8_jetmettau == 0) { cout << "delta_phi_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_p8_jetmet == 0) { cout << "delta_phi_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_p6_jetmettau, delta_phi_p6_jetmet, delta_phi_p8_jetmettau, delta_phi_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"delta_phi", "hadron-detector", "#Delta#phi", resultado, 6, detail);


   //delta eta
   if (detail) { cout << "Delta Eta..." << endl; }
   TH1D *delta_eta_p6_jetmettau = 0;
   TH1D *delta_eta_p6_jetmet = 0;
   TH1D *delta_eta_p8_jetmettau = 0;
   TH1D *delta_eta_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_eta",delta_eta_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_eta",delta_eta_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_eta",delta_eta_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_eta",delta_eta_p8_jetmet);


   if (delta_eta_p6_jetmettau == 0) { cout << "delta_eta_p6_jetmettau not found!" << endl; return; }
   if (delta_eta_p6_jetmet == 0) { cout << "delta_eta_p6_jetmet not found!" << endl; return; }
   if (delta_eta_p8_jetmettau == 0) { cout << "delta_eta_p8_jetmettau not found!" << endl; return; }
   if (delta_eta_p8_jetmet == 0) { cout << "delta_eta_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_eta_p6_jetmettau, delta_eta_p6_jetmet, delta_eta_p8_jetmettau, delta_eta_p8_jetmet, -0.10, 0.09, path_plots, plot_prefix+"delta_eta", "hadron-detector", "#Delta#eta", resultado, 7, detail);


   //delta phi deta1
   if (detail) { cout << "Delta Phi Deta1..." << endl; }
   TH1D *delta_phi_deta1_p6_jetmettau = 0;
   TH1D *delta_phi_deta1_p6_jetmet = 0;
   TH1D *delta_phi_deta1_p8_jetmettau = 0;
   TH1D *delta_phi_deta1_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta1",delta_phi_deta1_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta1",delta_phi_deta1_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta1",delta_phi_deta1_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta1",delta_phi_deta1_p8_jetmet);


   if (delta_phi_deta1_p6_jetmettau == 0) { cout << "delta_phi_deta1_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta1_p6_jetmet == 0) { cout << "delta_phi_deta1_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta1_p8_jetmettau == 0) { cout << "delta_phi_deta1_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta1_p8_jetmet == 0) { cout << "delta_phi_deta1_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta1_p6_jetmettau, delta_phi_deta1_p6_jetmet, delta_phi_deta1_p8_jetmettau, delta_phi_deta1_p8_jetmet, -0.08, 0.09, path_plots, plot_prefix+"delta_phi_deta1", "hadron-detector", "#Delta#phi", resultado, 8, detail);


   //delta phi deta2
   if (detail) { cout << "Delta Phi Deta2..." << endl; }
   TH1D *delta_phi_deta2_p6_jetmettau = 0;
   TH1D *delta_phi_deta2_p6_jetmet = 0;
   TH1D *delta_phi_deta2_p8_jetmettau = 0;
   TH1D *delta_phi_deta2_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta2",delta_phi_deta2_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta2",delta_phi_deta2_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta2",delta_phi_deta2_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta2",delta_phi_deta2_p8_jetmet);


   if (delta_phi_deta2_p6_jetmettau == 0) { cout << "delta_phi_deta2_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta2_p6_jetmet == 0) { cout << "delta_phi_deta2_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta2_p8_jetmettau == 0) { cout << "delta_phi_deta2_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta2_p8_jetmet == 0) { cout << "delta_phi_deta2_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta2_p6_jetmettau, delta_phi_deta2_p6_jetmet, delta_phi_deta2_p8_jetmettau, delta_phi_deta2_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"delta_phi_deta2", "hadron-detector", "#Delta#phi", resultado, 9, detail);


   //delta phi deta3
   if (detail) { cout << "Delta Phi Deta3..." << endl; }
   TH1D *delta_phi_deta3_p6_jetmettau = 0;
   TH1D *delta_phi_deta3_p6_jetmet = 0;
   TH1D *delta_phi_deta3_p8_jetmettau = 0;
   TH1D *delta_phi_deta3_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta3",delta_phi_deta3_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta3",delta_phi_deta3_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta3",delta_phi_deta3_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta3",delta_phi_deta3_p8_jetmet);


   if (delta_phi_deta3_p6_jetmettau == 0) { cout << "delta_phi_deta3_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta3_p6_jetmet == 0) { cout << "delta_phi_deta3_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta3_p8_jetmettau == 0) { cout << "delta_phi_deta3_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta3_p8_jetmet == 0) { cout << "delta_phi_deta3_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta3_p6_jetmettau, delta_phi_deta3_p6_jetmet, delta_phi_deta3_p8_jetmettau, delta_phi_deta3_p8_jetmet, -0.11, 0.11, path_plots, plot_prefix+"delta_phi_deta3", "hadron-detector", "#Delta#phi", resultado, 10, detail);


   //delta phi deta4
   if (detail) { cout << "Delta Phi Deta4..." << endl; }
   TH1D *delta_phi_deta4_p6_jetmettau = 0;
   TH1D *delta_phi_deta4_p6_jetmet = 0;
   TH1D *delta_phi_deta4_p8_jetmettau = 0;
   TH1D *delta_phi_deta4_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta4",delta_phi_deta4_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta4",delta_phi_deta4_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta4",delta_phi_deta4_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta4",delta_phi_deta4_p8_jetmet);


   if (delta_phi_deta4_p6_jetmettau == 0) { cout << "delta_phi_deta4_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta4_p6_jetmet == 0) { cout << "delta_phi_deta4_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta4_p8_jetmettau == 0) { cout << "delta_phi_deta4_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta4_p8_jetmet == 0) { cout << "delta_phi_deta4_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta4_p6_jetmettau, delta_phi_deta4_p6_jetmet, delta_phi_deta4_p8_jetmettau, delta_phi_deta4_p8_jetmet, -0.12, 0.12, path_plots, plot_prefix+"delta_phi_deta4", "hadron-detector", "#Delta#phi", resultado, 11, detail);


   //delta phi gap
   if (detail) { cout << "Delta Phi Gap ..." << endl; }
   TH1D *delta_phi_gap_p6_jetmettau = 0;
   TH1D *delta_phi_gap_p6_jetmet = 0;
   TH1D *delta_phi_gap_p8_jetmettau = 0;
   TH1D *delta_phi_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_gap",delta_phi_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_gap",delta_phi_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_gap",delta_phi_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_gap",delta_phi_gap_p8_jetmet);


   if (delta_phi_gap_p6_jetmettau == 0) { cout << "delta_phi_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_gap_p6_jetmet == 0) { cout << "delta_phi_gap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_gap_p8_jetmettau == 0) { cout << "delta_phi_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_gap_p8_jetmet == 0) { cout << "delta_phi_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_gap_p6_jetmettau, delta_phi_gap_p6_jetmet, delta_phi_gap_p8_jetmettau, delta_phi_gap_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"delta_phi_gap", "hadron-detector", "#Delta#phi", resultado, 12, detail);


   //delta eta gap
   if (detail) { cout << "Delta Eta Gap..." << endl; }
   TH1D *delta_eta_gap_p6_jetmettau = 0;
   TH1D *delta_eta_gap_p6_jetmet = 0;
   TH1D *delta_eta_gap_p8_jetmettau = 0;
   TH1D *delta_eta_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_eta_gap",delta_eta_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_eta_gap",delta_eta_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_eta_gap",delta_eta_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_eta_gap",delta_eta_gap_p8_jetmet);


   if (delta_eta_gap_p6_jetmettau == 0) { cout << "delta_eta_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_eta_gap_p6_jetmet == 0) { cout << "delta_eta_gap_p6_jetmet not found!" << endl; return; }
   if (delta_eta_gap_p8_jetmettau == 0) { cout << "delta_eta_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_eta_gap_p8_jetmet == 0) { cout << "delta_eta_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_eta_gap_p6_jetmettau, delta_eta_gap_p6_jetmet, delta_eta_gap_p8_jetmettau, delta_eta_gap_p8_jetmet, -0.09, 0.08, path_plots, plot_prefix+"delta_eta_gap", "hadron-detector", "#Delta#eta", resultado, 13, detail);


   //delta phi deta1 gap
   if (detail) { cout << "Delta Phi Deta1 Gap..." << endl; }
   TH1D *delta_phi_deta1_gap_p6_jetmettau = 0;
   TH1D *delta_phi_deta1_gap_p6_jetmet = 0;
   TH1D *delta_phi_deta1_gap_p8_jetmettau = 0;
   TH1D *delta_phi_deta1_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta1_gap",delta_phi_deta1_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta1_gap",delta_phi_deta1_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta1_gap",delta_phi_deta1_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta1_gap",delta_phi_deta1_gap_p8_jetmet);


   if (delta_phi_deta1_gap_p6_jetmettau == 0) { cout << "delta_phi_deta1_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta1_gap_p6_jetmet == 0) { cout << "delta_phi_deta1_gap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta1_gap_p8_jetmettau == 0) { cout << "delta_phi_deta1_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta1_gap_p8_jetmet == 0) { cout << "delta_phi_deta1_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta1_gap_p6_jetmettau, delta_phi_deta1_gap_p6_jetmet, delta_phi_deta1_gap_p8_jetmettau, delta_phi_deta1_gap_p8_jetmet, -0.08, 0.09, path_plots, plot_prefix+"delta_phi_deta1_gap", "hadron-detector", "#Delta#phi", resultado, 14, detail);


   //delta phi deta2 gap
   if (detail) { cout << "Delta Phi Deta2 Gap..." << endl; }
   TH1D *delta_phi_deta2_gap_p6_jetmettau = 0;
   TH1D *delta_phi_deta2_gap_p6_jetmet = 0;
   TH1D *delta_phi_deta2_gap_p8_jetmettau = 0;
   TH1D *delta_phi_deta2_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta2_gap",delta_phi_deta2_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta2_gap",delta_phi_deta2_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta2_gap",delta_phi_deta2_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta2_gap",delta_phi_deta2_gap_p8_jetmet);


   if (delta_phi_deta2_gap_p6_jetmettau == 0) { cout << "delta_phi_deta2_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta2_gap_p6_jetmet == 0) { cout << "delta_phi_deta2_gap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta2_gap_p8_jetmettau == 0) { cout << "delta_phi_deta2_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta2_gap_p8_jetmet == 0) { cout << "delta_phi_deta2_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta2_gap_p6_jetmettau, delta_phi_deta2_gap_p6_jetmet, delta_phi_deta2_gap_p8_jetmettau, delta_phi_deta2_gap_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"delta_phi_deta2_gap", "hadron-detector", "#Delta#phi", resultado, 15, detail);


   //delta phi deta3 gap
   if (detail) { cout << "Delta Phi Deta3 Gap..." << endl; }
   TH1D *delta_phi_deta3_gap_p6_jetmettau = 0;
   TH1D *delta_phi_deta3_gap_p6_jetmet = 0;
   TH1D *delta_phi_deta3_gap_p8_jetmettau = 0;
   TH1D *delta_phi_deta3_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta3_gap",delta_phi_deta3_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta3_gap",delta_phi_deta3_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta3_gap",delta_phi_deta3_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta3_gap",delta_phi_deta3_gap_p8_jetmet);


   if (delta_phi_deta3_gap_p6_jetmettau == 0) { cout << "delta_phi_deta3_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta3_gap_p6_jetmet == 0) { cout << "delta_phi_deta3_gap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta3_gap_p8_jetmettau == 0) { cout << "delta_phi_deta3_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta3_gap_p8_jetmet == 0) { cout << "delta_phi_deta3_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta3_gap_p6_jetmettau, delta_phi_deta3_gap_p6_jetmet, delta_phi_deta3_gap_p8_jetmettau, delta_phi_deta3_gap_p8_jetmet, -0.11, 0.11, path_plots, plot_prefix+"delta_phi_deta3_gap", "hadron-detector", "#Delta#phi", resultado, 16, detail);


   //delta phi deta4 gap
   if (detail) { cout << "Delta Phi Deta4 Gap..." << endl; }
   TH1D *delta_phi_deta4_gap_p6_jetmettau = 0;
   TH1D *delta_phi_deta4_gap_p6_jetmet = 0;
   TH1D *delta_phi_deta4_gap_p8_jetmettau = 0;
   TH1D *delta_phi_deta4_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta4_gap",delta_phi_deta4_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta4_gap",delta_phi_deta4_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta4_gap",delta_phi_deta4_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta4_gap",delta_phi_deta4_gap_p8_jetmet);


   if (delta_phi_deta4_gap_p6_jetmettau == 0) { cout << "delta_phi_deta4_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta4_gap_p6_jetmet == 0) { cout << "delta_phi_deta4_gap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta4_gap_p8_jetmettau == 0) { cout << "delta_phi_deta4_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta4_gap_p8_jetmet == 0) { cout << "delta_phi_deta4_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta4_gap_p6_jetmettau, delta_phi_deta4_gap_p6_jetmet, delta_phi_deta4_gap_p8_jetmettau, delta_phi_deta4_gap_p8_jetmet, -0.12, 0.12, path_plots, plot_prefix+"delta_phi_deta4_gap", "hadron-detector", "#Delta#phi", resultado, 17, detail);


   //delta phi nogap
   if (detail) { cout << "Delta Phi Nogap ..." << endl; }
   TH1D *delta_phi_nogap_p6_jetmettau = 0;
   TH1D *delta_phi_nogap_p6_jetmet = 0;
   TH1D *delta_phi_nogap_p8_jetmettau = 0;
   TH1D *delta_phi_nogap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_nogap",delta_phi_nogap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_nogap",delta_phi_nogap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_nogap",delta_phi_nogap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_nogap",delta_phi_nogap_p8_jetmet);


   if (delta_phi_nogap_p6_jetmettau == 0) { cout << "delta_phi_nogap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_nogap_p6_jetmet == 0) { cout << "delta_phi_nogap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_nogap_p8_jetmettau == 0) { cout << "delta_phi_nogap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_nogap_p8_jetmet == 0) { cout << "delta_phi_nogap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_nogap_p6_jetmettau, delta_phi_nogap_p6_jetmet, delta_phi_nogap_p8_jetmettau, delta_phi_nogap_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"delta_phi_nogap", "hadron-detector", "#Delta#phi", resultado, 18, detail);


   //delta eta nogap
   if (detail) { cout << "Delta Eta Nogap..." << endl; }
   TH1D *delta_eta_nogap_p6_jetmettau = 0;
   TH1D *delta_eta_nogap_p6_jetmet = 0;
   TH1D *delta_eta_nogap_p8_jetmettau = 0;
   TH1D *delta_eta_nogap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_eta_nogap",delta_eta_nogap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_eta_nogap",delta_eta_nogap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_eta_nogap",delta_eta_nogap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_eta_nogap",delta_eta_nogap_p8_jetmet);


   if (delta_eta_nogap_p6_jetmettau == 0) { cout << "delta_eta_nogap_p6_jetmettau not found!" << endl; return; }
   if (delta_eta_nogap_p6_jetmet == 0) { cout << "delta_eta_nogap_p6_jetmet not found!" << endl; return; }
   if (delta_eta_nogap_p8_jetmettau == 0) { cout << "delta_eta_nogap_p8_jetmettau not found!" << endl; return; }
   if (delta_eta_nogap_p8_jetmet == 0) { cout << "delta_eta_nogap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_eta_nogap_p6_jetmettau, delta_eta_nogap_p6_jetmet, delta_eta_nogap_p8_jetmettau, delta_eta_nogap_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"delta_eta_nogap", "hadron-detector", "#Delta#eta", resultado, 19, detail);


   //delta phi deta1 nogap
   if (detail) { cout << "Delta Phi Deta1 Nogap..." << endl; }
   TH1D *delta_phi_deta1_nogap_p6_jetmettau = 0;
   TH1D *delta_phi_deta1_nogap_p6_jetmet = 0;
   TH1D *delta_phi_deta1_nogap_p8_jetmettau = 0;
   TH1D *delta_phi_deta1_nogap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta1_nogap",delta_phi_deta1_nogap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta1_nogap",delta_phi_deta1_nogap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta1_nogap",delta_phi_deta1_nogap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta1_nogap",delta_phi_deta1_nogap_p8_jetmet);


   if (delta_phi_deta1_nogap_p6_jetmettau == 0) { cout << "delta_phi_deta1_nogap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta1_nogap_p6_jetmet == 0) { cout << "delta_phi_deta1_nogap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta1_nogap_p8_jetmettau == 0) { cout << "delta_phi_deta1_nogap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta1_nogap_p8_jetmet == 0) { cout << "delta_phi_deta1_nogap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta1_nogap_p6_jetmettau, delta_phi_deta1_nogap_p6_jetmet, delta_phi_deta1_nogap_p8_jetmettau, delta_phi_deta1_nogap_p8_jetmet, -0.11, 0.11, path_plots, plot_prefix+"delta_phi_deta1_nogap", "hadron-detector", "#Delta#phi", resultado, 20, detail);


   //delta phi deta2 nogap
   if (detail) { cout << "Delta Phi Deta2 Nogap..." << endl; }
   TH1D *delta_phi_deta2_nogap_p6_jetmettau = 0;
   TH1D *delta_phi_deta2_nogap_p6_jetmet = 0;
   TH1D *delta_phi_deta2_nogap_p8_jetmettau = 0;
   TH1D *delta_phi_deta2_nogap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta2_nogap",delta_phi_deta2_nogap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta2_nogap",delta_phi_deta2_nogap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta2_nogap",delta_phi_deta2_nogap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta2_nogap",delta_phi_deta2_nogap_p8_jetmet);


   if (delta_phi_deta2_nogap_p6_jetmettau == 0) { cout << "delta_phi_deta2_nogap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta2_nogap_p6_jetmet == 0) { cout << "delta_phi_deta2_nogap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta2_nogap_p8_jetmettau == 0) { cout << "delta_phi_deta2_nogap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta2_nogap_p8_jetmet == 0) { cout << "delta_phi_deta2_nogap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta2_nogap_p6_jetmettau, delta_phi_deta2_nogap_p6_jetmet, delta_phi_deta2_nogap_p8_jetmettau, delta_phi_deta2_nogap_p8_jetmet, -0.11, 0.11, path_plots, plot_prefix+"delta_phi_deta2_nogap", "hadron-detector", "#Delta#phi", resultado, 21, detail);


   //delta phi deta3 nogap
   if (detail) { cout << "Delta Phi Deta3 Nogap..." << endl; }
   TH1D *delta_phi_deta3_nogap_p6_jetmettau = 0;
   TH1D *delta_phi_deta3_nogap_p6_jetmet = 0;
   TH1D *delta_phi_deta3_nogap_p8_jetmettau = 0;
   TH1D *delta_phi_deta3_nogap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta3_nogap",delta_phi_deta3_nogap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta3_nogap",delta_phi_deta3_nogap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta3_nogap",delta_phi_deta3_nogap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta3_nogap",delta_phi_deta3_nogap_p8_jetmet);


   if (delta_phi_deta3_nogap_p6_jetmettau == 0) { cout << "delta_phi_deta3_nogap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta3_nogap_p6_jetmet == 0) { cout << "delta_phi_deta3_nogap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta3_nogap_p8_jetmettau == 0) { cout << "delta_phi_deta3_nogap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta3_nogap_p8_jetmet == 0) { cout << "delta_phi_deta3_nogap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta3_nogap_p6_jetmettau, delta_phi_deta3_nogap_p6_jetmet, delta_phi_deta3_nogap_p8_jetmettau, delta_phi_deta3_nogap_p8_jetmet, -0.11, 0.11, path_plots, plot_prefix+"delta_phi_deta3_nogap", "hadron-detector", "#Delta#phi", resultado, 22, detail);


   //delta phi deta4 nogap
   if (detail) { cout << "Delta Phi Deta4 Nogap..." << endl; }
   TH1D *delta_phi_deta4_nogap_p6_jetmettau = 0;
   TH1D *delta_phi_deta4_nogap_p6_jetmet = 0;
   TH1D *delta_phi_deta4_nogap_p8_jetmettau = 0;
   TH1D *delta_phi_deta4_nogap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_phi_deta4_nogap",delta_phi_deta4_nogap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_phi_deta4_nogap",delta_phi_deta4_nogap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_phi_deta4_nogap",delta_phi_deta4_nogap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_phi_deta4_nogap",delta_phi_deta4_nogap_p8_jetmet);


   if (delta_phi_deta4_nogap_p6_jetmettau == 0) { cout << "delta_phi_deta4_nogap_p6_jetmettau not found!" << endl; return; }
   if (delta_phi_deta4_nogap_p6_jetmet == 0) { cout << "delta_phi_deta4_nogap_p6_jetmet not found!" << endl; return; }
   if (delta_phi_deta4_nogap_p8_jetmettau == 0) { cout << "delta_phi_deta4_nogap_p8_jetmettau not found!" << endl; return; }
   if (delta_phi_deta4_nogap_p8_jetmet == 0) { cout << "delta_phi_deta4_nogap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(delta_phi_deta4_nogap_p6_jetmettau, delta_phi_deta4_nogap_p6_jetmet, delta_phi_deta4_nogap_p8_jetmettau, delta_phi_deta4_nogap_p8_jetmet, -0.09, 0.09, path_plots, plot_prefix+"delta_phi_deta4_nogap", "hadron-detector", "#Delta#phi", resultado, 23, detail);


   //leading pt inside gap
   if (detail) { cout << "Leading pT Inside Gap..." << endl; }
   TH1D *leading_pt_inside_gap_p6_jetmettau = 0;
   TH1D *leading_pt_inside_gap_p6_jetmet = 0;
   TH1D *leading_pt_inside_gap_p8_jetmettau = 0;
   TH1D *leading_pt_inside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_pt_inside_gap",leading_pt_inside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_pt_inside_gap",leading_pt_inside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_pt_inside_gap",leading_pt_inside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_pt_inside_gap",leading_pt_inside_gap_p8_jetmet);


   if (leading_pt_inside_gap_p6_jetmettau == 0) { cout << "leading_pt_inside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_pt_inside_gap_p6_jetmet == 0) { cout << "leading_pt_inside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_pt_inside_gap_p8_jetmettau == 0) { cout << "leading_pt_inside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_pt_inside_gap_p8_jetmet == 0) { cout << "leading_pt_inside_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_pt_inside_gap_p6_jetmettau, leading_pt_inside_gap_p6_jetmet, leading_pt_inside_gap_p8_jetmettau, leading_pt_inside_gap_p8_jetmet, -0.15, 0.25, path_plots, plot_prefix+"leading_pt_inside_gap", "hadron-detector/hadron", "p_{T}^{inside}", resultado, 24, detail);


   //leading eta inside gap
   if (detail) { cout << "Leading Eta Inside Gap..." << endl; }
   TH1D *leading_eta_inside_gap_p6_jetmettau = 0;
   TH1D *leading_eta_inside_gap_p6_jetmet = 0;
   TH1D *leading_eta_inside_gap_p8_jetmettau = 0;
   TH1D *leading_eta_inside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_eta_inside_gap",leading_eta_inside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_eta_inside_gap",leading_eta_inside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_eta_inside_gap",leading_eta_inside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_eta_inside_gap",leading_eta_inside_gap_p8_jetmet);


   if (leading_eta_inside_gap_p6_jetmettau == 0) { cout << "leading_eta_inside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_eta_inside_gap_p6_jetmet == 0) { cout << "leading_eta_inside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_eta_inside_gap_p8_jetmettau == 0) { cout << "leading_eta_inside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_eta_inside_gap_p8_jetmet == 0) { cout << "leading_eta_inside_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_eta_inside_gap_p6_jetmettau, leading_eta_inside_gap_p6_jetmet, leading_eta_inside_gap_p8_jetmettau, leading_eta_inside_gap_p8_jetmet, -0.07, 0.07, path_plots, plot_prefix+"leading_eta_inside_gap", "hadron-detector", "#eta^{inside}", resultado, 25, detail);


   //leading phi inside gap
   if (detail) { cout << "Leading Phi Inside Gap..." << endl; }
   TH1D *leading_phi_inside_gap_p6_jetmettau = 0;
   TH1D *leading_phi_inside_gap_p6_jetmet = 0;
   TH1D *leading_phi_inside_gap_p8_jetmettau = 0;
   TH1D *leading_phi_inside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_phi_inside_gap",leading_phi_inside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_phi_inside_gap",leading_phi_inside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_phi_inside_gap",leading_phi_inside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_phi_inside_gap",leading_phi_inside_gap_p8_jetmet);


   if (leading_phi_inside_gap_p6_jetmettau == 0) { cout << "leading_phi_inside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_phi_inside_gap_p6_jetmet == 0) { cout << "leading_phi_inside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_phi_inside_gap_p8_jetmettau == 0) { cout << "leading_phi_inside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_phi_inside_gap_p8_jetmet == 0) { cout << "leading_phi_inside_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_phi_inside_gap_p6_jetmettau, leading_phi_inside_gap_p6_jetmet, leading_phi_inside_gap_p8_jetmettau, leading_phi_inside_gap_p8_jetmet, -0.09, 0.09, path_plots, plot_prefix+"leading_phi_inside_gap", "hadron-detector", "#phi^{inside}", resultado, 26, detail);


   //leading eta star inside gap
   if (detail) { cout << "Leading Eta Star Inside Gap..." << endl; }
   TH1D *leading_eta_star_inside_gap_p6_jetmettau = 0;
   TH1D *leading_eta_star_inside_gap_p6_jetmet = 0;
   TH1D *leading_eta_star_inside_gap_p8_jetmettau = 0;
   TH1D *leading_eta_star_inside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_eta_star_inside_gap",leading_eta_star_inside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_eta_star_inside_gap",leading_eta_star_inside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_eta_star_inside_gap",leading_eta_star_inside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_eta_star_inside_gap",leading_eta_star_inside_gap_p8_jetmet);


   if (leading_eta_star_inside_gap_p6_jetmettau == 0) { cout << "leading_eta_star_inside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_eta_star_inside_gap_p6_jetmet == 0) { cout << "leading_eta_star_inside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_eta_star_inside_gap_p8_jetmettau == 0) { cout << "leading_eta_star_inside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_eta_star_inside_gap_p8_jetmet == 0) { cout << "leading_eta_star_inside_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_eta_star_inside_gap_p6_jetmettau, leading_eta_star_inside_gap_p6_jetmet, leading_eta_star_inside_gap_p8_jetmettau, leading_eta_star_inside_gap_p8_jetmet, -0.10, 0.10, path_plots, plot_prefix+"leading_eta_star_inside_gap", "hadron-detector", "#eta*^{inside}", resultado, 27, detail);


   //leading pt outside gap
   if (detail) { cout << "Leading pT Outside Gap..." << endl; }
   TH1D *leading_pt_outside_gap_p6_jetmettau = 0;
   TH1D *leading_pt_outside_gap_p6_jetmet = 0;
   TH1D *leading_pt_outside_gap_p8_jetmettau = 0;
   TH1D *leading_pt_outside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_pt_outside_gap",leading_pt_outside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_pt_outside_gap",leading_pt_outside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_pt_outside_gap",leading_pt_outside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_pt_outside_gap",leading_pt_outside_gap_p8_jetmet);


   if (leading_pt_outside_gap_p6_jetmettau == 0) { cout << "leading_pt_outside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_pt_outside_gap_p6_jetmet == 0) { cout << "leading_pt_outside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_pt_outside_gap_p8_jetmettau == 0) { cout << "leading_pt_outside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_pt_outside_gap_p8_jetmet == 0) { cout << "leading_pt_outside_gap_p8_jetmet not found!" << endl; return; }

   ////fit_and_plot(leading_pt_outside_gap_p6_jetmettau, leading_pt_outside_gap_p6_jetmet, leading_pt_outside_gap_p8_jetmettau, leading_pt_outside_gap_p8_jetmet, -0.15, 0.25, path_plots, plot_prefix+"leading_pt_outside_gap", "hadron-detector/hadron", "p_{T}^{outside}", resultado, 28, detail);


   //leading eta outside gap
   if (detail) { cout << "Leading Eta Outside Gap..." << endl; }
   TH1D *leading_eta_outside_gap_p6_jetmettau = 0;
   TH1D *leading_eta_outside_gap_p6_jetmet = 0;
   TH1D *leading_eta_outside_gap_p8_jetmettau = 0;
   TH1D *leading_eta_outside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_eta_outside_gap",leading_eta_outside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_eta_outside_gap",leading_eta_outside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_eta_outside_gap",leading_eta_outside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_eta_outside_gap",leading_eta_outside_gap_p8_jetmet);


   if (leading_eta_outside_gap_p6_jetmettau == 0) { cout << "leading_eta_outside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_eta_outside_gap_p6_jetmet == 0) { cout << "leading_eta_outside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_eta_outside_gap_p8_jetmettau == 0) { cout << "leading_eta_outside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_eta_outside_gap_p8_jetmet == 0) { cout << "leading_eta_outside_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_eta_outside_gap_p6_jetmettau, leading_eta_outside_gap_p6_jetmet, leading_eta_outside_gap_p8_jetmettau, leading_eta_outside_gap_p8_jetmet, -0.06, 0.06, path_plots, plot_prefix+"leading_eta_outside_gap", "hadron-detector", "#eta^{outside}", resultado, 29, detail);


   //leading phi outside gap
   if (detail) { cout << "Leading Phi Outside Gap..." << endl; }
   TH1D *leading_phi_outside_gap_p6_jetmettau = 0;
   TH1D *leading_phi_outside_gap_p6_jetmet = 0;
   TH1D *leading_phi_outside_gap_p8_jetmettau = 0;
   TH1D *leading_phi_outside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_leading_phi_outside_gap",leading_phi_outside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_leading_phi_outside_gap",leading_phi_outside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_leading_phi_outside_gap",leading_phi_outside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_leading_phi_outside_gap",leading_phi_outside_gap_p8_jetmet);


   if (leading_phi_outside_gap_p6_jetmettau == 0) { cout << "leading_phi_outside_gap_p6_jetmettau not found!" << endl; return; }
   if (leading_phi_outside_gap_p6_jetmet == 0) { cout << "leading_phi_outside_gap_p6_jetmet not found!" << endl; return; }
   if (leading_phi_outside_gap_p8_jetmettau == 0) { cout << "leading_phi_outside_gap_p8_jetmettau not found!" << endl; return; }
   if (leading_phi_outside_gap_p8_jetmet == 0) { cout << "leading_phi_outside_gap_p8_jetmet not found!" << endl; return; }

   ///fit_and_plot(leading_phi_outside_gap_p6_jetmettau, leading_phi_outside_gap_p6_jetmet, leading_phi_outside_gap_p8_jetmettau, leading_phi_outside_gap_p8_jetmet, -0.08, 0.08, path_plots, plot_prefix+"leading_phi_outside_gap", "hadron-detector", "#phi^{outside}", resultado, 30, detail);


   //leading delta eta outside gap
   if (detail) { cout << "Leading Delta Eta Outside Gap..." << endl; }
   TH1D *delta_eta_outside_gap_p6_jetmettau = 0;
   TH1D *delta_eta_outside_gap_p6_jetmet = 0;
   TH1D *delta_eta_outside_gap_p8_jetmettau = 0;
   TH1D *delta_eta_outside_gap_p8_jetmet = 0;

   p6_z2_jetmettau->GetObject("res_delta_eta_outside_gap",delta_eta_outside_gap_p6_jetmettau);
   p6_z2_jetmet->GetObject("res_delta_eta_outside_gap",delta_eta_outside_gap_p6_jetmet);
   p8_4c_jetmettau->GetObject("res_delta_eta_outside_gap",delta_eta_outside_gap_p8_jetmettau);
   p8_4c_jetmet->GetObject("res_delta_eta_outside_gap",delta_eta_outside_gap_p8_jetmet);


   if (delta_eta_outside_gap_p6_jetmettau == 0) { cout << "delta_eta_outside_gap_p6_jetmettau not found!" << endl; return; }
   if (delta_eta_outside_gap_p6_jetmet == 0) { cout << "delta_eta_outside_gap_p6_jetmet not found!" << endl; return; }
   if (delta_eta_outside_gap_p8_jetmettau == 0) { cout << "delta_eta_outside_gap_p8_jetmettau not found!" << endl; return; }
   if (delta_eta_outside_gap_p8_jetmet == 0) { cout << "delta_eta_outside_gap_p8_jetmet not found!" << endl; return; }

   fit_and_plot(delta_eta_outside_gap_p6_jetmettau, delta_eta_outside_gap_p6_jetmet, delta_eta_outside_gap_p8_jetmettau, delta_eta_outside_gap_p8_jetmet, -0.05, 0.08, path_plots, plot_prefix+"delta_eta_outside_gap", "hadron-detector", "#Delta#eta^{outside}", resultado, 31, detail);

//  TH1D *res_delta_eta_outside_gap;

   if (show_result)
	{
   	for (int i = 0; i < 4*dist; i++)
		{
		resultado[i] = resultado[i]*1000;
		}
	show_results(resultado);
	}

}
