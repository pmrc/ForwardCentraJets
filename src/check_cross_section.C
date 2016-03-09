// Pedro Cipriano, Feb 2016
// FESB, CMS
// Last Update: 15 Feb 2016
//

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "common_methods.h"

void check_cross_section(string input_file1, string label1, bool isgen1, string input_file2, string label2 = "Label 2", bool isgen2 = false, string input_file3 = "dummy3.root", string label3 = "Label 3", bool isgen3 = false, bool detail = true)
{

   if (detail) { cout<<"Input file1 :             "<<input_file1<<endl; }
   if (detail) { cout<<"Label for file1 :         "<<label1<<endl; }
   if (detail) { cout<<"File1 is Gen?:            "<<isgen1<<endl; }
   if (detail) { cout<<"Input file2 :             "<<input_file2<<endl; }
   if (detail) { cout<<"Label for file2 :         "<<label2<<endl; }
   if (detail) { cout<<"File2 is Gen?:            "<<isgen2<<endl; }
   if (detail) { cout<<"Input file3 :             "<<input_file3<<endl; }
   if (detail) { cout<<"Label for file3 :         "<<label3<<endl; }
   if (detail) { cout<<"File3 is Gen?:            "<<isgen3<<endl; }
   if (detail) { cout<<"Detail :                  "<<detail<<endl; }

//opening the input and output data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *file1 = new TFile( input_file1.c_str() );
    TFile *file2 = new TFile( input_file2.c_str() );
    TFile *file3 = new TFile( input_file3.c_str() );

//reference values
    float central[7] = {10000., 5200., 1900., 590., 135., 28., 6.};
    float forward[7] = {21000., 9200., 2900., 690., 110., 10., 0.4};

//bin width values
    int bins[7] = {10, 12, 15, 18, 30, 30, 50};


//get the histograms for the leading central jet
    if (detail) { cout<<"Leading Central Jet"<<endl; }

    TH1D *leading_central_pt_file1 = 0;
    TH1D *leading_central_pt_file2 = 0;
    TH1D *leading_central_pt_file3 = 0;
    TString leading_central_pt_file1_name = "ak5PF_leading_central_pt";
	if (isgen1) { leading_central_pt_file1_name = "ak5Gen_leading_central_pt"; }
    TString leading_central_pt_file2_name = "ak5PF_leading_central_pt";
	if (isgen2) { leading_central_pt_file2_name = "ak5Gen_leading_central_pt"; }
    TString leading_central_pt_file3_name = "ak5PF_leading_central_pt";
	if (isgen3) { leading_central_pt_file3_name = "ak5Gen_leading_central_pt"; }

    file1->GetObject(leading_central_pt_file1_name,leading_central_pt_file1);
    if (leading_central_pt_file1 == 0) { cout << leading_central_pt_file1_name << " not found!" << endl; return; }
    file2->GetObject(leading_central_pt_file2_name,leading_central_pt_file2);
    if (leading_central_pt_file2 == 0) { cout << leading_central_pt_file2_name << " not found!" << endl; return; }
    file3->GetObject(leading_central_pt_file3_name,leading_central_pt_file3);
    if (leading_central_pt_file3 == 0) { cout << leading_central_pt_file3_name << " not found!" << endl; return; }

	for (int a=1; a <= leading_central_pt_file1->GetNbinsX(); a++ )
		{
		if (detail)
			{
			cout << "Bin: " << a << " -> " << central[a-1] << endl;
			double val1 = leading_central_pt_file1->GetBinContent(a)/(5.6*bins[a-1]);
			double val2 = leading_central_pt_file2->GetBinContent(a)/(5.6*bins[a-1]);
			double val3 = leading_central_pt_file3->GetBinContent(a)/(5.6*bins[a-1]);
			cout << val1 << " / " << val2  << " / " << val3  << endl;
			double dif1 = 100*(central[a-1] - val1)/central[a-1];
			double dif2 = 100*(central[a-1] - val2)/central[a-1];
			double dif3 = 100*(central[a-1] - val3)/central[a-1];
			cout << dif1 << " / " << dif2  << " / " << dif3  << endl;
			} 
		}


//get the histograms for the leading forward jet
    if (detail) { cout<<"Leading Forward Jet"<<endl; }

    TH1D *leading_forward_pt_file1 = 0;
    TH1D *leading_forward_pt_file2 = 0;
    TH1D *leading_forward_pt_file3 = 0;
    TString leading_forward_pt_file1_name = "ak5PF_leading_forward_pt";
	if (isgen1) { leading_forward_pt_file1_name = "ak5Gen_leading_forward_pt"; }
    TString leading_forward_pt_file2_name = "ak5PF_leading_forward_pt";
	if (isgen2) { leading_forward_pt_file2_name = "ak5Gen_leading_forward_pt"; }
    TString leading_forward_pt_file3_name = "ak5PF_leading_forward_pt";
	if (isgen3) { leading_forward_pt_file3_name = "ak5Gen_leading_forward_pt"; }

    file1->GetObject(leading_forward_pt_file1_name,leading_forward_pt_file1);
    if (leading_forward_pt_file1 == 0) { cout << leading_forward_pt_file1_name << " not found!" << endl; return; }
    file2->GetObject(leading_forward_pt_file2_name,leading_forward_pt_file2);
    if (leading_forward_pt_file2 == 0) { cout << leading_forward_pt_file2_name << " not found!" << endl; return; }
    file3->GetObject(leading_forward_pt_file3_name,leading_forward_pt_file3);
    if (leading_forward_pt_file3 == 0) { cout << leading_forward_pt_file3_name << " not found!" << endl; return; }

	for (int a=1; a <= leading_forward_pt_file1->GetNbinsX(); a++ )
		{
		if (detail)
			{
			cout << "Bin: " << a << " -> " << forward[a-1] << endl;
			double val1 = leading_forward_pt_file1->GetBinContent(a)/(3.0*bins[a-1]);
			double val2 = leading_forward_pt_file2->GetBinContent(a)/(3.0*bins[a-1]);
			double val3 = leading_forward_pt_file3->GetBinContent(a)/(3.0*bins[a-1]);
			cout << val1 << " / " << val2  << " / " << val3  << endl;
			double dif1 = 100*(forward[a-1] - val1)/forward[a-1];
			double dif2 = 100*(forward[a-1] - val2)/forward[a-1];
			double dif3 = 100*(forward[a-1] - val3)/forward[a-1];
			cout << dif1 << " / " << dif2  << " / " << dif3  << endl;
			} 
		}



}
