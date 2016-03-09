// Pedro Cipriano, Set 2013
// DESY, CMS
// Last Update: 26 Set 2013
//
//piotr_analysis(string *data_in, string data_out, double *data_lumi, int n_files, string data_type = "DATA", string sel_mode = "1vertex", bool detail = false)
// reads the Ntuples and created a root file with cross-section histograms

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <string>

#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "common_methods.h"

using namespace std;

void piotr_analysis(string *data_in, string data_out, double *data_lumi, int n_files, string data_type = "DATA", string sel_mode = "1vertex", string vertex_weights = "", TString vertex_sufix = "_v0", bool detail = false, bool test = false)
{ 

//output the configuration
   if (detail) { cout<<"Read NTuple Configuration"<<endl; }
   if (detail) { cout<<"Data Type :            "<<data_type<<endl; }
   if (detail) { cout<<"Selection mode :       "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files :      "<<n_files<<endl; }
   if (detail) { cout<<"Vertex Weights File :  "<<vertex_weights<<endl; }
   if (detail) { cout<<"Vertex Weights Sufix : "<<vertex_sufix<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :          "<<data_out<<endl; }
   if (detail) { cout<<"Detail level :         "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :            "<<test<<endl; }

//opening the vertex weight file
    bool vertex_weighted = false;
    TH1D *vertex_hist = 0;
    if ((vertex_weights != "") and (data_type == "MC_DET" or data_type == "MC_GEN"))
    {
    if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
    TFile *vertex_file = new TFile( vertex_weights.c_str() );
    vertex_file->GetObject("vertex_weights"+vertex_sufix,vertex_hist);
    if (vertex_hist != 0) { vertex_weighted = true; }
    }

//setting the prefix
   TString prefix;
   if (data_type == "DATA" or data_type == "MC_DET") { prefix = "ak5PF_"; }
   if (data_type == "MC_GEN") { prefix = "ak5Gen_"; }
   if (detail) { cout<<"Prefix :          "<<prefix<<endl; }

//main cuts
   double pt_min = 20.0;
 
   double pt, eta, phi, unc, pt_up, pt_down, gen_pt, gen_eta, gen_phi;
   double prescale, prescale_prov, vertex_factor, leading_pt;
   unsigned int number_of_jets;
   bool pass_tight, pv_pass, jet_selected, event_selected;

   double pt1, pt2, pt3, pt4, eta1, eta2, eta3, phi1, phi2, phi3, delta_phi13;

//declare the main binning


//declare the histograms
   TH1D *hist_delta_phi;

   hist_delta_phi =  new TH1D(prefix+"delta_phi","#Delta#phi;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", 10, 0, 2*3.15);

   hist_delta_phi->Sumw2();


  int counter_hlt = 0, counter_pv = 0, counter_jet = 0, counter_selected = 0, counter_entries = 0;
  int counter_jet15u = 0, counter_jet30u = 0, counter_jet50u = 0;
  double counter_weigth = 0.0, counter_check = 0.0;
  int nentries = 0;
  double pu_scale = 0.0;

// multiple triggers
  int number_triggers = 7;
  string trigger_list[number_triggers];
  int hlt_list[number_triggers];
//  double probs[number_triggers];
//  double prescales[number_triggers];
  bool hltPass(false);
//  bool trig_pass[number_triggers];
  trigger_list[0] = "HLT_Jet15U";
  trigger_list[1] = "HLT_Jet30U";
  trigger_list[2] = "HLT_Jet50U";
  trigger_list[3] = "HLT_DiJetAve15U_8E29";
  trigger_list[4] = "HLT_DiJetAve30U_8E29";
  trigger_list[5] = "HLT_DoubleJet15U_ForwardBackward";
  trigger_list[6] = "HLT_FwdJet20U";

//loop over the files
for (int z = 0; z < n_files; z++)
{
//open the file
  string file = data_in[z];
  cout<<z+1<<"/"<<n_files<<" -> "<<file<<endl;
  TFile *inf = TFile::Open( file.c_str() );

//define the tree branch
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);

//get the trigger information
  if (data_type == "DATA")
  {
  TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
  cout<<"Finding trigger mapping: "<<endl;
  //----------- loop over the X-axis labels -----------------
  for (int trg = 0; trg < number_triggers; trg++)
	{
	for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
		TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
		if (ss == trigger_list[trg]) {
			hlt_list[trg] = ibin;
			continue;
			}
		}
  	if (hlt_list[trg] == -1) {
		cout<<"The requested trigger ("<<trigger_list[trg]<<") is not found." << endl;
    		//break;
  		}
	else {
		cout<<trigger_list[trg]<<" --> "<<hlt_list[trg]<<endl;
		}
	}
  }

  prescale = 0.0;
  prescale_prov = 0.0;
  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  if (data_type == "MC_DET" or data_type == "MC_GEN")
  {
  prescale = data_lumi[z];
  cout<<"Cross-section per event : "<<prescale<<endl;
  }
  if (data_type == "DATA")
  {
  prescale_prov =  1.0 / data_lumi[z];
  cout<<"Luminosity for " << trigger_list[z] << " = " << data_lumi[z] << endl;
  cout<<"Cross-section per event : "<<prescale_prov<<endl;
  }

  int decade = 0;

  if (test) { nentries = 100000; } //reduced number of read entries, usefull for testing

   for (Int_t i=0;i<nentries;i++) {
     counter_entries++; 
     double progress = 10.0*i/(1.0*nentries);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
      cout<<k*10<<" % "<<endl;
     decade = k; 
     
     tr->GetEntry(i);

      //-------- check if the primary vertex is good ----  && Event->evtHdr().nVtxGood() == 1
      pv_pass = false;
      if (data_type == "MC_GEN")
	{
	if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)	{ pv_pass = true; }
	if (sel_mode == "allvertex") 					{ pv_pass = true; }
	}
      if (data_type == "MC_DET")
	{
	if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)						{ pv_pass = true; }
	if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) 					{ pv_pass = true; }
	if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)	{ pv_pass = true; }
	}
      if (data_type == "DATA")
	{
	if ((sel_mode == "up" || sel_mode == "down" || sel_mode == "allvertex") && Event->evtHdr().isPVgood() == 1)	{ pv_pass = true; }
	if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)    	  	{ pv_pass = true; }
	}

      if (pv_pass == true)
      {
        counter_pv++;
	vertex_factor = 1.0;
        if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1);} 

     if (data_type == "DATA")
     {
     // multiple triggers
     hltPass = false;
     prescale = 0.0;
	/*
	for (int trg = 0; trg < 3; trg++)
		{
		prescales[trg] = 1.0;
		trig_pass[trg] = false;
		if (hlt_list[trg] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[trg]) > 0) {
				hltPass = true;
				trig_pass[trg] = true;
				prescales[trg] = Event->preL1(hlt_list[trg]) * Event->preHLT(hlt_list[trg]);
				if (trg == 0) { counter_jet15u++; }
				if (trg == 1) { counter_jet30u++; }
				if (trg == 2) { counter_jet50u++; }
				}
			}
		}
	*/
	if (Event->fired(hlt_list[z]) > 0) { hltPass = true; }
     }

     if (hltPass and data_type == "DATA")
     {
     
     leading_pt = 0.0;

     for(unsigned int j=0; j<Event->nPFJets(); j++)
        {
	pt = Event->pfjet(j).ptCor();
	if (leading_pt < pt && Event->pfjet(j).tightID()) { leading_pt = pt; }
	}

	if (z == 0 and leading_pt > 35.0 and leading_pt < 70.0) { prescale = prescale_prov; }
	if (z == 1 and leading_pt > 70.0 and leading_pt < 110.0) { prescale = prescale_prov; }
	if (z == 2 and leading_pt > 110.0) { prescale = prescale_prov; }

     }

     number_of_jets = 0;
     pt = 0.0;
     eta = -10.0;
     phi = -10.0;
     pt1 = 0.0;
     pt2 = 0.0;
     pt3 = 0.0;
     pt4 = 0.0;
     eta1 = -10.0;
     eta2 = -10.0;
     eta3 = -10.0;
     phi1 = -10.0;
     phi2 = -10.0;
     phi3 = -10.0;
     delta_phi13 = -10.0;
     pass_tight = false;
     event_selected = false;

     if (prescale > 0.0) { //or prescale > 0 hltPass_monitor
      counter_hlt++;
      //cout<<"Prescale = "<<prescale<<endl;

     if (data_type == "DATA" or data_type == "MC_DET") { number_of_jets = Event->nPFJets(); }
     if (data_type == "MC_GEN") { number_of_jets = Event->nGenJets(); }

    for(unsigned int j=0; j<number_of_jets; j++) {
	if (data_type == "DATA" or data_type == "MC_DET")
	{
        pass_tight = Event->pfjet(j).tightID();
	pt = Event->pfjet(j).ptCor();
	eta = Event->pfjet(j).eta();
	phi = Event->pfjet(j).phi();
	}
	if (data_type == "MC_GEN")
	{
	pass_tight = true;
	pt = Event->genjet(j).pt();
	eta = Event->genjet(j).eta();
	phi = Event->genjet(j).phi();
	}
	if (data_type == "MC_DET")
		{
		gen_pt = Event->genjet(j).pt();
		gen_eta = Event->genjet(j).eta();
		gen_phi = Event->genjet(j).phi();
		pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,test);
		}
	if (data_type == "DATA")
	{
	unc = Event->pfjet(j).unc();
	pt_up = pt * (1+unc);
	pt_down = pt * (1-unc);
	if (sel_mode == "up") { pt = pt_up;}
	if (sel_mode == "down") { pt = pt_down;}
	}
     if (pt >= pt_min && pass_tight)
	{
	jet_selected = false;
     	counter_jet++;
	//if (pt > pt1 and abs(eta) <= 2.8) { pt1 = pt; eta1 = eta; phi1 = phi; jet_selected = true; }
	//if (pt > pt2 and abs(eta) <= 2.8 and !jet_selected)  { pt2 = pt; eta2 = eta; phi2 = phi; jet_selected = true; }
	if (pt > pt1 and abs(eta) >= 3.2 and abs(eta) <= 4.7)
		{
		if (pt1 > pt2 and pt1 > 0)
			{
			if (pt2 > 0 and pt2 > pt3) {pt3 = pt2; eta3 = eta2; phi3 = phi2;}
			pt2 = pt1; eta2 = eta1; phi2 = phi1;
			}
		pt1 = pt; eta1 = eta; phi1 = phi; jet_selected = true;
		}

	if (pt > pt2 and abs(eta) >= 3.2 and abs(eta) <= 4.7 and !jet_selected)
		{
		if (pt2 > pt3 and pt1 > 0)
			{ pt3 = pt2; eta3 = eta2; phi3 = phi2;}
		pt2 = pt; eta2 = eta; phi2 = phi; jet_selected = true;
		}
	if (pt > pt3 and abs(eta) >= 3.2 and abs(eta) <= 4.7 and !jet_selected)
		{ pt3 = pt; eta3 = eta; phi3 = phi; jet_selected = true; }
        if (pt > pt4 and !jet_selected) { pt4 = pt; }
  	} 

	}

if (pt1 > 0.0 and pt2 > 0.0 and pt3 > 0.0 and pt3 < pt2 and pt4 == 0.0)
	{
	event_selected = true;
	counter_selected++;
	counter_check = counter_check + prescale * vertex_factor;
	delta_phi13 = abs(phi1 - phi3);
	if (test) { cout << "Delta Phi = " << delta_phi13 << endl; }
	hist_delta_phi->Fill(delta_phi13, prescale * vertex_factor);
	}

if (event_selected and test) { cout << "Jets pT:        " << pt1 << " " << pt2 << " " << pt3 << " " << pt4 << endl; }
if (event_selected and test) { cout << "Jets Eta:       " << eta1 << " " << eta2 << " " << eta3 << endl; }
if (event_selected and test) { cout << "Jets Phi:       " << phi1 << " " << phi2 << " " << phi3 << endl; }
if (event_selected and test) { cout << "Jets Delta Phi: " << delta_phi13 << endl; }

	}
  }  } 
  
    cout<<"Events read:                      "<<nentries<<endl; 
 } 
  
  cout<<endl<<endl;
  cout<<"Events read:                      "<<counter_entries<<endl;
  cout<<"Events after the PV cut:          "<<counter_pv<<endl;
  if (data_type == "DATA")
	{
	cout<<"Events after the trigger cut:     "<<counter_hlt<<endl;
	cout<<"Events which fired HLT_Jet15U:    "<<counter_jet15u<<endl;
	cout<<"Events which fired HLT_Jet30U:    "<<counter_jet30u<<endl;
	cout<<"Events which fired HLT_Jet50U:    "<<counter_jet50u<<endl;
	}
  cout<<"Events Main Selection:            "<<counter_selected<<endl;
  cout<<"Events weight count:              "<<counter_weigth<<endl;
  cout<<"Total number of Jets:             "<<counter_jet<<endl;
  cout<<"Cross-section Check:              "<<counter_check<<endl;

     pu_scale = (double) nentries/ (double) counter_pv;
     cout<<endl;
     cout<<"Pile up correction =  "<<pu_scale<<endl;
     cout<<endl;

     //normalize the histograms
     normalize_histogram(hist_delta_phi, "Delta Phi");
    
     //Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     //save histograms on the file
     hist_delta_phi->Write();

     //close the output file
     data_output->Close();

     //delete the histograms to avoid memory leak
     delete(hist_delta_phi);

}
