// Pedro Cipriano, Jul 2013
// DESY, CMS
// Last Update: 30 Jul 2013
//
// 

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
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

void compute_resolution(string *data_in, string data_out, double *data_lumi, int n_files, string sel_mode = "1vertex", string vertex_weights = "", TString vertex_sufix = "_v0", bool detail = false, bool test = false)
{

//output the configuration
   if (detail) { cout<<"Read NTuple Configuration"<<endl; }
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
    if (vertex_weights != "")
    {
    if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
    TFile *vertex_file = new TFile( vertex_weights.c_str() );
    vertex_file->GetObject("vertex_weights"+vertex_sufix,vertex_hist);
    if (vertex_hist != 0) { vertex_weighted = true; }
    }

//main cuts
   double pt_min = 35.0;
   double gap_req = 20.0;

//declaring variables
   int nentries;
   double prescale, vertex_factor, weightMC;
   double pt, eta, phi, gen_pt, gen_eta, gen_phi;
   double deta_out1, deta_out2;
   double gen_central_pt, det_central_pt, gen_forward_pt, det_forward_pt;
   double gen_inside_pt, det_inside_pt, gen_outside_pt, det_outside_pt;
   double gen_central_eta, det_central_eta, gen_forward_eta, det_forward_eta;
   double gen_inside_eta, det_inside_eta, gen_outside_eta, det_outside_eta;
   double gen_central_phi, det_central_phi, gen_forward_phi, det_forward_phi;
   double gen_inside_phi, det_inside_phi, gen_outside_phi, det_outside_phi;
   double gen_inside_eta_star, det_inside_eta_star, gen_delta_eta_outside, det_delta_eta_outside;
   double gen_delta_phi, det_delta_phi, gen_delta_eta, det_delta_eta;
   double res_central_pt_val, res_forward_pt_val, res_inside_pt_val, res_outside_pt_val;;
   double res_central_eta_val, res_forward_eta_val, res_inside_eta_val, res_outside_eta_val;
   double res_central_phi_val, res_forward_phi_val, res_inside_phi_val, res_outside_phi_val;
   double res_delta_phi_val, res_delta_eta_val, res_inside_eta_star_val, res_delta_eta_outside_val;

   bool pv_pass, pass_tight;
   bool gen_selected, det_selected;
   bool gen_gap, det_gap, gen_nogap, det_nogap, gen_out, det_out;
   bool gen_deta1, gen_deta2, gen_deta3, gen_deta4;
   bool det_deta1, det_deta2, det_deta3, det_deta4;

   int counter_entries = 0, counter_pv = 0, gen_counter_jet = 0, det_counter_jet = 0;
   int det_counter_selected = 0, gen_counter_selected = 0;
   int det_counter_gap = 0, gen_counter_gap = 0;
   int det_counter_nogap = 0, gen_counter_nogap = 0;
   int det_counter_outside = 0, gen_counter_outside = 0;
   double det_counter_weigth = 0.0, gen_counter_weigth = 0.0;

//declaring histograms
  TH1D *res_central_pt;
  TH1D *res_central_eta;
  TH1D *res_central_phi;
  TH1D *res_forward_pt;
  TH1D *res_forward_eta;
  TH1D *res_forward_phi;
  TH1D *res_delta_phi;
  TH1D *res_delta_eta;
  TH1D *res_delta_phi_deta1;
  TH1D *res_delta_phi_deta2;
  TH1D *res_delta_phi_deta3;
  TH1D *res_delta_phi_deta4;
  TH1D *res_delta_phi_gap;
  TH1D *res_delta_eta_gap;
  TH1D *res_delta_phi_deta1_gap;
  TH1D *res_delta_phi_deta2_gap;
  TH1D *res_delta_phi_deta3_gap;
  TH1D *res_delta_phi_deta4_gap;
  TH1D *res_delta_phi_nogap;
  TH1D *res_delta_eta_nogap;
  TH1D *res_delta_phi_deta1_nogap;
  TH1D *res_delta_phi_deta2_nogap;
  TH1D *res_delta_phi_deta3_nogap;
  TH1D *res_delta_phi_deta4_nogap;
  TH1D *res_leading_pt_inside_gap;
  TH1D *res_leading_eta_inside_gap;
  TH1D *res_leading_phi_inside_gap;
  TH1D *res_leading_eta_star_inside_gap;
  TH1D *res_leading_pt_outside_gap;
  TH1D *res_leading_eta_outside_gap;
  TH1D *res_leading_phi_outside_gap;
  TH1D *res_delta_eta_outside_gap;

  res_central_pt =  new TH1D("res_central_pt","Central pt;p_T;Events", 800, -1.1, 1.0);
  res_central_eta =  new TH1D("res_central_eta","Central eta;#eta;Events", 4000, -2.0, 2.0);
  res_central_phi =  new TH1D("res_central_phi","Central phi;#phi;Events", 4000, -2.0, 2.0);
  res_forward_pt =  new TH1D("res_forward_pt","Forward pt;p_T;Events", 800, -1.2, 1.0);
  res_forward_eta =  new TH1D("res_forward_eta","Forward eta;#eta;Events", 4000, -2.0, 2.0);
  res_forward_phi =  new TH1D("res_forward_phi","Forward phi;#phi;Events", 4000, -2.0, 2.0);
  res_delta_phi =  new TH1D("res_delta_phi","#Delta#phi main scenario;#Delta#phi;Events", 4000, -2.0, 2.0);
  res_delta_eta =  new TH1D("res_delta_eta","#Delta#eta main scenario;#Delta#eta;Events", 4000, -2.0, 2.0);
  res_delta_phi_deta1 =  new TH1D("res_delta_phi_deta1","#Delta#phi deta1 main scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta2 =  new TH1D("res_delta_phi_deta2","#Delta#phi deta2 main scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta3 =  new TH1D("res_delta_phi_deta3","#Delta#phi deta3 main scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta4 =  new TH1D("res_delta_phi_deta4","#Delta#phi deta4 main scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_gap =  new TH1D("res_delta_phi_gap","#Delta#phi inside-jet veto scenario;#Delta#phi;Events", 4000, -2.0, 2.0);
  res_delta_eta_gap =  new TH1D("res_delta_eta_gap","#Delta#eta inside-jet veto scenario;#Delta#eta;Events", 4000, -2.0, 2.0);
  res_delta_phi_deta1_gap =  new TH1D("res_delta_phi_deta1_gap","#Delta#phi deta1 inside-jet veto scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta2_gap =  new TH1D("res_delta_phi_deta2_gap","#Delta#phi deta2 inside-jet veto scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta3_gap =  new TH1D("res_delta_phi_deta3_gap","#Delta#phi deta3 inside-jet veto scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta4_gap =  new TH1D("res_delta_phi_deta4_gap","#Delta#phi deta4 inside-jet veto scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_nogap =  new TH1D("res_delta_phi_nogap","#Delta#phi inside-jet tag scenario;#Delta#phi;Events", 4000, -2.0, 2.0);
  res_delta_eta_nogap =  new TH1D("res_delta_eta_nogap","#Delta#eta inside-jet tag scenario;#Delta#eta;Events", 4000, -2.0, 2.0);
  res_delta_phi_deta1_nogap =  new TH1D("res_delta_phi_deta1_nogap","#Delta#phi deta1 inside-jet tag scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta2_nogap =  new TH1D("res_delta_phi_deta2_nogap","#Delta#phi deta2 inside-jet tag scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta3_nogap =  new TH1D("res_delta_phi_deta3_nogap","#Delta#phi deta3 inside-jet tag scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_delta_phi_deta4_nogap =  new TH1D("res_delta_phi_deta4_nogap","#Delta#phi deta4 inside-jet tag scenario;#Delta#phi;Events", 800, -2.0, 2.0);
  res_leading_pt_inside_gap =  new TH1D("res_leading_pt_inside_gap","Leading pT Inside Gap;p_T;Events", 400, -1.1, 1.0);
  res_leading_eta_inside_gap =  new TH1D("res_leading_eta_inside_gap","Leading #eta Inside Gap;#eta;Events", 400, -2.0, 2.0);
  res_leading_phi_inside_gap =  new TH1D("res_leading_phi_inside_gap","Leading #phi Inside Gap;#phi;Events", 400, -1.0, 1.0);
  res_leading_eta_star_inside_gap =  new TH1D("res_leading_eta_star_inside_gap","Leading #eta* Inside Gap;#eta*;Events", 400, -2.0, 2.0);
  res_leading_pt_outside_gap =  new TH1D("res_leading_pt_outside_gap","Leading pT Outside Gap;p_T;Events", 400, -1.1, 1.0);
  res_leading_eta_outside_gap =  new TH1D("res_leading_eta_outside_gap","Leading #eta Outside Gap;#eta;Events", 400, -2.0, 2.0);
  res_leading_phi_outside_gap =  new TH1D("res_leading_phi_outside_gap","Leading #phi Outside Gap;#phi;Events", 400, -2.0, 2.0);
  res_delta_eta_outside_gap =  new TH1D("res_delta_eta_outside_gap","#Delta#eta* Outside Gap;#Delta#eta*;Events", 400, -2.0, 2.0);

  res_central_pt->Sumw2();
  res_central_eta->Sumw2();
  res_central_phi->Sumw2();
  res_forward_pt->Sumw2();
  res_forward_eta->Sumw2();
  res_forward_phi->Sumw2();
  res_delta_phi->Sumw2();
  res_delta_eta->Sumw2();
  res_delta_phi_deta1->Sumw2();
  res_delta_phi_deta2->Sumw2();
  res_delta_phi_deta3->Sumw2();
  res_delta_phi_deta4->Sumw2();
  res_delta_phi_gap->Sumw2();
  res_delta_eta_gap->Sumw2();
  res_delta_phi_deta1_gap->Sumw2();
  res_delta_phi_deta2_gap->Sumw2();
  res_delta_phi_deta3_gap->Sumw2();
  res_delta_phi_deta4_gap->Sumw2();
  res_delta_phi_nogap->Sumw2();
  res_delta_eta_nogap->Sumw2();
  res_delta_phi_deta1_nogap->Sumw2();
  res_delta_phi_deta2_nogap->Sumw2();
  res_delta_phi_deta3_nogap->Sumw2();
  res_delta_phi_deta4_nogap->Sumw2();
  res_leading_pt_inside_gap->Sumw2();
  res_leading_eta_inside_gap->Sumw2();
  res_leading_phi_inside_gap->Sumw2();
  res_leading_eta_star_inside_gap->Sumw2();
  res_leading_pt_outside_gap->Sumw2();
  res_leading_eta_outside_gap->Sumw2();
  res_leading_phi_outside_gap->Sumw2();
  res_delta_eta_outside_gap->Sumw2();

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

  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  prescale = data_lumi[z];
  cout<<"Cross-section per event : "<<prescale<<endl;

  int decade = 0;

  for (Int_t i=0;i<nentries;i++)
    {
    counter_entries++; 
    double progress = 10.0*i/(1.0*nentries);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
    cout<<k*10<<" % "<<endl;
    decade = k;

    tr->GetEntry(i);
     if (Event->evtHdr().weight() > 0.0)
	{ weightMC = Event->evtHdr().weight(); }
     else { weightMC = 1.0; }

    //-------- check if the primary vertex is good ----
    pv_pass = false;

    if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)						{ pv_pass = true; }
    if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) 					{ pv_pass = true; }
    if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)	{ pv_pass = true; }

    if (pv_pass == true)
    {
    	counter_pv++;
    	if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1);}

   	gen_selected = false;
	det_selected = false;
	gen_gap = false;
	det_gap = false;
	gen_nogap = false;
	det_nogap = false;
	gen_out = false;
	det_out = false;
	gen_deta1 = false;
	gen_deta2 = false;
	gen_deta3 = false;
	gen_deta4 = false;
	det_deta1 = false;
	det_deta2 = false;
	det_deta3 = false;
	det_deta4 = false;

	gen_central_pt = -1.0;
	gen_central_eta = -10.0;
	gen_central_phi = -1.0;
	gen_forward_pt = -1.0;
	gen_forward_eta = -10.0;
	gen_forward_phi = -1.0;
	gen_inside_pt = -1.0;
	gen_inside_eta = -10.0;
	gen_inside_phi = -1.0;
	gen_inside_eta_star = -10.0;
	gen_outside_pt = -1.0;
	gen_outside_eta = -10.0;
	gen_outside_phi = -1.0;
	gen_delta_eta_outside = -10.0;
	gen_delta_phi = -1.0;
	gen_delta_eta = -1.0;

	det_central_pt = -1.0;
	det_central_eta = -10.0;
	det_central_phi = -1.0;
	det_forward_pt = -1.0;
	det_forward_eta = -10.0;
	det_forward_phi = -1.0;
	det_inside_pt = -1.0;
	det_inside_eta = -10.0;
	det_inside_phi = -1.0;
	det_inside_eta_star = 10.0;
	det_outside_pt = -1.0;
	det_outside_eta = -1.0;
	det_outside_phi = -1.0;
	det_delta_eta_outside = -10.0;
	det_delta_phi = -1.0;
	det_delta_eta = -1.0;

	res_central_pt_val = -1.0;
	res_central_eta_val = -10.0;
	res_central_phi_val = -10.0;
	res_forward_pt_val = -1.0;
	res_forward_eta_val = -10.0;
	res_forward_phi_val = -10.0;
	res_inside_pt_val = -1.0;
	res_inside_eta_val = -10.0;
	res_inside_phi_val = -10.0;
	res_inside_eta_star_val = -10.0;
	res_outside_pt_val = -1.0;
	res_outside_eta_val = -10.0;
	res_outside_phi_val = -10.0;
	res_delta_eta_outside_val = -10.0;
	res_delta_phi_val = -10.0;
	res_delta_eta_val = -10.0;

	//Gen Jets
    	for(unsigned int j=0; j<Event->nGenJets(); j++)
		{
		pt = Event->genjet(j).pt();
		eta = Event->genjet(j).eta();
		phi = Event->genjet(j).phi();
     		if (pt >= pt_min)
			{
     			gen_counter_jet++;
     			if (eta <= 2.8 && eta >= -2.8 && pt > gen_central_pt)
     			{ gen_central_pt = pt; gen_central_eta = eta; gen_central_phi = phi; }
     			if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > gen_forward_pt)
     			{ gen_forward_pt = pt; gen_forward_eta = eta; gen_forward_phi = phi; }
			}   
		}

	if (gen_central_pt > pt_min and gen_forward_pt > pt_min)
		{
   		gen_selected = true;
     		gen_counter_selected++;
     		gen_counter_weigth = gen_counter_weigth + prescale*weightMC*vertex_weighted;
     		gen_delta_eta = gen_forward_eta - gen_central_eta;
     		if (gen_delta_eta < 0) { gen_delta_eta = -gen_delta_eta; }
     		gen_delta_phi = calc_delta_phi(gen_forward_phi, gen_central_phi);
		if (gen_delta_eta >= 0.4 and gen_delta_eta < 2.5) { gen_deta1 = true; }
		if (gen_delta_eta >= 2.5 and gen_delta_eta < 3.5) { gen_deta2 = true; }
		if (gen_delta_eta >= 3.5 and gen_delta_eta < 4.5) { gen_deta3 = true; }
		if (gen_delta_eta >= 4.5 and gen_delta_eta < 7.5) { gen_deta4 = true; }

    		for(unsigned int j=0; j<Event->nGenJets(); j++)
			{
			pt = Event->genjet(j).pt();
			eta = Event->genjet(j).eta();
			phi = Event->genjet(j).phi();
     			if (pt >= gap_req)
				{
     				if (gen_central_eta > gen_forward_eta && eta > gen_forward_eta && eta < gen_central_eta && pt > gen_inside_pt)
					{gen_inside_pt = pt; gen_inside_eta = eta; gen_inside_phi = phi;}
     				if (gen_central_eta < gen_forward_eta && eta < gen_forward_eta && eta > gen_central_eta && pt > gen_inside_pt)
 					{gen_inside_pt = pt; gen_inside_eta = eta; gen_inside_phi = phi;}
     				if (gen_central_eta > gen_forward_eta && (eta < gen_forward_eta || eta > gen_central_eta) && pt > gen_outside_pt)
					{gen_outside_pt = pt; gen_outside_eta = eta; gen_outside_phi = phi; }
     				if (gen_central_eta < gen_forward_eta && (eta > gen_forward_eta || eta < gen_central_eta) && pt > gen_outside_pt) 
					{gen_outside_pt = pt; gen_outside_eta = eta; gen_outside_phi = phi; }
				}   
			}

     		if (gen_inside_pt > gap_req )
     			{
			gen_nogap = true;
			gen_counter_nogap++;
			gen_inside_eta_star = gen_inside_eta - (gen_forward_eta + gen_central_eta)/2;
			}
		else
			{
			gen_counter_gap++;
			gen_gap = true;
			}

     		if (gen_outside_pt > gap_req )
     			{
			gen_out = true;
			gen_counter_outside++;
     			deta_out1 = gen_central_eta - gen_outside_eta;
     			if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     			deta_out2 = gen_forward_eta - gen_outside_eta;
     			if (deta_out2 < 0) { deta_out2 = -deta_out2; }
     			if (deta_out1 < deta_out2) { gen_delta_eta_outside = deta_out1; }
     			if (deta_out2 < deta_out1) { gen_delta_eta_outside = deta_out2; }
			}

		}


	//Det Jets
    	for(unsigned int j=0; j<Event->nPFJets(); j++)
		{
        	pass_tight = Event->pfjet(j).tightID();
		pt = Event->pfjet(j).ptCor();
		eta = Event->pfjet(j).eta();
		phi = Event->pfjet(j).phi();
		gen_pt = Event->genjet(j).pt();
		gen_eta = Event->genjet(j).eta();
		gen_phi = Event->genjet(j).phi();
		pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt);
		if (pt >= pt_min && pass_tight)
			{
     			det_counter_jet++;
     			if (eta <= 2.8 && eta >= -2.8 && pt > det_central_pt)
     			{ det_central_pt = pt; det_central_eta = eta; det_central_phi = phi; }
     			if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > det_forward_pt)
     			{ det_forward_pt = pt; det_forward_eta = eta; det_forward_phi = phi; }
			}
		}

	if (det_central_pt > pt_min and det_forward_pt > pt_min)
		{
		det_selected = true;
     		det_counter_selected++;
     		det_counter_weigth = det_counter_weigth + prescale*weightMC*vertex_weighted;
     		det_delta_eta = det_forward_eta - det_central_eta;
     		if (det_delta_eta < 0) { det_delta_eta = -det_delta_eta; }
     		det_delta_phi = calc_delta_phi(det_forward_phi, det_central_phi);
		if (det_delta_eta >= 0.4 and det_delta_eta < 2.5) { det_deta1 = true; }
		if (det_delta_eta >= 2.5 and det_delta_eta < 3.5) { det_deta2 = true; }
		if (det_delta_eta >= 3.5 and det_delta_eta < 4.5) { det_deta3 = true; }
		if (det_delta_eta >= 4.5 and det_delta_eta < 7.5) { det_deta4 = true; }
    		for(unsigned int j=0; j<Event->nPFJets(); j++)
			{
        		pass_tight = Event->pfjet(j).tightID();
			pt = Event->pfjet(j).ptCor();
			eta = Event->pfjet(j).eta();
			phi = Event->pfjet(j).phi();
			gen_pt = Event->genjet(j).pt();
			gen_eta = Event->genjet(j).eta();
			gen_phi = Event->genjet(j).phi();
			pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt);
			if (pt >= gap_req && pass_tight)
				{
     				if (det_central_eta > det_forward_eta && eta > det_forward_eta && eta < det_central_eta && pt > det_inside_pt)
					{det_inside_pt = pt; det_inside_eta = eta; det_inside_phi = phi;}
     				if (det_central_eta < det_forward_eta && eta < det_forward_eta && eta > det_central_eta && pt > det_inside_pt)
 					{det_inside_pt = pt; det_inside_eta = eta; det_inside_phi = phi;}
     				if (det_central_eta > det_forward_eta && (eta < det_forward_eta || eta > det_central_eta) && pt > det_outside_pt)
					{det_outside_pt = pt; det_outside_eta = eta; det_outside_phi = phi; }
     				if (det_central_eta < det_forward_eta && (eta > det_forward_eta || eta < det_central_eta) && pt > det_outside_pt) 
					{det_outside_pt = pt; det_outside_eta = eta; det_outside_phi = phi; }
				}
			}
     		if (det_inside_pt > gap_req )
     			{
			det_nogap = true;
			det_counter_nogap++;
			det_inside_eta_star = det_inside_eta - (det_forward_eta + det_central_eta)/2;
			}
		else
			{
			det_gap = true;
			det_counter_gap++;
			}

     		if (det_outside_pt > gap_req )
     			{
			det_out = true;
			det_counter_outside++;
     			deta_out1 = det_central_eta - det_outside_eta;
     			if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     			deta_out2 = det_forward_eta - det_outside_eta;
     			if (deta_out2 < 0) { deta_out2 = -deta_out2; }
     			if (deta_out1 < deta_out2) { det_delta_eta_outside = deta_out1; }
     			if (deta_out2 < deta_out1) { det_delta_eta_outside = deta_out2; }
			}
		}

	if (gen_central_pt > 0 and det_central_pt > 0) { res_central_pt_val = (gen_central_pt - det_central_pt) / gen_central_pt; }
	if (gen_central_phi > 0 and det_central_phi > 0) { res_central_phi_val = gen_central_phi - det_central_phi; }
	if (gen_central_eta > -10 and det_central_eta > -10) { res_central_eta_val = gen_central_eta - det_central_eta; }
	if (gen_forward_pt > 0 and det_forward_pt > 0) { res_forward_pt_val = (gen_forward_pt - det_forward_pt) / gen_forward_pt; }
	if (gen_forward_phi > 0 and det_forward_phi > 0) { res_forward_phi_val = gen_forward_phi - det_forward_phi; }
	if (gen_forward_eta > -10 and det_forward_eta > -10) { res_forward_eta_val = gen_forward_eta - det_forward_eta; }
	if (gen_inside_pt > 0 and det_inside_pt > 0) { res_inside_pt_val = (gen_inside_pt - det_inside_pt) / gen_inside_pt; }
	if (gen_inside_phi > 0 and det_inside_phi > 0) { res_inside_phi_val = gen_inside_phi - det_inside_phi; }
	if (gen_inside_eta > -10 and det_inside_eta > -10) { res_inside_eta_val = gen_inside_eta - det_inside_eta; }
	if (gen_inside_eta_star > -10 and det_inside_eta_star > -10)
		{ res_inside_eta_star_val = gen_inside_eta_star - det_inside_eta_star; }
	if (gen_outside_pt > 0 and det_outside_pt > 0) { res_outside_pt_val = (gen_outside_pt - det_outside_pt) / gen_outside_pt; }
	if (gen_outside_phi > 0 and det_outside_phi > 0) { res_outside_phi_val = gen_outside_phi - det_outside_phi; }
	if (gen_outside_eta > -10 and det_outside_eta > -10) { res_outside_eta_val = gen_outside_eta - det_outside_eta; }
	if (gen_delta_eta_outside > -10 and det_delta_eta_outside > -10)
		{ res_delta_eta_outside_val = gen_delta_eta_outside - det_delta_eta_outside; }
	if (gen_selected and det_selected) { res_delta_phi_val = gen_delta_phi - det_delta_phi; }
	if (gen_selected and det_selected) { res_delta_eta_val = gen_delta_eta - det_delta_eta; }

	if (res_central_pt_val > -1.0) { res_central_pt->Fill(res_central_pt_val, prescale*weightMC*vertex_weighted); }
	if (res_central_phi_val > -10.0) { res_central_phi->Fill(res_central_phi_val, prescale*weightMC*vertex_weighted); }
	if (res_central_eta_val > -10.0) { res_central_eta->Fill(res_central_eta_val, prescale*weightMC*vertex_weighted); }
	if (res_forward_pt_val > -1.0) { res_forward_pt->Fill(res_forward_pt_val, prescale*weightMC*vertex_weighted); }
	if (res_forward_phi_val > -10.0) { res_forward_phi->Fill(res_forward_phi_val, prescale*weightMC*vertex_weighted); }
	if (res_forward_eta_val > -10.0) { res_forward_eta->Fill(res_forward_eta_val, prescale*weightMC*vertex_weighted); }
	if (res_inside_pt_val > -1.0) { res_leading_pt_inside_gap->Fill(res_inside_pt_val, prescale*weightMC*vertex_weighted); }
	if (res_inside_phi_val > -10.0) { res_leading_phi_inside_gap->Fill(res_inside_phi_val, prescale*weightMC*vertex_weighted); }
	if (res_inside_eta_val > -10.0) { res_leading_eta_inside_gap->Fill(res_inside_eta_val, prescale*weightMC*vertex_weighted); }
	if (res_inside_eta_star_val > -10.0) { res_leading_eta_star_inside_gap->Fill(res_inside_eta_star_val, prescale*weightMC*vertex_weighted); }
	if (res_outside_pt_val > -1.0) { res_leading_pt_outside_gap->Fill(res_outside_pt_val, prescale*weightMC*vertex_weighted); }
	if (res_outside_phi_val > -10.0) { res_leading_phi_outside_gap->Fill(res_outside_phi_val, prescale*weightMC*vertex_weighted); }
	if (res_outside_eta_val > -10.0) { res_leading_eta_outside_gap->Fill(res_outside_eta_val, prescale*weightMC*vertex_weighted); }
	if (res_delta_eta_outside_val > -10.0) { res_delta_eta_outside_gap->Fill(res_delta_eta_outside_val, prescale*weightMC*vertex_weighted); }
	if (res_delta_phi_val > -10.0) { res_delta_phi->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (res_delta_eta_val > -10.0) { res_delta_eta->Fill(res_delta_eta_val, prescale*weightMC*vertex_weighted); }
	if (det_deta1 and gen_deta1 and res_delta_phi_val > -10.0) { res_delta_phi_deta1->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta2 and gen_deta2 and res_delta_phi_val > -10.0) { res_delta_phi_deta2->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta3 and gen_deta3 and res_delta_phi_val > -10.0) { res_delta_phi_deta3->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta4 and gen_deta4 and res_delta_phi_val > -10.0) { res_delta_phi_deta4->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_gap and gen_gap and res_delta_phi_val > -10.0) { res_delta_phi_gap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_gap and gen_gap and res_delta_eta_val > -10.0) { res_delta_eta_gap->Fill(res_delta_eta_val, prescale*weightMC*vertex_weighted); }
	if (det_deta1 and gen_deta1 and det_gap and gen_gap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta1_gap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta2 and gen_deta2 and det_gap and gen_gap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta2_gap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta3 and gen_deta3 and det_gap and gen_gap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta3_gap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta4 and gen_deta4 and det_gap and gen_gap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta4_gap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_nogap and gen_nogap and res_delta_phi_val > -10.0) { res_delta_phi_nogap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_nogap and gen_nogap and res_delta_eta_val > -10.0) { res_delta_eta_nogap->Fill(res_delta_eta_val, prescale*weightMC*vertex_weighted); }
	if (det_deta1 and gen_deta1 and det_nogap and gen_nogap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta1_nogap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta2 and gen_deta2 and det_nogap and gen_nogap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta2_nogap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta3 and gen_deta3 and det_nogap and gen_nogap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta3_nogap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
	if (det_deta4 and gen_deta4 and det_nogap and gen_nogap and res_delta_phi_val > -10.0)
		{ res_delta_phi_deta4_nogap->Fill(res_delta_phi_val, prescale*weightMC*vertex_weighted); }
    } 

    } 

  cout<<"Events read:                      "<<nentries<<endl; 

}


  
     cout<<endl<<endl;
     cout<<"Events read:                      "<<counter_entries<<endl;
     cout<<"Events after the PV cut:          "<<counter_pv<<endl;
     cout<<" "<<endl;
     cout<<"Generator Level Jets:             "<<gen_counter_jet<<endl;
     cout<<"Generator Selected Events:        "<<gen_counter_selected<<endl;
     cout<<"Generator Cross-Section:          "<<gen_counter_weigth<<endl;
     cout<<"Generator Gap Events:             "<<gen_counter_gap<<endl;
     cout<<"Generator Nogap Events:           "<<gen_counter_nogap<<endl;
     cout<<"Generator Outside Events:         "<<gen_counter_outside<<endl;
     cout<<" "<<endl;
     cout<<"Detector Level Jets:              "<<det_counter_jet<<endl;
     cout<<"Detector Selected Events:         "<<det_counter_selected<<endl;
     cout<<"Detector Cross-Section:           "<<det_counter_weigth<<endl;
     cout<<"Detector Gap Events:              "<<det_counter_gap<<endl;
     cout<<"Detector Nogap Events:            "<<det_counter_nogap<<endl;
     cout<<"Detector Outside Events:          "<<det_counter_outside<<endl;
     cout<<" "<<endl;

     //Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     //save histograms on the file
     res_central_pt->Write();
     res_central_eta->Write();
     res_central_phi->Write();
     res_forward_pt->Write();
     res_forward_eta->Write();
     res_forward_phi->Write();
     res_delta_phi->Write();
     res_delta_eta->Write();
     res_delta_phi_deta1->Write();
     res_delta_phi_deta2->Write();
     res_delta_phi_deta3->Write();
     res_delta_phi_deta4->Write();
     res_delta_phi_gap->Write();
     res_delta_eta_gap->Write();
     res_delta_phi_deta1_gap->Write();
     res_delta_phi_deta2_gap->Write();
     res_delta_phi_deta3_gap->Write();
     res_delta_phi_deta4_gap->Write();
     res_delta_phi_nogap->Write();
     res_delta_eta_nogap->Write();
     res_delta_phi_deta1_nogap->Write();
     res_delta_phi_deta2_nogap->Write();
     res_delta_phi_deta3_nogap->Write();
     res_delta_phi_deta4_nogap->Write();
     res_leading_pt_inside_gap->Write();
     res_leading_eta_inside_gap->Write();
     res_leading_phi_inside_gap->Write();
     res_leading_eta_star_inside_gap->Write();
     res_leading_pt_outside_gap->Write();
     res_leading_eta_outside_gap->Write();
     res_leading_phi_outside_gap->Write();
     res_delta_eta_outside_gap->Write();

     //close the output file
     data_output->Close();

     //delete the histograms to avoid memory leak
     delete(res_central_pt);
     delete(res_central_eta);
     delete(res_central_phi);
     delete(res_forward_pt);
     delete(res_forward_eta);
     delete(res_forward_phi);
     delete(res_delta_phi);
     delete(res_delta_eta);
     delete(res_delta_phi_deta1);
     delete(res_delta_phi_deta2);
     delete(res_delta_phi_deta3);
     delete(res_delta_phi_deta4);
     delete(res_delta_phi_gap);
     delete(res_delta_eta_gap);
     delete(res_delta_phi_deta1_gap);
     delete(res_delta_phi_deta2_gap);
     delete(res_delta_phi_deta3_gap);
     delete(res_delta_phi_deta4_gap);
     delete(res_delta_phi_nogap);
     delete(res_delta_eta_nogap);
     delete(res_delta_phi_deta1_nogap);
     delete(res_delta_phi_deta2_nogap);
     delete(res_delta_phi_deta3_nogap);
     delete(res_delta_phi_deta4_nogap);
     delete(res_leading_pt_inside_gap);
     delete(res_leading_eta_inside_gap);
     delete(res_leading_phi_inside_gap);
     delete(res_leading_eta_star_inside_gap);
     delete(res_leading_pt_outside_gap);
     delete(res_leading_eta_outside_gap);
     delete(res_leading_phi_outside_gap);
     delete(res_delta_eta_outside_gap);

}
