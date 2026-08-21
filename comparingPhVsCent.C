// This is needed for calling standalone classes (not needed on RACF)
#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <ctime>
#include <cmath>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"

TGraphErrors* MakeAverage( TGraphErrors *g1, TGraphErrors *g2 );

void comparingPhVsCent(){
    TFile *inputFirst = new TFile("../output/output_PhVsCent_SinInclude_14p5.root", "read");
    TFile *inputSecond = new TFile("../output/output_PhVsCent_ConstSine_14p5.root", "read");
    TFile *inputThird = new TFile("../output/output_PhModulationsVsCent_forPsi2_ConstSine_14p5.root", "read");
    

    TGraphErrors *PxVsCentLamFirst = (TGraphErrors*)inputFirst->Get("PhVsCent_Lam");
    TGraphErrors *PxVsCentLamBarFirst = (TGraphErrors*)inputFirst->Get("PhVsCent_LamBar");
    TGraphErrors *PxVsCentLamSecond = (TGraphErrors*)inputSecond->Get("PhVsCent_Lam");    
    TGraphErrors *PxVsCentLamBarSecond = (TGraphErrors*)inputSecond->Get("PhVsCent_LamBar");
    TGraphErrors *PxVsCentLamThird = (TGraphErrors*)inputThird->Get("PhVsCent_Lam_1");
    TGraphErrors *PxVsCentLamBarThird = (TGraphErrors*)inputThird->Get("PhVsCent_LamBar_1");

    TGraphErrors *PxVsCentFirst = MakeAverage(PxVsCentLamFirst, PxVsCentLamBarFirst);
    TGraphErrors *PxVsCentSecond = MakeAverage(PxVsCentLamSecond, PxVsCentLamBarSecond);
    TGraphErrors *PxVsCentThird = MakeAverage(PxVsCentLamThird, PxVsCentLamBarThird);

    PxVsCentFirst->SetMarkerSize(3);
    PxVsCentFirst->SetLineWidth(2);
    PxVsCentFirst->SetMarkerStyle(29);
    PxVsCentFirst->SetMarkerColor(634);
    PxVsCentFirst->SetLineColor(634);

    PxVsCentSecond->SetMarkerSize(3);
    PxVsCentSecond->SetLineWidth(2);
    PxVsCentSecond->SetMarkerStyle(30);
    PxVsCentSecond->SetMarkerColor(602);
    PxVsCentSecond->SetLineColor(602);

    PxVsCentLamSecond->SetMarkerSize(3);
    PxVsCentLamSecond->SetLineWidth(2);
    PxVsCentLamSecond->SetMarkerStyle(29);
    PxVsCentLamSecond->SetMarkerColor(634);
    PxVsCentLamSecond->SetLineColor(634);

    PxVsCentLamBarSecond->SetMarkerSize(3);
    PxVsCentLamBarSecond->SetLineWidth(2);
    PxVsCentLamBarSecond->SetMarkerStyle(76);
    PxVsCentLamBarSecond->SetMarkerColor(602);
    PxVsCentLamBarSecond->SetLineColor(602);

    PxVsCentLamThird->SetMarkerSize(3);
    PxVsCentLamThird->SetLineWidth(2);
    PxVsCentLamThird->SetMarkerStyle(34);
    PxVsCentLamThird->SetMarkerColor(635);
    PxVsCentLamThird->SetLineColor(635);

    PxVsCentLamBarThird->SetMarkerSize(3);
    PxVsCentLamBarThird->SetLineWidth(2);
    PxVsCentLamBarThird->SetMarkerStyle(75);
    PxVsCentLamBarThird->SetMarkerColor(603);
    PxVsCentLamBarThird->SetLineColor(603);

    for(int iCent = 0; iCent!=PxVsCentSecond->GetN(); iCent++){
        double x = PxVsCentSecond->GetPointX(iCent);
        PxVsCentSecond->SetPointX(iCent, x+0.5);
    }

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);

    TCanvas *c_sinconst_sininc = new TCanvas("c_sinconst_sininc", "", 0, 1024, 1200, 1024);
    
    c_sinconst_sininc->cd();
    TMultiGraph *PxVsCent_First_Second = new TMultiGraph();

    PxVsCent_First_Second->Add(PxVsCentFirst);
    PxVsCent_First_Second->Add(PxVsCentSecond);

    PxVsCent_First_Second->GetXaxis()->SetTitle("centrality %");
    PxVsCent_First_Second->GetYaxis()->SetTitle("<P_{x}sin(2#phi - 2#Psi_{2})> %");
    PxVsCent_First_Second->SetMinimum(-0.4);
    PxVsCent_First_Second->SetMaximum(2.4);

    TLegend *leg_first_second = new TLegend();

    leg_first_second->AddEntry((TObject*) 0, "P_{h} of #Lambda + #bar{#Lambda}, #Psi_{2,comb} used", "");
    leg_first_second->AddEntry(PxVsCentFirst, "<cos(#Psi_{1} - #phi)> extracted", "p");
    leg_first_second->AddEntry(PxVsCentSecond, "<cos(#Psi_{1} - #phi)sin(#theta)> extracted", "p");
    leg_first_second->SetBorderSize(0);

    PxVsCent_First_Second->Draw("AP");
    funcOfZero->Draw("same");
    leg_first_second->Draw("same");

    //
    TCanvas *c_Second_Third = new TCanvas("c_Second_Third", "", 0, 1024, 1200, 1024);
    
    c_Second_Third->cd();
    TMultiGraph *PxVsCent_Second_Third = new TMultiGraph();

    PxVsCent_Second_Third->Add(PxVsCentLamSecond);
    PxVsCent_Second_Third->Add(PxVsCentLamBarSecond);
    PxVsCent_Second_Third->Add(PxVsCentLamThird);
    PxVsCent_Second_Third->Add(PxVsCentLamBarThird);

    PxVsCent_Second_Third->GetXaxis()->SetTitle("centrality %");
    PxVsCent_Second_Third->GetYaxis()->SetTitle("<P_{h}> %");
    PxVsCent_Second_Third->SetMinimum(-0.4);
    PxVsCent_Second_Third->SetMaximum(2.4);

    TLegend *leg_Second_Third = new TLegend();

    leg_Second_Third->AddEntry((TObject*) 0, "P_{h} of #Lambda + #bar{#Lambda}, #Psi_{2,comb} used", "");
    leg_Second_Third->AddEntry(PxVsCentLamSecond, "#Lambda, Ph average", "p");
    leg_Second_Third->AddEntry(PxVsCentLamBarSecond, "#bar{#Lambda} Ph average", "p");
    leg_Second_Third->AddEntry(PxVsCentLamThird, "#Lambda, Ph extracter", "p");
    leg_Second_Third->AddEntry(PxVsCentLamBarThird, "#bar{#Lambda} Ph extracted", "p");
    leg_Second_Third->SetBorderSize(0);

    PxVsCent_Second_Third->Draw("AP");
    funcOfZero->Draw("same");
    leg_Second_Third->Draw("same");

    //

    TMultiGraph *MgPolar_HVsCent = new TMultiGraph();
    
    MgPolar_HVsCent->Add(P_HVsCent);
    MgPolar_HVsCent->Add(P_HVsCentEgorik);
    MgPolar_HVsCent->Add(P_HVsCentChinese);
    
    MgPolar_HVsCent->SetTitle("<P_{H}> of #Lambda Vs Centrality");
    MgPolar_HVsCent->GetXaxis()->SetTitle("Centrality, %");
    MgPolar_HVsCent->GetYaxis()->SetTitle("<P_{H}>, %");
    

    MgPolar_HVsCent->Draw("AP");
    funcOfZero->Draw("same");

    legend->AddEntry(P_HVsCent, "This research", "p");
    legend->AddEntry(P_HVsCentEgorik, "Egor research", "p");
    legend->AddEntry(P_HVsCentChinese, "Chinese research", "p");
    legend->Draw("same");

    TCanvas *can1 = new TCanvas("can1", "can1", 0, 1024, 1800, 1024);
    can1->cd();

    TLegend *legendBar = new TLegend();
    TMultiGraph *MgPolar_HBarVsCent = new TMultiGraph();
    
    MgPolar_HBarVsCent->Add(P_HVsCentBar);
    MgPolar_HBarVsCent->Add(P_HVsCentBarEgorik);
    MgPolar_HBarVsCent->Add(P_HVsCentBarChinese);
    
    MgPolar_HBarVsCent->SetTitle("<P_{H}> of #bar{#Lambda} Vs Centrality");
    MgPolar_HBarVsCent->GetXaxis()->SetTitle("Centrality, %");
    MgPolar_HBarVsCent->GetYaxis()->SetTitle("<P_{H}>, %");
    

    MgPolar_HBarVsCent->Draw("AP");
    funcOfZero->Draw("same");

    legendBar->AddEntry(P_HVsCentBar, "This research", "p");
    legendBar->AddEntry(P_HVsCentBarEgorik, "Egor research", "p");
    legendBar->AddEntry(P_HVsCentBarChinese, "Chinese research", "p");
    legendBar->Draw("same");


    
}

TGraphErrors* MakeAverage( TGraphErrors *g1, TGraphErrors *g2 ){

	if( !g1 || !g2 ){ cout <<"MakeAverage(): Input data is NULL!!"<< endl; return NULL; }
	//if( g1->GetN()!=g2->GetN() || g1->GetN()!=gOut->GetN() ){ cout <<"Avergae(): Input graphs have different bins!!"<< endl; return; }

	TGraphErrors *gOut = new TGraphErrors( g1->GetN() );
	for( int ix=0; ix<g1->GetN(); ix++ ){

		float er1 = g1->GetEY()[ix];
		float er2 = g2->GetEY()[ix];
		float w1 = fabs(1.0 / er1 / er1);
		float w2 = fabs(1.0 / er2 / er2);

		float val =  ( w1*g1->GetY()[ix] + w2*g2->GetY()[ix]) / (w1+w2);
		float err = sqrt( w1*w1*er1*er1 + w2*w2*er2*er2 ) / (w1+w2);
		gOut->SetPoint( ix, g1->GetX()[ix], val );
		gOut->SetPointError( ix, 0.0, err );
	}
    gOut->SetMarkerSize(3);
    gOut->SetLineWidth(2);
	return gOut;
}