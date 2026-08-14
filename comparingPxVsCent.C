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

void comparingPxVsCent(){
    TFile *inputFirst = new TFile("../output/output_PxVsCent_forPsi2_ConstSine_14p5.root", "read");
    TFile *inputSecond = new TFile("../output/output_PxVsCent_forPsi2_SinInclude_14p5.root", "read");
    TFile *inputThird = new TFile("../output/output_PxVsCent_forPsi3_ConstSine_14p5.root", "read");
    TFile *inputFourth = new TFile("../output/output_PxVsCent_forPsi3_SinInclude_14p5.root", "read");

    TGraphErrors *PxVsCentFirst = (TGraphErrors*)inputFirst->Get("PxVsCent_Tog_1");
    TGraphErrors *PxVsCentSecond = (TGraphErrors*)inputSecond->Get("PxVsCent_Tog_1");
    TGraphErrors *PxVsCentThird = (TGraphErrors*)inputThird->Get("PxVsCent_Tog_1");
    TGraphErrors *PxVsCentFourth = (TGraphErrors*)inputFourth->Get("PxVsCent_Tog_1");

    PxVsCentFirst->SetMarkerSize(3);
    PxVsCentFirst->SetLineWidth(2);
    PxVsCentFirst->SetMarkerStyle(34);
    PxVsCentFirst->SetMarkerColor(1);
    PxVsCentFirst->SetLineColor(1);

    PxVsCentSecond->SetMarkerSize(3);
    PxVsCentSecond->SetLineWidth(2);
    PxVsCentSecond->SetMarkerStyle(75);
    PxVsCentSecond->SetMarkerColor(804);
    PxVsCentSecond->SetLineColor(804);

    for(int iCent = 0; iCent!=PxVsCentSecond->GetN(); iCent++){
        double x = PxVsCentSecond->GetPointX(iCent);
        PxVsCentSecond->SetPointX(iCent, x+0.5);
    }

    PxVsCentThird->SetMarkerSize(3);
    PxVsCentThird->SetLineWidth(2);
    PxVsCentThird->SetMarkerStyle(30);
    PxVsCentThird->SetMarkerColor(602);
    PxVsCentThird->SetLineColor(602);

    for(int iCent = 0; iCent!=PxVsCentSecond->GetN(); iCent++){
        double x = PxVsCentThird->GetPointX(iCent);
        PxVsCentThird->SetPointX(iCent, x+0.5);
    }

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);

    TCanvas *c_sinconst_sininc = new TCanvas("c_sinconst_sininc", "", 0, 1024, 1200, 1024);
    TCanvas *c_psi2_psi3 = new TCanvas("c_psi2_psi3", "", 0, 1024, 1200, 1024);

    c_sinconst_sininc->cd();
    TMultiGraph *PxVsCent_First_Second = new TMultiGraph();

    PxVsCent_First_Second->Add(PxVsCentFirst);
    PxVsCent_First_Second->Add(PxVsCentSecond);

    PxVsCent_First_Second->GetXaxis()->SetTitle("centrality %");
    PxVsCent_First_Second->GetYaxis()->SetTitle("<P_{x}sin(2#phi - 2#Psi_{2})> %");

    TLegend *leg_first_second = new TLegend();

    leg_first_second->AddEntry((TObject*) 0, "P_{x} of #Lambda + #bar{#Lambda}, #Psi_{2,comb} used", "");
    leg_first_second->AddEntry(PxVsCentFirst, "<cos(#Psi_{1} - #phi)> extracted", "p");
    leg_first_second->AddEntry(PxVsCentSecond, "<cos(#Psi_{1} - #phi)sin(#theta)> extracted", "p");
    leg_first_second->SetBorderSize(0);

    PxVsCent_First_Second->Draw("AP");
    funcOfZero->Draw("same");
    leg_first_second->Draw("same");

    //
    c_psi2_psi3->cd();

    TMultiGraph *PxVsCent_First_Third = new TMultiGraph();

    PxVsCent_First_Third->Add(PxVsCentFirst);
    PxVsCent_First_Third->Add(PxVsCentThird);

    PxVsCent_First_Third->GetXaxis()->SetTitle("centrality %");
    PxVsCent_First_Third->GetYaxis()->SetTitle("<P_{x}sin(n(#phi - #Psi_{n}))> %");
    PxVsCent_First_Third->GetXaxis()->SetLimits(0, 65);
    PxVsCent_First_Third->SetMinimum(-1.5);
    PxVsCent_First_Third->SetMaximum(1.5);

    TLegend *leg_first_Third = new TLegend();

    leg_first_Third->AddEntry((TObject*) 0, "P_{x} of #Lambda + #bar{#Lambda}, #Psi_{n,comb} used", "");
    leg_first_Third->AddEntry(PxVsCentFirst, "n=2", "p");
    leg_first_Third->AddEntry(PxVsCentThird, "n=3", "p");
    leg_first_Third->SetBorderSize(0);

    PxVsCent_First_Third->Draw("AP");
    funcOfZero->Draw("same");
    leg_first_Third->Draw("same");
}