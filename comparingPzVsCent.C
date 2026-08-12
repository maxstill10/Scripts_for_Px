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

void comparingPzVsCent(){
    TFile *inputFirst = new TFile("../output/output_PzVsCent_forPsi2_14p5.root", "read");
    TFile *inputThird = new TFile("../output/output_PzVsCent_forPsi3_14p5.root", "read");

    TGraphErrors *PxVsCentFirst = (TGraphErrors*)inputFirst->Get("PzVsCent_Tog_0");
    TGraphErrors *PxVsCentThird = (TGraphErrors*)inputThird->Get("PzVsCent_Tog_0");

    PxVsCentFirst->SetMarkerSize(3);
    PxVsCentFirst->SetLineWidth(2);
    PxVsCentFirst->SetMarkerStyle(29);
    PxVsCentFirst->SetMarkerColor(634);
    PxVsCentFirst->SetLineColor(634);


    PxVsCentThird->SetMarkerSize(3);
    PxVsCentThird->SetLineWidth(2);
    PxVsCentThird->SetMarkerStyle(30);
    PxVsCentThird->SetMarkerColor(602);
    PxVsCentThird->SetLineColor(602);

    for(int iCent = 0; iCent!=PxVsCentThird->GetN(); iCent++){
        double x = PxVsCentThird->GetPointX(iCent);
        PxVsCentThird->SetPointX(iCent, x+0.5);
    }

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);

    TCanvas *c_psi2_psi3 = new TCanvas("c_psi2_psi3", "", 0, 1024, 1200, 1024);

    //
    c_psi2_psi3->cd();

    TMultiGraph *PxVsCent_First_Third = new TMultiGraph();

    PxVsCent_First_Third->Add(PxVsCentFirst);
    PxVsCent_First_Third->Add(PxVsCentThird);

    PxVsCent_First_Third->GetXaxis()->SetTitle("centrality %");
    PxVsCent_First_Third->GetYaxis()->SetTitle("<P_{z}sin(n(#phi - #Psi_{n}))> %");
    PxVsCent_First_Third->GetXaxis()->SetLimits(0, 60);
    PxVsCent_First_Third->SetMinimum(-1.4);
    PxVsCent_First_Third->SetMaximum(1.4);

    TLegend *leg_first_Third = new TLegend();

    leg_first_Third->AddEntry((TObject*) 0, "P_{z} of #Lambda + #bar{#Lambda}, #Psi_{n,ew} used", "");
    leg_first_Third->AddEntry(PxVsCentFirst, "n=2", "p");
    leg_first_Third->AddEntry(PxVsCentThird, "n=3", "p");
    leg_first_Third->SetBorderSize(0);

    PxVsCent_First_Third->Draw("AP");
    funcOfZero->Draw("same");
    leg_first_Third->Draw("same");
}