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

void comp_cossin_cos(){
    TFile *input = new TFile("../output/output_PxVsCent_ConstSine_14p5.root", "read");
    TFile *input_withsin = new TFile("../output/output_PxVsCent_SinInclude_14p5.root", "read");

    TGraphErrors *PhVsCent_Lam = (TGraphErrors*)input->Get("PxVsCent_Lam_0");
    TGraphErrors *PhVsCent_LamBar = (TGraphErrors*)input->Get("PxVsCent_LamBar_0");

    TGraphErrors *PhVsCent_Lam_sin = (TGraphErrors*)input_withsin->Get("PxVsCent_Lam_0");
    TGraphErrors *PhVsCent_LamBar_sin = (TGraphErrors*)input_withsin->Get("PxVsCent_LamBar_0");

    PhVsCent_Lam_sin->SetMarkerStyle(30);
    PhVsCent_Lam_sin->SetMarkerColor(602);
    PhVsCent_Lam_sin->SetMarkerSize(3);
    PhVsCent_Lam_sin->SetLineWidth(2);
    PhVsCent_Lam->SetMarkerStyle(29);
    PhVsCent_Lam->SetMarkerColor(634);
    PhVsCent_Lam->SetMarkerSize(3);
    PhVsCent_Lam->SetLineWidth(2);

    PhVsCent_LamBar_sin->SetMarkerStyle(30);
    PhVsCent_LamBar_sin->SetMarkerColor(602);
    PhVsCent_LamBar_sin->SetMarkerSize(3);
    PhVsCent_LamBar_sin->SetLineWidth(2);
    PhVsCent_LamBar->SetMarkerStyle(29);
    PhVsCent_LamBar->SetMarkerColor(634);
    PhVsCent_LamBar->SetMarkerSize(3);
    PhVsCent_LamBar->SetLineWidth(2);


    TMultiGraph *PzVswCent_Lam = new TMultiGraph();

    PzVswCent_Lam->Add(PhVsCent_Lam);
    PzVswCent_Lam->Add(PhVsCent_Lam_sin);

    TMultiGraph *PzVswCent_LamBar = new TMultiGraph();

    PzVswCent_LamBar->Add(PhVsCent_LamBar);
    PzVswCent_LamBar->Add(PhVsCent_LamBar_sin);

    TLegend *leg = new TLegend();

    leg->AddEntry(PhVsCent_Lam, "Sin is const", "p");
    leg->AddEntry(PhVsCent_Lam_sin, "Sin included", "p");

    TCanvas *c_lam = new TCanvas("c_lam", "c_lam", 0, 1024, 1200, 1024);
    TCanvas *c_lamBar = new TCanvas("c_lamBar", "c_lamBar", 0, 1024, 1200, 1024);

    c_lam->cd();
    PzVswCent_Lam->Draw("AP");
    leg->Draw("same");

    c_lamBar->cd();
    PzVswCent_LamBar->Draw("AP");
    leg->Draw("same");
}