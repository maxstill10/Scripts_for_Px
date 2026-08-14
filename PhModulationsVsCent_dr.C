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

TGraphErrors* MakeAverageBin(TGraphErrors *g, double n_i, double n_f);

void PhModulationsVsCent_dr(){
    TFile *inputFirst = new TFile("../output/../output/output_PhModulationsVsCent_forPsi2_ConstSine_14p5.root", "read");
    TFile *inputSecond = new TFile("../output/../output/output_PhModulationsVsCent_forPsi3_ConstSine_14p5.root", "read");

    TGraphErrors *PxVsCentFirst = (TGraphErrors*)inputFirst->Get("PhCosVsCent_Tog_1");
    TGraphErrors *PxVsCentSecond = (TGraphErrors*)inputSecond->Get("PhCosVsCent_Tog_1");
    TGraphErrors *PxVsCentThird = (TGraphErrors*)inputFirst->Get("PhSinVsCent_Tog_1");
    TGraphErrors *PxVsCentFourth = (TGraphErrors*)inputSecond->Get("PhSinVsCent_Tog_1");

    //PxVsCentFourth = MakeAverageBin(PxVsCentFourth, 7, 8);

    for(int iCent = 0; iCent!=PxVsCentSecond->GetN(); iCent++){
        double x = PxVsCentSecond->GetPointX(iCent);
        PxVsCentSecond->SetPointX(iCent, x+0.5);
    }

    for(int iCent = 0; iCent!=PxVsCentFourth->GetN(); iCent++){
        double x = PxVsCentFourth->GetPointX(iCent);
        PxVsCentFourth->SetPointX(iCent, x+0.5);
    }

    PxVsCentFirst->SetMarkerSize(2);
    PxVsCentFirst->SetLineWidth(2);
    PxVsCentFirst->SetMarkerStyle(34);
    PxVsCentFirst->SetMarkerColor(1);
    PxVsCentFirst->SetLineColor(1);

    PxVsCentSecond->SetMarkerSize(2);
    PxVsCentSecond->SetLineWidth(2);
    PxVsCentSecond->SetMarkerStyle(28);
    PxVsCentSecond->SetMarkerColor(804);
    PxVsCentSecond->SetLineColor(804);

    PxVsCentThird->SetMarkerSize(2);
    PxVsCentThird->SetLineWidth(2);
    PxVsCentThird->SetMarkerStyle(34);
    PxVsCentThird->SetMarkerColor(1);
    PxVsCentThird->SetLineColor(1);

    PxVsCentFourth->SetMarkerSize(2);
    PxVsCentFourth->SetLineWidth(2);
    PxVsCentFourth->SetMarkerStyle(28);
    PxVsCentFourth->SetMarkerColor(804);
    PxVsCentFourth->SetLineColor(804);

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);

    //average polar calc
    TGraphErrors *PxFirst = MakeAverageBin(PxVsCentFirst, 4, 6);
    TGraphErrors *PxSecond = MakeAverageBin(PxVsCentSecond, 4, 6);
    TGraphErrors *PxThird = MakeAverageBin(PxVsCentThird, 4, 6);
    TGraphErrors *PxFourth = MakeAverageBin(PxVsCentFourth, 4, 6);

    //Start drawing
    TCanvas *c_px_pz = new TCanvas("c_px_pz", "", 0, 1024, 800, 1024);
    c_px_pz->Divide(1, 2, 0, 0);

    c_px_pz->cd(1);
    gPad->SetBottomMargin(0);
    gPad->SetLeftMargin(0.2);
    TMultiGraph *PxVsCent_Third_Fourth = new TMultiGraph();

    PxVsCent_Third_Fourth->Add(PxVsCentThird);
    PxVsCent_Third_Fourth->Add(PxVsCentFourth);

    PxVsCent_Third_Fourth->GetXaxis()->SetTitle("centrality %");
    PxVsCent_Third_Fourth->GetYaxis()->SetTitle("<P_{h}sin(n(#phi - #Psi_{n}))> %");
    PxVsCent_Third_Fourth->GetYaxis()->SetTitleSize(0.06);
    PxVsCent_Third_Fourth->GetYaxis()->SetLabelSize(0.06);
    PxVsCent_Third_Fourth->GetYaxis()->SetTitleOffset(1);
    PxVsCent_Third_Fourth->SetMaximum(1.7);
    PxVsCent_Third_Fourth->SetMinimum(-1.7);

    TLegend *leg_Third_Fourth = new TLegend();

    leg_Third_Fourth->AddEntry((TObject*) 0, "P_{h} of #Lambda + #bar{#Lambda}, #Psi_{n,comb} used", "");
    leg_Third_Fourth->AddEntry(PxVsCentThird, "n = 2", "p");
    leg_Third_Fourth->AddEntry(PxVsCentFourth, "n = 3", "p");
    leg_Third_Fourth->SetBorderSize(0);

    PxVsCent_Third_Fourth->Draw("AP");
    funcOfZero->Draw("same");
    leg_Third_Fourth->Draw("same");
    
    //
    c_px_pz->cd(2);
    gPad->SetTopMargin(0);
    gPad->SetLeftMargin(0.2);
    TMultiGraph *PxVsCent_First_Second = new TMultiGraph();

    PxVsCent_First_Second->Add(PxVsCentFirst);
    PxVsCent_First_Second->Add(PxVsCentSecond);

    PxVsCent_First_Second->GetXaxis()->SetTitle("centrality %");
    PxVsCent_First_Second->GetYaxis()->SetTitle("<P_{h}cos(n(#phi - #Psi_{n}))> %");
    PxVsCent_First_Second->GetYaxis()->SetTitleSize(0.06);
    PxVsCent_First_Second->GetYaxis()->SetLabelSize(0.06);
    PxVsCent_First_Second->GetYaxis()->SetTitleOffset(1);
    PxVsCent_First_Second->SetMaximum(1.7);
    PxVsCent_First_Second->SetMinimum(-1.7);

    TLegend *leg_first_second = new TLegend();

    leg_first_second->AddEntry((TObject*) 0, "P_{h} of #Lambda + #bar{#Lambda}, #Psi_{n,comb} used", "");
    leg_first_second->AddEntry(PxVsCentFirst, "n = 2", "p");
    leg_first_second->AddEntry(PxVsCentSecond, "n = 3", "p");
    leg_first_second->SetBorderSize(0);

    PxVsCent_First_Second->Draw("AP");
    funcOfZero->Draw("same");
    leg_first_second->Draw("same");
}

//n_i = number of the first point
//n_f = number of the last point
//number starts with 1
TGraphErrors* MakeAverageBin(TGraphErrors *g, double n_i, double n_f){
    TGraphErrors* gOut = new TGraphErrors(g->GetN() - n_f + n_i);

    double x, ex, y, ey, w, sum_w = 0., sum_x = 0., sum_y = 0., sum_y_err = 0.;
    for(int i = 0; i!=n_i - 1; i++){
        x = g->GetPointX(i);
        y = g->GetPointY(i);
        ex = g->GetErrorX(i);
        ey = g->GetErrorY(i);

        gOut->SetPoint(i, x, y);
        gOut->SetPointError(i, ex, ey);
    }

    for(int i = n_f; i!=g->GetN(); i++){
        x = g->GetPointX(i);
        y = g->GetPointY(i);
        ex = g->GetErrorX(i);
        ey = g->GetErrorY(i);

        gOut->SetPoint(i - n_f + n_i, x, y);
        gOut->SetPointError(i - n_f + n_i, ex, ey);
    }

    for(int i = n_i - 1; i!=n_f; i++){
        w = 1/pow(g->GetErrorY(i), 2);
        sum_w += w;
        sum_y += w*g->GetPointY(i);
        sum_y_err += w*w*pow(g->GetErrorY(i), 2);
        sum_x += g->GetPointX(i);
    }

    gOut->SetPoint(n_i - 1, sum_x/(n_f - n_i + 1), sum_y/sum_w);
    gOut->SetPointError(n_i - 1, 0., sqrt(sum_y_err)/sum_w);

    return gOut;
}