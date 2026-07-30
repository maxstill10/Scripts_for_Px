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

const Double_t alpha_Lam = 0.732;
const Double_t alpha_LamBar = -0.732;

TGraphErrors* ConvertProfToGraph(TProfile* prof);

void comparingCosSinDep(){
    TFile *input = new TFile("../Files/14p5GeV/output_Polar_14p5.root", "read");

    TH1F *CosOfDiffFlatt_1 = (TH1F*)input->Get("CosOfDiff_1");

    TGraphErrors *Res_1VsCent[2];

    for(int i=0; i!=2; i++){
        Res_1VsCent[i] = new TGraphErrors(9);
    }

    double resFlatt, resErr;
    double cent[9] = {1.5, 6.5, 14, 24, 34, 44, 54, 64, 74};
    
    for(int iCen=0; iCen!=9;iCen++){
        
        resFlatt = pow(2*TMath::Abs(CosOfDiffFlatt_1->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(resFlatt))*CosOfDiffFlatt_1->GetBinError(iCen+1);
        
        Res_1VsCent[0]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_1VsCent[0]->SetPointError(8-iCen, 0., resErr);

        resFlatt = pow(TMath::Abs(CosOfDiffFlatt_1->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(2*resFlatt))*CosOfDiffFlatt_1->GetBinError(iCen+1);
        
        Res_1VsCent[1]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_1VsCent[1]->SetPointError(8-iCen, 0., resErr);
        
    } 

    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(1);

    TProfile *prCos_diffPhiPsi1[9][2];
    TProfile *prCos_diffPhiPsi1Sin_theta[9][2];

    TProfile *prCos_diffPhiPsi1_LamBar[9][2];
    TProfile *prCos_diffPhiPsi1Sin_theta_LamBar[9][2];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iSub=0; iSub!=2; iSub++){
            prCos_diffPhiPsi1[iCent][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1_%i_%i", iCent, iSub));

            prCos_diffPhiPsi1_LamBar[iCent][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1_LamBar_%i_%i", iCent, iSub));
        }
    }

    TGraphErrors *grCos_diffPhiPsi1[9][2];

    TGraphErrors *grCos_diffPhiPsi1_LamBar[9][2];

    /*TCanvas *c[24];
    for(int i=0; i!=24; i++){
        c[i] = new TCanvas(Form("c_%i", i), Form("c_%i", i), 0, 1024, 1200, 1024);
        c[i]->Divide(3,3);
    }*/

    TGaxis::SetMaxDigits(1);


    int cent_lim[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

    /*for(int iCent=0; iCent!=9; iCent++){
        for(int iSub=0; iSub!=2; iSub++){
            
            grSin_diffPhiPsi1[iCent][iSub] = ConvertProfToGraph(prSin_diffPhiPsi1[iCent][iSub]);
            grCos_diffPhiPsi1[iCent][iSub] = ConvertProfToGraph(prCos_diffPhiPsi1[iCent][iSub]);
            grCos_theta[iCent][iSub] = ConvertProfToGraph(prCos_theta[iCent][iSub]);
            grSin_theta[iCent][iSub] = ConvertProfToGraph(prSin_theta[iCent][iSub]);
            grSin_diffPhiPsi1Sin_theta[iCent][iSub] = ConvertProfToGraph(prSin_diffPhiPsi1Sin_theta[iCent][iSub]);
            grCos_diffPhiPsi1Sin_theta[iCent][iSub] = ConvertProfToGraph(prCos_diffPhiPsi1Sin_theta[iCent][iSub]);

            grSin_diffPhiPsi1_LamBar[iCent][iSub] = ConvertProfToGraph(prSin_diffPhiPsi1_LamBar[iCent][iSub]);
            grCos_diffPhiPsi1_LamBar[iCent][iSub] = ConvertProfToGraph(prCos_diffPhiPsi1_LamBar[iCent][iSub]);
            grCos_theta_LamBar[iCent][iSub] = ConvertProfToGraph(prCos_theta_LamBar[iCent][iSub]);
            grSin_theta_LamBar[iCent][iSub] = ConvertProfToGraph(prSin_theta_LamBar[iCent][iSub]);
            grSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub] = ConvertProfToGraph(prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]);
            grCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub] = ConvertProfToGraph(prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]);

            //Set Y titles
            grSin_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitle("<sin(#phi_{#Lambda}-#Psi_{1})>");
            grSin_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            grCos_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})>");
            grCos_theta[iCent][iSub]->GetYaxis()->SetTitle("<cos(#theta)>");
            grSin_theta[iCent][iSub]->GetYaxis()->SetTitle("<sin(#theta)>");
            grSin_diffPhiPsi1Sin_theta[iCent][iSub]->GetYaxis()->SetTitle("<sin(#phi_{#Lambda}-#Psi_{1})sin(#theta)>");
            grCos_diffPhiPsi1Sin_theta[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})sin(#theta)>");

            grSin_diffPhiPsi1_LamBar[iCent][iSub]->GetYaxis()->SetTitle("<sin(#phi_{#bar{#Lambda}}-#Psi_{1})>");
            grCos_diffPhiPsi1_LamBar[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#bar{#Lambda}}-#Psi_{1})>");
            grCos_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitle("<cos(#theta)>");
            grSin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitle("<sin(#theta)>");
            grSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitle("<sin(#phi_{#bar{#Lambda}}-#Psi_{1})sin(#theta)>");
            grCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#bar{#Lambda}}-#Psi_{1})sin(#theta)>");

            //Set graph titles
            grSin_diffPhiPsi1[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_diffPhiPsi1[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_theta[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grSin_theta[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grSin_diffPhiPsi1Sin_theta[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_diffPhiPsi1Sin_theta[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));

            grSin_diffPhiPsi1_LamBar[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_diffPhiPsi1_LamBar[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_theta_LamBar[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grSin_theta_LamBar[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            

            //Drawing
            c[iSub]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grSin_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grSin_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grSin_diffPhiPsi1[iCent][iSub]->Draw("AP");
            

            c[iSub+2]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_diffPhiPsi1[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_diffPhiPsi1[iCent][iSub]->Draw("AP");
            

            c[iSub+4]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_theta[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_theta[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_theta[iCent][iSub]->Draw("AP");
            
            
            c[iSub+6]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grSin_theta[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grSin_theta[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grSin_theta[iCent][iSub]->Draw("AP");
            

            c[iSub+8]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grSin_diffPhiPsi1Sin_theta[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grSin_diffPhiPsi1Sin_theta[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grSin_diffPhiPsi1Sin_theta[iCent][iSub]->Draw("AP");
            

            c[iSub+10]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_diffPhiPsi1Sin_theta[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_diffPhiPsi1Sin_theta[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_diffPhiPsi1Sin_theta[iCent][iSub]->Draw("AP");
            

            //LamBar
            c[iSub+12]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grSin_diffPhiPsi1_LamBar[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grSin_diffPhiPsi1_LamBar[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grSin_diffPhiPsi1_LamBar[iCent][iSub]->Draw("AP");
            

            c[iSub+14]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_diffPhiPsi1_LamBar[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_diffPhiPsi1_LamBar[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_diffPhiPsi1_LamBar[iCent][iSub]->Draw("AP");
            

            c[iSub+16]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_theta_LamBar[iCent][iSub]->Draw("AP");
            

            c[iSub+18]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grSin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grSin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grSin_theta_LamBar[iCent][iSub]->Draw("AP");
            

            c[iSub+20]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grSin_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->Draw("AP");

            c[iSub+22]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_diffPhiPsi1Sin_theta_LamBar[iCent][iSub]->Draw("AP");
        }
    }*/

    //20-50% range centrality drawing
    TProfile *prCos_diffPhiPsi1_20_50[2];

    TProfile *prCos_diffPhiPsi1_LamBar_20_50[2];

    for(int iSub=0; iSub!=2; iSub++){
        prCos_diffPhiPsi1_20_50[iSub] = (TProfile*)prCos_diffPhiPsi1[0][iSub]->Clone();
        prCos_diffPhiPsi1_20_50[iSub]->Reset();

        prCos_diffPhiPsi1_LamBar_20_50[iSub] = (TProfile*)prCos_diffPhiPsi1_LamBar[0][iSub]->Clone();
        prCos_diffPhiPsi1_LamBar_20_50[iSub]->Reset();
    }

    for(int iCent = 3; iCent<=5; iCent++){
        for(int iSub=0; iSub!=2; iSub++){
            prCos_diffPhiPsi1_20_50[iSub]->Add(prCos_diffPhiPsi1[iCent][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));

            prCos_diffPhiPsi1_LamBar_20_50[iSub]->Add(prCos_diffPhiPsi1_LamBar[iCent][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
        }
        
    }
    
    TGraphErrors *grCos_diffPhiPsi1_20_50[2];

    TGraphErrors *grCos_diffPhiPsi1_LamBar_20_50[2];

    for(int iSub=0; iSub!=2; iSub++){
        grCos_diffPhiPsi1_20_50[iSub] = ConvertProfToGraph(prCos_diffPhiPsi1_20_50[iSub]);

        grCos_diffPhiPsi1_LamBar_20_50[iSub] = ConvertProfToGraph(prCos_diffPhiPsi1_LamBar_20_50[iSub]);

    }


    //Set marker styles

    grCos_diffPhiPsi1_20_50[0]->SetMarkerColor(kRed);
    grCos_diffPhiPsi1_20_50[0]->SetMarkerStyle(29);
    grCos_diffPhiPsi1_20_50[0]->SetLineColor(kRed);

    //LamBar
    grCos_diffPhiPsi1_LamBar_20_50[0]->SetMarkerColor(kRed);
    grCos_diffPhiPsi1_LamBar_20_50[0]->SetMarkerStyle(34);
    grCos_diffPhiPsi1_LamBar_20_50[0]->SetLineColor(kRed);

    grCos_diffPhiPsi1_20_50[1]->SetMarkerColor(kBlue);
    grCos_diffPhiPsi1_20_50[1]->SetMarkerStyle(30);
    grCos_diffPhiPsi1_20_50[1]->SetLineColor(kBlue);

    //LamBar
    grCos_diffPhiPsi1_LamBar_20_50[1]->SetMarkerColor(kBlue);
    grCos_diffPhiPsi1_LamBar_20_50[1]->SetMarkerStyle(28);
    grCos_diffPhiPsi1_LamBar_20_50[1]->SetLineColor(kBlue);

    TCanvas *can = new TCanvas("can", "can", 0, 1024, 1200, 1000);

    
    TMultiGraph *mgCos_diffPhiPsi1_20_50 = new TMultiGraph();


    TLegend *leg;

    leg = new TLegend();
    leg->AddEntry(grCos_diffPhiPsi1_20_50[0], "#Lambda Psi_e/w was used", "p");
    leg->AddEntry(grCos_diffPhiPsi1_LamBar_20_50[0], "#bar{#Lambda} Psi_e/w was used", "p");
    leg->AddEntry(grCos_diffPhiPsi1_20_50[1], "#Lambda Psi_comb was used", "p");
    leg->AddEntry(grCos_diffPhiPsi1_LamBar_20_50[1], "#bar{#Lambda} Psi_comb was used", "p");
    leg->AddEntry((TObject*)0, "20-50 \% centrality", "");
    leg->SetBorderSize(0);

    
    for(int iSub=0; iSub!=2; iSub++){
        //
        mgCos_diffPhiPsi1_20_50->Add(grCos_diffPhiPsi1_20_50[iSub]);
        mgCos_diffPhiPsi1_20_50->Add(grCos_diffPhiPsi1_LamBar_20_50[iSub]);
    }        
    //Set Titles for tgraph

    mgCos_diffPhiPsi1_20_50->GetXaxis()->SetTitle("#phi - #Psi_{2}");
    mgCos_diffPhiPsi1_20_50->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})>");

    //Drawing

    can->cd();
    mgCos_diffPhiPsi1_20_50->Draw("AP");
    funcOfZero->Draw("same");
    leg->Draw("same");
}

TGraphErrors* ConvertProfToGraph(TProfile* prof){
    TGraphErrors *grout = new TGraphErrors(prof->GetNbinsX());
    double x, ex, y, ey;
    for(int iBin = 0; iBin!=prof->GetNbinsX(); iBin++){
        x = prof->GetBinCenter(iBin + 1);
        ex = 0.;
        y = prof->GetBinContent(iBin + 1);
        ey = prof->GetBinError(iBin+1);

        grout->SetPoint(iBin, x, y);
        grout->SetPointError(iBin, ex, ey);
    }

    grout->SetMarkerColor(1);
    grout->SetMarkerStyle(20);
    grout->SetMarkerSize(3);
    
    grout->GetXaxis()->SetTitle("#phi - #Psi_{2}");
    //grout->GetYaxis()->SetTitle("<cos(#theta)>");
    grout->GetXaxis()->SetLabelSize(0.1);
    grout->GetYaxis()->SetLabelSize(0.1);
    
    return grout;
}