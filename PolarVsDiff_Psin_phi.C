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
double GetCoef(double x, TF1 *Sig, TF1 *Bg);
double MyFitFunc(double *x, double *par);
TGraphErrors* GetCosTetaSub(TProfile* prof, TGraphErrors* graph);

TF1 *FuncGausKFP = new TF1("FuncGausKFP", "gaus(0)+gaus(3)", 1.1, 1.13);
TF1 *LineKFP = new TF1("LineKFP", "[0] + [1]*x", 1.1, 1.13);

void PolarVsDiff_Psin_phi(){
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

    TProfile *prSin_diffPhiPsi1[9][6][2];
    TProfile *prCos_diffPhiPsi1[9][6][2];
    TProfile *prCos_theta[9][6][2];
    TProfile *prSin_theta[9][6][2];
    TProfile *prSin_diffPhiPsi1Sin_theta[9][6][2];
    TProfile *prCos_diffPhiPsi1Sin_theta[9][6][2];

    TProfile *prSin_diffPhiPsi1_LamBar[9][6][2];
    TProfile *prCos_diffPhiPsi1_LamBar[9][6][2];
    TProfile *prCos_theta_LamBar[9][6][2];
    TProfile *prSin_theta_LamBar[9][6][2];
    TProfile *prSin_diffPhiPsi1Sin_theta_LamBar[9][6][2];
    TProfile *prCos_diffPhiPsi1Sin_theta_LamBar[9][6][2];

    TH1F *InvMLamDist[9][6][2];
    TH1F *InvMLamBarDist[9][6][2];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iphi=0; iphi!=6; iphi++){
            for(int iSub=0; iSub!=2; iSub++){
                prSin_diffPhiPsi1[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_%i_%i%i", iCent, iphi, iSub));
                prCos_diffPhiPsi1[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1_%i_%i_%i", iCent, iphi, iSub));
                prCos_theta[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_theta_%i_%i_%i", iCent, iphi, iSub));
                prSin_theta[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_theta_%i_%i_%i", iCent, iphi, iSub));
                prSin_diffPhiPsi1Sin_theta[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, iSub));
                prCos_diffPhiPsi1Sin_theta[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, iSub));

                prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, iSub));
                prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, iSub));
                prCos_theta_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_theta_LamBar_%i_%i_%i", iCent, iphi, iSub));
                prSin_theta_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub));
                prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub));
                prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub));

                //lam invM dist
                InvMLamDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub));
                InvMLamBarDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub));
            }
        }
        
    }

    TGraphErrors *grSin_diffPhiPsi1[9][2];
    TGraphErrors *grCos_diffPhiPsi1[9][2];
    TGraphErrors *grCos_theta[9][2];
    TGraphErrors *grSin_theta[9][2];
    TGraphErrors *grSin_diffPhiPsi1Sin_theta[9][2];
    TGraphErrors *grCos_diffPhiPsi1Sin_theta[9][2];

    TGraphErrors *grSin_diffPhiPsi1_LamBar[9][2];
    TGraphErrors *grCos_diffPhiPsi1_LamBar[9][2];
    TGraphErrors *grCos_theta_LamBar[9][2];
    TGraphErrors *grSin_theta_LamBar[9][2];
    TGraphErrors *grSin_diffPhiPsi1Sin_theta_LamBar[9][2];
    TGraphErrors *grCos_diffPhiPsi1Sin_theta_LamBar[9][2];

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

    //........................20-50% range centrality drawing.......................

    //20-50% range cetrality merging
    TProfile *prSin_diffPhiPsi1_20_50[6][2];
    TProfile *prCos_diffPhiPsi1_20_50[6][2];
    TProfile *prCos_theta_20_50[6][2];
    TProfile *prSin_theta_20_50[6][2];
    TProfile *prSin_diffPhiPsi1Sin_theta_20_50[6][2];
    TProfile *prCos_diffPhiPsi1Sin_theta_20_50[6][2];

    TProfile *prSin_diffPhiPsi1_LamBar_20_50[6][2];
    TProfile *prCos_diffPhiPsi1_LamBar_20_50[6][2];
    TProfile *prCos_theta_LamBar_20_50[6][2];
    TProfile *prSin_theta_LamBar_20_50[6][2];
    TProfile *prSin_diffPhiPsi1Sin_theta_LamBar_20_50[6][2];
    TProfile *prCos_diffPhiPsi1Sin_theta_LamBar_20_50[6][2];

    //lam invM dist
    TH1F *InvMLamDist_20_50[6][2];
    TH1F *InvMLamBarDist_20_50[6][2];

    for(int iSub=0; iSub!=2; iSub++){
        for(int iphi=0; iphi!=6; iphi++){
            prSin_diffPhiPsi1_20_50[iphi][iSub] = (TProfile*)prSin_diffPhiPsi1[0][iphi][iSub]->Clone();
            prSin_diffPhiPsi1_20_50[iphi][iSub]->Reset();
            prCos_diffPhiPsi1_20_50[iphi][iSub] = (TProfile*)prCos_diffPhiPsi1[0][iphi][iSub]->Clone();
            prCos_diffPhiPsi1_20_50[iphi][iSub]->Reset();
            prCos_theta_20_50[iphi][iSub] = (TProfile*)prCos_theta[0][iphi][iSub]->Clone();
            prCos_theta_20_50[iphi][iSub]->Reset();
            prSin_theta_20_50[iphi][iSub] = (TProfile*)prSin_theta[0][iphi][iSub]->Clone();
            prSin_theta_20_50[iphi][iSub]->Reset();
            prSin_diffPhiPsi1Sin_theta_20_50[iphi][iSub] = (TProfile*)prSin_diffPhiPsi1Sin_theta[0][iphi][iSub]->Clone();
            prSin_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Reset();
            prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub] = (TProfile*)prCos_diffPhiPsi1Sin_theta[0][iphi][iSub]->Clone();
            prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Reset();

            prSin_diffPhiPsi1_LamBar_20_50[iphi][iSub] = (TProfile*)prSin_diffPhiPsi1_LamBar[0][iphi][iSub]->Clone();
            prSin_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Reset();
            prCos_diffPhiPsi1_LamBar_20_50[iphi][iSub] = (TProfile*)prCos_diffPhiPsi1_LamBar[0][iphi][iSub]->Clone();
            prCos_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Reset();
            prCos_theta_LamBar_20_50[iphi][iSub] = (TProfile*)prCos_theta_LamBar[0][iphi][iSub]->Clone();
            prCos_theta_LamBar_20_50[iphi][iSub]->Reset();
            prSin_theta_LamBar_20_50[iphi][iSub] = (TProfile*)prSin_theta_LamBar[0][iphi][iSub]->Clone();
            prSin_theta_LamBar_20_50[iphi][iSub]->Reset();
            prSin_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub] = (TProfile*)prSin_diffPhiPsi1Sin_theta_LamBar[0][iphi][iSub]->Clone();
            prSin_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Reset();
            prCos_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub] = (TProfile*)prCos_diffPhiPsi1Sin_theta_LamBar[0][iphi][iSub]->Clone();
            prCos_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Reset();

            //lam invM dist
            InvMLamDist_20_50[iphi][iSub] = (TH1F*)InvMLamDist[0][iphi][iSub]->Clone();
            InvMLamDist_20_50[iphi][iSub]->Reset();
            InvMLamBarDist_20_50[iphi][iSub] = (TH1F*)InvMLamBarDist[0][iphi][iSub]->Clone();
            InvMLamBarDist_20_50[iphi][iSub]->Reset();

            for(int iCent = 3; iCent!=6; iCent++){
                prSin_diffPhiPsi1_20_50[iphi][iSub]->Add(prSin_diffPhiPsi1[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
                prCos_diffPhiPsi1_20_50[iphi][iSub]->Add(prCos_diffPhiPsi1[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
                prCos_theta_20_50[iphi][iSub]->Add(prCos_theta[iCent][iphi][iSub]);
                prSin_theta_20_50[iphi][iSub]->Add(prSin_theta[iCent][iphi][iSub]);
                prSin_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Add(prSin_diffPhiPsi1Sin_theta[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
                prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Add(prCos_diffPhiPsi1Sin_theta[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));

                prSin_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Add(prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
                prCos_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Add(prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
                prCos_theta_LamBar_20_50[iphi][iSub]->Add(prCos_theta_LamBar[iCent][iphi][iSub]);
                prSin_theta_LamBar_20_50[iphi][iSub]->Add(prSin_theta_LamBar[iCent][iphi][iSub]);
                prSin_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Add(prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));
                prCos_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Add(prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][iSub], Res_1VsCent[iSub]->GetPointY(8-iCent));

                //lam invM dist
                InvMLamDist_20_50[iphi][iSub]->Add(InvMLamDist[iCent][iphi][iSub]);
                InvMLamBarDist_20_50[iphi][iSub]->Add(InvMLamBarDist[iCent][iphi][iSub]);
            }        
        }       
    }

    //...................................Fitting by invM..................................
    double Pz, Pz_err;

    //Canvas init
    TCanvas *c_lam_20_50[2];
    TCanvas *c_lamBar_20_50[2];

    TCanvas *c_Sin_diffPhiPsi1_20_50[2];
    TCanvas *c_Sin_diffPhiPsi1_LamBar_20_50[2];

    TCanvas *c_Cos_diffPhiPsi1_20_50[2];
    TCanvas *c_Cos_diffPhiPsi1_LamBar_20_50[2];

    TCanvas *c_Cos_theta_20_50[2];
    TCanvas *c_Cos_theta_LamBar_20_50[2];

    TCanvas *c_Sin_theta_20_50[2];
    TCanvas *c_Sin_theta_LamBar_20_50[2];

    TCanvas *c_Sin_diffPhiPsi1Sin_theta_20_50[2];
    TCanvas *c_Sin_diffPhiPsi1Sin_theta_LamBar_20_50[2];

    TCanvas *c_Cos_diffPhiPsi1Sin_theta_20_50[2];
    TCanvas *c_Cos_diffPhiPsi1Sin_theta_LamBar_20_50[2];

    for(int iSub=0; iSub!=2; iSub++){
        c_lam_20_50[iSub] = new TCanvas(Form("c_lam_20_50_%i", iSub), Form("c_lam_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_lamBar_20_50[iSub] = new TCanvas(Form("c_lamBar_20_50_%i", iSub), Form("c_lamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Sin_diffPhiPsi1_20_50[iSub] = new TCanvas(Form("c_Sin_diffPhiPsi1_20_50_%i", iSub), Form("c_Sin_diffPhiPsi1_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_Sin_diffPhiPsi1_LamBar_20_50[iSub] = new TCanvas(Form("c_Sin_diffPhiPsi1_LamBar_20_50_%i", iSub), Form("c_Sin_diffPhiPsi1_LamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Cos_diffPhiPsi1_20_50[iSub] = new TCanvas(Form("c_Cos_diffPhiPsi1_20_50_%i", iSub), Form("c_Cos_diffPhiPsi1_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_Cos_diffPhiPsi1_LamBar_20_50[iSub] = new TCanvas(Form("c_Cos_diffPhiPsi1_LamBar_20_50_%i", iSub), Form("c_Cos_diffPhiPsi1_LamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Cos_theta_20_50[iSub] = new TCanvas(Form("c_Cos_theta_20_50_%i", iSub), Form("c_Cos_theta_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_Cos_theta_LamBar_20_50[iSub] = new TCanvas(Form("c_Cos_theta_LamBar_20_50_%i", iSub), Form("c_Cos_theta_LamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Sin_theta_20_50[iSub] = new TCanvas(Form("c_Sin_theta_20_50_%i", iSub), Form("c_Sin_theta_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_Sin_theta_LamBar_20_50[iSub] = new TCanvas(Form("c_Sin_theta_LamBar_20_50_%i", iSub), Form("c_Sin_theta_LamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Sin_diffPhiPsi1Sin_theta_20_50[iSub] = new TCanvas(Form("c_Sin_diffPhiPsi1Sin_theta_20_50_%i", iSub), Form("c_Sin_diffPhiPsi1Sin_theta_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_Sin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub] = new TCanvas(Form("c_Sin_diffPhiPsi1Sin_theta_LamBar_20_50_%i", iSub), Form("c_Sin_diffPhiPsi1Sin_theta_LamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Cos_diffPhiPsi1Sin_theta_20_50[iSub] = new TCanvas(Form("c_Cos_diffPhiPsi1Sin_theta_20_50_%i", iSub), Form("c_Cos_diffPhiPsi1Sin_theta_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_Cos_diffPhiPsi1Sin_theta_LamBar_20_50[iSub] = new TCanvas(Form("c_Cos_diffPhiPsi1Sin_theta_LamBar_20_50_%i", iSub), Form("c_Cos_diffPhiPsi1Sin_theta_LamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);
        
        
        c_lam_20_50[iSub]->Divide(3, 2);
        c_lamBar_20_50[iSub]->Divide(3, 2);

        c_Sin_diffPhiPsi1_20_50[iSub]->Divide(3, 2);
        c_Sin_diffPhiPsi1_LamBar_20_50[iSub]->Divide(3, 2);

        c_Cos_diffPhiPsi1_20_50[iSub]->Divide(3, 2);
        c_Cos_diffPhiPsi1_LamBar_20_50[iSub]->Divide(3, 2);

        c_Cos_theta_20_50[iSub]->Divide(3, 2);
        c_Cos_theta_LamBar_20_50[iSub]->Divide(3, 2);

        c_Sin_theta_20_50[iSub]->Divide(3, 2);
        c_Sin_theta_LamBar_20_50[iSub]->Divide(3, 2);

        c_Sin_diffPhiPsi1Sin_theta_20_50[iSub]->Divide(3, 2);
        c_Sin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->Divide(3, 2);

        c_Cos_diffPhiPsi1Sin_theta_20_50[iSub]->Divide(3, 2);
        c_Cos_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->Divide(3, 2);
    }
    
    //Fit function init
    TF1 *FitCosTetaFunc;

    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.1, 1.13);
	FitFuncKFP->SetParLimits(2, 0.0005, 0.001);    
    FitFuncKFP->SetParLimits(4, 1.113, 1.12);
	FitFuncKFP->SetParLimits(5, 0.001, 0.005);

    //graphs init
    TGraphErrors *grSin_diffPhiPsi1_20_50[2];
    TGraphErrors *grCos_diffPhiPsi1_20_50[2];
    TGraphErrors *grCos_theta_20_50[2];
    TGraphErrors *grSin_theta_20_50[2];
    TGraphErrors *grSin_diffPhiPsi1Sin_theta_20_50[2];
    TGraphErrors *grCos_diffPhiPsi1Sin_theta_20_50[2];

    TGraphErrors *grSin_diffPhiPsi1_LamBar_20_50[2];
    TGraphErrors *grCos_diffPhiPsi1_LamBar_20_50[2];
    TGraphErrors *grCos_theta_LamBar_20_50[2];
    TGraphErrors *grSin_theta_LamBar_20_50[2];
    TGraphErrors *grSin_diffPhiPsi1Sin_theta_LamBar_20_50[2];
    TGraphErrors *grCos_diffPhiPsi1Sin_theta_LamBar_20_50[2];

    //This need to subtract distributions
    TProfile *prCos_diffPhiPsi1_20_50_aver[2];
    TProfile *prCos_theta_20_50_aver[2];
    TProfile *prCos_diffPhiPsi1Sin_theta_20_50_aver[2];

    TProfile *prCos_diffPhiPsi1_LamBar_20_50_aver[2];
    TProfile *prCos_theta_LamBar_20_50_aver[2];
    TProfile *prCos_diffPhiPsi1Sin_theta_LamBar_20_50_aver[2];

    for(int iSub=0; iSub!=2; iSub++){
        prCos_diffPhiPsi1_20_50_aver[iSub] = (TProfile*)prCos_diffPhiPsi1_20_50[0][iSub]->Clone();
        prCos_diffPhiPsi1_20_50_aver[iSub]->Reset();
        prCos_theta_20_50_aver[iSub] = (TProfile*)prCos_theta_20_50[0][iSub]->Clone();
        prCos_theta_20_50_aver[iSub]->Reset();
        prCos_diffPhiPsi1Sin_theta_20_50_aver[iSub] = (TProfile*)prCos_diffPhiPsi1Sin_theta_20_50[0][iSub]->Clone();
        prCos_diffPhiPsi1Sin_theta_20_50_aver[iSub]->Reset();

        prCos_diffPhiPsi1_LamBar_20_50_aver[iSub] = (TProfile*)prCos_diffPhiPsi1_LamBar_20_50[0][iSub]->Clone();
        prCos_diffPhiPsi1_LamBar_20_50_aver[iSub]->Reset();
        prCos_theta_LamBar_20_50_aver[iSub] = (TProfile*)prCos_theta_LamBar_20_50[0][iSub]->Clone();
        prCos_theta_LamBar_20_50_aver[iSub]->Reset();
        prCos_diffPhiPsi1Sin_theta_LamBar_20_50_aver[iSub] = (TProfile*)prCos_diffPhiPsi1Sin_theta_LamBar_20_50[0][iSub]->Clone();
        prCos_diffPhiPsi1Sin_theta_LamBar_20_50_aver[iSub]->Reset();

        for(int iphi=0; iphi!=6; iphi++){
            //....................Lam invM..................
            c_lam_20_50[iSub]->cd(iphi+1);
            InvMLamDist_20_50[iphi][iSub]->GetXaxis()->SetTitle("M_{inv}");
            InvMLamDist_20_50[iphi][iSub]->GetXaxis()->SetLabelSize(0.06);
            InvMLamDist_20_50[iphi][iSub]->GetXaxis()->SetTitleSize(0.06);
            InvMLamDist_20_50[iphi][iSub]->Fit(FitFuncKFP, "R");        
            InvMLamDist_20_50[iphi][iSub]->Draw();

            FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
            FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
            FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
            FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
            FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
            FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

            LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
            LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));
            
            //<Sin(phi-Psi1)>
            c_Sin_diffPhiPsi1_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prSin_diffPhiPsi1_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prSin_diffPhiPsi1_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            double phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grSin_diffPhiPsi1_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grSin_diffPhiPsi1_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Cos(phi-Psi1)>
            prCos_diffPhiPsi1_20_50_aver[iSub]->Add(prCos_diffPhiPsi1_20_50[iphi][iSub]);

            c_Cos_diffPhiPsi1_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prCos_diffPhiPsi1_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prCos_diffPhiPsi1_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grCos_diffPhiPsi1_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_diffPhiPsi1_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Cos(Theta)>
            prCos_theta_20_50_aver[iSub]->Add(prCos_diffPhiPsi1_20_50[iphi][iSub]);

            c_Cos_theta_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prCos_theta_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prCos_theta_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grCos_theta_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_theta_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Sin(Theta)>
            c_Sin_theta_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prSin_theta_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prSin_theta_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grSin_theta_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grSin_theta_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Sin(Psi1-phi)Sin(Theta)>
            c_Sin_diffPhiPsi1Sin_theta_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prSin_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prSin_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Cos(Psi1-phi)Sin(Theta)>
            prCos_diffPhiPsi1Sin_theta_20_50_aver[iSub]->Add(prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub]);

            c_Cos_diffPhiPsi1Sin_theta_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grCos_diffPhiPsi1Sin_theta_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_diffPhiPsi1Sin_theta_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //....................LamBar invM..................
            c_lamBar_20_50[iSub]->cd(iphi+1);
            InvMLamBarDist_20_50[iphi][iSub]->GetXaxis()->SetTitle("M_{inv}");
            InvMLamBarDist_20_50[iphi][iSub]->GetXaxis()->SetLabelSize(0.06);
            InvMLamBarDist_20_50[iphi][iSub]->GetXaxis()->SetTitleSize(0.06);
            InvMLamBarDist_20_50[iphi][iSub]->Fit(FitFuncKFP, "R");        
            InvMLamBarDist_20_50[iphi][iSub]->Draw();

            FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
            FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
            FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
            FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
            FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
            FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

            LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
            LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));
            
            //<Sin(phi-Psi1)>
            c_Sin_diffPhiPsi1_LamBar_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prSin_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prSin_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Cos(phi-Psi1)>
            prCos_diffPhiPsi1_LamBar_20_50_aver[iSub]->Add(prCos_diffPhiPsi1_20_50[iphi][iSub]);

            c_Cos_diffPhiPsi1_LamBar_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prCos_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prCos_diffPhiPsi1_LamBar_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grCos_diffPhiPsi1_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_diffPhiPsi1_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Cos(Theta)>
            prCos_theta_LamBar_20_50_aver[iSub]->Add(prCos_diffPhiPsi1_20_50[iphi][iSub]);

            c_Cos_theta_LamBar_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prCos_theta_LamBar_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prCos_theta_LamBar_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grCos_theta_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_theta_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Sin(Theta)>
            c_Sin_theta_LamBar_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prSin_theta_LamBar_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prSin_theta_LamBar_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grSin_theta_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grSin_theta_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Sin(Psi1-phi)Sin(Theta)>
            c_Sin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prSin_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prSin_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

            //<Cos(Psi1-phi)Sin(Theta)>
            prCos_diffPhiPsi1Sin_theta_LamBar_20_50_aver[iSub]->Add(prCos_diffPhiPsi1Sin_theta_20_50[iphi][iSub]);

            c_Cos_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->cd(iphi+1);        

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            prCos_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            prCos_diffPhiPsi1Sin_theta_LamBar_20_50[iphi][iSub]->Draw();

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            grCos_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

        }//for(int iphi=0; iphi!=6; iphi++)
    }//for(int iSub=0; iSub!=2; iSub++)

    //...................................end of Fitting by invM..................................

    //...................................subtracting.............................................
    
    TGraphErrors *grCos_diffPhiPsi1_sub_20_50[2];
    TGraphErrors *grCos_theta_sub_20_50[2];
    TGraphErrors *grCos_diffPhiPsi1Sin_theta_sub_20_50[2];

    TGraphErrors *grCos_diffPhiPsi1_LamBar_sub_20_50[2];
    TGraphErrors *grCos_theta_LamBar_sub_20_50[2];
    TGraphErrors *grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[2];

    for(int iSub=0; iSub!=2; iSub++){
        grCos_diffPhiPsi1_sub_20_50[iSub] = GetCosTetaSub(prCos_diffPhiPsi1_20_50_aver[iSub], grCos_diffPhiPsi1_20_50[iSub]);
        grCos_theta_sub_20_50[iSub] = GetCosTetaSub(prCos_theta_20_50_aver[iSub], grCos_theta_20_50[iSub]);
        grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub] = GetCosTetaSub(prCos_diffPhiPsi1Sin_theta_20_50_aver[iSub], grCos_diffPhiPsi1Sin_theta_20_50[iSub]);

        grCos_diffPhiPsi1_LamBar_sub_20_50[iSub] = GetCosTetaSub(prCos_diffPhiPsi1_LamBar_20_50_aver[iSub], grCos_diffPhiPsi1_LamBar_20_50[iSub]);
        grCos_theta_LamBar_sub_20_50[iSub] = GetCosTetaSub(prCos_theta_LamBar_20_50_aver[iSub], grCos_theta_LamBar_20_50[iSub]);
        grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub] = GetCosTetaSub(prCos_diffPhiPsi1Sin_theta_LamBar_20_50_aver[iSub], grCos_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]);

        //Set Marker Styles
        //<Sin(Psi1-phi)>
        grSin_diffPhiPsi1_20_50[iSub]->SetMarkerSize(3);
        grSin_diffPhiPsi1_20_50[iSub]->SetLineWidth(2);
        grSin_diffPhiPsi1_20_50[iSub]->SetMarkerStyle(29);
        grSin_diffPhiPsi1_20_50[iSub]->SetMarkerColor(634);
        grSin_diffPhiPsi1_20_50[iSub]->SetLineColor(634);

        grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetMarkerSize(3);
        grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetLineWidth(2);
        grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetMarkerStyle(30);
        grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetMarkerColor(602);
        grSin_diffPhiPsi1_LamBar_20_50[iSub]->SetLineColor(602);

        //<Cos(Psi1-phi)>
        grCos_diffPhiPsi1_sub_20_50[iSub]->SetMarkerSize(3);
        grCos_diffPhiPsi1_sub_20_50[iSub]->SetLineWidth(2);
        grCos_diffPhiPsi1_sub_20_50[iSub]->SetMarkerStyle(29);
        grCos_diffPhiPsi1_sub_20_50[iSub]->SetMarkerColor(634);
        grCos_diffPhiPsi1_sub_20_50[iSub]->SetLineColor(634);

        grCos_diffPhiPsi1_LamBar_sub_20_50[iSub]->SetMarkerSize(3);
        grCos_diffPhiPsi1_LamBar_sub_20_50[iSub]->SetLineWidth(2);
        grCos_diffPhiPsi1_LamBar_sub_20_50[iSub]->SetMarkerStyle(30);
        grCos_diffPhiPsi1_LamBar_sub_20_50[iSub]->SetMarkerColor(602);
        grCos_diffPhiPsi1_LamBar_sub_20_50[iSub]->SetLineColor(602);

        //<Cos(Theta)>
        grCos_theta_sub_20_50[iSub]->SetMarkerSize(3);
        grCos_theta_sub_20_50[iSub]->SetLineWidth(2);
        grCos_theta_sub_20_50[iSub]->SetMarkerStyle(29);
        grCos_theta_sub_20_50[iSub]->SetMarkerColor(634);
        grCos_theta_sub_20_50[iSub]->SetLineColor(634);

        grCos_theta_LamBar_sub_20_50[iSub]->SetMarkerSize(3);
        grCos_theta_LamBar_sub_20_50[iSub]->SetLineWidth(2);
        grCos_theta_LamBar_sub_20_50[iSub]->SetMarkerStyle(30);
        grCos_theta_LamBar_sub_20_50[iSub]->SetMarkerColor(602);
        grCos_theta_LamBar_sub_20_50[iSub]->SetLineColor(602);
        
        //<Sin(Theta)>
        grSin_theta_20_50[iSub]->SetMarkerSize(3);
        grSin_theta_20_50[iSub]->SetLineWidth(2);
        grSin_theta_20_50[iSub]->SetMarkerStyle(29);
        grSin_theta_20_50[iSub]->SetMarkerColor(634);
        grSin_theta_20_50[iSub]->SetLineColor(634);

        grSin_theta_LamBar_20_50[iSub]->SetMarkerSize(3);
        grSin_theta_LamBar_20_50[iSub]->SetLineWidth(2);
        grSin_theta_LamBar_20_50[iSub]->SetMarkerStyle(30);
        grSin_theta_LamBar_20_50[iSub]->SetMarkerColor(602);
        grSin_theta_LamBar_20_50[iSub]->SetLineColor(602);

        //<Sin(Psi1-phi)Sin(Theta)>
        grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetMarkerSize(3);
        grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetLineWidth(2);
        grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetMarkerStyle(29);
        grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetMarkerColor(634);
        grSin_diffPhiPsi1Sin_theta_20_50[iSub]->SetLineColor(634);

        grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetMarkerSize(3);
        grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetLineWidth(2);
        grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetMarkerStyle(30);
        grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetMarkerColor(602);
        grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]->SetLineColor(602);

        //<Cos(Psi1-phi)Sin(Theta)>
        grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub]->SetMarkerSize(3);
        grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub]->SetLineWidth(2);
        grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub]->SetMarkerStyle(29);
        grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub]->SetMarkerColor(634);
        grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub]->SetLineColor(634);

        grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub]->SetMarkerSize(3);
        grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub]->SetLineWidth(2);
        grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub]->SetMarkerStyle(30);
        grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub]->SetMarkerColor(602);
        grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub]->SetLineColor(602);
        
        
    }

    //...................................end of subtracting.............................................

    //...................................Drawing........................................................
    TCanvas *can[12];

    for(int i=0; i!=12; i++){
        can[i] = new TCanvas(Form("can_%i", i), Form("can_%i", i), 0, 1024, 1200, 1000);
    }
    
    TMultiGraph *mgSin_diffPhiPsi1_20_50[2];
    TMultiGraph *mgCos_diffPhiPsi1_20_50[2];
    TMultiGraph *mgCos_theta_20_50[2];
    TMultiGraph *mgSin_theta_20_50[2];
    TMultiGraph *mgSin_diffPhiPsi1Sin_theta_20_50[2];
    TMultiGraph *mgCos_diffPhiPsi1Sin_theta_20_50[2];

    TLegend *leg[2];
    for(int iSub=0; iSub!=2; iSub++){
        leg[iSub] = new TLegend();

        leg[iSub]->AddEntry(grCos_theta_20_50[0], "#Lambda", "p");
        leg[iSub]->AddEntry(grCos_theta_LamBar_20_50[0], "#bar{#Lambda}", "p");
        leg[iSub]->AddEntry((TObject*)0, "20-50 \% centrality", "");
        leg[iSub]->SetBorderSize(0);
    }


    leg[0]->AddEntry((TObject*)0, "Psi_comb was used", "");
    leg[1]->AddEntry((TObject*)0, "Psi_e/w was used", "");

    for(int iSub=0; iSub!=2; iSub++){
        mgSin_diffPhiPsi1_20_50[iSub] = new TMultiGraph();
        mgCos_diffPhiPsi1_20_50[iSub] = new TMultiGraph();
        mgCos_theta_20_50[iSub] = new TMultiGraph();
        mgSin_theta_20_50[iSub] = new TMultiGraph();
        mgSin_diffPhiPsi1Sin_theta_20_50[iSub] = new TMultiGraph();
        mgCos_diffPhiPsi1Sin_theta_20_50[iSub] = new TMultiGraph();

        //
        mgSin_diffPhiPsi1_20_50[iSub]->Add(grSin_diffPhiPsi1_20_50[iSub]);
        mgSin_diffPhiPsi1_20_50[iSub]->Add(grSin_diffPhiPsi1_LamBar_20_50[iSub]);

        //
        mgCos_diffPhiPsi1_20_50[iSub]->Add(grCos_diffPhiPsi1_sub_20_50[iSub]);
        mgCos_diffPhiPsi1_20_50[iSub]->Add(grCos_diffPhiPsi1_LamBar_sub_20_50[iSub]);

        //
        mgCos_theta_20_50[iSub]->Add(grCos_theta_sub_20_50[iSub]);
        mgCos_theta_20_50[iSub]->Add(grCos_theta_LamBar_sub_20_50[iSub]);

        //
        mgSin_theta_20_50[iSub]->Add(grSin_theta_20_50[iSub]);
        mgSin_theta_20_50[iSub]->Add(grSin_theta_LamBar_20_50[iSub]);

        //
        mgSin_diffPhiPsi1Sin_theta_20_50[iSub]->Add(grSin_diffPhiPsi1Sin_theta_20_50[iSub]);
        mgSin_diffPhiPsi1Sin_theta_20_50[iSub]->Add(grSin_diffPhiPsi1Sin_theta_LamBar_20_50[iSub]);

        //
        mgCos_diffPhiPsi1Sin_theta_20_50[iSub]->Add(grCos_diffPhiPsi1Sin_theta_sub_20_50[iSub]);
        mgCos_diffPhiPsi1Sin_theta_20_50[iSub]->Add(grCos_diffPhiPsi1Sin_theta_LamBar_sub_20_50[iSub]);

        //Set Titles for tgraph
        mgSin_diffPhiPsi1_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgSin_diffPhiPsi1_20_50[iSub]->GetYaxis()->SetTitle("<sin(#phi_{#Lambda}-#Psi_{1})>");

        mgCos_diffPhiPsi1_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgCos_diffPhiPsi1_20_50[iSub]->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})>");

        mgCos_theta_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgCos_theta_20_50[iSub]->GetYaxis()->SetTitle("<cos(#theta)>");

        mgSin_theta_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgSin_theta_20_50[iSub]->GetYaxis()->SetTitle("<sin(#theta)>");

        mgSin_diffPhiPsi1Sin_theta_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgSin_diffPhiPsi1Sin_theta_20_50[iSub]->GetYaxis()->SetTitle("<sin(#phi_{#Lambda}-#Psi_{1})sin(#theta)>");

        mgCos_diffPhiPsi1Sin_theta_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgCos_diffPhiPsi1Sin_theta_20_50[iSub]->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})sin(#theta)>");

        //Drawing
        can[iSub]->cd();
        mgSin_diffPhiPsi1_20_50[iSub]->Draw("AP");
        leg[iSub]->Draw("same");
        funcOfZero->Draw("same");

        can[iSub+2]->cd();
        mgCos_diffPhiPsi1_20_50[iSub]->Draw("AP");
        leg[iSub]->Draw("same");
        funcOfZero->Draw("same");

        can[iSub+4]->cd();
        mgCos_theta_20_50[iSub]->Draw("AP");
        leg[iSub]->Draw("same");
        funcOfZero->Draw("same");

        can[iSub+6]->cd();
        mgSin_theta_20_50[iSub]->Draw("AP");
        leg[iSub]->Draw("same");
        funcOfZero->Draw("same");

        can[iSub+8]->cd();
        mgSin_diffPhiPsi1Sin_theta_20_50[iSub]->Draw("AP");
        leg[iSub]->Draw("same");
        funcOfZero->Draw("same");

        can[iSub+10]->cd();
        mgCos_diffPhiPsi1Sin_theta_20_50[iSub]->Draw("AP");
        leg[iSub]->Draw("same");
        funcOfZero->Draw("same");

    }
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

double GetCoef(double x, TF1 *Sig, TF1 *Bg){
    double coef =0.;
    coef = Bg->Eval(x)/(Bg->Eval(x) + Sig->Eval(x));
    return coef;
}

double MyFitFunc(double *x, double *par){
    double coef = GetCoef(x[0], FuncGausKFP, LineKFP);
    double f = (1-coef)*par[0] + coef*(par[1] + par[2]*x[0]);
    return f;
}

TGraphErrors* GetCosTetaSub(TProfile* prof, TGraphErrors* graph){
    TProfile *prof_new = (TProfile*)prof->Clone();
    TH1 *hist_average = prof_new->Rebin(prof->GetNbinsX());
    TGraphErrors *graph_cos_sub = new TGraphErrors(6);

    double x, ex, y, ey;

    for(int iBin=0; iBin!=6; iBin++){
        x = graph->GetPointX(iBin);
        ex = 0.;
        y = graph->GetPointY(iBin) - hist_average->GetBinContent(1);
        ey = sqrt(pow(graph->GetErrorY(iBin), 2) + pow(hist_average->GetBinError(1), 2));
        
        graph_cos_sub->SetPoint(iBin, x, y);
        graph_cos_sub->SetPointError(iBin, ex, ey);
    }

    graph_cos_sub->SetMarkerColor(1);
    graph_cos_sub->SetMarkerStyle(20);
    graph_cos_sub->SetMarkerSize(3);
    
    graph_cos_sub->GetXaxis()->SetTitle("#phi - #Psi_{2}");
    graph_cos_sub->GetYaxis()->SetTitle("<cos(#theta)>");
    graph_cos_sub->GetXaxis()->SetLabelSize(0.1);
    graph_cos_sub->GetYaxis()->SetLabelSize(0.1);
    
    return graph_cos_sub;
}