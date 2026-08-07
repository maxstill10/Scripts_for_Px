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

//0-Psi_e/w
//1-Psi_comb

const Double_t alpha_Lam = 0.732;
const Double_t alpha_LamBar = -0.732;

//true, if <cos(Psi_1-phi)sin(theta)> calculated
//false, if <cos(Psi_1-phi)> calculated
const bool sinUsage = false;

TGraphErrors* ConvertProfToGraph(TProfile* prof);
TGraphErrors* GetCosTetaSub(TProfile* prof);
TGraphErrors* GetCosTetaSub(TProfile* prof, TGraphErrors* graph);
TGraphErrors* MultiplyLamFactor(TGraphErrors *g, double coef);
TGraphErrors* MakeAverage( TGraphErrors *g1, TGraphErrors *g2 );
TGraphErrors* ResFactCorr(TGraphErrors* cosSubVsCent, TGraphErrors* resVsCent);
TGraphErrors* MakeAverageBin(TGraphErrors *g, double n_i, double n_f);

double GetCoef(double x, TF1 *Sig, TF1 *Bg);
double MyFitFunc(double *x, double *par);

TF1 *FuncGausKFP = new TF1("FuncGausKFP", "gaus(0)+gaus(3)", 1.1, 1.13);
TF1 *LineKFP = new TF1("LineKFP", "[0] + [1]*x", 1.1, 1.13);

void PxVsCent_calc(){
    TFile *input = new TFile("../Files/14p5GeV/output_Polar_14p5.root", "read");
    TFile *output;

    if(sinUsage) output = TFile::Open("../output/output_PxVsCent_SinInclude_14p5.root", "recreate");
    if(!sinUsage) output = TFile::Open("../output/output_PxVsCent_ConstSine_14p5.root", "recreate");

    TH1F *CosOfDiffFlatt_1 = (TH1F*)input->Get("CosOfDiff_1");
    TH1F *CosOfDiffFlatt_2 = (TH1F*)input->Get("CosOfDiff_2");

    TGraphErrors *Res_1VsCent[2];
    TGraphErrors *Res_2VsCent[2];

    for(int i=0; i!=2; i++){
        Res_1VsCent[i] = new TGraphErrors(9);
        Res_2VsCent[i] = new TGraphErrors(9);
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

        resFlatt = pow(2*TMath::Abs(CosOfDiffFlatt_2->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(resFlatt))*CosOfDiffFlatt_2->GetBinError(iCen+1);
        
        Res_2VsCent[0]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_2VsCent[0]->SetPointError(8-iCen, 0., resErr);

        resFlatt = pow(TMath::Abs(CosOfDiffFlatt_2->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(2*resFlatt))*CosOfDiffFlatt_2->GetBinError(iCen+1);
        
        Res_2VsCent[1]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_2VsCent[1]->SetPointError(8-iCen, 0., resErr);
        
    }

    TCanvas *c_res = new TCanvas("c_res", "c_res", 0, 1024, 1200, 1024);
    c_res->cd();

    Res_1VsCent[0]->SetMarkerStyle(23);
    Res_1VsCent[0]->SetMarkerColor(634);

    Res_1VsCent[1]->SetMarkerStyle(29);
    Res_1VsCent[1]->SetMarkerColor(633);

    Res_2VsCent[0]->SetMarkerStyle(34);
    Res_2VsCent[0]->SetMarkerColor(1);

    Res_2VsCent[1]->SetMarkerStyle(20);
    Res_2VsCent[1]->SetMarkerColor(922);

    TMultiGraph *mgResVsCent = new TMultiGraph();

    for(int iSub=0; iSub!=2; iSub++){
        Res_1VsCent[iSub]->SetMarkerSize(3);
        Res_2VsCent[iSub]->SetMarkerSize(3);

        mgResVsCent->Add(Res_1VsCent[iSub]);
        mgResVsCent->Add(Res_2VsCent[iSub]);
    }

    mgResVsCent->GetXaxis()->SetTitle("Centrality, %");
    mgResVsCent->GetYaxis()->SetTitle("Res(#Psi)");
    mgResVsCent->SetMinimum(0);
    mgResVsCent->SetMaximum(0.8);

    TLegend *leg_res = new TLegend();

    leg_res->AddEntry(Res_1VsCent[0], "Full Res(Psi1)", "p");
    leg_res->AddEntry(Res_1VsCent[1], "Half Res(Psi1)", "p");
    leg_res->AddEntry(Res_2VsCent[0], "Full Res(Psi2)", "p");
    leg_res->AddEntry(Res_2VsCent[1], "Half Res(Psi2)", "p");
    leg_res->SetBorderSize(0);

    mgResVsCent->Draw("AP");
    leg_res->Draw("same");

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);


    //..................................Px calculation starts here.................................
    TProfile *prCos_diffPhiPsi1[9][6][2];
    TProfile *prCos_diffPhiPsi1_LamBar[9][6][2];

    TH1F *InvMLamDist[9][6][2];
    TH1F *InvMLamBarDist[9][6][2];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iphi=0; iphi!=6; iphi++){
            for(int iSub=0; iSub!=2; iSub++){
                if(sinUsage){
                    prCos_diffPhiPsi1[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, iSub));
                    prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, iSub));
                }
                
                if(!sinUsage){
                    prCos_diffPhiPsi1[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1_%i_%i_%i", iCent, iphi, iSub));
                    prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub));
                }

                //lam invM dist
                InvMLamDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub));
                InvMLamBarDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub));
            }
        }        
    }

    //...................................Fitting by invM..................................
    TGraphErrors *grCos_diffPhiPsi1[9][2];
    TGraphErrors *grCos_diffPhiPsi1_LamBar[9][2];

    TProfile *prCos_diffPhiPsi1_aver[9][2];
    TProfile *prCos_diffPhiPsi1_LamBar_aver[9][2];

    for(int iCent = 0; iCent!=9; iCent++){
        for(int iSub = 0; iSub!=2; iSub++){           
            grCos_diffPhiPsi1[iCent][iSub] = new TGraphErrors(6);
            grCos_diffPhiPsi1_LamBar[iCent][iSub] = new TGraphErrors(6);
        }
    }
    
    double Pz, Pz_err;
    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.113, 1.12);
	FitFuncKFP->SetParLimits(2, 0.0005, 0.001);    
    FitFuncKFP->SetParLimits(4, 1.113, 1.12);
	FitFuncKFP->SetParLimits(5, 0.001, 0.005);

    TF1 *FitCosTetaFunc;

    TCanvas *canv[6][2];
    TCanvas *canv1[6][2];
    TCanvas *canv2[6][2];
    TCanvas *canv3[6][2];
    

    for(int iSub=0; iSub!=2; iSub++){
        for(int iphi = 0; iphi !=6; iphi++){
            canv[iphi][iSub] = new TCanvas(Form("canv_%i_%i", iphi, iSub), Form("canv_%i sub = %i, #Lambda", iphi, iSub), 0, 1024, 1800, 1024);
            canv1[iphi][iSub] = new TCanvas(Form("canv1_%i_%i", iphi, iSub), Form("canv1_%i sub = %i, #Lambda", iphi, iSub), 0, 1024, 1800, 1024);
            canv2[iphi][iSub] = new TCanvas(Form("canv2_%i_%i", iphi, iSub), Form("canv2_%i sub = %i, #bar{#Lambda}", iphi, iSub), 0, 1024, 1800, 1024);
            canv3[iphi][iSub] = new TCanvas(Form("canv3_%i_%i", iphi, iSub), Form("canv3_%i sub = %i, #bar{#Lambda}", iphi, iSub), 0, 1024, 1800, 1024);

            canv[iphi][iSub]->Divide(3,3);
            canv1[iphi][iSub]->Divide(3,3);
            canv2[iphi][iSub]->Divide(3,3);
            canv3[iphi][iSub]->Divide(3,3);
        }
    }
    

    for(int iSub=0; iSub!=2; iSub++){
        for(int iCent = 0; iCent!=9; iCent++){
            prCos_diffPhiPsi1_aver[iCent][iSub] = (TProfile*)prCos_diffPhiPsi1[iCent][0][iSub]->Clone();
            prCos_diffPhiPsi1_aver[iCent][iSub]->Reset();
            prCos_diffPhiPsi1_LamBar_aver[iCent][iSub] = (TProfile*)prCos_diffPhiPsi1_LamBar[iCent][0][iSub]->Clone();
            prCos_diffPhiPsi1_LamBar_aver[iCent][iSub]->Reset();
            for(int iphi = 0; iphi!=6; iphi++){
                
                prCos_diffPhiPsi1_aver[iCent][iSub]->Add(prCos_diffPhiPsi1[iCent][iphi][iSub]);            
                prCos_diffPhiPsi1_LamBar_aver[iCent][iSub]->Add(prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub]);

                canv[iphi][iSub]->cd(iCent+1);
            
                InvMLamDist[iCent][iphi][iSub]->GetXaxis()->SetTitle("M_{inv}");
                InvMLamDist[iCent][iphi][iSub]->GetXaxis()->SetLabelSize(0.06);
                InvMLamDist[iCent][iphi][iSub]->GetXaxis()->SetTitleSize(0.06);
                InvMLamDist[iCent][iphi][iSub]->Fit(FitFuncKFP, "R");        
                InvMLamDist[iCent][iphi][iSub]->Draw();

                FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
                FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
                FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
                FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
                FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
                FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

                LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
                LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));

                canv1[iphi][iSub]->cd(iCent+1);
                

                FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);

                //FitCosTetaFunc->FixParameter(1, 0.);
                //FitCosTetaFunc->FixParameter(2, 0.);
                
                prCos_diffPhiPsi1[iCent][iphi][iSub]->Fit(FitCosTetaFunc);
                prCos_diffPhiPsi1[iCent][iphi][iSub]->Draw();

                Pz = FitCosTetaFunc->GetParameter(0);
                Pz_err = FitCosTetaFunc->GetParError(0);

                double phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
                grCos_diffPhiPsi1[iCent][iSub]->SetPoint(iphi, phi, Pz);
                grCos_diffPhiPsi1[iCent][iSub]->SetPointError(iphi, 0., Pz_err);
                
                //AntiLambda
                canv2[iphi][iSub]->cd(iCent+1);
                
                InvMLamBarDist[iCent][iphi][iSub]->GetXaxis()->SetTitle("M_{inv}");
                InvMLamBarDist[iCent][iphi][iSub]->GetXaxis()->SetLabelSize(0.05);
                InvMLamBarDist[iCent][iphi][iSub]->GetXaxis()->SetTitleSize(0.05);
                InvMLamBarDist[iCent][iphi][iSub]->Fit(FitFuncKFP, "R");        
                InvMLamBarDist[iCent][iphi][iSub]->Draw();

                FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
                FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
                FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
                FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
                FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
                FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

                LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
                LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));

                canv3[iphi][iSub]->cd(iCent+1);
                

                FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);

                //FitCosTetaFunc->FixParameter(1, 0.);
                //FitCosTetaFunc->FixParameter(2, 0.);
                
                prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub]->Fit(FitCosTetaFunc);
                prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub]->Draw();

                Pz = FitCosTetaFunc->GetParameter(0);
                Pz_err = FitCosTetaFunc->GetParError(0);

                phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
                grCos_diffPhiPsi1_LamBar[iCent][iSub]->SetPoint(iphi, phi, Pz);
                grCos_diffPhiPsi1_LamBar[iCent][iSub]->SetPointError(iphi, 0., Pz_err);
            }//for(int iphi = 0; iphi!=6; iphi++)    
        }//for(int iCent = 0; iCent!=9; iCent++)
    }//for(int iSub=0; iSub!=2; iSub++)

    //...................................end of Fitting by invM..................................

    //...................................Subtracting.............................................
    TGraphErrors *grCos_diffPhiPsi1_sub[9][2];
    TGraphErrors *grCos_diffPhiPsi1_LamBar_sub[9][2];

    for(int iSub=0; iSub!=2; iSub++){
        for(int iCent=0; iCent!=9; iCent++){
            grCos_diffPhiPsi1_sub[iCent][iSub] = GetCosTetaSub(prCos_diffPhiPsi1_aver[iCent][iSub], grCos_diffPhiPsi1[iCent][iSub]);
            grCos_diffPhiPsi1_LamBar_sub[iCent][iSub] = GetCosTetaSub(prCos_diffPhiPsi1_LamBar_aver[iCent][iSub], grCos_diffPhiPsi1_LamBar[iCent][iSub]);
        }
    }
    //...................................end of Subtracting........................................

    //...................................Fitting by phi............................................

    TGraphErrors *PxVsCent_Lam[2];
    TGraphErrors *PxVsCent_LamBar[2];

    for(int iSub=0; iSub!=2; iSub++){
        PxVsCent_Lam[iSub] = new TGraphErrors(9);
        PxVsCent_LamBar[iSub] = new TGraphErrors(9);
    }

    TCanvas *c[4];
    for(int i=0; i!=4; i++){
        c[i] = new TCanvas(Form("c_%i", i), "", 0, 1024, 1200, 1024);
        c[i]->Divide(3,3);
    }
    c[0]->SetTitle("Lambda, Psi_e/w");
    c[1]->SetTitle("Lambda, Psi_comb");
    c[2]->SetTitle("LambdaBar, Psi_e/w");
    c[3]->SetTitle("LambdaBar, Psi_comb");

    //TGaxis::SetMaxDigits(1);
    gStyle->SetOptFit(0001);


    int cent_lim[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

    TF1 *fitFuncOfCos = new TF1("fitFuncOfCos", "[0] + 2*[1]*sin(2*x) + 2*[2]*cos(2*x)", 0, TMath::Pi());

    fitFuncOfCos->FixParameter(0, 0);
    fitFuncOfCos->SetParLimits(1, -1, 1);
    fitFuncOfCos->FixParameter(2, 0);

    double Px, Px_err;
    
    for(int iCent=0; iCent!=9; iCent++){
        for(int iSub=0; iSub!=2; iSub++){

            //Set Y titles
            grCos_diffPhiPsi1_sub[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})>");
            grCos_diffPhiPsi1_LamBar_sub[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#bar{#Lambda}}-#Psi_{1})>");

            //Set graph titles
            grCos_diffPhiPsi1_sub[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_diffPhiPsi1_LamBar_sub[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            

            //Drawing           

            c[iSub]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_diffPhiPsi1_sub[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_diffPhiPsi1_sub[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_diffPhiPsi1_sub[iCent][iSub]->Fit("fitFuncOfCos");
            grCos_diffPhiPsi1_sub[iCent][iSub]->Draw("AP");
            
            Px = fitFuncOfCos->GetParameter(1);
            Px_err = fitFuncOfCos->GetParError(1);
            PxVsCent_Lam[iSub]->SetPoint(8-iCent, cent[8-iCent], Px);
            PxVsCent_Lam[iSub]->SetPointError(8-iCent, 0, Px_err);


            c[iSub+2]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_diffPhiPsi1_LamBar_sub[iCent][iSub]->GetYaxis()->SetTitleOffset(0.7);
                grCos_diffPhiPsi1_LamBar_sub[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
            }
            grCos_diffPhiPsi1_LamBar_sub[iCent][iSub]->Draw("AP");
            grCos_diffPhiPsi1_LamBar_sub[iCent][iSub]->Fit("fitFuncOfCos");

            Px = fitFuncOfCos->GetParameter(1);
            Px_err = fitFuncOfCos->GetParError(1);
            PxVsCent_LamBar[iSub]->SetPoint(8-iCent, cent[8-iCent]+0.5, Px);
            PxVsCent_LamBar[iSub]->SetPointError(8-iCent, 0, Px_err);
        }
    }
    //...................................end of Fitting by phi........................................

    //Res factor correction
    for(int iSub=0; iSub!=2; iSub++){
        PxVsCent_Lam[iSub] = ResFactCorr(PxVsCent_Lam[iSub], Res_1VsCent[iSub]);
        PxVsCent_LamBar[iSub] = ResFactCorr(PxVsCent_LamBar[iSub], Res_1VsCent[iSub]);

        PxVsCent_Lam[iSub] = ResFactCorr(PxVsCent_Lam[iSub], Res_2VsCent[iSub]);
        PxVsCent_LamBar[iSub] = ResFactCorr(PxVsCent_LamBar[iSub], Res_2VsCent[iSub]);

        if(sinUsage){
            PxVsCent_Lam[iSub] = MultiplyLamFactor(PxVsCent_Lam[iSub], 3/alpha_Lam);
            PxVsCent_LamBar[iSub] = MultiplyLamFactor(PxVsCent_LamBar[iSub], 3/alpha_LamBar);
        }

        if(!sinUsage){
            PxVsCent_Lam[iSub] = MultiplyLamFactor(PxVsCent_Lam[iSub], 8/alpha_Lam / TMath::Pi());
            PxVsCent_LamBar[iSub] = MultiplyLamFactor(PxVsCent_LamBar[iSub], 8/alpha_LamBar / TMath::Pi());
        }
        

        PxVsCent_Lam[iSub] = MultiplyLamFactor(PxVsCent_Lam[iSub], 100);
        PxVsCent_LamBar[iSub] = MultiplyLamFactor(PxVsCent_LamBar[iSub], 100);

        PxVsCent_Lam[iSub] = MakeAverageBin(PxVsCent_Lam[iSub], 8, 9);
        PxVsCent_LamBar[iSub] = MakeAverageBin(PxVsCent_LamBar[iSub], 8, 9);

    }

    //Writing histograms to file
    TGraphErrors *PxVsCent_Tog[2];
    
    for(int iSub=0; iSub!=2; iSub++){

        //Lam+LamBar
        PxVsCent_Tog[iSub] = MakeAverage(PxVsCent_Lam[iSub], PxVsCent_LamBar[iSub]);

        //writing histograms
        PxVsCent_Lam[iSub]->SetName(Form("PxVsCent_Lam_%i", iSub));
        PxVsCent_LamBar[iSub]->SetName(Form("PxVsCent_LamBar_%i", iSub));
        PxVsCent_Tog[iSub]->SetName(Form("PxVsCent_Tog_%i", iSub));
        if(iSub == 0){
            PxVsCent_Lam[iSub]->SetTitle("PxVsCent for #Lambda, if #Psi2_{e/w} used");
            PxVsCent_LamBar[iSub]->SetTitle("PxVsCent for #bar{#Lambda}, if #Psi2_{e/w} used");
            PxVsCent_Tog[iSub]->SetTitle("PxVsCent for #Lambda + #bar{#Lambda}, if #Psi2_{e/w} used");
        }
        if(iSub == 1){
            PxVsCent_Lam[iSub]->SetTitle("PxVsCent for #Lambda, if #Psi2_{e/w} used");
            PxVsCent_LamBar[iSub]->SetTitle("PxVsCent for #bar{#Lambda}, if #Psi2_{e/w} used");
            PxVsCent_Tog[iSub]->SetTitle("PxVsCent for #Lambda + #bar{#Lambda}, if #Psi2_{e/w} used");
        }
        PxVsCent_Lam[iSub]->Write();
        PxVsCent_LamBar[iSub]->Write();
        PxVsCent_Tog[iSub]->Write();
    }
    
    output->Close();

    
    //Drawing Px Vs Centrality
    TCanvas *c_px[2];
    TCanvas *c_px_full = new TCanvas("c_px_full", "", 0, 1024, 1800, 900);
    c_px_full->Divide(2, 1, 0.0, 0.0);
    TMultiGraph *mgPxVsCent[2];
    TLegend *leg_px[2];
    for(int iSub = 0; iSub!=2; iSub++){
        c_px[iSub] = new TCanvas(Form("c_px_%i", iSub), "", 0, 1024, 1200, 1024);

        c_px[iSub]->cd();
        PxVsCent_Lam[iSub]->SetMarkerSize(3);
        PxVsCent_Lam[iSub]->SetLineWidth(2);
        PxVsCent_Lam[iSub]->SetMarkerStyle(29);
        PxVsCent_Lam[iSub]->SetMarkerColor(634);
        PxVsCent_Lam[iSub]->SetLineColor(634);

        PxVsCent_LamBar[iSub]->SetMarkerSize(3);
        PxVsCent_LamBar[iSub]->SetLineWidth(2);
        PxVsCent_LamBar[iSub]->SetMarkerStyle(30);
        PxVsCent_LamBar[iSub]->SetMarkerColor(602);
        PxVsCent_LamBar[iSub]->SetLineColor(602);

        mgPxVsCent[iSub] = new TMultiGraph();
        mgPxVsCent[iSub]->Add(PxVsCent_Lam[iSub]);
        mgPxVsCent[iSub]->Add(PxVsCent_LamBar[iSub]);

        mgPxVsCent[iSub]->GetXaxis()->SetTitle("centrality \%");
        mgPxVsCent[iSub]->GetYaxis()->SetTitle("<P_{x}sin(2#phi - 2#Psi_{2})> \%");
        mgPxVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");

        leg_px[iSub] = new TLegend();
        leg_px[iSub]->SetBorderSize(0);
        leg_px[iSub]->AddEntry(PxVsCent_Lam[iSub], "#Lambda", "p");
        leg_px[iSub]->AddEntry(PxVsCent_LamBar[iSub], "#bar{#Lambda}", "p");
        if(iSub==0) leg_px[iSub]->AddEntry((TObject*) 0, "#Psi_{e/w} was used", "");
        if(iSub==1) leg_px[iSub]->AddEntry((TObject*) 0, "#Psi_{comb} was used", "");

        leg_px[iSub]->Draw("same");

        c_px_full->cd(iSub+1);
        if(iSub==0) {gPad->SetRightMargin(0); gPad->SetLeftMargin(0.2); }
        if(iSub==1) {gPad->SetLeftMargin(0); gPad->SetRightMargin(0.2);}
        mgPxVsCent[iSub]->SetMaximum(2);
        mgPxVsCent[iSub]->SetMinimum(-2.2);
        mgPxVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        leg_px[iSub]->Draw("same");
        
    }

    
    TCanvas *c_pxtog = new TCanvas("c_pxtog", "", 0, 1024, 1200, 1024);
    c_pxtog->cd();

    
    TMultiGraph *mgPxVsCent_Tog = new TMultiGraph();

    
    for(int iCent=0; iCent!=PxVsCent_Tog[0]->GetN(); iCent++){
        PxVsCent_Tog[1]->SetPointX(iCent, cent[iCent] + 0.5);
    }

    PxVsCent_Tog[0]->SetMarkerStyle(34);
    PxVsCent_Tog[0]->SetMarkerColor(1);
    PxVsCent_Tog[0]->SetLineColor(1);

    PxVsCent_Tog[1]->SetMarkerStyle(75);
    PxVsCent_Tog[1]->SetMarkerColor(634);
    PxVsCent_Tog[1]->SetLineColor(634);

    mgPxVsCent_Tog->Add(PxVsCent_Tog[0]);
    mgPxVsCent_Tog->Add(PxVsCent_Tog[1]);

    mgPxVsCent_Tog->GetXaxis()->SetTitle("centrality \%");
    mgPxVsCent_Tog->GetYaxis()->SetTitle("<P_{x}sin(2#phi - 2#Psi)> \%");
    mgPxVsCent_Tog->SetMinimum(-1.6);
    mgPxVsCent_Tog->SetMaximum(1.4);
    mgPxVsCent_Tog->Draw("AP");
    funcOfZero->Draw("same");

    TLegend *leg_pxtog = new TLegend();
    leg_pxtog->SetBorderSize(0);

    leg_pxtog->AddEntry(PxVsCent_Tog[0], "#Lambda + #bar{#Lambda}, if #Psi_{e/w} used", "p");
    leg_pxtog->AddEntry(PxVsCent_Tog[1], "#Lambda + #bar{#Lambda}, if #Psi_{comb} used", "p");
    leg_pxtog->AddEntry((TObject*) 0, "Average in 20-50\% range centrality:", "");
    

    TGraphErrors *Px_20_50[2];

    for(int iSub = 0; iSub!=2; iSub++){
        Px_20_50[iSub] = MakeAverageBin(PxVsCent_Tog[iSub], 4, 6);
        double Px_aver = Px_20_50[iSub]->GetPointY(3);
        double Px_aver_err = Px_20_50[iSub]->GetErrorY(3);
        leg_pxtog->AddEntry(PxVsCent_Tog[iSub], Form("<P_{x}> = %.2f #pm %.2f", Px_aver, Px_aver_err), "p");
    }

    
    
    leg_pxtog->Draw("same");

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

TGraphErrors* GetCosTetaSub(TProfile* prof){
    TProfile *prof_new = (TProfile*)prof->Clone();
    TH1 *hist_average = prof_new->Rebin(prof->GetNbinsX());
    TGraphErrors *graph_cos_sub = new TGraphErrors(6);

    
    double x, ex, y, ey;

    for(int iBin=0; iBin!=6; iBin++){
        x = prof->GetBinCenter(iBin + 1);
        ex = 0.;
        y = prof->GetBinContent(iBin + 1) - hist_average->GetBinContent(1);
        ey = sqrt(pow(prof->GetBinError(iBin+1), 2) + pow(hist_average->GetBinError(1), 2));
        
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

TGraphErrors* MultiplyLamFactor(TGraphErrors *g, double coef){
    TGraphErrors *g_out = new TGraphErrors(g->GetN());

    double x, y, ex, ey;
    for(int i=0; i!=g->GetN(); i++){
        x = g->GetPointX(i);
        ex = g->GetErrorX(i);

        y = coef * g->GetPointY(i);
        ey = abs(coef) * g->GetErrorY(i);

        g_out->SetPoint(i, x, y);
        g_out->SetPointError(i, ex, ey);
    }

    g_out->SetMarkerSize(3);

    return g_out;
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

TGraphErrors* ResFactCorr(TGraphErrors* cosSubVsCent, TGraphErrors* resVsCent){
    TGraphErrors *gr = new TGraphErrors(9);

    double x, ex, y_cos, y_res, ey_cos, ey_res, ey, y;
    for(int i=0; i!=9; i++){
        x = cosSubVsCent->GetPointX(i);
        ex = 0.;
        y_cos = cosSubVsCent->GetPointY(i);
        y_res = resVsCent->GetPointY(i);
        y = y_cos/y_res;
        ey_cos = cosSubVsCent->GetErrorY(i);
        ey_res = resVsCent->GetErrorY(i);

        ey = sqrt(pow(ey_cos/y_res, 2) + pow(y_cos*ey_res/y_res/y_res, 2));

        gr->SetPoint(i, x, y);
        gr->SetPointError(i, ex, ey);
    }

    return gr;
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