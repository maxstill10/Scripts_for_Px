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

//on/off fit drawing
const bool mfitDr = false;

//order of EP
//const int nOrd = 2;
const int nOrd = 3;

TGraphErrors* ConvertProfToGraph(TProfile* prof);
TGraphErrors* GetCosTetaSub(TProfile* prof);
TGraphErrors* GetCosTetaSub(TProfile* prof, TGraphErrors* graph);
TGraphErrors* GetCosTetaSub(TProfile* prof, TProfile* prof_aver);
TGraphErrors* MultiplyLamFactor(TGraphErrors *g, double coef);
TGraphErrors* MakeAverage( TGraphErrors *g1, TGraphErrors *g2 );
TGraphErrors* ResFactCorr(TGraphErrors* cosSubVsCent, TGraphErrors* resVsCent);
TGraphErrors* MakeAverageBin(TGraphErrors *g, double n_i, double n_f);

double GetCoef(double x, TF1 *Sig, TF1 *Bg);
double MyFitFunc(double *x, double *par);

TF1 *FuncGausKFP = new TF1("FuncGausKFP", "gaus(0)+gaus(3)", 1.1, 1.13);
TF1 *LineKFP = new TF1("LineKFP", "[0] + [1]*x", 1.1, 1.13);

void PhModulationsVsCent_calc(){
    TFile *input = new TFile("../Files/14p5GeV/outputPolar_Psi1_2_3Dep_14p5.root", "read");
    TFile *output;

    if(sinUsage) output = TFile::Open(Form("../output/output_PhModulationsVsCent_forPsi%i_SinInclude_14p5.root", nOrd), "recreate");
    if(!sinUsage) output = TFile::Open(Form("../output/output_PhModulationsVsCent_forPsi%i_ConstSine_14p5.root", nOrd), "recreate");

    TH1F *CosOfDiffFlatt_1 = (TH1F*)input->Get("CosOfDiff_1");
    TH1F *CosOfDiffFlatt_2 = (TH1F*)input->Get("CosOfDiff_2");
    TH1F *CosOfDiffFlatt_3 = (TH1F*)input->Get("CosOfDiff_3");

    TGraphErrors *Res_1VsCent[2];
    TGraphErrors *Res_2VsCent[2];
    TGraphErrors *Res_3VsCent[2];

    for(int i=0; i!=2; i++){
        Res_1VsCent[i] = new TGraphErrors(9);
        Res_2VsCent[i] = new TGraphErrors(9);
        Res_3VsCent[i] = new TGraphErrors(9);
    }

    double resFlatt, resErr;
    double cent[9] = {1.5, 6.5, 14, 24, 34, 44, 54, 64, 74};
    
    for(int iCen=0; iCen!=9;iCen++){
        
        //Psi1
        resFlatt = pow(TMath::Abs(CosOfDiffFlatt_1->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(2*resFlatt))*CosOfDiffFlatt_1->GetBinError(iCen+1);
        
        Res_1VsCent[0]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_1VsCent[0]->SetPointError(8-iCen, 0., resErr);

        resFlatt = pow(2*TMath::Abs(CosOfDiffFlatt_1->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(resFlatt))*CosOfDiffFlatt_1->GetBinError(iCen+1);
        
        Res_1VsCent[1]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_1VsCent[1]->SetPointError(8-iCen, 0., resErr);

        //Psi2
        resFlatt = pow(TMath::Abs(CosOfDiffFlatt_2->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(2*resFlatt))*CosOfDiffFlatt_2->GetBinError(iCen+1);
        
        Res_2VsCent[0]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_2VsCent[0]->SetPointError(8-iCen, 0., resErr);

        resFlatt = pow(2*TMath::Abs(CosOfDiffFlatt_2->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(resFlatt))*CosOfDiffFlatt_2->GetBinError(iCen+1);
        
        Res_2VsCent[1]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_2VsCent[1]->SetPointError(8-iCen, 0., resErr);

        //Psi3
        resFlatt = pow(TMath::Abs(CosOfDiffFlatt_3->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(2*resFlatt))*CosOfDiffFlatt_3->GetBinError(iCen+1);
        
        Res_3VsCent[0]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_3VsCent[0]->SetPointError(8-iCen, 0., resErr);

        resFlatt = pow(2*TMath::Abs(CosOfDiffFlatt_3->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(resFlatt))*CosOfDiffFlatt_3->GetBinError(iCen+1);
        
        Res_3VsCent[1]->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_3VsCent[1]->SetPointError(8-iCen, 0., resErr);
        
    }

    TCanvas *c_res = new TCanvas("c_res", "c_res", 0, 1024, 1200, 1024);
    c_res->cd();

    Res_1VsCent[1]->SetMarkerStyle(23);
    Res_1VsCent[1]->SetMarkerColor(634);

    Res_1VsCent[0]->SetMarkerStyle(29);
    Res_1VsCent[0]->SetMarkerColor(633);

    Res_2VsCent[1]->SetMarkerStyle(34);
    Res_2VsCent[1]->SetMarkerColor(1);

    Res_2VsCent[0]->SetMarkerStyle(20);
    Res_2VsCent[0]->SetMarkerColor(922);

    Res_3VsCent[1]->SetMarkerStyle(33);
    Res_3VsCent[1]->SetMarkerColor(603);

    Res_3VsCent[0]->SetMarkerStyle(43);
    Res_3VsCent[0]->SetMarkerColor(601);

    TMultiGraph *mgResVsCent = new TMultiGraph();

    for(int iSub=0; iSub!=2; iSub++){
        Res_1VsCent[iSub]->SetMarkerSize(3);
        Res_2VsCent[iSub]->SetMarkerSize(3);
        Res_3VsCent[iSub]->SetMarkerSize(3);

        mgResVsCent->Add(Res_1VsCent[iSub]);
        mgResVsCent->Add(Res_2VsCent[iSub]);
        mgResVsCent->Add(Res_3VsCent[iSub]);
    }

    mgResVsCent->GetXaxis()->SetTitle("Centrality, %");
    mgResVsCent->GetYaxis()->SetTitle("Res(#Psi)");
    mgResVsCent->SetMinimum(0);
    mgResVsCent->SetMaximum(0.8);

    TLegend *leg_res = new TLegend();

    leg_res->AddEntry(Res_1VsCent[0], "Half Res(Psi1)", "p");
    leg_res->AddEntry(Res_1VsCent[1], "Full Res(Psi1)", "p");
    leg_res->AddEntry(Res_2VsCent[0], "Half Res(Psi2)", "p");
    leg_res->AddEntry(Res_2VsCent[1], "Full Res(Psi2)", "p");
    leg_res->AddEntry(Res_3VsCent[0], "Half Res(Psi3)", "p");
    leg_res->AddEntry(Res_3VsCent[1], "Full Res(Psi3)", "p");
    leg_res->SetBorderSize(0);

    mgResVsCent->Draw("AP");
    leg_res->Draw("same");

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);


    //..................................Ph calculation starts here.................................
    TProfile *prSin_diffPhiPsi1[9][6][2];
    TProfile *prSin_diffPhiPsi1_LamBar[9][6][2];

    TH1F *InvMLamDist[9][6][2];
    TH1F *InvMLamBarDist[9][6][2];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iphi=0; iphi!=6; iphi++){
            for(int iSub=0; iSub!=2; iSub++){
                if(sinUsage){
                    prSin_diffPhiPsi1[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, iSub + (nOrd-2)*2));
                    prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub + (nOrd-2)*2));
                }
                
                if(!sinUsage){
                    prSin_diffPhiPsi1[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_%i_%i_%i", iCent, iphi, iSub + (nOrd-2)*2));
                    prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, iSub + (nOrd-2)*2));
                }

                //lam invM dist
                InvMLamDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub));
                InvMLamBarDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub));
            }
        }        
    }

    //...................................Fitting by invM..................................
    TGraphErrors *grSin_diffPhiPsi1Vsdphi[9][2];
    TGraphErrors *grSin_diffPhiPsi1_LamBarVsdphi[9][2];

    for(int iCent = 0; iCent!=9; iCent++){
        for(int iSub = 0; iSub!=2; iSub++){           
            grSin_diffPhiPsi1Vsdphi[iCent][iSub] = new TGraphErrors(6);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub] = new TGraphErrors(6);
        }
    }

    TGraphErrors *grSin_diffPhiPsi1_Sub[9][6][2];
    TGraphErrors *grSin_diffPhiPsi1_LamBar_Sub[9][6][2];

    TProfile *prSin_diffPhiPsi1_aver[9][2];
    TProfile *prSin_diffPhiPsi1_LamBar_aver[9][2];

    for(int iSub=0; iSub!=2; iSub++){
        for(int iCent = 0; iCent!=9; iCent++){
            prSin_diffPhiPsi1_aver[iCent][iSub] = (TProfile*)prSin_diffPhiPsi1[iCent][0][iSub]->Clone();
            prSin_diffPhiPsi1_aver[iCent][iSub]->Reset();
            prSin_diffPhiPsi1_LamBar_aver[iCent][iSub] = (TProfile*)prSin_diffPhiPsi1_LamBar[iCent][0][iSub]->Clone();
            prSin_diffPhiPsi1_LamBar_aver[iCent][iSub]->Reset();
            for(int iphi = 0; iphi!=6; iphi++){
                
                prSin_diffPhiPsi1_aver[iCent][iSub]->Add(prSin_diffPhiPsi1[iCent][iphi][iSub]);            
                prSin_diffPhiPsi1_LamBar_aver[iCent][iSub]->Add(prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub]);
            }
        }
    }



    

    
    double Pz, Pz_err;
    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.1145, 1.1165);
	FitFuncKFP->SetParLimits(2, 0.0005, 0.0011);    
    FitFuncKFP->SetParLimits(4, 1.114, 1.117);
	FitFuncKFP->SetParLimits(5, 0.0015, 0.005);

    TF1 *FitCosTetaFunc;

    if(!mfitDr){
        gROOT->SetBatch(kTRUE);
    }
    
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
            for(int iphi = 0; iphi!=6; iphi++){

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

                FitCosTetaFunc->FixParameter(1, 0.);
                FitCosTetaFunc->FixParameter(2, 0.);

                grSin_diffPhiPsi1_Sub[iCent][iphi][iSub] = ConvertProfToGraph(prSin_diffPhiPsi1[iCent][iphi][iSub]);
                
                grSin_diffPhiPsi1_Sub[iCent][iphi][iSub]->Fit(FitCosTetaFunc);
                grSin_diffPhiPsi1_Sub[iCent][iphi][iSub]->Draw("AP");

                Pz = FitCosTetaFunc->GetParameter(0);
                Pz_err = FitCosTetaFunc->GetParError(0);

                double phi = 2*TMath::Pi() / nOrd /12 + iphi*2*TMath::Pi() / nOrd /6;
                grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetPoint(iphi, phi, Pz);
                grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetPointError(iphi, 0., Pz_err);
                
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

                FitCosTetaFunc->FixParameter(1, 0.);
                FitCosTetaFunc->FixParameter(2, 0.);

                //Subtracting
                grSin_diffPhiPsi1_LamBar_Sub[iCent][iphi][iSub] = ConvertProfToGraph(prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub]);
                
                grSin_diffPhiPsi1_LamBar_Sub[iCent][iphi][iSub]->Fit(FitCosTetaFunc);
                grSin_diffPhiPsi1_LamBar_Sub[iCent][iphi][iSub]->Draw("AP");

                Pz = FitCosTetaFunc->GetParameter(0);
                Pz_err = FitCosTetaFunc->GetParError(0);

                phi = 2*TMath::Pi() / nOrd /12 + iphi*2*TMath::Pi() / nOrd /6;
                grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetPoint(iphi, phi, Pz);
                grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetPointError(iphi, 0., Pz_err);
            }//for(int iphi = 0; iphi!=6; iphi++)
            
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMarkerSize(2);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetLineWidth(2);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMarkerStyle(20);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMarkerColor(1);

            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMarkerSize(2);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetLineWidth(2);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMarkerStyle(20);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMarkerColor(1);
        }//for(int iCent = 0; iCent!=9; iCent++)
    }//for(int iSub=0; iSub!=2; iSub++)

    //...................................end of Fitting by invM..................................

    if(!mfitDr){
        gROOT->SetBatch(kFALSE);
    }

    //...................................Fitting by phi............................................

    TGraphErrors *PhVsCent_Lam[2];
    TGraphErrors *PhVsCent_LamBar[2];
    TGraphErrors *PhSinVsCent_Lam[2];
    TGraphErrors *PhSinVsCent_LamBar[2];
    TGraphErrors *PhCosVsCent_Lam[2];
    TGraphErrors *PhCosVsCent_LamBar[2];

    for(int iSub=0; iSub!=2; iSub++){
        PhVsCent_Lam[iSub] = new TGraphErrors(9);
        PhVsCent_LamBar[iSub] = new TGraphErrors(9);
        PhSinVsCent_Lam[iSub] = new TGraphErrors(9);
        PhSinVsCent_LamBar[iSub] = new TGraphErrors(9);
        PhCosVsCent_Lam[iSub] = new TGraphErrors(9);
        PhCosVsCent_LamBar[iSub] = new TGraphErrors(9);
    }

    TCanvas *c[4];
    for(int i=0; i!=4; i++){
        c[i] = new TCanvas(Form("c_%i", i), "", 0, 1024, 1100, 1024);
        c[i]->Divide(3,3,0,0);
    }
    c[0]->SetTitle("Lambda, Psi_e/w");
    c[1]->SetTitle("Lambda, Psi_comb");
    c[2]->SetTitle("LambdaBar, Psi_e/w");
    c[3]->SetTitle("LambdaBar, Psi_comb");

    TGaxis::SetMaxDigits(1);
    //gStyle->SetOptFit(0001);


    int cent_lim[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

    TF1 *fitFuncOfCos = new TF1("fitFuncOfCos", Form("[0] + 2*[1]*sin(%i*x) + 2*[2]*cos(%i*x)", nOrd, nOrd), 0, 2*TMath::Pi()/nOrd);

    //fitFuncOfCos->FixParameter(0, 0);
    fitFuncOfCos->SetParLimits(1, -1, 1);
    //fitFuncOfCos->FixParameter(2, 0);

    fitFuncOfCos->SetLineWidth(2);

    double Ph, Ph_err;

    TText tFit;
    tFit.SetTextFont(42);

    TLatex latFit;
    latFit.SetTextFont(42);
    latFit.SetTextSize(0.1);
    
    for(int iCent=0; iCent!=9; iCent++){
        for(int iSub=0; iSub!=2; iSub++){

            //Set Y titles
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetYaxis()->SetTitle("<sin(#Psi_{1} - #phi^{*}_{p})>");
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetTitle("<sin(#Psi_{1} - #phi^{*}_{#bar{p}})>");
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetXaxis()->SetTitle(Form("#phi_{#Lambda} - #Psi_{%i}", nOrd));
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetXaxis()->SetTitle(Form("#phi_{#bar{#Lambda}} - #Psi_{%i}", nOrd));
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetXaxis()->SetTitleSize(0.08);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetXaxis()->SetTitleSize(0.08);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetXaxis()->SetLabelSize(0.08);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetXaxis()->SetLabelSize(0.08);
            if(nOrd == 3){
                grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetXaxis()->SetNdivisions(505); 
                grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetXaxis()->SetNdivisions(505);
            }

            if(nOrd == 2){
                if(iSub == 0){grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMinimum(-0.0005); grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMaximum(0.004);}
                if(iSub == 1){grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMinimum(-0.0001); grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMaximum(0.004);}
                if(iSub == 0){grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMinimum(-0.002); grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMaximum(0.007);}
                if(iSub == 1){grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMinimum(-0.005); grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMaximum(0.007);}
            }

            if(nOrd == 3){
                if(iSub == 0){grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMinimum(-0.0005); grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMaximum(0.004);}
                if(iSub == 1){grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMinimum(-0.0003); grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetMaximum(0.005);}
                if(iSub == 0){grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMinimum(-0.002); grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMaximum(0.007);}
                if(iSub == 1){grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMinimum(-0.005); grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetMaximum(0.008);}
            }

            //Set graph titles
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->SetTitle("");
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->SetTitle("");
            

            //Drawing           

            c[iSub]->cd(9-iCent);
            gPad->SetBottomMargin(0.25);
            gPad->SetLeftMargin(0.25);
            gPad->SetRightMargin(0.02);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetYaxis()->SetTitleOffset(1);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetYaxis()->SetTitleSize(0.08);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->GetYaxis()->SetLabelSize(0.08);
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->Fit("fitFuncOfCos");
            grSin_diffPhiPsi1Vsdphi[iCent][iSub]->Draw("AP");

            if(nOrd==2){
                if(iCent == 4 && iSub == 0) latFit.DrawLatex(1, 0, Form("#bar{#Lambda}, #Psi_{%i,ew}", nOrd));
                if(iCent == 4 && iSub == 1) latFit.DrawLatex(0.5, 0.0002, Form("#bar{#Lambda}, #Psi_{%i,comb}", nOrd));
                if(iSub == 0) tFit.DrawText(0.2, 0.0035, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
                if(iSub == 1) tFit.DrawText(0.2, 0.0035, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
            }

            if(nOrd==3){
                if(iCent == 4 && iSub == 0) latFit.DrawLatex(0.6, 0., Form("#bar{#Lambda}, #Psi_{%i,ew}", nOrd));
                if(iCent == 4 && iSub == 1) latFit.DrawLatex(0.6, 0.0002, Form("#bar{#Lambda}, #Psi_{%i,comb}", nOrd));
                if(iSub == 0) tFit.DrawText(0.8, 0.0035, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
                if(iSub == 1) tFit.DrawText(1, 0.0045, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
            }
            
            Ph = fitFuncOfCos->GetParameter(0);
            Ph_err = fitFuncOfCos->GetParError(0);
            PhVsCent_Lam[iSub]->SetPoint(8-iCent, cent[8-iCent], Ph);
            PhVsCent_Lam[iSub]->SetPointError(8-iCent, 0, Ph_err);
            PhSinVsCent_Lam[iSub]->SetPoint(8-iCent, cent[8-iCent], fitFuncOfCos->GetParameter(1));
            PhSinVsCent_Lam[iSub]->SetPointError(8-iCent, 0, fitFuncOfCos->GetParError(1));
            PhCosVsCent_Lam[iSub]->SetPoint(8-iCent, cent[8-iCent], fitFuncOfCos->GetParameter(2));
            PhCosVsCent_Lam[iSub]->SetPointError(8-iCent, 0, fitFuncOfCos->GetParError(2));


            c[iSub+2]->cd(9-iCent);
            
            gPad->SetBottomMargin(0.25);
            gPad->SetLeftMargin(0.25);
            gPad->SetRightMargin(0.02);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetTitleOffset(1);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetTitleSize(0.08);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetLabelSize(0.08);
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->Draw("AP");
            grSin_diffPhiPsi1_LamBarVsdphi[iCent][iSub]->Fit("fitFuncOfCos");

            if(nOrd==2){
                if(iCent == 4 && iSub == 0) latFit.DrawLatex(1, -0.001, Form("#bar{#Lambda}, #Psi_{%i,ew}", nOrd));
                if(iCent == 4 && iSub == 1) latFit.DrawLatex(1, -0.001, Form("#bar{#Lambda}, #Psi_{%i,comb}", nOrd));
                if(iSub == 0) tFit.DrawText(1.5, 0.0063, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
                if(iSub == 1) tFit.DrawText(1.6, -0.004, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
            }

            if(nOrd==3){
                if(iCent == 4 && iSub == 0) latFit.DrawLatex(0.6, -0.001, Form("#bar{#Lambda}, #Psi_{%i,ew}", nOrd));
                if(iCent == 4 && iSub == 1) latFit.DrawLatex(0.6, -0.002, Form("#bar{#Lambda}, #Psi_{%i,comb}", nOrd));
                if(iSub == 0) tFit.DrawText(1, 0.006, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
                if(iSub == 1) tFit.DrawText(1, 0.007, Form("Centrality, %i-%i %%", cent_lim[8-iCent], cent_lim[9-iCent]));
            }

            Ph = fitFuncOfCos->GetParameter(0);
            Ph_err = fitFuncOfCos->GetParError(0);

            PhVsCent_LamBar[iSub]->SetPoint(8-iCent, cent[8-iCent]+0.5, Ph);
            PhVsCent_LamBar[iSub]->SetPointError(8-iCent, 0, Ph_err);
            PhSinVsCent_LamBar[iSub]->SetPoint(8-iCent, cent[8-iCent]+0.5, fitFuncOfCos->GetParameter(1));
            PhSinVsCent_LamBar[iSub]->SetPointError(8-iCent, 0, fitFuncOfCos->GetParError(1));
            PhCosVsCent_LamBar[iSub]->SetPoint(8-iCent, cent[8-iCent]+0.5, fitFuncOfCos->GetParameter(2));
            PhCosVsCent_LamBar[iSub]->SetPointError(8-iCent, 0, fitFuncOfCos->GetParError(2));
        }
    }
    //...................................end of Fitting by phi........................................

    //Res factor correction
    for(int iSub=0; iSub!=2; iSub++){
        PhVsCent_Lam[iSub] = ResFactCorr(PhVsCent_Lam[iSub], Res_1VsCent[1]);
        PhVsCent_LamBar[iSub] = ResFactCorr(PhVsCent_LamBar[iSub], Res_1VsCent[1]);
        PhSinVsCent_Lam[iSub] = ResFactCorr(PhSinVsCent_Lam[iSub], Res_1VsCent[1]);
        PhSinVsCent_LamBar[iSub] = ResFactCorr(PhSinVsCent_LamBar[iSub], Res_1VsCent[1]);
        PhCosVsCent_Lam[iSub] = ResFactCorr(PhCosVsCent_Lam[iSub], Res_1VsCent[1]);
        PhCosVsCent_LamBar[iSub] = ResFactCorr(PhCosVsCent_LamBar[iSub], Res_1VsCent[1]);

        if(nOrd == 2)
        {
            PhSinVsCent_Lam[iSub] = ResFactCorr(PhSinVsCent_Lam[iSub], Res_2VsCent[iSub]);
            PhSinVsCent_LamBar[iSub] = ResFactCorr(PhSinVsCent_LamBar[iSub], Res_2VsCent[iSub]);
            PhCosVsCent_Lam[iSub] = ResFactCorr(PhCosVsCent_Lam[iSub], Res_2VsCent[iSub]);
            PhCosVsCent_LamBar[iSub] = ResFactCorr(PhCosVsCent_LamBar[iSub], Res_2VsCent[iSub]);
        }

        if(nOrd == 3)
        {
            PhSinVsCent_Lam[iSub] = ResFactCorr(PhSinVsCent_Lam[iSub], Res_3VsCent[iSub]);
            PhSinVsCent_LamBar[iSub] = ResFactCorr(PhSinVsCent_LamBar[iSub], Res_3VsCent[iSub]);
            PhCosVsCent_Lam[iSub] = ResFactCorr(PhCosVsCent_Lam[iSub], Res_3VsCent[iSub]);
            PhCosVsCent_LamBar[iSub] = ResFactCorr(PhCosVsCent_LamBar[iSub], Res_3VsCent[iSub]);
        }

        if(sinUsage){
            PhVsCent_Lam[iSub] = MultiplyLamFactor(PhVsCent_Lam[iSub], 3/alpha_Lam);
            PhVsCent_LamBar[iSub] = MultiplyLamFactor(PhVsCent_LamBar[iSub], 3/alpha_Lam);
            PhSinVsCent_Lam[iSub] = MultiplyLamFactor(PhSinVsCent_Lam[iSub], 3/alpha_Lam);
            PhSinVsCent_LamBar[iSub] = MultiplyLamFactor(PhSinVsCent_LamBar[iSub], 3/alpha_Lam);
            PhCosVsCent_Lam[iSub] = MultiplyLamFactor(PhCosVsCent_Lam[iSub], 3/alpha_Lam);
            PhCosVsCent_LamBar[iSub] = MultiplyLamFactor(PhCosVsCent_LamBar[iSub], 3/alpha_Lam);
        }

        if(!sinUsage){
            PhVsCent_Lam[iSub] = MultiplyLamFactor(PhVsCent_Lam[iSub], 8/alpha_Lam / TMath::Pi());
            PhVsCent_LamBar[iSub] = MultiplyLamFactor(PhVsCent_LamBar[iSub], 8/alpha_Lam / TMath::Pi());
            PhSinVsCent_Lam[iSub] = MultiplyLamFactor(PhSinVsCent_Lam[iSub], 8/alpha_Lam / TMath::Pi());
            PhSinVsCent_LamBar[iSub] = MultiplyLamFactor(PhSinVsCent_LamBar[iSub], 8/alpha_Lam / TMath::Pi());
            PhCosVsCent_Lam[iSub] = MultiplyLamFactor(PhCosVsCent_Lam[iSub], 8/alpha_Lam / TMath::Pi());
            PhCosVsCent_LamBar[iSub] = MultiplyLamFactor(PhCosVsCent_LamBar[iSub], 8/alpha_Lam / TMath::Pi());
        }
        

        PhVsCent_Lam[iSub] = MultiplyLamFactor(PhVsCent_Lam[iSub], 100);
        PhVsCent_LamBar[iSub] = MultiplyLamFactor(PhVsCent_LamBar[iSub], 100);
        PhSinVsCent_Lam[iSub] = MultiplyLamFactor(PhSinVsCent_Lam[iSub], 100);
        PhSinVsCent_LamBar[iSub] = MultiplyLamFactor(PhSinVsCent_LamBar[iSub], 100);
        PhCosVsCent_Lam[iSub] = MultiplyLamFactor(PhCosVsCent_Lam[iSub], 100);
        PhCosVsCent_LamBar[iSub] = MultiplyLamFactor(PhCosVsCent_LamBar[iSub], 100);

        PhSinVsCent_Lam[iSub] = MakeAverageBin(PhSinVsCent_Lam[iSub], 8, 9);
        PhSinVsCent_LamBar[iSub] = MakeAverageBin(PhSinVsCent_LamBar[iSub], 8, 9);
        PhCosVsCent_Lam[iSub] = MakeAverageBin(PhCosVsCent_Lam[iSub], 8, 9);
        PhCosVsCent_LamBar[iSub] = MakeAverageBin(PhCosVsCent_LamBar[iSub], 8, 9);

    }

    //Writing histograms to file
    TGraphErrors *PhVsCent_Tog[2];
    TGraphErrors *PhSinVsCent_Tog[2];
    TGraphErrors *PhCosVsCent_Tog[2];
    
    for(int iSub=0; iSub!=2; iSub++){

        //PhVsCent_Lam[iSub]=MakeAverageBin(PhVsCent_Lam[iSub], 7, 8);
        //PhVsCent_LamBar[iSub]=MakeAverageBin(PhVsCent_LamBar[iSub], 7, 8);

        //Lam+LamBar
        PhVsCent_Tog[iSub] = MakeAverage(PhVsCent_Lam[iSub], PhVsCent_LamBar[iSub]);
        PhSinVsCent_Tog[iSub] = MakeAverage(PhSinVsCent_Lam[iSub], PhSinVsCent_LamBar[iSub]);
        PhCosVsCent_Tog[iSub] = MakeAverage(PhCosVsCent_Lam[iSub], PhCosVsCent_LamBar[iSub]);

        //writing histograms
        PhVsCent_Lam[iSub]->SetName(Form("PhVsCent_Lam_%i", iSub));
        PhVsCent_LamBar[iSub]->SetName(Form("PhVsCent_LamBar_%i", iSub));
        PhVsCent_Tog[iSub]->SetName(Form("PhVsCent_Tog_%i", iSub));
        PhSinVsCent_Lam[iSub]->SetName(Form("PhSinVsCent_Lam_%i", iSub));
        PhSinVsCent_LamBar[iSub]->SetName(Form("PhSinVsCent_LamBar_%i", iSub));
        PhSinVsCent_Tog[iSub]->SetName(Form("PhSinVsCent_Tog_%i", iSub));
        PhCosVsCent_Lam[iSub]->SetName(Form("PhCosVsCent_Lam_%i", iSub));
        PhCosVsCent_LamBar[iSub]->SetName(Form("PhCosVsCent_LamBar_%i", iSub));
        PhCosVsCent_Tog[iSub]->SetName(Form("PhCosVsCent_Tog_%i", iSub));

        if(iSub == 0){
            PhVsCent_Lam[iSub]->SetTitle(Form("PhVsCent for #Lambda, if #Psi%i_{e/w} used", nOrd));
            PhVsCent_LamBar[iSub]->SetTitle(Form("PhVsCent for #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd));
            PhVsCent_Tog[iSub]->SetTitle(Form("PhVsCent for #Lambda + #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd));

            PhSinVsCent_Lam[iSub]->SetTitle(Form("P_{h}sin(%i#Delta#phi) Vs Cent for #Lambda, if #Psi%i_{e/w} used", nOrd, nOrd));
            PhSinVsCent_LamBar[iSub]->SetTitle(Form("P_{h}sin(%i#Delta#phi) Vs Cent for #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd, nOrd));
            PhSinVsCent_Tog[iSub]->SetTitle(Form("P_{h}sin(%i#Delta#phi) Vs Cent for #Lambda + #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd, nOrd));

            PhCosVsCent_Lam[iSub]->SetTitle(Form("P_{h}cos(%i#Delta#phi) Vs Cent for #Lambda, if #Psi%i_{e/w} used", nOrd, nOrd));
            PhCosVsCent_LamBar[iSub]->SetTitle(Form("P_{h}cos(%i#Delta#phi) Vs Cent for #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd, nOrd));
            PhCosVsCent_Tog[iSub]->SetTitle(Form("P_{h}cos(%i#Delta#phi) Vs Cent for #Lambda + #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd, nOrd));
        }
        if(iSub == 1){
            PhVsCent_Lam[iSub]->SetTitle(Form("PhVsCent for #Lambda, if #Psi%i_{comb} used", nOrd));
            PhVsCent_LamBar[iSub]->SetTitle(Form("PhVsCent for #bar{#Lambda}, if #Psi%i_{comb} used", nOrd));
            PhVsCent_Tog[iSub]->SetTitle(Form("PhVsCent for #Lambda + #bar{#Lambda}, if #Psi%i_{comb} used", nOrd));

            PhSinVsCent_Lam[iSub]->SetTitle(Form("P_{h}sin(%i#Delta#phi) Vs Cent for #Lambda, if #Psi%i_{comb} used", nOrd, nOrd));
            PhSinVsCent_LamBar[iSub]->SetTitle(Form("P_{h}sin(%i#Delta#phi) Vs Cent for #bar{#Lambda}, if #Psi%i_{comb} used", nOrd, nOrd));
            PhSinVsCent_Tog[iSub]->SetTitle(Form("P_{h}sin(%i#Delta#phi) Vs Cent for #Lambda + #bar{#Lambda}, if #Psi%i_{comb} used", nOrd, nOrd));

            PhCosVsCent_Lam[iSub]->SetTitle(Form("P_{h}cos(%i#Delta#phi) Vs Cent for #Lambda, if #Psi%i_{comb} used", nOrd, nOrd));
            PhCosVsCent_LamBar[iSub]->SetTitle(Form("P_{h}cos(%i#Delta#phi) Vs Cent for #bar{#Lambda}, if #Psi%i_{comb} used", nOrd, nOrd));
            PhCosVsCent_Tog[iSub]->SetTitle(Form("P_{h}cos(%i#Delta#phi) Vs Cent for #Lambda + #bar{#Lambda}, if #Psi%i_{comb} used", nOrd, nOrd));
        }
        PhVsCent_Lam[iSub]->Write();
        PhVsCent_LamBar[iSub]->Write();
        PhVsCent_Tog[iSub]->Write();
        PhSinVsCent_Lam[iSub]->Write();
        PhSinVsCent_LamBar[iSub]->Write();
        PhSinVsCent_Tog[iSub]->Write();
        PhCosVsCent_Lam[iSub]->Write();
        PhCosVsCent_LamBar[iSub]->Write();
        PhCosVsCent_Tog[iSub]->Write();
    }
    
    output->Close();

    

   
    //Drawing Ph Vs Centrality
    TCanvas *c_Ph[2];
    TCanvas *c_Ph_full = new TCanvas("c_Ph_full", "", 0, 1024, 1800, 900);
    c_Ph_full->Divide(2, 1, 0.0, 0.0);
    TMultiGraph *mgPhVsCent[2];
    TLegend *leg_Ph[2];
    for(int iSub = 0; iSub!=2; iSub++){
        c_Ph[iSub] = new TCanvas(Form("c_Ph_%i", iSub), "", 0, 1024, 1200, 1024);

        c_Ph[iSub]->cd();
        PhVsCent_Lam[iSub]->SetMarkerSize(3);
        PhVsCent_Lam[iSub]->SetLineWidth(2);
        PhVsCent_Lam[iSub]->SetMarkerStyle(29);
        PhVsCent_Lam[iSub]->SetMarkerColor(634);
        PhVsCent_Lam[iSub]->SetLineColor(634);

        PhVsCent_LamBar[iSub]->SetMarkerSize(3);
        PhVsCent_LamBar[iSub]->SetLineWidth(2);
        PhVsCent_LamBar[iSub]->SetMarkerStyle(30);
        PhVsCent_LamBar[iSub]->SetMarkerColor(602);
        PhVsCent_LamBar[iSub]->SetLineColor(602);

        mgPhVsCent[iSub] = new TMultiGraph();
        mgPhVsCent[iSub]->Add(PhVsCent_Lam[iSub]);
        mgPhVsCent[iSub]->Add(PhVsCent_LamBar[iSub]);

        mgPhVsCent[iSub]->GetXaxis()->SetTitle("centrality %");
        mgPhVsCent[iSub]->GetYaxis()->SetTitle("<P_{h}> %");
        mgPhVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");

        leg_Ph[iSub] = new TLegend();
        leg_Ph[iSub]->SetBorderSize(0);
        leg_Ph[iSub]->AddEntry(PhVsCent_Lam[iSub], "#Lambda", "p");
        leg_Ph[iSub]->AddEntry(PhVsCent_LamBar[iSub], "#bar{#Lambda}", "p");
        if(iSub==0) leg_Ph[iSub]->AddEntry((TObject*) 0, "#Psi_{e/w} was used", "");
        if(iSub==1) leg_Ph[iSub]->AddEntry((TObject*) 0, "#Psi_{comb} was used", "");

        leg_Ph[iSub]->Draw("same");

        c_Ph_full->cd(iSub+1);
        if(iSub==0) {gPad->SetRightMargin(0); gPad->SetLeftMargin(0.2); }
        if(iSub==1) {gPad->SetLeftMargin(0); gPad->SetRightMargin(0.2);}
        //mgPhVsCent[iSub]->SetMaximum(4);
        //mgPhVsCent[iSub]->SetMinimum(-4);
        mgPhVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        leg_Ph[iSub]->Draw("same");
        
    }

    
    TCanvas *c_Phtog = new TCanvas("c_Phtog", "", 0, 1024, 1200, 1024);
    c_Phtog->cd();

    
    TMultiGraph *mgPhVsCent_Tog = new TMultiGraph();

    
    for(int iCent=0; iCent!=PhVsCent_Tog[0]->GetN(); iCent++){
        PhVsCent_Tog[1]->SetPointX(iCent, PhVsCent_Tog[0]->GetPointX(iCent) + 0.5);
    }

    PhVsCent_Tog[0]->SetMarkerStyle(34);
    PhVsCent_Tog[0]->SetMarkerColor(1);
    PhVsCent_Tog[0]->SetLineColor(1);

    PhVsCent_Tog[1]->SetMarkerStyle(75);
    PhVsCent_Tog[1]->SetMarkerColor(634);
    PhVsCent_Tog[1]->SetLineColor(634);

    mgPhVsCent_Tog->Add(PhVsCent_Tog[0]);
    mgPhVsCent_Tog->Add(PhVsCent_Tog[1]);

    mgPhVsCent_Tog->GetXaxis()->SetTitle("centrality %");
    mgPhVsCent_Tog->GetYaxis()->SetTitle("<P_{h}> %");
    //mgPhVsCent_Tog->SetMinimum(-2);
    //mgPhVsCent_Tog->SetMaximum(2);
    mgPhVsCent_Tog->Draw("AP");
    funcOfZero->Draw("same");

    TLegend *leg_Phtog = new TLegend();
    leg_Phtog->SetBorderSize(0);

    leg_Phtog->AddEntry(PhVsCent_Tog[0], "#Lambda + #bar{#Lambda}, if #Psi_{e/w} used", "p");
    leg_Phtog->AddEntry(PhVsCent_Tog[1], "#Lambda + #bar{#Lambda}, if #Psi_{comb} used", "p");
    leg_Phtog->AddEntry((TObject*) 0, "Average in 20-50\% range centrality:", "");
    

    TGraphErrors *Ph_20_50[2];

    for(int iSub = 0; iSub!=2; iSub++){
        Ph_20_50[iSub] = MakeAverageBin(PhVsCent_Tog[iSub], 4, 6);
        double Ph_aver = Ph_20_50[iSub]->GetPointY(3);
        double Ph_aver_err = Ph_20_50[iSub]->GetErrorY(3);
        leg_Phtog->AddEntry(PhVsCent_Tog[iSub], Form("<P_{x}> = %.2f #pm %.2f", Ph_aver, Ph_aver_err), "p");
    }    
    
    leg_Phtog->Draw("same");


    //Drawing PhSin Vs Centrality
    TCanvas *c_PhSin[2];
    TCanvas *c_PhSin_full = new TCanvas("c_PhSin_full", "", 0, 1024, 1800, 900);
    c_PhSin_full->Divide(2, 1, 0.0, 0.0);
    TMultiGraph *mgPhSinVsCent[2];
    TLegend *leg_PhSin[2];
    for(int iSub = 0; iSub!=2; iSub++){
        c_PhSin[iSub] = new TCanvas(Form("c_PhSin_%i", iSub), "", 0, 1024, 1200, 1024);

        c_PhSin[iSub]->cd();
        PhSinVsCent_Lam[iSub]->SetMarkerSize(3);
        PhSinVsCent_Lam[iSub]->SetLineWidth(2);
        PhSinVsCent_Lam[iSub]->SetMarkerStyle(29);
        PhSinVsCent_Lam[iSub]->SetMarkerColor(634);
        PhSinVsCent_Lam[iSub]->SetLineColor(634);

        PhSinVsCent_LamBar[iSub]->SetMarkerSize(3);
        PhSinVsCent_LamBar[iSub]->SetLineWidth(2);
        PhSinVsCent_LamBar[iSub]->SetMarkerStyle(30);
        PhSinVsCent_LamBar[iSub]->SetMarkerColor(602);
        PhSinVsCent_LamBar[iSub]->SetLineColor(602);

        mgPhSinVsCent[iSub] = new TMultiGraph();
        mgPhSinVsCent[iSub]->Add(PhSinVsCent_Lam[iSub]);
        mgPhSinVsCent[iSub]->Add(PhSinVsCent_LamBar[iSub]);

        mgPhSinVsCent[iSub]->GetXaxis()->SetTitle("centrality %");
        mgPhSinVsCent[iSub]->GetYaxis()->SetTitle(Form("<P_{h}sin(%i(#phi - #Psi_{%i}))> %%", nOrd, nOrd));
        mgPhSinVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");

        leg_PhSin[iSub] = new TLegend();
        leg_PhSin[iSub]->SetBorderSize(0);
        leg_PhSin[iSub]->AddEntry(PhSinVsCent_Lam[iSub], "#Lambda", "p");
        leg_PhSin[iSub]->AddEntry(PhSinVsCent_LamBar[iSub], "#bar{#Lambda}", "p");
        if(iSub==0) leg_PhSin[iSub]->AddEntry((TObject*) 0, "#Psi_{e/w} was used", "");
        if(iSub==1) leg_PhSin[iSub]->AddEntry((TObject*) 0, "#Psi_{comb} was used", "");

        leg_PhSin[iSub]->Draw("same");

        c_PhSin_full->cd(iSub+1);
        if(iSub==0) {gPad->SetRightMargin(0); gPad->SetLeftMargin(0.2); }
        if(iSub==1) {gPad->SetLeftMargin(0); gPad->SetRightMargin(0.2);}
        //mgPhSinVsCent[iSub]->SetMaximum(7.5);
        //mgPhSinVsCent[iSub]->SetMinimum(-6);
        mgPhSinVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        leg_PhSin[iSub]->Draw("same");
        
    }

    
    TCanvas *c_PhSintog = new TCanvas("c_PhSintog", "", 0, 1024, 1200, 1024);
    c_PhSintog->cd();

    
    TMultiGraph *mgPhSinVsCent_Tog = new TMultiGraph();

    
    for(int iCent=0; iCent!=PhSinVsCent_Tog[0]->GetN(); iCent++){
        PhSinVsCent_Tog[1]->SetPointX(iCent, PhSinVsCent_Tog[0]->GetPointX(iCent) + 0.5);
    }

    PhSinVsCent_Tog[0]->SetMarkerStyle(34);
    PhSinVsCent_Tog[0]->SetMarkerColor(1);
    PhSinVsCent_Tog[0]->SetLineColor(1);

    PhSinVsCent_Tog[1]->SetMarkerStyle(75);
    PhSinVsCent_Tog[1]->SetMarkerColor(634);
    PhSinVsCent_Tog[1]->SetLineColor(634);

    mgPhSinVsCent_Tog->Add(PhSinVsCent_Tog[0]);
    mgPhSinVsCent_Tog->Add(PhSinVsCent_Tog[1]);

    mgPhSinVsCent_Tog->GetXaxis()->SetTitle("centrality %");
    mgPhSinVsCent_Tog->GetYaxis()->SetTitle(Form("<P_{h}sin(%i(#phi - #Psi_{%i}))> %%", nOrd, nOrd));
    //mgPhSinVsCent_Tog->SetMinimum(-2);
    //mgPhSinVsCent_Tog->SetMaximum(2);
    mgPhSinVsCent_Tog->Draw("AP");
    funcOfZero->Draw("same");

    TLegend *leg_PhSintog = new TLegend();
    leg_PhSintog->SetBorderSize(0);

    leg_PhSintog->AddEntry(PhVsCent_Tog[0], "#Lambda + #bar{#Lambda}, if #Psi_{e/w} used", "p");
    leg_PhSintog->AddEntry(PhVsCent_Tog[1], "#Lambda + #bar{#Lambda}, if #Psi_{comb} used", "p");
    leg_PhSintog->AddEntry((TObject*) 0, "Average in 20-50\% range centrality:", "");
    

    TGraphErrors *PhSin_20_50[2];

    for(int iSub = 0; iSub!=2; iSub++){
        PhSin_20_50[iSub] = MakeAverageBin(PhSinVsCent_Tog[iSub], 4, 6);
        double PhSin_aver = PhSin_20_50[iSub]->GetPointY(3);
        double PhSin_aver_err = PhSin_20_50[iSub]->GetErrorY(3);
        leg_PhSintog->AddEntry(PhSinVsCent_Tog[iSub], Form("<P_{x}> = %.2f #pm %.2f", PhSin_aver, PhSin_aver_err), "p");
    }    
    
    leg_PhSintog->Draw("same");


    //Drawing PhCos Vs Centrality
    TCanvas *c_PhCos[2];
    TCanvas *c_PhCos_full = new TCanvas("c_PhCos_full", "", 0, 1024, 1800, 900);
    c_PhCos_full->Divide(2, 1, 0.0, 0.0);
    TMultiGraph *mgPhCosVsCent[2];
    TLegend *leg_PhCos[2];
    for(int iSub = 0; iSub!=2; iSub++){
        c_PhCos[iSub] = new TCanvas(Form("c_PhCos_%i", iSub), "", 0, 1024, 1200, 1024);

        c_PhCos[iSub]->cd();
        PhCosVsCent_Lam[iSub]->SetMarkerSize(3);
        PhCosVsCent_Lam[iSub]->SetLineWidth(2);
        PhCosVsCent_Lam[iSub]->SetMarkerStyle(29);
        PhCosVsCent_Lam[iSub]->SetMarkerColor(634);
        PhCosVsCent_Lam[iSub]->SetLineColor(634);

        PhCosVsCent_LamBar[iSub]->SetMarkerSize(3);
        PhCosVsCent_LamBar[iSub]->SetLineWidth(2);
        PhCosVsCent_LamBar[iSub]->SetMarkerStyle(30);
        PhCosVsCent_LamBar[iSub]->SetMarkerColor(602);
        PhCosVsCent_LamBar[iSub]->SetLineColor(602);

        mgPhCosVsCent[iSub] = new TMultiGraph();
        mgPhCosVsCent[iSub]->Add(PhCosVsCent_Lam[iSub]);
        mgPhCosVsCent[iSub]->Add(PhCosVsCent_LamBar[iSub]);

        mgPhCosVsCent[iSub]->GetXaxis()->SetTitle("centrality %");
        mgPhCosVsCent[iSub]->GetYaxis()->SetTitle(Form("<P_{h}cos(%i(#phi - #Psi_{%i}))> %%", nOrd, nOrd));
        mgPhCosVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");

        leg_PhCos[iSub] = new TLegend();
        leg_PhCos[iSub]->SetBorderSize(0);
        leg_PhCos[iSub]->AddEntry(PhCosVsCent_Lam[iSub], "#Lambda", "p");
        leg_PhCos[iSub]->AddEntry(PhCosVsCent_LamBar[iSub], "#bar{#Lambda}", "p");
        if(iSub==0) leg_PhCos[iSub]->AddEntry((TObject*) 0, "#Psi_{e/w} was used", "");
        if(iSub==1) leg_PhCos[iSub]->AddEntry((TObject*) 0, "#Psi_{comb} was used", "");

        leg_PhCos[iSub]->Draw("same");

        c_PhCos_full->cd(iSub+1);
        if(iSub==0) {gPad->SetRightMargin(0); gPad->SetLeftMargin(0.2); }
        if(iSub==1) {gPad->SetLeftMargin(0); gPad->SetRightMargin(0.2);}
        //mgPhCosVsCent[iSub]->SetMaximum(2.5);
        //mgPhCosVsCent[iSub]->SetMinimum(-6);
        mgPhCosVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        leg_PhCos[iSub]->Draw("same");
        
    }

    
    TCanvas *c_PhCostog = new TCanvas("c_PhCostog", "", 0, 1024, 1200, 1024);
    c_PhCostog->cd();

    
    TMultiGraph *mgPhCosVsCent_Tog = new TMultiGraph();

    
    for(int iCent=0; iCent!=PhCosVsCent_Tog[0]->GetN(); iCent++){
        PhCosVsCent_Tog[1]->SetPointX(iCent, PhCosVsCent_Tog[0]->GetPointX(iCent) + 0.5);
    }

    PhCosVsCent_Tog[0]->SetMarkerStyle(34);
    PhCosVsCent_Tog[0]->SetMarkerColor(1);
    PhCosVsCent_Tog[0]->SetLineColor(1);

    PhCosVsCent_Tog[1]->SetMarkerStyle(75);
    PhCosVsCent_Tog[1]->SetMarkerColor(634);
    PhCosVsCent_Tog[1]->SetLineColor(634);

    mgPhCosVsCent_Tog->Add(PhCosVsCent_Tog[0]);
    mgPhCosVsCent_Tog->Add(PhCosVsCent_Tog[1]);

    mgPhCosVsCent_Tog->GetXaxis()->SetTitle("centrality %");
    mgPhCosVsCent_Tog->GetYaxis()->SetTitle(Form("<P_{h}cos(%i(#phi - #Psi_{%i}))> %%", nOrd, nOrd));
    //mgPhCosVsCent_Tog->SetMinimum(-6);
    //mgPhCosVsCent_Tog->SetMaximum(2);
    mgPhCosVsCent_Tog->Draw("AP");
    funcOfZero->Draw("same");

    TLegend *leg_PhCostog = new TLegend();
    leg_PhCostog->SetBorderSize(0);

    leg_PhCostog->AddEntry(PhCosVsCent_Tog[0], "#Lambda + #bar{#Lambda}, if #Psi_{e/w} used", "p");
    leg_PhCostog->AddEntry(PhCosVsCent_Tog[1], "#Lambda + #bar{#Lambda}, if #Psi_{comb} used", "p");
    leg_PhCostog->AddEntry((TObject*) 0, "Average in 20-50\% range centrality:", "");
    

    TGraphErrors *PhCos_20_50[2];

    for(int iSub = 0; iSub!=2; iSub++){
        PhCos_20_50[iSub] = MakeAverageBin(PhCosVsCent_Tog[iSub], 4, 6);
        double PhCos_aver = PhCos_20_50[iSub]->GetPointY(3);
        double PhCos_aver_err = PhCos_20_50[iSub]->GetErrorY(3);
        leg_PhCostog->AddEntry(PhCosVsCent_Tog[iSub], Form("<P_{x}> = %.2f #pm %.2f", PhCos_aver, PhCos_aver_err), "p");
    }    
    
    leg_PhCostog->Draw("same");

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
    grout->SetMarkerSize(1);
    
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

TGraphErrors* GetCosTetaSub(TProfile* prof, TProfile* prof_aver){
    TGraphErrors *gr_out = new TGraphErrors(prof->GetNbinsX());

    double y, ey, ey_1, ey_2;

    for(int iBin=0; iBin!=prof->GetNbinsX(); iBin++){
        y = prof->GetBinContent(iBin+1) - prof_aver->GetBinContent(iBin+1);
        ey_1 = prof->GetBinError(iBin+1);
        ey_2 = prof_aver->GetBinError(iBin+1);
        ey = pow(ey_1*ey_1 + ey_2*ey_2, 0.5);
        gr_out->SetPoint(iBin, prof->GetBinCenter(iBin+1), y);
        gr_out->SetPointError(iBin, 0., ey);
    }

    gr_out->SetMarkerColor(1);
    gr_out->SetMarkerStyle(20);
    //gr_out->SetMarkerSize(3);
    
    gr_out->GetXaxis()->SetTitle("#phi - #Psi_{2}");
    gr_out->GetYaxis()->SetTitle("<cos(#theta)>");
    gr_out->GetXaxis()->SetLabelSize(0.1);
    gr_out->GetYaxis()->SetLabelSize(0.1);
    
    return gr_out;
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