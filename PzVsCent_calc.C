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


//on/off fit drawing
const bool mfitDr = false;

//order of EP
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

void PzVsCent_calc(){
    TFile *input = new TFile("../Files/14p5GeV/outputPolar_Psi1_2_3Dep_14p5.root", "read");
    TFile *output;

    output = TFile::Open(Form("../output/output_PzVsCent_forPsi%i_14p5.root", nOrd), "recreate");

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


    //..................................Pz calculation starts here.................................
    TProfile *prCos_theta[9][6][2];
    TProfile *prCos_theta_LamBar[9][6][2];

    TH1F *InvMLamDist[9][6][2];
    TH1F *InvMLamBarDist[9][6][2];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iphi=0; iphi!=6; iphi++){
            for(int iSub=0; iSub!=2; iSub++){
                
                prCos_theta[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_theta_%i_%i_%i", iCent, iphi, iSub + (nOrd-2)*2));
                prCos_theta_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_theta_LamBar_%i_%i_%i", iCent, iphi, iSub + (nOrd-2)*2));
                

                //lam invM dist
                InvMLamDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub));
                InvMLamBarDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub));
            }
        }        
    }

    //...................................Fitting by invM..................................
    TGraphErrors *grCos_thetaVsdphi[9][2];
    TGraphErrors *grCos_theta_LamBarVsdphi[9][2];

    for(int iCent = 0; iCent!=9; iCent++){
        for(int iSub = 0; iSub!=2; iSub++){           
            grCos_thetaVsdphi[iCent][iSub] = new TGraphErrors(6);
            grCos_theta_LamBarVsdphi[iCent][iSub] = new TGraphErrors(6);
        }
    }

    TGraphErrors *grCos_theta_Sub[9][6][2];
    TGraphErrors *grCos_theta_LamBar_Sub[9][6][2];

    TProfile *prCos_theta_aver[9][2];
    TProfile *prCos_theta_LamBar_aver[9][2];

    for(int iSub=0; iSub!=2; iSub++){
        for(int iCent = 0; iCent!=9; iCent++){
            prCos_theta_aver[iCent][iSub] = (TProfile*)prCos_theta[iCent][0][iSub]->Clone();
            prCos_theta_aver[iCent][iSub]->Reset();
            prCos_theta_LamBar_aver[iCent][iSub] = (TProfile*)prCos_theta_LamBar[iCent][0][iSub]->Clone();
            prCos_theta_LamBar_aver[iCent][iSub]->Reset();
            for(int iphi = 0; iphi!=6; iphi++){
                
                prCos_theta_aver[iCent][iSub]->Add(prCos_theta[iCent][iphi][iSub]);            
                prCos_theta_LamBar_aver[iCent][iSub]->Add(prCos_theta_LamBar[iCent][iphi][iSub]);
            }
        }
    }



    

    
    double Pz, Pz_err;
    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.114, 1.117);
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

                //FitCosTetaFunc->FixParameter(1, 0.);
                //FitCosTetaFunc->FixParameter(2, 0.);

                grCos_theta_Sub[iCent][iphi][iSub] = GetCosTetaSub(prCos_theta[iCent][iphi][iSub], prCos_theta_aver[iCent][iSub]);
                
                grCos_theta_Sub[iCent][iphi][iSub]->Fit(FitCosTetaFunc);
                grCos_theta_Sub[iCent][iphi][iSub]->Draw("AP");

                Pz = FitCosTetaFunc->GetParameter(0);
                Pz_err = FitCosTetaFunc->GetParError(0);

                double phi = 2*TMath::Pi() / nOrd /12 + iphi*2*TMath::Pi() / nOrd /6;
                grCos_thetaVsdphi[iCent][iSub]->SetPoint(iphi, phi, Pz);
                grCos_thetaVsdphi[iCent][iSub]->SetPointError(iphi, 0., Pz_err);
                
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

                //Subtracting
                grCos_theta_LamBar_Sub[iCent][iphi][iSub] = GetCosTetaSub(prCos_theta_LamBar[iCent][iphi][iSub], prCos_theta_LamBar_aver[iCent][iSub]);
                
                grCos_theta_LamBar_Sub[iCent][iphi][iSub]->Fit(FitCosTetaFunc);
                grCos_theta_LamBar_Sub[iCent][iphi][iSub]->Draw("AP");

                Pz = FitCosTetaFunc->GetParameter(0);
                Pz_err = FitCosTetaFunc->GetParError(0);

                phi = 2*TMath::Pi() / nOrd /12 + iphi*2*TMath::Pi() / nOrd /6;
                grCos_theta_LamBarVsdphi[iCent][iSub]->SetPoint(iphi, phi, Pz);
                grCos_theta_LamBarVsdphi[iCent][iSub]->SetPointError(iphi, 0., Pz_err);
            }//for(int iphi = 0; iphi!=6; iphi++)
            
            grCos_thetaVsdphi[iCent][iSub]->SetMarkerSize(3);
            grCos_thetaVsdphi[iCent][iSub]->SetLineWidth(2);
            grCos_thetaVsdphi[iCent][iSub]->SetMarkerStyle(20);
            grCos_thetaVsdphi[iCent][iSub]->SetMarkerColor(1);

            grCos_theta_LamBarVsdphi[iCent][iSub]->SetMarkerSize(3);
            grCos_theta_LamBarVsdphi[iCent][iSub]->SetLineWidth(2);
            grCos_theta_LamBarVsdphi[iCent][iSub]->SetMarkerStyle(20);
            grCos_theta_LamBarVsdphi[iCent][iSub]->SetMarkerColor(1);
        }//for(int iCent = 0; iCent!=9; iCent++)
    }//for(int iSub=0; iSub!=2; iSub++)

    //...................................end of Fitting by invM..................................

    if(!mfitDr){
        gROOT->SetBatch(kFALSE);
    }

    //...................................Fitting by phi............................................

    TGraphErrors *PzVsCent_Lam[2];
    TGraphErrors *PzVsCent_LamBar[2];

    for(int iSub=0; iSub!=2; iSub++){
        PzVsCent_Lam[iSub] = new TGraphErrors(9);
        PzVsCent_LamBar[iSub] = new TGraphErrors(9);
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

    TGaxis::SetMaxDigits(1);
    //gStyle->SetOptFit(0001);


    int cent_lim[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

    TF1 *fitFuncOfCos = new TF1("fitFuncOfCos", Form("[0] + 2*[1]*sin(%i*x) + 2*[2]*cos(%i*x)", nOrd, nOrd), 0, 2*TMath::Pi()/nOrd);

    //fitFuncOfCos->FixParameter(0, 0);
    fitFuncOfCos->SetParLimits(1, -1, 1);
    //fitFuncOfCos->FixParameter(2, 0);

    
    for(int iCent=0; iCent!=9; iCent++){
        for(int iSub=0; iSub!=2; iSub++){

            //Set Y titles
            grCos_thetaVsdphi[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#Lambda}-#Psi_{1})>");
            grCos_theta_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetTitle("<cos(#phi_{#bar{#Lambda}}-#Psi_{1})>");

            if(iSub == 0){grCos_thetaVsdphi[iCent][iSub]->SetMinimum(-0.002); grCos_thetaVsdphi[iCent][iSub]->SetMaximum(0.003);}
            if(iSub == 1){grCos_thetaVsdphi[iCent][iSub]->SetMinimum(-0.006); grCos_thetaVsdphi[iCent][iSub]->SetMaximum(0.006);}
            if(iSub == 0){grCos_theta_LamBarVsdphi[iCent][iSub]->SetMinimum(-0.003); grCos_theta_LamBarVsdphi[iCent][iSub]->SetMaximum(0.005);}
            if(iSub == 1){grCos_theta_LamBarVsdphi[iCent][iSub]->SetMinimum(-0.008); grCos_theta_LamBarVsdphi[iCent][iSub]->SetMaximum(0.008);}

            //Set graph titles
            grCos_thetaVsdphi[iCent][iSub]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            grCos_theta_LamBarVsdphi[iCent][iSub]->SetTitle(Form("#bar{#Lambda} centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
            

            //Drawing           

            c[iSub]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_thetaVsdphi[iCent][iSub]->GetYaxis()->SetTitleOffset(1);
                grCos_thetaVsdphi[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
                grCos_thetaVsdphi[iCent][iSub]->GetYaxis()->SetLabelSize(0.1);
            }
            grCos_thetaVsdphi[iCent][iSub]->Fit("fitFuncOfCos");
            grCos_thetaVsdphi[iCent][iSub]->Draw("AP");
            
            Pz = fitFuncOfCos->GetParameter(1);
            Pz_err = fitFuncOfCos->GetParError(1);
            PzVsCent_Lam[iSub]->SetPoint(8-iCent, cent[8-iCent], Pz);
            PzVsCent_Lam[iSub]->SetPointError(8-iCent, 0, Pz_err);


            c[iSub+2]->cd(9-iCent);
            if(iCent == 2 || iCent == 5 || iCent == 8){
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.02);
                grCos_theta_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetTitleOffset(1);
                grCos_theta_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetTitleSize(0.1);
                grCos_theta_LamBarVsdphi[iCent][iSub]->GetYaxis()->SetLabelSize(0.1);
            }
            grCos_theta_LamBarVsdphi[iCent][iSub]->Draw("AP");
            grCos_theta_LamBarVsdphi[iCent][iSub]->Fit("fitFuncOfCos");

            Pz = fitFuncOfCos->GetParameter(1);
            Pz_err = fitFuncOfCos->GetParError(1);
            PzVsCent_LamBar[iSub]->SetPoint(8-iCent, cent[8-iCent]+0.5, Pz);
            PzVsCent_LamBar[iSub]->SetPointError(8-iCent, 0, Pz_err);
        }
    }
    //...................................end of Fitting by phi........................................

    //Res factor correction
    for(int iSub=0; iSub!=2; iSub++){

        if(nOrd == 2)
        {
            PzVsCent_Lam[iSub] = ResFactCorr(PzVsCent_Lam[iSub], Res_2VsCent[iSub]);
            PzVsCent_LamBar[iSub] = ResFactCorr(PzVsCent_LamBar[iSub], Res_2VsCent[iSub]);
        }

        if(nOrd == 3)
        {
            PzVsCent_Lam[iSub] = ResFactCorr(PzVsCent_Lam[iSub], Res_3VsCent[iSub]);
            PzVsCent_LamBar[iSub] = ResFactCorr(PzVsCent_LamBar[iSub], Res_3VsCent[iSub]);
        }

        PzVsCent_Lam[iSub] = MultiplyLamFactor(PzVsCent_Lam[iSub], 3/alpha_Lam);
        PzVsCent_LamBar[iSub] = MultiplyLamFactor(PzVsCent_LamBar[iSub], 3/alpha_LamBar);        

        PzVsCent_Lam[iSub] = MultiplyLamFactor(PzVsCent_Lam[iSub], 100);
        PzVsCent_LamBar[iSub] = MultiplyLamFactor(PzVsCent_LamBar[iSub], 100);

        PzVsCent_Lam[iSub] = MakeAverageBin(PzVsCent_Lam[iSub], 8, 9);
        PzVsCent_LamBar[iSub] = MakeAverageBin(PzVsCent_LamBar[iSub], 8, 9);

    }

    //Writing histograms to file
    TGraphErrors *PzVsCent_Tog[2];
    
    for(int iSub=0; iSub!=2; iSub++){

        //Lam+LamBar
        PzVsCent_Tog[iSub] = MakeAverage(PzVsCent_Lam[iSub], PzVsCent_LamBar[iSub]);

        //writing histograms
        PzVsCent_Lam[iSub]->SetName(Form("PzVsCent_Lam_%i", iSub));
        PzVsCent_LamBar[iSub]->SetName(Form("PzVsCent_LamBar_%i", iSub));
        PzVsCent_Tog[iSub]->SetName(Form("PzVsCent_Tog_%i", iSub));
        if(iSub == 0){
            PzVsCent_Lam[iSub]->SetTitle(Form("PzVsCent for #Lambda, if #Psi%i_{e/w} used", nOrd));
            PzVsCent_LamBar[iSub]->SetTitle(Form("PzVsCent for #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd));
            PzVsCent_Tog[iSub]->SetTitle(Form("PzVsCent for #Lambda + #bar{#Lambda}, if #Psi%i_{e/w} used", nOrd));
        }
        if(iSub == 1){
            PzVsCent_Lam[iSub]->SetTitle(Form("PzVsCent for #Lambda, if #Psi%i_{comb} used", nOrd));
            PzVsCent_LamBar[iSub]->SetTitle(Form("PzVsCent for #bar{#Lambda}, if #Psi%i_{comb} used", nOrd));
            PzVsCent_Tog[iSub]->SetTitle(Form("PzVsCent for #Lambda + #bar{#Lambda}, if #Psi%i_{comb} used", nOrd));
        }
        PzVsCent_Lam[iSub]->Write();
        PzVsCent_LamBar[iSub]->Write();
        PzVsCent_Tog[iSub]->Write();
    }
    
    output->Close();  

    
    //Drawing Pz Vs Centrality
    TCanvas *c_Pz[2];
    TCanvas *c_Pz_full = new TCanvas("c_Pz_full", "", 0, 1024, 1800, 900);
    c_Pz_full->Divide(2, 1, 0.0, 0.0);
    TMultiGraph *mgPzVsCent[2];
    TLegend *leg_Pz[2];
    for(int iSub = 0; iSub!=2; iSub++){
        c_Pz[iSub] = new TCanvas(Form("c_Pz_%i", iSub), "", 0, 1024, 1200, 1024);

        c_Pz[iSub]->cd();
        PzVsCent_Lam[iSub]->SetMarkerSize(3);
        PzVsCent_Lam[iSub]->SetLineWidth(2);
        PzVsCent_Lam[iSub]->SetMarkerStyle(29);
        PzVsCent_Lam[iSub]->SetMarkerColor(634);
        PzVsCent_Lam[iSub]->SetLineColor(634);

        PzVsCent_LamBar[iSub]->SetMarkerSize(3);
        PzVsCent_LamBar[iSub]->SetLineWidth(2);
        PzVsCent_LamBar[iSub]->SetMarkerStyle(30);
        PzVsCent_LamBar[iSub]->SetMarkerColor(602);
        PzVsCent_LamBar[iSub]->SetLineColor(602);

        mgPzVsCent[iSub] = new TMultiGraph();
        mgPzVsCent[iSub]->Add(PzVsCent_Lam[iSub]);
        mgPzVsCent[iSub]->Add(PzVsCent_LamBar[iSub]);

        mgPzVsCent[iSub]->GetXaxis()->SetTitle("centrality %");
        mgPzVsCent[iSub]->GetYaxis()->SetTitle(Form("<P_{z}sin(%i(#phi - #Psi_{%i}))> %%", nOrd, nOrd));
        mgPzVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");

        leg_Pz[iSub] = new TLegend();
        leg_Pz[iSub]->SetBorderSize(0);
        leg_Pz[iSub]->AddEntry(PzVsCent_Lam[iSub], "#Lambda", "p");
        leg_Pz[iSub]->AddEntry(PzVsCent_LamBar[iSub], "#bar{#Lambda}", "p");
        if(iSub==0) leg_Pz[iSub]->AddEntry((TObject*) 0, "#Psi_{e/w} was used", "");
        if(iSub==1) leg_Pz[iSub]->AddEntry((TObject*) 0, "#Psi_{comb} was used", "");

        leg_Pz[iSub]->Draw("same");

        c_Pz_full->cd(iSub+1);
        if(iSub==0) {gPad->SetRightMargin(0); gPad->SetLeftMargin(0.2); }
        if(iSub==1) {gPad->SetLeftMargin(0); gPad->SetRightMargin(0.2);}
        //mgPzVsCent[iSub]->SetMaximum(1);
        //mgPzVsCent[iSub]->SetMinimum(-1.5);
        mgPzVsCent[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        leg_Pz[iSub]->Draw("same");
        
    }

    
    TCanvas *c_Pztog = new TCanvas("c_Pztog", "", 0, 1024, 1200, 1024);
    c_Pztog->cd();

    
    TMultiGraph *mgPzVsCent_Tog = new TMultiGraph();

    
    for(int iCent=0; iCent!=PzVsCent_Tog[0]->GetN(); iCent++){
        PzVsCent_Tog[1]->SetPointX(iCent, cent[iCent] + 0.5);
    }

    PzVsCent_Tog[0]->SetMarkerStyle(34);
    PzVsCent_Tog[0]->SetMarkerColor(1);
    PzVsCent_Tog[0]->SetLineColor(1);

    PzVsCent_Tog[1]->SetMarkerStyle(75);
    PzVsCent_Tog[1]->SetMarkerColor(634);
    PzVsCent_Tog[1]->SetLineColor(634);

    mgPzVsCent_Tog->Add(PzVsCent_Tog[0]);
    mgPzVsCent_Tog->Add(PzVsCent_Tog[1]);

    mgPzVsCent_Tog->GetXaxis()->SetTitle("centrality %");
    mgPzVsCent_Tog->GetYaxis()->SetTitle(Form("<P_{z}sin(%i(#phi - #Psi_{%i}))> %%", nOrd, nOrd));
    //mgPzVsCent_Tog->SetMinimum(-1);
    //mgPzVsCent_Tog->SetMaximum(1);
    mgPzVsCent_Tog->Draw("AP");
    funcOfZero->Draw("same");

    TLegend *leg_Pztog = new TLegend();
    leg_Pztog->SetBorderSize(0);

    leg_Pztog->AddEntry(PzVsCent_Tog[0], "#Lambda + #bar{#Lambda}, if #Psi_{e/w} used", "p");
    leg_Pztog->AddEntry(PzVsCent_Tog[1], "#Lambda + #bar{#Lambda}, if #Psi_{comb} used", "p");
    leg_Pztog->AddEntry((TObject*) 0, "Average in 20-50\% range centrality:", "");
    

    TGraphErrors *Pz_20_50[2];

    for(int iSub = 0; iSub!=2; iSub++){
        Pz_20_50[iSub] = MakeAverageBin(PzVsCent_Tog[iSub], 4, 6);
        double Pz_aver = Pz_20_50[iSub]->GetPointY(3);
        double Pz_aver_err = Pz_20_50[iSub]->GetErrorY(3);
        leg_Pztog->AddEntry(PzVsCent_Tog[iSub], Form("<P_{x}> = %.2f #pm %.2f", Pz_aver, Pz_aver_err), "p");
    }

    
    
    leg_Pztog->Draw("same");

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