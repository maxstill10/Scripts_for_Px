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

//true, if <cos(Psi_1-phi)sin(theta)> calculated
//false, if <cos(Psi_1-phi)> calculated

//const int nOrd = 2;
const int nOrd = 4;

TGraphErrors* ConvertProfToGraph(TProfile* prof);
TGraphErrors* MultiplyLamFactor(TGraphErrors *g, double coef);
TGraphErrors* GetCosTetaSub(TProfile* prof);
TGraphErrors* GetCosTetaSub(TProfile* prof, TGraphErrors* graph);
TGraphErrors* GetCosTetaSub(TProfile* prof, TProfile* prof_aver);
double GetCoef(double x, TF1 *Sig, TF1 *Bg);
double MyFitFunc(double *x, double *par);

TF1 *FuncGausKFP = new TF1("FuncGausKFP", "gaus(0)+gaus(3)", 1.1, 1.13);
TF1 *LineKFP = new TF1("LineKFP", "[0] + [1]*x", 1.1, 1.13);

void PzVsDiff_Psi_phi(){
    TFile *input = new TFile("../Files/14p5GeV/outputPolar_Psi1_2_3Dep_14p5.root", "read");


    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);


    //..............................................Start Pz vs Phi calculation................................
    TProfile *prCos_theta[9][6][nOrd];
    TProfile *prCos_theta_LamBar[9][6][nOrd];

    TH1F *InvMLamDist[9][6][nOrd];
    TH1F *InvMLamBarDist[9][6][nOrd];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iphi=0; iphi!=6; iphi++){
            for(int iSub=0; iSub!=nOrd; iSub++){      
                
                prCos_theta[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_theta_%i_%i_%i", iCent, iphi, iSub));
                prCos_theta_LamBar[iCent][iphi][iSub] = (TProfile*)input->Get(Form("prCos_theta_LamBar_%i_%i_%i", iCent, iphi, iSub));                

                //lam invM dist
                InvMLamDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub));
                InvMLamBarDist[iCent][iphi][iSub] = (TH1F*)input->Get(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub));
            }
        }        
    }

    //20-50% range cetrality merging
    TProfile *prCos_theta_20_50[6][nOrd];
    TProfile *prCos_theta_LamBar_20_50[6][nOrd];
    

    //lam invM dist
    TH1F *InvMLamDist_20_50[6][nOrd];
    TH1F *InvMLamBarDist_20_50[6][nOrd];

    for(int iSub=0; iSub!=nOrd; iSub++){
        for(int iphi=0; iphi!=6; iphi++){
            prCos_theta_20_50[iphi][iSub] = (TProfile*)prCos_theta[0][iphi][iSub]->Clone();
            prCos_theta_20_50[iphi][iSub]->Reset();
            prCos_theta_LamBar_20_50[iphi][iSub] = (TProfile*)prCos_theta_LamBar[0][iphi][iSub]->Clone();
            prCos_theta_LamBar_20_50[iphi][iSub]->Reset();

            //lam invM dist
            InvMLamDist_20_50[iphi][iSub] = (TH1F*)InvMLamDist[0][iphi][iSub]->Clone();
            InvMLamDist_20_50[iphi][iSub]->Reset();
            InvMLamBarDist_20_50[iphi][iSub] = (TH1F*)InvMLamBarDist[0][iphi][iSub]->Clone();
            InvMLamBarDist_20_50[iphi][iSub]->Reset();

            for(int iCent = 3; iCent!=6; iCent++){
                prCos_theta_20_50[iphi][iSub]->Add(prCos_theta[iCent][iphi][iSub]);
                prCos_theta_LamBar_20_50[iphi][iSub]->Add(prCos_theta_LamBar[iCent][iphi][iSub]);

                //lam invM dist
                InvMLamDist_20_50[iphi][iSub]->Add(InvMLamDist[iCent][iphi][iSub]);
                InvMLamBarDist_20_50[iphi][iSub]->Add(InvMLamBarDist[iCent][iphi][iSub]);
            }        
        }       
    }

    //...................................Fitting by invM..................................
    double Pz, Pz_err, phi;

    //Canvas init
    TCanvas *c_lam_20_50[nOrd];
    TCanvas *c_lamBar_20_50[nOrd];

    TCanvas *c_Cos_theta_20_50[nOrd];
    TCanvas *c_Cos_theta_LamBar_20_50[nOrd];

    for(int iSub=0; iSub!=nOrd; iSub++){
        c_lam_20_50[iSub] = new TCanvas(Form("c_lam_20_50_%i", iSub), Form("c_lam_20_50_%i", iSub), 0, 1024, 1800, 1024);
        c_lamBar_20_50[iSub] = new TCanvas(Form("c_lamBar_20_50_%i", iSub), Form("c_lamBar_20_50_%i", iSub), 0, 1024, 1800, 1024);

        c_Cos_theta_20_50[iSub] = new TCanvas(Form("c_Cos_theta_20_50_%i", iSub), Form("c_Cos_theta_20_50_%i", iSub), 0, 1024, 1600, 1024);
        c_Cos_theta_LamBar_20_50[iSub] = new TCanvas(Form("c_Cos_theta_LamBar_20_50_%i", iSub), Form("c_Cos_theta_LamBar_20_50_%i", iSub), 0, 1024, 1600, 1024);

        
        c_lam_20_50[iSub]->Divide(3, 2);
        c_lamBar_20_50[iSub]->Divide(3, 2);

        c_Cos_theta_20_50[iSub]->Divide(3, 2, 0, 0);
        c_Cos_theta_LamBar_20_50[iSub]->Divide(3, 2, 0, 0);
    }
    
    //Fit function init
    TF1 *FitCosTetaFunc;

    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.1, 1.13);
	FitFuncKFP->SetParLimits(2, 0.0005, 0.001);    
    FitFuncKFP->SetParLimits(4, 1.113, 1.12);
	FitFuncKFP->SetParLimits(5, 0.001, 0.005);

    //graphs init
    TGraphErrors *grCos_thetaVsdphi_20_50[nOrd];
    TGraphErrors *grCos_thetaVsdphi_LamBar_20_50[nOrd];

    for(int iSub=0; iSub!=nOrd; iSub++){
        grCos_thetaVsdphi_20_50[iSub] = new TGraphErrors(6);
        grCos_thetaVsdphi_LamBar_20_50[iSub] = new TGraphErrors(6);
    }

    TGraphErrors *grCos_theta_Sub_20_50[6][nOrd];
    TGraphErrors *grCos_theta_LamBar_Sub_20_50[6][nOrd];

    //This need to subtract distributions
    TProfile *prCos_theta_20_50_aver[nOrd];
    TProfile *prCos_theta_LamBar_20_50_aver[nOrd];

    for(int iSub=0; iSub!=nOrd; iSub++){
        prCos_theta_20_50_aver[iSub] = (TProfile*)prCos_theta_20_50[0][iSub]->Clone();
        prCos_theta_20_50_aver[iSub]->Reset();
        prCos_theta_LamBar_20_50_aver[iSub] = (TProfile*)prCos_theta_LamBar_20_50[0][iSub]->Clone();
        prCos_theta_LamBar_20_50_aver[iSub]->Reset();

        for(int iphi=0; iphi!=6; iphi++){
            prCos_theta_20_50_aver[iSub]->Add(prCos_theta_20_50[iphi][iSub]);
            prCos_theta_LamBar_20_50_aver[iSub]->Add(prCos_theta_LamBar_20_50[iphi][iSub]);
        }
    }

    
    for(int iSub=0; iSub!=nOrd; iSub++){
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

            //<Cos(phi-Psi1)>
            

            c_Cos_theta_20_50[iSub]->cd(iphi+1);
            if(iphi<3) gPad->SetBottomMargin(0);
            else gPad->SetTopMargin(0);
            if(iphi!=2 && iphi!=5) gPad->SetRightMargin(0);
            if(iphi!=0 && iphi!=3) gPad->SetLeftMargin(0);
            
            //subtracting
            grCos_theta_Sub_20_50[iphi][iSub] = GetCosTetaSub(prCos_theta_20_50[iphi][iSub], prCos_theta_20_50_aver[iSub]);

            if(iSub == 0) {grCos_theta_Sub_20_50[iphi][iSub]->SetMaximum(0.01); grCos_theta_Sub_20_50[iphi][iSub]->SetMinimum(-0.015);}
            if(iSub == 1) {grCos_theta_Sub_20_50[iphi][iSub]->SetMaximum(0.01); grCos_theta_Sub_20_50[iphi][iSub]->SetMinimum(-0.015);}
            if(iSub == 2) {grCos_theta_Sub_20_50[iphi][iSub]->SetMaximum(0.01); grCos_theta_Sub_20_50[iphi][iSub]->SetMinimum(-0.01);}
            if(iSub == 3) {grCos_theta_Sub_20_50[iphi][iSub]->SetMaximum(0.025); grCos_theta_Sub_20_50[iphi][iSub]->SetMinimum(-0.025);}

            if(iSub<2) grCos_theta_Sub_20_50[iphi][iSub]->SetTitle(Form("#Lambda, #Delta#phi = (%i#pi/6, %i#pi/6)", iphi, iphi+1));
            else if(iSub<4) grCos_theta_Sub_20_50[iphi][iSub]->SetTitle(Form("#Lambda, #Delta#phi = (%i#pi/9, %i#pi/9)", iphi, iphi+1));
            //grCos_theta_Sub_20_50[iphi][iSub] = ConvertProfToGraph(prCos_theta_20_50[iphi][iSub]);

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            grCos_theta_Sub_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            grCos_theta_Sub_20_50[iphi][iSub]->Draw("AP");

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            if( iSub < 2 ) phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            else if( iSub < 4 ) phi = TMath::Pi()/18 + iphi*TMath::Pi()/9.;
            grCos_thetaVsdphi_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_thetaVsdphi_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

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

            //<Cos(phi-Psi1)>
            c_Cos_theta_LamBar_20_50[iSub]->cd(iphi+1);
            if(iphi<3) gPad->SetBottomMargin(0);
            else gPad->SetTopMargin(0);
            if(iphi!=2 && iphi!=5) gPad->SetRightMargin(0);
            if(iphi!=0 && iphi!=3) gPad->SetLeftMargin(0);
            
            //Subtracting
            grCos_theta_LamBar_Sub_20_50[iphi][iSub] = GetCosTetaSub(prCos_theta_LamBar_20_50[iphi][iSub], prCos_theta_LamBar_20_50_aver[iSub]);

            if(iSub == 0) {grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMaximum(0.02); grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMinimum(-0.025);}
            if(iSub == 1) {grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMaximum(0.04); grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMinimum(-0.025);}
            if(iSub == 2) {grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMaximum(0.03); grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMinimum(-0.02);}
            if(iSub == 3) {grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMaximum(0.06); grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetMinimum(-0.07);}

            //grCos_theta_LamBar_Sub_20_50[iphi][iSub] = ConvertProfToGraph(prCos_theta_LamBar_20_50[iphi][iSub]);
            if(iSub<2) grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetTitle(Form("#bar{#Lambda}, #Delta#phi = (%i#pi/6, %i#pi/6)", iphi, iphi+1));
            else if(iSub<4) grCos_theta_LamBar_Sub_20_50[iphi][iSub]->SetTitle(Form("#bar{#Lambda}, #Delta#phi = (%i#pi/9, %i#pi/9)", iphi, iphi+1));

            FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
            
            grCos_theta_LamBar_Sub_20_50[iphi][iSub]->Fit(FitCosTetaFunc);
            grCos_theta_LamBar_Sub_20_50[iphi][iSub]->Draw("AP");

            Pz = FitCosTetaFunc->GetParameter(0);
            Pz_err = FitCosTetaFunc->GetParError(0);

            if( iSub < 2 ) phi = TMath::Pi()/12 + iphi*TMath::Pi()/6;
            else if( iSub < 4 ) phi = TMath::Pi()/18 + iphi*TMath::Pi()/9.;
            grCos_thetaVsdphi_LamBar_20_50[iSub]->SetPoint(iphi, phi, Pz);
            grCos_thetaVsdphi_LamBar_20_50[iSub]->SetPointError(iphi, 0., Pz_err);

        }//for(int iphi=0; iphi!=6; iphi++)

        grCos_thetaVsdphi_20_50[iSub]->SetMarkerSize(3);
        grCos_thetaVsdphi_20_50[iSub]->SetLineWidth(2);
        grCos_thetaVsdphi_20_50[iSub]->SetMarkerStyle(29);
        grCos_thetaVsdphi_20_50[iSub]->SetMarkerColor(634);
        grCos_thetaVsdphi_20_50[iSub]->SetLineColor(634);

        grCos_thetaVsdphi_LamBar_20_50[iSub] = MultiplyLamFactor(grCos_thetaVsdphi_LamBar_20_50[iSub], -1);
        grCos_thetaVsdphi_LamBar_20_50[iSub]->SetMarkerSize(3);
        grCos_thetaVsdphi_LamBar_20_50[iSub]->SetLineWidth(2);
        grCos_thetaVsdphi_LamBar_20_50[iSub]->SetMarkerStyle(30);
        grCos_thetaVsdphi_LamBar_20_50[iSub]->SetMarkerColor(602);
        grCos_thetaVsdphi_LamBar_20_50[iSub]->SetLineColor(602);

    }//for(int iSub=0; iSub!=2; iSub++)

    //...................................end of Fitting by invM..................................


    TGaxis::SetMaxDigits(1);

    TF1 *fitFuncOfCos = new TF1("fitFuncOfCos", "[0] + 2*[1]*sin(2*x) + 2*[2]*cos(2*x)", 0, TMath::Pi());

    fitFuncOfCos->FixParameter(0, 0);
    fitFuncOfCos->SetParLimits(1, -1, 1);
    //fitFuncOfCos->FixParameter(2, 0);

    TF1 *FuncOfFourier[2];
    TF1 *FuncOfFourier_Bar[2];
    TFitResultPtr r[2];
    TFitResultPtr r_bar[2];
    
    for(int iSub=0; iSub!=2; iSub++){
        grCos_thetaVsdphi_20_50[iSub]->Fit("fitFuncOfCos");
        FuncOfFourier[iSub] = new TF1(Form("FuncOfFourier_%i", iSub), "[0] + 2*[1]*sin(2*x) + 2*[2]*cos(2*x)", 0, TMath::Pi());
        FuncOfFourier[iSub]->SetParameter(0, fitFuncOfCos->GetParameter(0));
        FuncOfFourier[iSub]->SetParameter(1, fitFuncOfCos->GetParameter(1));
        FuncOfFourier[iSub]->SetParameter(2, fitFuncOfCos->GetParameter(2));
        FuncOfFourier[iSub]->SetLineColor(634);
        FuncOfFourier[iSub]->SetLineWidth(3);
        r[iSub] = grCos_thetaVsdphi_20_50[iSub]->Fit(fitFuncOfCos,"S");

        grCos_thetaVsdphi_LamBar_20_50[iSub]->Fit("fitFuncOfCos");
        FuncOfFourier_Bar[iSub] = new TF1(Form("FuncOfFourier_Bar_%i", iSub), "[0] + 2*[1]*sin(2*x) + 2*[2]*cos(2*x)", 0, TMath::Pi());
        FuncOfFourier_Bar[iSub]->SetParameter(0, fitFuncOfCos->GetParameter(0));
        FuncOfFourier_Bar[iSub]->SetParameter(1, fitFuncOfCos->GetParameter(1));
        FuncOfFourier_Bar[iSub]->SetParameter(2, fitFuncOfCos->GetParameter(2));
        FuncOfFourier_Bar[iSub]->SetLineColor(602);
        FuncOfFourier_Bar[iSub]->SetLineWidth(3);
        r_bar[iSub] = grCos_thetaVsdphi_LamBar_20_50[iSub]->Fit(fitFuncOfCos,"S");
    }
    
    TCanvas *can[2];
    
    for(int iSub=0; iSub!=2; iSub++){
        can[iSub] = new TCanvas(Form("can_%i", iSub), Form("can_%i", iSub), 0, 1024, 1200, 1000);
    }

    TCanvas *can_tog = new TCanvas("can_tog", "", 0, 1024, 1800, 900);
    can_tog->Divide(2, 1, 0.0, 0.0);
    

    
    TMultiGraph *mgCos_theta_20_50[2];


    TLegend *leg[2];

    for(int iSub=0; iSub!=2; iSub++){
        leg[iSub] = new TLegend();
        leg[iSub]->AddEntry(grCos_thetaVsdphi_20_50[iSub], "#Lambda", "p");
        leg[iSub]->AddEntry((TObject*)0, Form("p_{1} = %.3f #pm %.3f", r[iSub]->Value(1)*300/alpha_Lam/0.6 , r[iSub]->ParError(1)*300/alpha_Lam/0.6 ), "");
        leg[iSub]->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", r[iSub]->Chi2()/r[iSub]->Ndf()), "");
        leg[iSub]->AddEntry(grCos_thetaVsdphi_LamBar_20_50[iSub], "#bar{#Lambda}", "p");
        leg[iSub]->AddEntry((TObject*)0, Form("p_{1} = %.3f #pm %.3f", r_bar[iSub]->Value(1)*300/alpha_Lam/0.6, r_bar[iSub]->ParError(1)*300/alpha_Lam/0.6), "");
        leg[iSub]->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", r_bar[iSub]->Chi2()/r_bar[iSub]->Ndf()), "");
        leg[iSub]->SetBorderSize(0);
    }

    

    TText t;
    t.SetTextFont(42); 

    
    for(int iSub=0; iSub!=2; iSub++){
        //
        mgCos_theta_20_50[iSub] = new TMultiGraph();
        mgCos_theta_20_50[iSub]->Add(grCos_thetaVsdphi_20_50[iSub]);
        mgCos_theta_20_50[iSub]->Add(grCos_thetaVsdphi_LamBar_20_50[iSub]);


        mgCos_theta_20_50[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{2}");
        mgCos_theta_20_50[iSub]->GetYaxis()->SetTitle("<cos(#theta)>");

        mgCos_theta_20_50[iSub]->SetMinimum(-0.0015);
        mgCos_theta_20_50[iSub]->SetMaximum(0.0015);
        mgCos_theta_20_50[iSub]->GetXaxis()->SetLimits(0, 3.15);

        can[iSub]->cd();
        mgCos_theta_20_50[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        FuncOfFourier[iSub]->Draw("same");
        FuncOfFourier_Bar[iSub]->Draw("same");
        leg[iSub]->Draw("same");
        if(iSub == 1){
            t.DrawText(1.5, 0.0013, "Centrality 20 - 50%");
            t.DrawText(1.5, 0.0011, "Psi_comb was used");
        }
        if(iSub == 0){
            t.DrawText(1.5, 0.0013, "Centrality 20 - 50%");
            t.DrawText(1.5, 0.0011, "Psi_e/w was used");
        }

        can_tog->cd(iSub+1);
        mgCos_theta_20_50[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        FuncOfFourier[iSub]->Draw("same");
        FuncOfFourier_Bar[iSub]->Draw("same");
        leg[iSub]->Draw("same");
        if(iSub == 1){
            gPad->SetLeftMargin(0);
            gPad->SetRightMargin(0.2);
            t.DrawText(1.5, 0.0013, "Centrality 20 - 50%");
            t.DrawText(1.5, 0.0011, "Psi_comb was used");
        }
        if(iSub == 0){
            gPad->SetRightMargin(0);
            gPad->SetLeftMargin(0.2);
            t.DrawText(1.5, 0.0013, "Centrality 20 - 50%");
            t.DrawText(1.5, 0.0011, "Psi_e/w was used");
        }
    }
    
    //...............................................For Psi3 calculations....................................

    if(nOrd <= 2) return;


    TF1 *fitFuncOfCos_Psi3 = new TF1("fitFuncOfCos_Psi3", "[0] + 2*[1]*sin(3*x) + 2*[2]*cos(3*x)", 0, 2*TMath::Pi()/3);

    fitFuncOfCos_Psi3->FixParameter(0, 0);
    fitFuncOfCos_Psi3->SetParLimits(1, -1, 1);
    fitFuncOfCos_Psi3->FixParameter(2, 0);

    TF1 *FuncOfFourier_Psi3[1];
    TF1 *FuncOfFourier_Bar_Psi3[1];
    TFitResultPtr r_Psi3[1];
    TFitResultPtr r_bar_Psi3[1];
    
    for(int iSub=0; iSub!=1; iSub++){
        grCos_thetaVsdphi_20_50[2+iSub]->Fit("fitFuncOfCos_Psi3");
        FuncOfFourier_Psi3[iSub] = new TF1(Form("FuncOfFourier_%i", iSub), "[0] + 2*[1]*sin(3*x) + 2*[2]*cos(3*x)", 0, 2*TMath::Pi()/3.);
        FuncOfFourier_Psi3[iSub]->SetParameter(0, fitFuncOfCos_Psi3->GetParameter(0));
        FuncOfFourier_Psi3[iSub]->SetParameter(1, fitFuncOfCos_Psi3->GetParameter(1));
        FuncOfFourier_Psi3[iSub]->SetParameter(2, fitFuncOfCos_Psi3->GetParameter(2));
        FuncOfFourier_Psi3[iSub]->SetLineColor(634);
        FuncOfFourier_Psi3[iSub]->SetLineWidth(3);
        r_Psi3[iSub] = grCos_thetaVsdphi_20_50[2+iSub]->Fit(fitFuncOfCos_Psi3,"S");

        grCos_thetaVsdphi_LamBar_20_50[2+iSub]->Fit("fitFuncOfCos_Psi3");
        FuncOfFourier_Bar_Psi3[iSub] = new TF1(Form("FuncOfFourier_Bar_%i", iSub), "[0] + 2*[1]*sin(3*x) + 2*[2]*cos(3*x)", 0, 2*TMath::Pi()/3.);
        FuncOfFourier_Bar_Psi3[iSub]->SetParameter(0, fitFuncOfCos_Psi3->GetParameter(0));
        FuncOfFourier_Bar_Psi3[iSub]->SetParameter(1, fitFuncOfCos_Psi3->GetParameter(1));
        FuncOfFourier_Bar_Psi3[iSub]->SetParameter(2, fitFuncOfCos_Psi3->GetParameter(2));
        FuncOfFourier_Bar_Psi3[iSub]->SetLineColor(602);
        FuncOfFourier_Bar_Psi3[iSub]->SetLineWidth(3);
        r_bar_Psi3[iSub] = grCos_thetaVsdphi_LamBar_20_50[2+iSub]->Fit(fitFuncOfCos_Psi3,"S");
    }
    
    TCanvas *can_Psi3[2];
    
    for(int iSub=0; iSub!=2; iSub++){
        can_Psi3[iSub] = new TCanvas(Form("can_Psi3_%i", iSub), Form("can_Psi3_%i", iSub), 0, 1024, 1200, 1000);
    }

    TCanvas *can_tog_Psi3 = new TCanvas("can_tog_Psi3", "", 0, 1024, 1800, 900);
    can_tog_Psi3->Divide(2, 1, 0.0, 0.0);
    

    
    TMultiGraph *mgCos_theta_20_50_Psi3[2];


    TLegend *leg_Psi3[2];

    for(int iSub=0; iSub!=2; iSub++){
        leg_Psi3[iSub] = new TLegend();
        leg_Psi3[iSub]->AddEntry(grCos_thetaVsdphi_20_50[2+iSub], "#Lambda", "p");
        leg_Psi3[iSub]->AddEntry((TObject*)0, Form("p_{1} = %.2f #pm %.2f", r_Psi3[iSub]->Value(1)*300/alpha_Lam/0.6 , r_Psi3[iSub]->ParError(1)*300/alpha_Lam/0.6 ), "");
        leg_Psi3[iSub]->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", r_Psi3[iSub]->Chi2()/r_Psi3[iSub]->Ndf()), "");
        leg_Psi3[iSub]->AddEntry(grCos_thetaVsdphi_LamBar_20_50[2+iSub], "#bar{#Lambda}", "p");
        leg_Psi3[iSub]->AddEntry((TObject*)0, Form("p_{1} = %.2f #pm %.2f", r_bar_Psi3[iSub]->Value(1)*300/alpha_Lam/0.6, r_bar_Psi3[iSub]->ParError(1)*300/alpha_Lam/0.6), "");
        leg_Psi3[iSub]->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", r_bar_Psi3[iSub]->Chi2()/r_bar_Psi3[iSub]->Ndf()), "");
        leg_Psi3[iSub]->SetBorderSize(0);
    }

    
    
    for(int iSub=0; iSub!=2; iSub++){
        //
        mgCos_theta_20_50_Psi3[iSub] = new TMultiGraph();
        mgCos_theta_20_50_Psi3[iSub]->Add(grCos_thetaVsdphi_20_50[2+iSub]);
        mgCos_theta_20_50_Psi3[iSub]->Add(grCos_thetaVsdphi_LamBar_20_50[2+iSub]);


        mgCos_theta_20_50_Psi3[iSub]->GetXaxis()->SetTitle("#phi - #Psi_{3}");
        mgCos_theta_20_50_Psi3[iSub]->GetYaxis()->SetTitle("<cos(#theta)>");

        //mgCos_theta_20_50_Psi3[iSub]->SetMinimum(-0.003);
        //mgCos_theta_20_50_Psi3[iSub]->SetMaximum(0.003);
        mgCos_theta_20_50_Psi3[iSub]->GetXaxis()->SetLimits(0, 2.1);

        can_Psi3[iSub]->cd();
        mgCos_theta_20_50_Psi3[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        if(iSub==0) FuncOfFourier_Psi3[iSub]->Draw("same");
        if(iSub==0) FuncOfFourier_Bar_Psi3[iSub]->Draw("same");
        leg_Psi3[iSub]->Draw("same");
        if(iSub == 1){
            t.DrawText(1, 0.0025, "Centrality 20 - 50%");
            t.DrawText(1, 0.002, "Psi_comb was used");
        }
        if(iSub == 0){
            t.DrawText(0.9, 0.0013, "Centrality 20 - 50%");
            t.DrawText(0.9, 0.0011, "Psi_e/w was used");
        }

        can_tog_Psi3->cd(iSub+1);
        mgCos_theta_20_50_Psi3[iSub]->Draw("AP");
        funcOfZero->Draw("same");
        if(iSub==0) FuncOfFourier_Psi3[iSub]->Draw("same");
        if(iSub==0) FuncOfFourier_Bar_Psi3[iSub]->Draw("same");
        leg_Psi3[iSub]->Draw("same");
        if(iSub == 1){
            gPad->SetLeftMargin(0);
            gPad->SetRightMargin(0.2);
            t.DrawText(1, 0.0025, "Centrality 20 - 50%");
            t.DrawText(1, 0.002, "Psi_comb was used");
        }
        if(iSub == 0){
            gPad->SetRightMargin(0);
            gPad->SetLeftMargin(0.2);
            t.DrawText(1, 0.0025, "Centrality 20 - 50%");
            t.DrawText(1, 0.002, "Psi_e/w was used");
        }
    }

    mgCos_theta_20_50_Psi3[0]->SetMinimum(-0.0015);
    mgCos_theta_20_50_Psi3[0]->SetMaximum(0.0015);

    TCanvas *can_psi2_psi3 = new TCanvas("can_psi2_psi3", "", 0, 1024, 1800, 900);
    can_psi2_psi3->Divide(2, 1, 0.0, 0.0);
    
    //
    can_psi2_psi3->cd(1);

    mgCos_theta_20_50[0]->Draw("AP");

    funcOfZero->Draw("same");
    FuncOfFourier[0]->Draw("same");
    FuncOfFourier_Bar[0]->Draw("same");
    leg[0]->Draw("same");

    gPad->SetRightMargin(0);
    gPad->SetLeftMargin(0.2);
    t.DrawText(1.5, 0.0013, "Centrality 20 - 50%");
    t.DrawText(1.5, 0.0011, "Psi_ew was used");

    //
    can_psi2_psi3->cd(2);

    mgCos_theta_20_50_Psi3[0]->Draw("AP");
    funcOfZero->Draw("same");
    FuncOfFourier_Psi3[0]->Draw("same");
    FuncOfFourier_Bar_Psi3[0]->Draw("same");
    leg_Psi3[0]->Draw("same");

    gPad->SetLeftMargin(0);
    gPad->SetRightMargin(0.2);
    t.DrawText(0.9, 0.0013, "Centrality 20 - 50%");
    t.DrawText(0.9, 0.0011, "Psi_ew was used");

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
    //graph_cos_sub->GetXaxis()->SetLabelSize(0.1);
    //graph_cos_sub->GetYaxis()->SetLabelSize(0.1);
    
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
    gr_out->GetXaxis()->SetLabelSize(0.05);
    gr_out->GetYaxis()->SetLabelSize(0.05);
    
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