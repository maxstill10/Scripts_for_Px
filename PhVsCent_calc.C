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

void PhVsCent_calc(){
    TFile *input = new TFile("../Files/14p5GeV/outputPolar_Psi1_2_3Dep_14p5.root", "read");
    TFile *output;

    if(sinUsage) output = TFile::Open("../output/output_PhVsCent_SinInclude_14p5.root", "recreate");
    if(!sinUsage) output = TFile::Open("../output/output_PhVsCent_ConstSine_14p5.root", "recreate");

    TH1F *CosOfDiffFlatt_1 = (TH1F*)input->Get("CosOfDiff_1");
    TH1F *CosOfDiffFlatt_2 = (TH1F*)input->Get("CosOfDiff_2");

    TGraphErrors *Res_1VsCent = new TGraphErrors(9);


    double resFlatt, resErr;
    double cent[9] = {1.5, 6.5, 14, 24, 34, 44, 54, 64, 74};
    int cent_lim[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
    
    for(int iCen=0; iCen!=9;iCen++){
        
        resFlatt = pow(2*TMath::Abs(CosOfDiffFlatt_1->GetBinContent(iCen+1)), 0.5);
        resErr = (1/(resFlatt))*CosOfDiffFlatt_1->GetBinError(iCen+1);
        
        Res_1VsCent->SetPoint(8-iCen, cent[8-iCen], resFlatt);
        Res_1VsCent->SetPointError(8-iCen, 0., resErr);
        
    }

    TCanvas *c_res = new TCanvas("c_res", "c_res", 0, 1024, 1200, 1024);
    c_res->cd();

    Res_1VsCent->SetMarkerStyle(23);
    Res_1VsCent->SetMarkerColor(634);
    Res_1VsCent->SetMarkerSize(3);    

    Res_1VsCent->GetXaxis()->SetTitle("Centrality, %");
    Res_1VsCent->GetYaxis()->SetTitle("Res(#Psi)");
    Res_1VsCent->SetMinimum(0);
    Res_1VsCent->SetMaximum(0.8);

    TLegend *leg_res = new TLegend();

    leg_res->AddEntry(Res_1VsCent, "Full Res(Psi1)", "p");
    leg_res->SetBorderSize(0);

    Res_1VsCent->Draw("AP");
    leg_res->Draw("same");

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);


    //..................................Ph calculation starts here.................................
    TProfile *prSin_diffPhiPsi1[9][6];
    TProfile *prSin_diffPhiPsi1_LamBar[9][6];

    TH1F *InvMLamDist[9][6];
    TH1F *InvMLamBarDist[9][6];

    for(int iCent=0; iCent!=9; iCent++){
        for(int iphi=0; iphi!=6; iphi++){
            if(sinUsage){
                prSin_diffPhiPsi1[iCent][iphi] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, 0));
                prSin_diffPhiPsi1_LamBar[iCent][iphi] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, 0));
            }
            
            if(!sinUsage){
                prSin_diffPhiPsi1[iCent][iphi] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_%i_%i_%i", iCent, iphi, 0));
                prSin_diffPhiPsi1_LamBar[iCent][iphi] = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, 0));
            }

            //lam invM dist
            InvMLamDist[iCent][iphi] = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", iCent, iphi, 0));
            InvMLamBarDist[iCent][iphi] = (TH1F*)input->Get(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, 0));
        }        
    }


    //...................................Fitting by invM..................................
    TProfile *prSin_diffPhiPsi1_aver[9];
    TProfile *prSin_diffPhiPsi1_LamBar_aver[9];
    TH1F *InvMLamDist_aver[9];
    TH1F *InvMLamBarDist_aver[9];

    for(int iCent=0; iCent!=9; iCent++){
        prSin_diffPhiPsi1_aver[iCent] = (TProfile*)prSin_diffPhiPsi1[iCent][0]->Clone();
        prSin_diffPhiPsi1_aver[iCent]->Reset();
        prSin_diffPhiPsi1_LamBar_aver[iCent] = (TProfile*)prSin_diffPhiPsi1_LamBar[iCent][0]->Clone();
        prSin_diffPhiPsi1_LamBar_aver[iCent]->Reset();

        InvMLamDist_aver[iCent] = (TH1F*)InvMLamDist[iCent][0]->Clone();
        InvMLamDist_aver[iCent]->Reset();
        InvMLamBarDist_aver[iCent] = (TH1F*)InvMLamBarDist[iCent][0]->Clone();
        InvMLamBarDist_aver[iCent]->Reset();

        for(int iphi = 0; iphi!=6; iphi++){
            
            prSin_diffPhiPsi1_aver[iCent]->Add(prSin_diffPhiPsi1[iCent][iphi]);            
            prSin_diffPhiPsi1_LamBar_aver[iCent]->Add(prSin_diffPhiPsi1_LamBar[iCent][iphi]);
            InvMLamDist_aver[iCent]->Add(InvMLamDist[iCent][iphi]);
            InvMLamBarDist_aver[iCent]->Add(InvMLamBarDist[iCent][iphi]);
            
        }//for(int iphi = 0; iphi!=6; iphi++)
    }

    TGraphErrors *PhVsCent_Lam = new TGraphErrors(9);
    TGraphErrors *PhVsCent_LamBar = new TGraphErrors(9);
    
    
    double Ph, Ph_err;
    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.113, 1.12);
	FitFuncKFP->SetParLimits(2, 0.0005, 0.001);    
    FitFuncKFP->SetParLimits(4, 1.113, 1.12);
	FitFuncKFP->SetParLimits(5, 0.001, 0.005);

    TF1 *FitCosTetaFunc;

    TCanvas *canv_0 = new TCanvas("canv_0", "canv #Lambda", 0, 1024, 1800, 1024);
    TCanvas *canv_1 = new TCanvas("canv_1", "canv1 #Lambda", 0, 1024, 1800, 1024);
    TCanvas *canv_2 = new TCanvas("canv_2", "canv2 #Lambda", 0, 1024, 1800, 1024);
    TCanvas *canv_3 = new TCanvas("canv_3", "canv3 #Lambda", 0, 1024, 1800, 1024);

    canv_0->Divide(3,3);
    canv_1->Divide(3,3);
    canv_2->Divide(3,3);
    canv_3->Divide(3,3);

    
    for(int iCent = 0; iCent!=9; iCent++){    
        
        canv_0->cd(9-iCent);
        
        InvMLamDist_aver[iCent]->GetXaxis()->SetTitle("M_{inv}");
        InvMLamDist_aver[iCent]->GetXaxis()->SetLabelSize(0.06);
        InvMLamDist_aver[iCent]->GetXaxis()->SetTitleSize(0.06);
        InvMLamDist_aver[iCent]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
        InvMLamDist_aver[iCent]->Fit(FitFuncKFP, "R");        
        InvMLamDist_aver[iCent]->Draw();

        FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
        FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
        FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
        FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
        FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
        FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

        LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
        LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));

        
        

        FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);

        canv_1->cd(9-iCent);

        //FitCosTetaFunc->FixParameter(1, 0.);
        //FitCosTetaFunc->FixParameter(2, 0.);
        
        
        prSin_diffPhiPsi1_aver[iCent]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
        prSin_diffPhiPsi1_aver[iCent]->Fit(FitCosTetaFunc);
        prSin_diffPhiPsi1_aver[iCent]->Draw();

        Ph = FitCosTetaFunc->GetParameter(0);
        Ph_err = FitCosTetaFunc->GetParError(0);

        PhVsCent_Lam->SetPoint(8-iCent, cent[8-iCent], Ph);
        PhVsCent_Lam->SetPointError(8-iCent, 0, Ph_err);
        
        //AntiLambda
        canv_2->cd(9-iCent);
        
        InvMLamBarDist_aver[iCent]->GetXaxis()->SetTitle("M_{inv}");
        InvMLamBarDist_aver[iCent]->GetXaxis()->SetLabelSize(0.05);
        InvMLamBarDist_aver[iCent]->GetXaxis()->SetTitleSize(0.05);
        InvMLamBarDist_aver[iCent]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
        InvMLamBarDist_aver[iCent]->Fit(FitFuncKFP, "R");        
        InvMLamBarDist_aver[iCent]->Draw();

        FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
        FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
        FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
        FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
        FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
        FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

        LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
        LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));

        
        

        FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);

        canv_3->cd(9-iCent);

        //FitCosTetaFunc->FixParameter(1, 0.);
        //FitCosTetaFunc->FixParameter(2, 0.);
        
        prSin_diffPhiPsi1_LamBar_aver[iCent]->SetTitle(Form("#Lambda centrality, %i-%i", cent_lim[8-iCent], cent_lim[9-iCent]));
        prSin_diffPhiPsi1_LamBar_aver[iCent]->Fit(FitCosTetaFunc);
        prSin_diffPhiPsi1_LamBar_aver[iCent]->Draw();

        Ph = FitCosTetaFunc->GetParameter(0);
        Ph_err = FitCosTetaFunc->GetParError(0);

        PhVsCent_LamBar->SetPoint(8-iCent, cent[8-iCent]+0.5, Ph);
        PhVsCent_LamBar->SetPointError(8-iCent, 0, Ph_err);
    }//for(int iCent = 0; iCent!=9; iCent++)

    //...................................end of Fitting by invM..................................

    //Res factor correction
    PhVsCent_Lam = ResFactCorr(PhVsCent_Lam, Res_1VsCent);
    PhVsCent_LamBar = ResFactCorr(PhVsCent_LamBar, Res_1VsCent);

    if(sinUsage){
        PhVsCent_Lam = MultiplyLamFactor(PhVsCent_Lam, 3/alpha_Lam);
        PhVsCent_LamBar = MultiplyLamFactor(PhVsCent_LamBar, 3/alpha_Lam);
    }
    if(!sinUsage){
        PhVsCent_Lam = MultiplyLamFactor(PhVsCent_Lam, 8/alpha_Lam / TMath::Pi());
        PhVsCent_LamBar = MultiplyLamFactor(PhVsCent_LamBar, 8/alpha_Lam / TMath::Pi());
    }
    

    PhVsCent_Lam = MultiplyLamFactor(PhVsCent_Lam, 100);
    PhVsCent_LamBar = MultiplyLamFactor(PhVsCent_LamBar, 100);

    /*PhVsCent_Lam = MakeAverageBin(PhVsCent_Lam, 8, 9);
    PhVsCent_LamBar = MakeAverageBin(PhVsCent_LamBar, 8, 9);*/

    //Writing histograms to file
    
    //writing histograms
    PhVsCent_Lam->SetName("PhVsCent_Lam");
    PhVsCent_LamBar->SetName("PhVsCent_LamBar");

    PhVsCent_Lam->SetTitle("PhVsCent for #Lambda");
    PhVsCent_LamBar->SetTitle("PhVsCent for #bar{#Lambda}");
  
    PhVsCent_Lam->Write();
    PhVsCent_LamBar->Write();

    output->Close();

    //Drawing Ph Vs Centrality
    TCanvas *c_ph = new TCanvas("c_ph", "", 0, 1024, 1200, 1024);
    TMultiGraph *mgPhVsCent = new TMultiGraph();
    TLegend *leg_ph;
    
    c_ph->cd();
    PhVsCent_Lam->SetMarkerSize(3);
    PhVsCent_Lam->SetLineWidth(2);
    PhVsCent_Lam->SetMarkerStyle(29);
    PhVsCent_Lam->SetMarkerColor(634);
    PhVsCent_Lam->SetLineColor(634);

    PhVsCent_LamBar->SetMarkerSize(3);
    PhVsCent_LamBar->SetLineWidth(2);
    PhVsCent_LamBar->SetMarkerStyle(30);
    PhVsCent_LamBar->SetMarkerColor(602);
    PhVsCent_LamBar->SetLineColor(602);

    mgPhVsCent->Add(PhVsCent_Lam);
    mgPhVsCent->Add(PhVsCent_LamBar);

    mgPhVsCent->GetXaxis()->SetTitle("centrality \%");
    mgPhVsCent->GetYaxis()->SetTitle("<P_{H}> \%");
    mgPhVsCent->Draw("AP");
    funcOfZero->Draw("same");

    leg_ph = new TLegend();
    leg_ph->SetBorderSize(0);
    leg_ph->AddEntry(PhVsCent_Lam, "#Lambda", "p");
    leg_ph->AddEntry(PhVsCent_LamBar, "#bar{#Lambda}", "p");

    leg_ph->Draw("same");
    

    TGraphErrors *Ph_Lam_20_50;
    TGraphErrors *Ph_LamBar_20_50;

    Ph_Lam_20_50 = MakeAverageBin(PhVsCent_Lam, 4, 6);
    Ph_LamBar_20_50 = MakeAverageBin(PhVsCent_LamBar, 4, 6);

    cout<<"PhLam = "<<Ph_Lam_20_50->GetPointY(3)<<" #pm "<<Ph_Lam_20_50->GetErrorY(3)<<endl;
    cout<<"PhLamBar = "<<Ph_LamBar_20_50->GetPointY(3)<<" #pm "<<Ph_LamBar_20_50->GetErrorY(3)<<endl;

    //.............................Known dependences global polarization vs centrality were get..............................
    double egorikP_HLam[8] = {0.336783, 0.671942, 0.822737, 0.945594, 1.32896, 1.94692, 2.33029, 3.07855};
    double egorikP_HLam_Err[8] = {0., 0., 0., 0., 0., 0., 0.26, 0.52};

    double egorikP_HLamBar[8] = {-0.00236701, 0.775628, 0.848169, 1.04928, 1.4328, 1.63368, 2.74675, 6.96085};
    double egorikP_HLamBar_Err[8] = {0.62, 0.41, 0.24, 0.18, 0.26, 0.34, 0.75, 1.54};

    TGraphErrors *P_HVsCentEgorik = new TGraphErrors(8);
    TGraphErrors *P_HVsCentBarEgorik = new TGraphErrors(8); 

    for(int i = 0; i!=8; i++){
        P_HVsCentEgorik->SetPoint(i, cent[i]+0.5, egorikP_HLam[i]);
        P_HVsCentEgorik->SetPointError(i, 0., egorikP_HLam_Err[i]);

        P_HVsCentBarEgorik->SetPoint(i, cent[i]+0.5, egorikP_HLamBar[i]);
        P_HVsCentBarEgorik->SetPointError(i, 0., egorikP_HLamBar_Err[i]);
    }

    P_HVsCentEgorik->SetMarkerColor(6);
    P_HVsCentEgorik->SetMarkerStyle(29);
    P_HVsCentEgorik->SetMarkerSize(3);

    P_HVsCentBarEgorik->SetMarkerColor(6);
    P_HVsCentBarEgorik->SetMarkerStyle(29);
    P_HVsCentBarEgorik->SetMarkerSize(3);

    double chineseP_HLam[8] = {0.323741, 0.651335, 0.907589, 1.16439, 1.52711, 1.96062, 1.96964, 1.51822};
    double chineseP_HLam_Err[8] = {0., 0., 0., 0., 0., 0.33, 0.54};

    double chineseP_HLamBar[8] = {0.176211, 0.493392, 0.845815, 1.40969, 1.26872, 1.40969, 2.36123, 2.43172};
    double chineseP_HLamBar_Err[8] = {0., 0., 0., 0., 0., 0.54, 0.84, 1.66};

    TGraphErrors *P_HVsCentChinese = new TGraphErrors(8);
    TGraphErrors *P_HVsCentBarChinese = new TGraphErrors(8);

    P_HVsCentChinese->SetPoint(0, 5., chineseP_HLam[0]);
    P_HVsCentChinese->SetPointError(0, 0., 0.);

    P_HVsCentBarChinese->SetPoint(0, 5., chineseP_HLamBar[0]);
    P_HVsCentBarChinese->SetPointError(0, 0., 0.);

    for(int i = 1; i!=8; i++){
        P_HVsCentChinese->SetPoint(i, cent[i+1]+1, chineseP_HLam[i]);
        P_HVsCentChinese->SetPointError(i, 0., chineseP_HLam_Err[i]);

        P_HVsCentBarChinese->SetPoint(i, cent[i+1]+1, chineseP_HLamBar[i]);
        P_HVsCentBarChinese->SetPointError(i, 0., chineseP_HLamBar_Err[i]);
    }

    P_HVsCentChinese->SetMarkerColor(1);
    P_HVsCentChinese->SetMarkerStyle(34);
    P_HVsCentChinese->SetMarkerSize(3);

    P_HVsCentBarChinese->SetMarkerColor(1);
    P_HVsCentBarChinese->SetMarkerStyle(34);
    P_HVsCentBarChinese->SetMarkerSize(3);

    TCanvas *can = new TCanvas("can", "can", 0, 1024, 1200, 1024);
    can->cd();

    TMultiGraph *MgPolar_HVsCent = new TMultiGraph();
    
    MgPolar_HVsCent->Add(PhVsCent_Lam);
    MgPolar_HVsCent->Add(P_HVsCentEgorik);
    MgPolar_HVsCent->Add(P_HVsCentChinese);
    
    MgPolar_HVsCent->SetTitle("<P_{H}> of #Lambda Vs Centrality");
    MgPolar_HVsCent->GetXaxis()->SetTitle("Centrality, %");
    MgPolar_HVsCent->GetYaxis()->SetTitle("<P_{H}>, %");
    

    MgPolar_HVsCent->Draw("AP");
    funcOfZero->Draw("same");

    TLegend *legend = new TLegend();

    legend->AddEntry(PhVsCent_Lam, "This research", "p");
    legend->AddEntry(P_HVsCentEgorik, "Egor research", "p");
    legend->AddEntry(P_HVsCentChinese, "Chinese research", "p");
    legend->Draw("same");

    TCanvas *can1 = new TCanvas("can1", "can1", 0, 1024, 1200, 1024);
    can1->cd();

    TLegend *legendBar = new TLegend();
    TMultiGraph *MgPolar_HBarVsCent = new TMultiGraph();
    
    MgPolar_HBarVsCent->Add(PhVsCent_LamBar);
    MgPolar_HBarVsCent->Add(P_HVsCentBarEgorik);
    MgPolar_HBarVsCent->Add(P_HVsCentBarChinese);
    
    MgPolar_HBarVsCent->SetTitle("<P_{H}> of #bar{#Lambda} Vs Centrality");
    MgPolar_HBarVsCent->GetXaxis()->SetTitle("Centrality, %");
    MgPolar_HBarVsCent->GetYaxis()->SetTitle("<P_{H}>, %");
    

    MgPolar_HBarVsCent->Draw("AP");
    funcOfZero->Draw("same");

    legendBar->AddEntry(PhVsCent_LamBar, "This research", "p");
    legendBar->AddEntry(P_HVsCentBarEgorik, "Egor research", "p");
    legendBar->AddEntry(P_HVsCentBarChinese, "Chinese research", "p");
    legendBar->Draw("same");



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