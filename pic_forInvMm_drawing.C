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

void pic_forInvMm_drawing(){
   TFile *input = new TFile("../Files/14p5GeV/outputPolar_Psi1_2_3Dep_14p5.root", "read");

   gStyle->SetOptStat(0);

    //Function of zero
    TF1 *funcOfZero = new TF1("funcOfZero", "[0]", -100, 300);
    funcOfZero->SetParameter(0, 0);
    funcOfZero->SetLineStyle(9);
    funcOfZero->SetLineColor(1);
    funcOfZero->SetLineWidth(2);


    //..............................................Start Ph vs Phi calculation................................
    TProfile *prSin_diffPhiPsi1;

    TH1F *InvMLamDist;

    prSin_diffPhiPsi1 = (TProfile*)input->Get(Form("prSin_diffPhiPsi1_%i_%i_%i", 3, 3, 1));                

    //lam invM dist
    InvMLamDist = (TH1F*)input->Get(Form("InvMLamDist_%i_%i_%i", 3, 3, 1));

    //...................................Fitting by invM..................................
    double Ph, Ph_err, phi;

    //Canvas init

    
    //Fit function init
    TF1 *FitCosTetaFunc;

    TF1 *FitFuncKFP = new TF1("FitFuncKFP", "gaus(0)+gaus(3) + [6] + [7]*x", 1.1, 1.13);
    FitFuncKFP->SetParLimits(1, 1.1, 1.13);
	FitFuncKFP->SetParLimits(2, 0.0005, 0.001);    
    FitFuncKFP->SetParLimits(4, 1.113, 1.12);
	FitFuncKFP->SetParLimits(5, 0.001, 0.005);
    FitFuncKFP->SetLineWidth(2);

    //graphs init
    TGraphErrors *grSin_diffPhiPsi1 = ConvertProfToGraph(prSin_diffPhiPsi1);

    TCanvas *can = new TCanvas("can", "", 0, 1024, 1800, 900);
    can->Divide(2, 1);
    
    
    //....................Lam invM..................
    can->cd(1);

    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    InvMLamDist->Fit(FitFuncKFP, "R");
    
    TF1 *Gaus1 = new TF1("Gaus1", "gaus(0)", 1.1, 1.13);
    TF1 *Gaus2 = new TF1("Gaus2", "gaus(0)", 1.1, 1.13);
    TF1 *BackGround = new TF1("BackGround", "[0] + [1]*x", 1.1, 1.13);
    Gaus1->SetParameter(0, FitFuncKFP->GetParameter(0));
    Gaus1->SetParameter(1, FitFuncKFP->GetParameter(1));
    Gaus1->SetParameter(2, FitFuncKFP->GetParameter(2));
    Gaus1->SetLineColor(4);
    Gaus2->SetParameter(0, FitFuncKFP->GetParameter(3));
    Gaus2->SetParameter(1, FitFuncKFP->GetParameter(4));
    Gaus2->SetParameter(2, FitFuncKFP->GetParameter(5));
    Gaus2->SetLineColor(1);
    BackGround->SetParameter(0, FitFuncKFP->GetParameter(6));
    BackGround->SetParameter(1, FitFuncKFP->GetParameter(7));
    BackGround->SetLineColor(6);
    
    
    InvMLamDist->SetTitle("Invariant mass of #Lambda");
    InvMLamDist->GetYaxis()->SetTitle("Entries");
    InvMLamDist->GetYaxis()->SetTitleFont(42);
    InvMLamDist->GetYaxis()->SetLabelFont(42);
    InvMLamDist->GetXaxis()->SetTitle("M_{inv} [GeV/c^{2}]");
    InvMLamDist->GetXaxis()->SetTitleFont(42);
    InvMLamDist->GetXaxis()->SetLabelFont(42);
    InvMLamDist->GetXaxis()->SetLabelSize(0.045);
	InvMLamDist->GetXaxis()->SetTitleSize(0.045);
    InvMLamDist->GetYaxis()->SetTitleOffset(1.5);
	InvMLamDist->GetYaxis()->SetLabelSize(0.045);
	InvMLamDist->GetYaxis()->SetTitleSize(0.045);
    InvMLamDist->SetAxisRange(-2e4, 2e5, "Y");
    InvMLamDist->Draw("e");
    
    double mu_1 = FitFuncKFP->GetParameter(1), mu_2 = FitFuncKFP->GetParameter(4);
    double mu_err_1 = FitFuncKFP->GetParError(1), mu_err_2 = FitFuncKFP->GetParError(4);
    double sigma_1 = FitFuncKFP->GetParameter(2), sigma_2 = FitFuncKFP->GetParameter(5);
    double sigma_err_1 = FitFuncKFP->GetParError(2), sigma_err_2 = FitFuncKFP->GetParError(5);
    
    double mu = (mu_1 + mu_2)*0.5, mu_err = pow(mu_err_1*mu_err_1 + mu_err_2*mu_err_2,0.5)/2.;
    double sigma = pow(2*sigma_1*sigma_1 + 2*sigma_2*sigma_2 + (mu_1-mu_2)*(mu_1-mu_2),0.5)/2;
    double sigma_err = pow(pow(sigma_1*sigma_err_1,2) + pow(sigma_2*sigma_err_2,2) + pow((mu_1-mu_2)*mu_err, 2),0.5) / (2*sigma);

    TLegend *leg_mass = new TLegend();
    leg_mass->AddEntry(FitFuncKFP, "N(#mu_{1}, #sigma_{1}) + N(#mu_{2}, #sigma_{2}) + kx + b");
    leg_mass->SetBorderSize(0.);
    leg_mass->Draw("same");

    TLatex tex;
    tex.SetTextSize(0.035);
    tex.DrawLatex(1.101, 1.65e5, Form("#mu = %.7f #pm (%.2f #upoint 10^{-6})", mu, mu_err*1e6));
    tex.DrawLatex(1.101, 1.5e5, Form("#sigma = %.6f #pm (%.1f #upoint 10^{-6})", sigma, sigma_err*1e6));

    
    FuncGausKFP->SetParameter(0, FitFuncKFP->GetParameter(0));
    FuncGausKFP->SetParameter(1, FitFuncKFP->GetParameter(1));
    FuncGausKFP->SetParameter(2, FitFuncKFP->GetParameter(2));
    FuncGausKFP->SetParameter(3, FitFuncKFP->GetParameter(3));
    FuncGausKFP->SetParameter(4, FitFuncKFP->GetParameter(4));
    FuncGausKFP->SetParameter(5, FitFuncKFP->GetParameter(5));

    LineKFP->SetParameter(0, FitFuncKFP->GetParameter(6));
    LineKFP->SetParameter(1, FitFuncKFP->GetParameter(7));

    can->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    //<Cos(phi-Psi1)>
    
    FitCosTetaFunc = new TF1("FitCosTetaFunc", MyFitFunc, 1.1, 1.13, 3);
    
    grSin_diffPhiPsi1->Fit(FitCosTetaFunc);
    //grSin_diffPhiPsi1->Draw("AP");

    TFitResultPtr r = grSin_diffPhiPsi1->Fit(FitCosTetaFunc, "S");

    grSin_diffPhiPsi1->SetTitle(" ");
    grSin_diffPhiPsi1->GetXaxis()->SetTitle("M_{inv}, [GeV/c^{2}]");
    grSin_diffPhiPsi1->GetXaxis()->SetTitleFont(42);
    grSin_diffPhiPsi1->GetXaxis()->SetLabelFont(42);
    grSin_diffPhiPsi1->GetYaxis()->SetTitle("<sin(#Psi_{1} - #phi)>");
    grSin_diffPhiPsi1->GetYaxis()->SetTitleFont(42);
    grSin_diffPhiPsi1->GetYaxis()->SetLabelFont(42);
    grSin_diffPhiPsi1->GetXaxis()->SetLabelSize(0.045);
	grSin_diffPhiPsi1->GetXaxis()->SetTitleSize(0.045);
    grSin_diffPhiPsi1->GetYaxis()->SetTitleOffset(1.5);
	grSin_diffPhiPsi1->GetYaxis()->SetLabelSize(0.045);
	grSin_diffPhiPsi1->GetYaxis()->SetTitleSize(0.045);
    grSin_diffPhiPsi1->SetMaximum(0.06);
    grSin_diffPhiPsi1->SetMinimum(-0.06);
    
    grSin_diffPhiPsi1->Draw("AP");

    TLegend *leg_Draw = new TLegend();
    
    leg_Draw->AddEntry(FitCosTetaFunc, "(1-f^{BG})p_{0} + f^{BG}(p_{1} + p_{2}M_{inv})", "l");
    leg_Draw->AddEntry((TObject*)0, Form("p_{0} = %.4f #pm %.4f", r->Value(0), r->ParError(0)), "");
    leg_Draw->AddEntry((TObject*)0, Form("#chi^{2} = %.1f/NDF = %i", r->Chi2(), r->Ndf()), "");   
    
    leg_Draw->SetBorderSize(0);
    leg_Draw->Draw("same");

       
    TText t_1;
    t_1.SetTextFont(42);
    t_1.DrawText(1.115, 0.04, "Centrality 30 - 40%");
    tex.DrawLatex(1.115, 0.03, "#Delta#phi = (3#pi/6, 4#pi/6)");



    //for(int iSub=0; iSub!=2; iSub++)

    //...................................end of Fitting by invM..................................
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