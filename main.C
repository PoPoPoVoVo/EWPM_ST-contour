#include "draw_ST_curve.h"
#include "draw_ST_curve_top.h"
// #include "log_Likelihood.h"
#include "log_Likelihood_p.h"
#include <iostream>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>


void draw_ST_graph(){
    TCanvas* c = new TCanvas("c", "S-T contour", 800, 600);
    gStyle->SetPalette(50);

    // Likelihood-Contour
    
    int N = 2000;

    double Smin = -6.0, Smax = 6.0;
    double Tmin = -4.0, Tmax = 4.0;

    TGraph2D *g2d = new TGraph2D();

    // 구역 내 최대 log-likelihood 추정
    double logL_max =-1.00e+9;
    double S_best = 0.00;   //best point : Maiximum Likelihood Point
    double T_best = 0.00;

    double S_cp = 0.00;     // 95% critical point
    double T_cp = 0.00;

    for (int i = 0; i < N; ++i) {
        double S = Smin + (Smax - Smin) * i / (N - 1);
        for (int j = 0; j < N; ++j) {
            double T = Tmin + (Tmax - Tmin) * j / (N - 1);
            double logL = log_likelihood_p(S, T);
            // std::cout << Form ("current_logL = %0.4f",logL) <<std::endl;
            if (logL > logL_max) {
                logL_max = logL;
                S_best = S;
                T_best = T;
            }
            //std::cout << Form ("logL_max = %0.4f",logL_max) <<std::endl;
        }
    }
    std::cout << Form("Best Fit Point: S=%4f, T= %4f",S_best, T_best) <<std::endl;


    // 모든 포인트 delta log-likelihood 추정

    int point = 0;
    for (int i = 0; i < N; ++i) {
        double S = Smin + (Smax - Smin) * i / (N - 1);
        for (int j = 0; j < N; ++j) {
            double T = Tmin + (Tmax - Tmin) * j / (N - 1);
            double logL = log_likelihood_p(S, T);
            double deltalogL= logL - logL_max;
            if (deltalogL < -20) deltalogL = -20;
            g2d->SetPoint(point++, S, T, deltalogL);

        }
    }

    g2d->SetTitle("S-T curve with CL;S;T;logL");
    g2d->GetXaxis()->SetNdivisions(510);
    g2d->GetYaxis()->SetNdivisions(510);

    const int nCL68 = 1;
    double CL68[nCL68] = {-1.15};
    g2d->GetHistogram()->SetContour(nCL68, CL68);
    g2d->Draw("CONT1Z");

    const int nCL95 = 1;
    double CL95[nCL95] = {-2.71};
    TGraph2D* g2d_95 = (TGraph2D*) g2d->Clone("g2d_95");
    g2d_95->GetHistogram()->SetContour(nCL95, CL95);
    g2d_95->Draw("CONT1Z SAME");
    gPad->GetListOfPrimitives()->Remove(gPad->FindObject("palette"));

    // Higgs mass dependence
    TGraph* curve1 = draw_ST_curve();
    curve1->Draw("LP SAME");

    // top quark mass dependence
    TGraph* curve2 = draw_ST_curve_top();
    curve2->Draw("LP SAME");

    // Reference Point
    TMarker* refPoint = new TMarker(0,0,29);
    refPoint->SetMarkerColor(kGreen+2);
    refPoint->SetMarkerSize(1.0);
    refPoint->Draw("SAME");

    // Best Point
    TMarker* bestPoint = new TMarker(S_best,T_best,29);
    bestPoint->SetMarkerColor(kRed);
    bestPoint->SetMarkerSize(1.0);
    bestPoint->Draw("SAME");

    // 범례추가
    TLegend* legend1 = new TLegend(0.15,0.75,0.55,0.88);
    legend1->SetTextSize(0.025);
    legend1->AddEntry(curve1, "SM prediction ( m_{h} )","lp");
    legend1->AddEntry(curve2, "SM prediction ( m_{t} )","lp");
    legend1->AddEntry(refPoint,"Reference (m_{h}=1000 GeV & m_{t}=150 GeV)","p");
    legend1->Draw();

    TLatex* label68 = new TLatex(-1.31, -1.13, "68% CL");
    label68->SetTextColor(kRed);     
    label68->SetTextSize(0.02);
    label68->SetTextAngle(35);      
    label68->Draw();

    TLatex* label90 = new TLatex(-1.18, -1.36, "90% CL");
    label90->SetTextColor(kBlue);    
    label90->SetTextSize(0.02);
    label90->SetTextAngle(35);
    label90->Draw();

    c->SaveAs("S-T contour.pdf");
}
