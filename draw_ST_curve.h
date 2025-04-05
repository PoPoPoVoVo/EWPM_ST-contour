#include <TMath.h>

TGraph* draw_ST_curve(){
    const double mH_ref = 1000.0;
    const double cos2w = 0.7689;
    const double mZ = 91.19;
    const double mW = 80.39;

    const double Bwh_ref = TMath::Power(mH_ref,2) - TMath::Power(mW,2);
    const double Bzh_ref = TMath::Power(mH_ref,2) - TMath::Power(mZ,2);
    const double Cwh_ref = TMath::Power(mW,2)/(TMath::Power(mH_ref,2)-TMath::Power(mW,2));
    const double Czh_ref = TMath::Power(mZ,2)/(TMath::Power(mH_ref,2)-TMath::Power(mZ,2));

    const double Slog_ref = -(1.0/3.0)*( (1.0/3.0)*(3.0*TMath::Power(Czh_ref,2)+3.0*Czh_ref+1.0)- (3.0/2.0)*(1.0+Czh_ref)*(2.0*Czh_ref+1.0) +3.0*TMath::Power((1.0+Czh_ref),2) - TMath::Power((1.0+Czh_ref),3)*log((1+Czh_ref)/Czh_ref));
    const double Slog2_ref = 3.0*Czh_ref*( (1.0/2.0)*(2.0*Czh_ref+1.0) -(2.0*Czh_ref+1) +(TMath::Power(Czh_ref,2)+Czh_ref )*log((1+Czh_ref)/Czh_ref));
    const double Tlog_ref = (1.0+Cwh_ref)*log(1.0+Cwh_ref)-Cwh_ref*log(Cwh_ref) -(1.0/cos2w)*((1.0+Czh_ref)*log(1.0+Czh_ref)-Czh_ref*log(Czh_ref));
    
    double S_cp = 0.00;
    double T_cp = 0.00;

    std::vector<double> mH_vals = {100.0,200.0,300.0,400.0,500.0,600.0};
    
    TGraph* graph = new TGraph();
    int i = 0;

    for (double mH =100.0 ; mH<= 5000.0; mH += 1.0){
        double Czh = (TMath::Power(mZ,2))/(TMath::Power(mH,2)-TMath::Power(mZ,2));
        double Cwh = (TMath::Power(mW,2))/(TMath::Power(mH,2)-TMath::Power(mW,2));
        double Bwh = TMath::Power(mH,2) - TMath::Power(mW,2);
        double Bzh = TMath::Power(mH,2) - TMath::Power(mZ,2);

        double Slog = -(1.0/3.0)*( (1.0/3.0)*(3.0*TMath::Power(Czh,2)+3.0*Czh+1.0)- (3.0/2.0)*(1+Czh)*(2.0*Czh+1.0) +3.0*TMath::Power((1.0+Czh),2) - TMath::Power((1.0+Czh),3)*log((1.0+Czh)/Czh));
        double Slog2 = 3.0*Czh*( (1.0/2.0)*(2.0*Czh+1.0) -(2.0*Czh+1.0) +(TMath::Power(Czh,2)+Czh)*log((1.0+Czh)/Czh));

        double Tlog = (1.0+Cwh)*log(1.0+Cwh)-Cwh*log(Cwh) -(1.0/cos2w)*((1.0+Czh)*log(1.0+Czh)-Czh*log(Czh));
        
        double S = (1.0/ (4.0*TMath::Pi())) * (Slog - Slog_ref+ Slog2 - Slog2_ref);
        double T = (3.0/ (16.0*TMath::Pi()* (1.0-cos2w) )) * (Tlog-Tlog_ref+ log(Bwh/Bwh_ref)- (1.0/cos2w)*log(Bzh/Bzh_ref));

        graph->SetPoint(i++,S,T);
        if (mH == 175.00) {
            S_cp = S;
            T_cp = T;
        }

        // std::cout << Form ("mH = %0.f GeV -> S= %.4f, T= %.4f", mH,S,T) <<std::endl;
    }

    graph->SetLineWidth(1);
    graph->SetLineColor(kBlack);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.1);
    graph->SetMarkerColor(kBlack);

    return graph;
}

