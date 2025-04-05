#include <TMath.h>


TGraph* draw_ST_curve_top(){
    const double mt_ref = 150.0;
    const double cos2w = 0.7689;
    const double mb = 4.18;
    const double mZ= 91.19;
    const double mW= 80.39;
    
    const double Slog_ref = log(TMath::Power(mt_ref,2));
    const double Tlog_ref = (1.0/2.0)*( (TMath::Power(mt_ref,4)/(TMath::Power(mt_ref,2)-TMath::Power(mb,2)))*log(TMath::Power(mt_ref,2)) - (TMath::Power(mb,4)/(TMath::Power(mt_ref,2)-TMath::Power(mb,2)))*log(TMath::Power(mb,2)) - (TMath::Power(mt_ref,2))*log((TMath::Power(mt_ref,2))) );

    // std::vector<double> mt_vals = {150,170,200,300,400};
    
    TGraph* graph = new TGraph();
    int i =0;

    for (double mt=100.0 ; mt<= 300.0; mt += 0.1 ){
        double Slog = log(TMath::Power(mt,2));
        double Tlog = (1.0/2.0)*( (TMath::Power(mt,4)/(TMath::Power(mt,2)-TMath::Power(mb,2)))*log(TMath::Power(mt,2)) - (TMath::Power(mb,4)/(TMath::Power(mt,2)-TMath::Power(mb,2)))*log(TMath::Power(mb,2)) - (TMath::Power(mt,2))*log((TMath::Power(mt,2))) );
        
        double S = (1.0/TMath::Pi())* (Slog - Slog_ref);
        double T = -(3.0/(4.0*TMath::Pi()*TMath::Power(mW,2)*(1-cos2w)))*( Tlog - Tlog_ref - (1.0/4.0)*(TMath::Power(mt,2)-TMath::Power(mt_ref,2)));

        graph->SetPoint(i++,S,T);
    }

    graph->SetLineWidth(1);
    graph->SetLineColor(kBlue);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.1);
    graph->SetMarkerColor(kBlue);

    return graph;
}