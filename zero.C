#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReadervalue.h"
#include "TGraph.h"
#include <TNtuple.h>
#include <TMath.h>
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

void zero() {

    // data vectors

    const Int_t set = 3; 
    const Int_t n = 9*set;
    Double_t x0[n] = {
        0.5, 2.5, 5.0, 8.7, 12.5, 17.5, 21.075, 25.0, 29.287
    };
    Double_t y0[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t ex0[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t ey0[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    Double_t y1[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t ey1[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    Double_t y2[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t ey2[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // read root files 

    for (Int_t i=0; i<n; i++){
        auto id = std::to_string(i);
        string fileName = "run"+id+".root";
        TFile f(fileName.c_str());
        TTree *t = (TTree*)f.Get("hits");
        t->Draw("Scnt");
        TH1D *h = (TH1D*)gPad->GetPrimitive("htemp");
        if(i<9){
            y0[i] = h->GetMean();
            ey0[i] = h->GetStdDev();
        }
        if(i>8){
            y1[i-9] = h->GetMean();
            ey1[i-9] = h->GetStdDev();
        }
        if(i>17){
            y2[i-18] = h->GetMean();
            ey2[i-18] = h->GetStdDev();
        }
    }

    // plot data
    auto c1 = new TCanvas("c1","c1",200,10,700,500);
    c1->Divide(2,3);
    c1->SetGrid();

    auto g0 = new TGraphErrors(n/set, x0, y0, ex0, ey0); // R=100%
    auto g1 = new TGraphErrors(n/set, x0, y1, ex0, ey1); // R=99%
    auto g2 = new TGraphErrors(n/set, x0, y2, ex0, ey2); // R=99.5%
    auto hig0 = new TGraph("hig.dat"); // R_exp
    auto hig1 = new TGraph("hig.dat");
    auto hig2 = new TGraph("hig.dat");
    g0->SetMarkerColor(4);
    g0->SetMarkerStyle(21);
    g1->SetMarkerColor(1);
    g1->SetMarkerStyle(21);
    g2->SetMarkerColor(2);
    g2->SetMarkerStyle(21);
    hig0->SetMarkerColor(3);
    hig0->SetMarkerStyle(21);
    hig1->SetMarkerColor(3);
    hig1->SetMarkerStyle(21);
    hig2->SetMarkerColor(3);
    hig2->SetMarkerStyle(21);

    Int_t n_hig = hig0->GetN();
    Double_t last_g0 = g0->GetPointY(1);
    Double_t last_g1 = g1->GetPointY(3);
    Double_t last_g2 = g1->GetPointY(3);
    Double_t last_hig0 = hig0->GetPointY(0); 
    Double_t last_hig1 = hig0->GetPointY(1); 
    Double_t last_hig2 = hig0->GetPointY(1); 
    Float_t scale0 = last_g0/last_hig0; // scaling factors
    Float_t scale1 = last_g1/last_hig1;
    Float_t scale2 = last_g2/last_hig2;
    hig0->Scale(scale0);
    hig1->Scale(scale1);
    hig2->Scale(scale2);

    std::vector<Double_t> v0; // analysis
    std::vector<Double_t> v1;
    std::vector<Double_t> v2;

    v0.insert( v0.end(), g0->GetX(), (g0->GetX() + g0->GetN()) );
    v0.insert( v0.end(), hig0->GetX(), (hig0->GetX() + hig0->GetN()) );
    std::sort( v0.begin(), v0.end() );
    v0.erase( std::unique( v0.begin(), v0.end() ), v0.end() );
    Int_t size = v0.size();
    auto diff0 = new TGraph(size-1);

    v1.insert( v1.end(), g1->GetX(), (g1->GetX() + g1->GetN()) );
    v1.insert( v1.end(), hig1->GetX(), (hig1->GetX() + hig1->GetN()) );
    std::sort( v1.begin(), v1.end() );
    v1.erase( std::unique( v1.begin(), v1.end() ), v1.end() );
    auto diff1 = new TGraph(size-1);

    v2.insert( v2.end(), g2->GetX(), (g2->GetX() + g2->GetN()) );
    v2.insert( v2.end(), hig2->GetX(), (hig2->GetX() + hig2->GetN()) );
    std::sort( v2.begin(), v2.end() );
    v2.erase( std::unique( v2.begin(), v2.end() ), v2.end() );
    auto diff2 = new TGraph(size-1);

    for (Int_t i = 0; i < size-1; i++){
        diff0->SetPoint(i, v0[i+1], (g0->Eval(v0[i+1]) - hig0->Eval(v0[i+1])) );
        diff0->SetTitle("diff R=1;d;#Delta N_{#gamma}");

        diff1->SetPoint(i, v1[i+1], (g1->Eval(v1[i+1]) - hig1->Eval(v1[i+1])) );
        diff1->SetTitle("diff R=0.99;d;#Delta N_{#gamma}");

        diff2->SetPoint(i, v2[i+1], (g2->Eval(v2[i+1]) - hig2->Eval(v2[i+1])) );
        diff2->SetTitle("diff R=0.995;d;#Delta N_{#gamma}");
    }

    auto g_err00 = new TGraph(n/set, x0, ey0);
    auto g_err01 = new TGraph(n/set, x0, ey0);
    g_err00->SetLineStyle(2);
    g_err01->SetLineStyle(2);
    g_err01->Scale(-1);
    auto Xaxis0 = diff0->GetXaxis();
    auto Yaxis0 = diff0->GetYaxis();
    Xaxis0->SetLimits(2.5, 30);
    Yaxis0->SetRangeUser(-30, 30);

    auto g_err10 = new TGraph(n/set, x0, ey1);
    auto g_err11 = new TGraph(n/set, x0, ey1);
    g_err10->SetLineStyle(2);
    g_err11->SetLineStyle(2);
    g_err11->Scale(-1);
    auto Xaxis1 = diff1->GetXaxis();
    auto Yaxis1 = diff1->GetYaxis();
    Xaxis1->SetLimits(2.5, 30);
    Yaxis1->SetRangeUser(-30, 30);

    auto g_err20 = new TGraph(n/set, x0, ey2);
    auto g_err21 = new TGraph(n/set, x0, ey2);
    g_err20->SetLineStyle(2);
    g_err21->SetLineStyle(2);
    g_err21->Scale(-1);
    auto Xaxis2 = diff2->GetXaxis();
    auto Yaxis2 = diff2->GetYaxis();
    Xaxis2->SetLimits(2.5, 30);
    Yaxis2->SetRangeUser(-30, 30);
    
    std::cout << "scaled up: " << scale0 << " ; " << scale1 << " ; " << scale2 << std::endl;
    c1->cd(1); g0->Draw("APL"); hig0->Draw("PL SAME");
    c1->cd(2); diff0->Draw("APL"); g_err00->Draw("PL SAME"); g_err01->Draw("PL SAME");
    c1->cd(3); g2->Draw("APL"); hig2->Draw("PL SAME");
    c1->cd(4); diff2->Draw("APL"); g_err20->Draw("PL SAME"); g_err21->Draw("PL SAME");
    c1->cd(5); g1->Draw("APL"); hig1->Draw("PL SAME");
    c1->cd(6); diff1->Draw("APL"); g_err10->Draw("PL SAME"); g_err11->Draw("PL SAME");
}
