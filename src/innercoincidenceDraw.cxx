#include "innercoincidence.h"
#include "/home/shiomi/RootUtils/RootUtils/TLegend_addfunc.h"
#include "/home/shiomi/RootUtils/RootUtils/TCanvas_opt.h"
#include "TStyle.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TGraphErrors.h"

using namespace std;

void innercoincidence::Draw(TString pdf)
{
    TCanvas_opt *c1 = new TCanvas_opt();
    gStyle->SetOptStat(0);
    c1->SetTopMargin(0.20);
    c1->Print(pdf + "[", "pdf");
    
    S_xy->Draw("colz");
    c1->Print(pdf,"pdf");

    N_etaphi->Draw("colz");
    c1->Print(pdf,"pdf");

    c1 -> Print( pdf + "]", "pdf" );
    delete c1;
}
