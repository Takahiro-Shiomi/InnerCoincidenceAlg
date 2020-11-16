#include "innercoincidence.h"
#include "/home/shiomi/RootUtils/src/rootlogon.C"
#include <iostream>
#include <stdlib.h>
#include <typeinfo>
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

int main(int argc, char *argv[])
{
    rootlogon();
    int palette[49];
    for(int i=0;i!=49;i++){palette[i] = 51 + i;}
    gStyle->SetPalette(49,palette);

    std::cout<<"<<Run START>>"<<std::endl;
    TChain *tree1 = new TChain("physics", "physics");

    tree1->Add(argv[1]);
    TString PdfLabel = argv[2];
    Int_t limit_entry = atoi(argv[3]);

    TFile *fout = new TFile(Form("./rootfile/%s.root", PdfLabel.Data()), "recreate");

    innercoincidence tls(tree1);
    std::cout<<"<<Loop START>>"<<std::endl;
    tls.Loop(limit_entry);
    std::cout<<"<<Draw START>>"<<std::endl;
    tls.Draw("./pdf/Draw_" + PdfLabel + ".pdf");
    fout -> Write();
    std::cout<<"<<Run END>>"<<std::endl;
    tls.End();

    return 0;
}
