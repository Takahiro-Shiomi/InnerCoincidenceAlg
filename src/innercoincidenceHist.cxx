#include "innercoincidence.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <typeinfo>
#include <TMath.h>

void innercoincidence::InitHist()
{
    S_xy = new TH2D("",";x^{NSW};y^{NSW}",1400,-7000,7000,1400,-7000,7000);
    N_etaphi = new TH2D("N_etaphi",";#eta^{NSW};#phi^{NSW}",100,-2.7,2.7,96,-TMath::Pi(),TMath::Pi());
}

void innercoincidence::End()
{
    if(S_xy!=0){delete S_xy;}
    if(N_etaphi!=0){delete N_etaphi;}
}
