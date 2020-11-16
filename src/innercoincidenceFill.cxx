#include "innercoincidence.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <typeinfo>

void innercoincidence::FillHist()
{
    Clear();
    for(int l1=0;l1<trig_L1_mu_n;l1++){
        if((*trig_L1_mu_thrNumber)[l1]!=6){continue;}
        if((*trig_L1_mu_sectorAddress)[l1]<64){continue;}
       
        float c_eta=-100,c_phi=-100,c_theta=-100,dR=100;
	    if(fabs((*trig_L1_mu_sectorAddress)[l1]) > 127){
            Endcap(l1, c_eta, c_phi, c_theta, dR);
            if(c_eta==-100){Tile(l1, c_eta, c_phi, c_theta, dR);}
        }else{
            Forward(l1, c_eta, c_phi, c_theta, dR);
        }
        if(c_eta!=-100){
            N_etaphi->Fill(c_eta,c_phi);
            /*
            deta = (*trig_L1_mu_eta)[l1] - c_eta;
            dphi = (*trig_L1_mu_phi)[l1] - c_phi;
            deta = 2*atan(exp(-(*trig_L1_mu_eta)[l1])) - c_theta;
            */
        }
    }
}

void innercoincidence::Clear()
{
}
