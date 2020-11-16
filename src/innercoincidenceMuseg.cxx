#include "innercoincidence.h"
#include <iostream>
#include <vector>
#include <TVector3.h>
#include <stdlib.h>
#include <typeinfo>
#include <TRandom.h>

void innercoincidence::Endcap(int l1, float &c_eta, float &c_phi, float &c_theta, float dR)
{
    float s_x=-100,s_y=-100;
    for(int seg=0;seg<museg_n;seg++){
        //select MDT EI segment stationName:7=EIS,8=EIL
        if((*museg_stationName)[seg]!=7 && (*museg_stationName)[seg]!=8){continue;}
        TVector3 vseg;
        vseg.SetXYZ((*museg_x)[seg], (*museg_y)[seg], (*museg_z)[seg]);
        float s_eta = vseg.Eta() + gRandom->Gaus(0,0.005);
        float s_theta = vseg.Theta() + gRandom->Gaus(0,0.001);

        float s_phi = TGCPRD(l1);

        //delta Eta matching roughly
	    if(fabs(s_eta - (*trig_L1_mu_eta)[l1])>0.1){continue;}
        //delta theta matching roughly
	    if(acos(cos(s_theta-2*atan(exp(-(*trig_L1_mu_eta)[l1]))))>0.016){continue;}
        //dphi matching roughly 
        if(acos(cos(s_phi-(*trig_L1_mu_phi)[l1]))>0.16){continue;}
	    if(hypot((*trig_L1_mu_eta)[l1]-s_eta, acos(cos((*trig_L1_mu_phi)[l1]-s_phi)))<dR){
            dR=hypot((*trig_L1_mu_eta)[l1]-s_eta, acos(cos((*trig_L1_mu_phi)[l1]-s_phi)));
            c_eta = s_eta;
            c_phi = s_phi;
            c_theta = s_theta;
            s_x = (*museg_x)[seg];
            s_y = (*museg_y)[seg];
        }
    }
    if(s_x!=-100){S_xy->Fill(s_x,s_y);}
}

void innercoincidence::Forward(int l1, float &c_eta, float &c_phi, float &c_theta, float dR)
{
    float s_x=-100,s_y=-100;
    for(int seg=0;seg<museg_n;seg++){
        //select CSC EI segment stationName:15=CSS,16=CSL
        if((*museg_stationName)[seg]!=15 && (*museg_stationName)[seg]!=16){continue;}
        TVector3 vseg;
        vseg.SetXYZ((*museg_x)[seg], (*museg_y)[seg], (*museg_z)[seg]);
        float s_eta = vseg.Eta() + gRandom->Gaus(0,0.005);
        float s_phi = vseg.Phi() + gRandom->Gaus(0,0.01);
        float s_theta = vseg.Theta() + gRandom->Gaus(0,0.001);

        //delta Eta matching roughly
	    if(fabs(s_eta - (*trig_L1_mu_eta)[l1])>0.1){continue;}
        //delta theta matching roughly
	    if(acos(cos(s_theta-2*atan(exp(-(*trig_L1_mu_eta)[l1]))))>0.016){continue;}
        //dphi matching roughly 
        if(acos(cos(s_phi-(*trig_L1_mu_phi)[l1]))>0.16){continue;}
	    if(hypot((*trig_L1_mu_eta)[l1]-s_eta, acos(cos((*trig_L1_mu_phi)[l1]-s_phi)))<dR){
            dR=hypot((*trig_L1_mu_eta)[l1]-s_eta, acos(cos((*trig_L1_mu_phi)[l1]-s_phi)));
            c_eta = s_eta;
            c_phi = s_phi;
            c_theta = s_theta;
            s_x = (*museg_x)[seg];
            s_y = (*museg_y)[seg];
        }
    }
    if(s_x!=-100){S_xy->Fill(s_x,s_y);}
}

Float_t innercoincidence::TGCPRD(int l1)
{
    float s_phi = 999.;
    for(int prd=0;prd<TGC_prd_n;prd++){
        if((*TGC_prd_station)[prd]!=47){continue;}
        if((*TGC_prd_eta)[prd] * (*trig_L1_mu_eta)[l1]<0){continue;}
        if((*TGC_prd_isStrip)[prd]!=1){continue;}
        if((*TGC_prd_bunch)[prd]==3){continue;}
        TVector3 vseg;
        vseg.SetXYZ((*TGC_prd_x)[prd], (*TGC_prd_y)[prd], (*TGC_prd_z)[prd]);
        float r_phi = acos(cos((*trig_L1_mu_phi)[l1] - vseg.Phi()));
        float n_phi = acos(cos((*trig_L1_mu_phi)[l1] - s_phi));
        if(r_phi<n_phi){s_phi = vseg.Phi();}
    }
    return s_phi+gRandom->Gaus(0,0.01);
}

void innercoincidence::Tile(int l1,float &c_eta, float &c_phi, float &c_theta, float dR)
{
    int r_sec = ((*trig_L1_mu_sectorAddress)[l1]/2) - 64;
    bool coin = false;
    for(int tile=0;tile<TILE_murcv_trig_n;tile++){
        if((*TILE_murcv_trig_part)[tile]==3 && (*trig_L1_mu_eta)[l1]<0){continue;}
        if((*TILE_murcv_trig_part)[tile]==4 && (*trig_L1_mu_eta)[l1]>0){continue;}
        int t_sec = (*TILE_murcv_trig_mod)[tile];

        if(r_sec>=2){
            int sec1 = (r_sec-2) + (int)(r_sec-2)/3;
            int sec2 = (r_sec-2) + 1 + (int)(r_sec-2)/3;
            if(t_sec==sec1 || t_sec==sec2){
                if((*TILE_murcv_trig_bit0)[tile]==1 || (*TILE_murcv_trig_bit2)[tile]==1){coin = true;}
            }
        }
        if(r_sec<2){
            int sec1 = r_sec + 61;
            int sec2 = r_sec + 62;
            if(t_sec==sec1 || t_sec==sec2){
                if((*TILE_murcv_trig_bit0)[tile]==1 || (*TILE_murcv_trig_bit2)[tile]==1){coin = true;}
            }
        }
    }
    if(coin){
        c_eta = (*trig_L1_mu_eta)[l1] + gRandom->Gaus(0,0.005);
        c_phi = (*trig_L1_mu_phi)[l1] + gRandom->Gaus(0,0.01);
        c_theta = 2*atan(exp(-(*trig_L1_mu_eta)[l1])) + gRandom->Gaus(0,0.001);
    }
}
