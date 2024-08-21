//------ Luda's version ----------------
// Correlation Functions CF(k*) with QS+FSI weights calculating by R.Lednicky's code 
//You have to compile it using MakeFile
// make
//./CFs_Led
//author of interface to R. Lednicky's fortran code Malinina L.V: malinina@lav01.sinp.msu.ru 


#include <TH3.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TH1D.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <iostream>
#include <fstream>
#include "WLedCOMMONS.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <cstring>
using namespace std;

#ifndef precision
#define precision 1e-8;
#endif
#ifndef euler_constant
#define euler_constant 0.5772156649015328606
#endif
#include <TFile.h>
#include <iostream>
#include <TRandom.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

extern "C" void   fsiini_();
extern "C" void   ltran12_();
extern "C" void   fsiw_();
extern FSICONSCommon FSICONS;
extern LEDWEIGHTCommon LEDWEIGHT;
extern MOMLABCommon MOMLAB;
extern FSIMOMCommon FSIMOM;
extern FSINSCommon FSINS;
extern FSICOORCommon FSICOOR;
extern FSIPOCCommon FSIPOC;
extern FSIPRFCommon FSIPRF;


const double HBARC=0.1973269602; //1/HBARC fm->GeV-1
double pi=3.14159265358979323844;


const int NPT = 1 ;
const int NEVPACK = 1, NBINS = 8 ;
 
int main()
{
                                                                                               

   TH1D *hCFnom;
   TH1D *hCFden;
   TH1D *hCF;
   TH1D *hCFP ;                                                                                                              
   TH1D *HNum2 ;     
   
   //pp
   //int nbins=150;
   //double qmin=0.0;
   //double qmax=1.5;
 
 //PbPb
   int nbins=100;
   double qmin=0.0;
   double qmax=0.9;
   
   
   
   //hCFnom  = new TH1D("hCFnom","hCFnom",150,0.,1.5);
   //hCFden  = new TH1D("hCFden","hCFden",150,0.,1.5);
   //hCF  = new TH1D("hCF","hCF",150,0.,1.5);
   //hCFP    = new TH1D("hCFP","hCFP",150,0.,1.5);  
   //HNum2=   new TH1D("HNum2","HNum2",150,0.,1.5);     
     hCFnom  = new TH1D("hCFnom","hCFnom",nbins,qmin,qmax);
     hCFden  = new TH1D("hCFden","hCFden",nbins,qmin,qmax);
     hCF  = new TH1D("hCF","hCF",nbins,qmin,qmax);
     hCFP    = new TH1D("hCFP","hCFP",nbins,qmin,qmax);  
     HNum2=   new TH1D("HNum2","HNum2",nbins,qmin,qmax);      


  const Int_t kMax = 1000000;


  int nqmax,icthetamin;
  const int nctheta=20;
  int iq,nsamples,isample,alpha,ictheta;
  double ctheta,stheta,ctheta_qr,r[4],rmag,q,psisquared;
  double Rgauss[4],zoffset=0.0;
  double psi2=0, psi2_5=0, psi2A=0, test=0, weight=0, k=0, kstar=0 ;

  double qMax = qmax;// 1.5; //in GeV/c
  nqmax=nbins;//150;//150;
  double qStep = qMax/nqmax;
  nsamples=1000;
  cout << "qmax "<<qMax << "qstep " << qStep <<"nsamples"<<nsamples<< endl;


//Source size in fm
  Rgauss[1]=5.0; //out
  Rgauss[2]=5.0; //side
  Rgauss[3]=5.0; //long
  zoffset=0.;



//-- Lednicky's weight calculation initialisation

  //initial parameters of model
  //  Bethe-Salpeter amplitude
  //   NS=1  Square well potential,
  //   NS=3  not used
  //   NS=4  scattered wave approximated by the spherical wave,
  //   NS=2  same as NS=4 but the approx. of equal emission times in PRF
  //         not required (t=0 approx. used in all other cases).

   FSINS.NS=4;
   
 // Switch for automatic setting of all parameters (0)
 //ITEST=1 any values of parameters ICH, IQS, ISI, I3C are allowed
//ITEST=0 physical values of these parameters are put automatically
   LEDWEIGHT.ITEST=1;
    //  Swith for Couloumb interaction in the pair
    FSINS.ICH =1;
    // Switches strong interactions
    FSINS.ISI=0;
    // Switch for quantum statistics
    FSINS.IQS=1;
    //Switches couloumb interaction with residual nucleus 
    FSINS.I3C=0;

  //initial parameters of model
 //if(fRandomPosition)
    LEDWEIGHT.IRANPOS=1;
 

//C----------------------------------------------------------------------
//C-   LL       1  2  3  4  5   6   7   8  9 10  11  12  13  14 15 16 17
//C-   part. 1: n  p  n  a  pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
//C-   part. 2: n  p  p  a  pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p
//C   NS=1 y/n: +  +  +  +  +   -   -   -  -  -   -   -   -  -  -  -  -
//C----------------------------------------------------------------------
//C-   LL       18 19 20 21 22 23  24 25 26 27 28 29 30 31  32
//C-   part. 1: d  d  t  t  K0 K0  d  p  p  p  n  /\ p  pi+ pi-
//C-   part. 2: d  a  t  a  K0 K0b t  t  a  /\ /\ /\ pb Xi- Xi-
//C   NS=1 y/n: -  -  -  -  -  -   -  -  -  +  +  +  -  -   -
//C----------------------------------------------------------------------
  
  
  
    FSINS.LL = 4; 

//pdg codes of particles pair

    Int_t pdg1= 321;//211;
    Int_t pdg2= 321;  //3312;//321;//2212;//321;//3312;

 
    TParticlePDG* tpart1 = TDatabasePDG::Instance()->GetParticle(pdg1);
    TParticlePDG* tpart2 = TDatabasePDG::Instance()->GetParticle(pdg2);

    double mass1=tpart1->Mass();
    double mass2=tpart2->Mass();

    
    FSIPOC.AM1=mass1;
    FSIPOC.C1=1.;  
    FSIPOC.AM2=mass2; //0.492;  //0.938;  //1.13213; 
    FSIPOC.C2=1.;  

    
    fsiini_();


    std::cout<<"charge1 "<<FSIPOC.C1<<" mass1 "<< FSIPOC.AM1<<std::endl;
    std::cout<<"charge2 "<<FSIPOC.C2<<" mass2 "<< FSIPOC.AM2<<std::endl;

//q --k*

  for(iq=0;iq<nqmax;iq++){

    q=qStep/2.+iq*qStep;

    for(ictheta=0;ictheta<nctheta;ictheta++){

      for(isample=0;isample<nsamples;isample++){
	r[0]=0.0;
	for(alpha=1;alpha<=3;alpha++){
	  r[alpha]=sqrt(2.0)*Rgauss[alpha] * (gRandom->Gaus());///HBARC; //fm --> 1/MeV
	  }
	r[3]+=zoffset;///HBARC;
	rmag=sqrt(r[1]*r[1]+r[2]*r[2]+r[3]*r[3]); //MeV-1

        ctheta=double(ictheta+0.5)/double(nctheta); // ctheta is q relative to zhat
        stheta=sqrt(1.0-ctheta*ctheta);
       
        ctheta_qr=(r[1]*stheta+r[3]*ctheta)/rmag;
        
         rmag=rmag/HBARC;
         
         //q- Qinv

	 if(q<qMax/*1.5*/){

                FSIPRF.PPX=0;  
                FSIPRF.PPY=0;  
                FSIPRF.PPZ=q*0.5;  
                FSIPRF.AK=q*0.5;  
                FSIPRF.AKS=q*q*0.25;  

                FSIPRF.X=r[1]/HBARC;  
                FSIPRF.Y=r[2]/HBARC;  
                FSIPRF.Z=r[3]/HBARC;
                FSIPRF.T=r[0]/HBARC;      
                FSIPRF.RP=rmag;  
                FSIPRF.RPS=rmag*rmag;  
                  


//      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! k*=(p1-p2)/2 and x1-x2
//     1               X,Y,Z,T,RP,RPS      ! in pair rest frame (PRF)

// Komponenty vektorov r* i k* nuzny.
// > vam dal'she PPX,PPY,PPZ ne nujni v fsiw.f ?
// > ja rabotaiu v peremennih k*, theta, r
//    Esli theta= ugol mezdu r* i k*, mozete polozit' naprimer
//       r*-vektor= r*{sin(theta),0,cos(theta)}
//          k*-vektor= k*(0,0,1)
          
                fsiw_();

   //get weight from Lednicky's commons		

                weight = float(LEDWEIGHT.WEIN);
	        hCFnom->Fill(q, weight);
                hCFden->Fill(q, 1.);
                HNum2->Fill(q,weight*weight);    
        
		 
      }
    }    
  } 

}



 TFile outputFile("alpha.root", "RECREATE");
 

  // for(int ipt=0; ipt<NPT; ipt++)
//		{
   hCFnom->Sumw2();
   hCFden->Sumw2();
   hCF->Sumw2();

   double err=0.001;

   for(int ii=1; ii<= hCFP->GetNbinsX(); ii++){
    hCFP->SetBinContent(ii, hCFnom->GetBinContent(ii)/hCFden->GetBinContent(ii));
    if(HNum2->GetBinContent(ii)>0 && hCFden->GetBinContent(ii)>0){err = sqrt((HNum2->GetBinContent(ii)/hCFden->GetBinContent(ii)-hCFP->GetBinContent(ii)*hCFP->GetBinContent(ii))/hCFden->GetBinContent(ii));
    //std::cout<<" ii "<<ii<<" "<<hCFnom->GetBinContent(ii) <<" err "<<err<< std::endl;
    std::cout<<" ii "<<ii<<" "<<ii*qStep<<" "<<hCFP->GetBinContent(ii) <<" err "<<err<< std::endl;
    if(abs(err)<=0.0001)err=0.001;
    hCFP->SetBinError(ii,err); 
    }
    else
    {
    hCFP->SetBinContent(ii,0.0);
    hCFP->SetBinError(ii,0.1); 
    }
    
//    HCFP.SetBinContent(ti, HNumP.GetBinContent(ti)/HDenP.GetBinContent(ti));
//    HCFP.SetBinError(ti, sqrt((HNum2P.GetBinContent(ti)/HDenP.GetBinContent(ti)-HCFP.GetBinContent(ti)*HCFP.GetBinContent(ti))/HDenP.GetBinContent(ti)));
//    HCFN.SetBinContent(ti, HNumN.GetBinContent(ti)/HDenN.GetBinContent(ti));
//    HCFN.SetBinError(ti, sqrt((HNum2N.GetBinContent(ti)/HDenN.GetBinContent(ti)-HCFN.GetBinContent(ti)*HCFN.GetBinContent(ti))/HDenN.GetBinContent(ti)));
  }
 

   hCF->Divide(hCFnom,hCFden, 1.0, 1.0);
   hCF->Write();
   hCFnom->Write();
   hCFden->Write();
    hCFP->Write();  

   delete hCF ;
   delete hCFnom ;
   delete hCFden ;

//		}


 
 outputFile.Write();


}
