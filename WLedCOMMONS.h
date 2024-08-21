#ifndef WLedCOMMON
#define WLedCOMMON

extern "C" {

#define f2cFortran
#include "cfortran.h"



//----------------------------------------------------------------
/*
 COMMON/CONS/PI,PI2,SPI,DR,W   
*/
	            
 typedef struct //FSI_CONS
   {
   Double_t PI;  // 3.141592654 
   Double_t PI2; // PI2=2*PI 
   Double_t SPI; // SPI=DSQRT(PI)  
   Double_t DR;  // DR=180.D0/PI  from radian to degree  
   Double_t W;   // W=1/.1973D0    from fm to 1/GeV 
  }FSICONSCommon;
 
#define FSICONS COMMON_BLOCK(FSICONS,fsicons)
COMMON_BLOCK_DEF(FSICONSCommon, FSICONS);
//----------------------------------------------------------------
//   COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS  
   
	            
 typedef struct //LEDWEIGHT
   {
   Double_t WEIF;  
   Double_t WEI;   
   Double_t WEIN;  
   Int_t ITEST;
   Int_t IRANPOS;
  }LEDWEIGHTCommon;
 
#define LEDWEIGHT COMMON_BLOCK(LEDWEIGHT,ledweight)
COMMON_BLOCK_DEF(LEDWEIGHTCommon, LEDWEIGHT);

//----------------------------------------------------------------- 
//---------------------------------------------------------------------------       
//             COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S 
//        INTEGER LL,NS,ICH,ISI,IQS,I3C,I3S

                                  
      typedef struct //FSI_NS                                                      
         {                                                                           
	Int_t LL; // [GeV/c]                                                     
	Int_t NS;                                                                
	Int_t ICH;                                                                
	Int_t ISI;                                                                 
	Int_t IQS;                                                                 
	Int_t I3C;                                                                
	Int_t I3S;
	}FSINSCommon;                                                       
			                                                                                      
#define FSINS COMMON_BLOCK(FSINS,fsins)                                  
COMMON_BLOCK_DEF(FSINSCommon, FSINS);                               


//---------------------------------------------------------------------------  
//     COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2                             
                                                                            
 typedef struct //FSI_POC                                                      
    {                                                                           
      Double_t AMN; //mass of the effective nucleus   [GeV/c**2]                   
      Double_t AM1;                                                                
      Double_t AM2;                                                                
      Double_t CN; //charge of the effective nucleus [elem. charge units]          
      Double_t C1;                                                                 
      Double_t C2;                                                                 
      Double_t AC1;                                                                
      Double_t AC2;                                                                
    }FSIPOCCommon;                                                       
		                                                                                     
#define FSIPOC COMMON_BLOCK(FSIPOC, fsipoc)                                  
COMMON_BLOCK_DEF(FSIPOCCommon, FSIPOC);                               
		                                                                                     

//----------------------------------------------------------------- 
/*       COMMON/MOMLAB/AM1,PXP1,PYP1,PZP1,AM2,PXP2,PYP2,PZP2                          
               REAL*8 AM1,PXP1,PYP1,PZP1,AM2,PXP2,PYP2,PZP2                                 
*/	                                                                                            
	            

 typedef struct //MOMLAB
   {
   Double_t AM1;
   Double_t PXP1;
   Double_t PYP1;
   Double_t PZP1;
   Double_t AM2;
   Double_t PXP2;
   Double_t PYP2;
   Double_t PZP2;
   
  }MOMLABCommon;
 
#define MOMLAB COMMON_BLOCK(MOMLAB,momlab)

COMMON_BLOCK_DEF(MOMLABCommon, MOMLAB);

//---------------------------------------------------------------------------  
/*           COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  ! particle momenta in the      
//                                               rest frame of effective nucleu
     1       P2X,P2Y,P2Z,E2,P2                                                 
     */                                                                             
      typedef struct //FSI_MOM                                                      
         {                                                                           
	Double_t P1X; // [GeV/c]                                                     
	Double_t P1Y;                                                                
	Double_t P1Z;                                                                
	Double_t E1;                                                                 
	Double_t P1;                                                                 
	Double_t P2X;                                                                
	Double_t P2Y;                                                                
	Double_t P2Z;                                                                
	Double_t E2;                                                                 
	Double_t P2;                                                                 
	}FSIMOMCommon;                                                       
			                                                                                      
#define FSIMOM COMMON_BLOCK(FSIMOM,fsimom)                                  
COMMON_BLOCK_DEF(FSIMOMCommon, FSIMOM);                               
			                                                                                      
//-----------------------------------------------------------------------                                                                                     
                                                                               
 typedef struct //FSI_COOR                                                     
    {                                                                           
    Double_t X1;                                                                 
    Double_t Y1;                                                                 
    Double_t Z1;                                                                 
    Double_t T1;                                                                 
    Double_t R1;                                                                 
    Double_t X2;                                                                 
    Double_t Y2;                                                                 
    Double_t Z2;                                                                 
    Double_t T2;                                                                 
    Double_t R2;                                                                 
    }FSICOORCommon;                                                      
			                                                                                 
#define FSICOOR COMMON_BLOCK(FSICOOR,fsicoor)                               
COMMON_BLOCK_DEF(FSICOORCommon, FSICOOR);                             


//---------------------------------------------------------------------------         
                                                                                      
/*     COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! k*=(p1-p2)/2 and x1-x2                    
                     X,Y,Z,T,RP,RPS      ! in pair rest frame (PRF)                   
		     */                                                                                    
		                                                                                           
typedef struct //FSI_PRF                                                             
 {                                                                                  
  Double_t PPX;                                                                       
  Double_t PPY;                                                                       
  Double_t PPZ;                                                                       
  Double_t AK;                                                                        
  Double_t AKS;                                                                       
  Double_t X;                                                                         
  Double_t Y;                                                                         
  Double_t Z;                                                                         
  Double_t T;                                                                         
  Double_t RP;                                                                        
  Double_t RPS;                                                                       
  }FSIPRFCommon;                                                              
						                                                                                       
#define FSIPRF COMMON_BLOCK(FSIPRF,fsiprf)                                         
COMMON_BLOCK_DEF(FSIPRFCommon, FSIPRF);                                      


			                                                                                      
 
}
#endif                     
