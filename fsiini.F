          SUBROUTINE FSIINI
C---Note:
C-- ICH= 0 (1) if the Coulomb interaction is absent (present);
C-- ISPIN= JJ= 1,2,..,MSPIN denote increasing values of the pair
C-- total spin S.
C-- To calculate the CF of two particles (with masses m1, m2 and
C-- charges C1, C2) the following information is required:
C-- AM= twice the reduced mass= 2*m1*m2/(m1+m2) in GeV/c^2,
C-- DM= (m1-m2)/(m1+m2), required if NS=2;
C-- SM= (m1+m2);
C-- AC= Bohr radius= 2*137.036*0.1973/(C1*C2*AMH) in fm;
C-- AC > 1.D9 if C1*C2= 0, AC < 0 if C1*C2 < 0;
C-- MSPIN= MSPINH(LL)= number of the values of the total pair spin S;
C-- FD= FDH(LL,JJ), RD= RDH(LL,JJ)= scattering length and effective
C-- radius for each value of the total pair spin S, JJ= 1,..,MSPIN;     ;
C-- the corresponding square well parameters EB= EBH(LL,JJ), RB=
C-- RBH(LL,JJ) (required if NS=1) may be calculated by sear.f;
C-- if the effective range approximation is not valid (as is the case,
C-- e.g., for two-pion system) a code for calculation of the scattering
C-- amplitude should be supplemented;
C-- RHO= RHOH(LL,JJ), SF= SFH(LL,JJ), SE= SEH(LL) are spin factors;
C-- RHO= the probability that the spins j1 and j2 of the two particles
C-- will combine in a total spin S;
C-- RHO= (2*S+1)/[(2j1+1)*(2j2+1)] for unpolarized particles;
C-- RHO= (1-P1*P2)/4 and (3+P1*P2)/4 correspond to S=0 and 1 in the
C-- case of spin-1/2 particles with polarizations P1 and P2;
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSIPOC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSIPRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSINS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSICONS/PI,PI2,SPI,DR,W
      COMPLEX*16 C
      COMMON/FSI_AA/AA,DW_INNER
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
     
c-ml     
      DIMENSION AAPIXI(20,2)
     
      DIMENSION FDH(60,10),RDH(60,10),EBH(60,10),RBH(60,10)
      DIMENSION RHOH(60,10)
      DIMENSION AM1H(60),AM2H(60),C1H(60),C2H(60),MSPINH(60)

C============= declarations pour l'appel de READ_FILE()============
      CHARACTER*10 KEY
      CHARACTER*8  CH8
      INTEGER*4    INT4
      REAL*8       REAL8
      INTEGER*4    IERR
C
C--- mass of the first and second particle
      DATA AM1H/.93956563D0,.93827231D0,.93956563D0,3.72737978D0,
     C          .13957D0,.13498D0,.13957D0, .93956563D0, .93827231D0,
     C          4*.13957D0,4*.493677D0,
     C          2*1.87561339D0,2*2.80892165D0,2*.497672D0,
     C          1.87561339D0,3*.93827231D0,.93956563D0,
     C          1.115684D0,.93827231D0,
     C          .13957D0,.13957D0,28*.0D0/
      DATA AM2H/.93956563D0,.93827231D0,.93827231D0,3.72737978D0,
     C          .13957D0,.13498D0,.13957D0, 2*1.87561339D0,
     C          2*.493677D0,2*.93827231D0,
     C          2*.493677D0,2*.93827231D0,
     C          1.87561339D0,3.72737978D0,2.80892165D0,3.72737978D0,
     C          2*.497672D0,2*2.80892165D0,3.72737978D0,
     C          3*1.115684D0,.93827231D0,
     C          1.32131D0,1.32131D0,28*.0D0/
c--------|---------|---------|---------|---------|---------|---------|----------
C---  charge of the first and second particle
      DATA C1H/0.D0,1.D0,0.D0,2.D0, 1.D0,0.D0,1.D0,0.D0,1.D0,
     C         3*1.D0,-1.D0,3*1.D0,-1.D0,
     C         4*1.D0,2*0.D0,4*1.D0,2*0.D0, 1.D0,
     C         1.D0,-1.D0,28*.0D0/
      DATA C2H/0.D0,1.D0,1.D0,2.D0,-1.D0,0.D0,3*1.D0,
     C         -1.D0,3*1.D0,-1.D0,3*1.D0,
     C         1.D0,2.D0,1.D0,2.D0,2*0.D0,2*1.D0,2.D0,3*0.D0,-1.D0,
     C         -1.D0,-1.D0,28*.0D0/
C---MSPIN vs (LL)
      DATA MSPINH/3*2,4*1,2*2,8*1,3,1,2,1,2*1,2*2,1,3*2, 2,
     C            2*1,28*.0D0/
C---Spin factors RHO vs (LL,ISPIN)
      DATA RHOH/3*.25D0, 4*1.D0, 2*.3333D0, 8*1.D0,
     1          .1111D0,1.D0,.25D0,1.D0,2*1.D0,
     1          .3333D0,.25D0,1.D0,3*.25D0, .25D0,
     1          2*1.D0,28*0.D0,
     2          3*.75D0, 4*0.D0, 2*.6667D0, 8*0.D0,
     2          .3333D0,.0D0,.75D0,.0D0,2*0.D0,
     2          .6667D0,.75D0,0.D0,3*.75D0, .75D0,
     2          2*0.D0,28*0.D0,
     3          17*.0D0,.5556D0,3*0.D0, 8*0.D0,1*0.D0,
     3          2*0.D0,28*0.D0, 420*0.D0/
C---Scattering length FD and effective radius RD in fm vs (LL,ISPIN)
      DATA FDH/17.0D0,7.8D0,23.7D0,2230.1218D0,
     1  .225D0,.081D0,-.063D0,
     1  -.65D0,-2.73D0,
     1  .137D0,-.071D0,-.148D0,.112D0,2*1.D-6,
     1  -.330D0,   
     1  2*1.D-6,1.344D0,6*1.D-6,-5.628D0,2.18D0,2.40D0,
     1  2.81D0, 
     1  1*0.001D0,
     1  2*.1D-6,28*0.D0,
     2  3*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2  1.93D0,1.84D0,
     2  0.50D0, 
     2  1*0.001D0,
     2  2*0.D0,28*0.D0,
     3  480*0.D0/
c--------|---------|---------|---------|---------|---------|---------|----------
      DATA RDH/2.7D0,2.8D0,2.7D0,1.12139906D0,-44.36D0,
     1  64.0D0,784.9D0,
     1  477.9D0, 2.27D0, 9*0.D0,-69.973D0, 6*0.D0,3.529D0,
     1  3.19D0,3.15D0,
     1  2.95D0, 
     1  1*0.D0,
     1  2*0.D0,28*0.D0,
     2  3*1.7D0,4*0.D0,2.0D0,2.63D0, 17*0.D0,3.35D0,3.37D0,
     2  2.95D0, 
     2  1*0.D0,
     2  2*0.D0,28*0.D0,
     3  480*0.D0/
     
C---Corresponding square well parameters RB (width in fm) and
C-- EB =SQRT(-AM*U) (in GeV/c); U is the well height
    
      DATA RBH/2.545739D0,   2.779789D0, 2.585795D0, 
     1  5.023544D0,
     1 .124673D0, .3925180D0,.09D0, 2.D0,
     1  4.058058D0, 17*0.D0,
     1  2.252623D0, 2.278575D0,
     1  2.234089D0, 
     1  3*0.001D0,28*0.D0,
     2  3*2.003144D0,
     2  4*0.D0, 2.D0, 4.132163D0, 17*0.D0,
     2  2.272703D0, 2.256355D0,
     2  2.234089D0, 
     2  3*0.001D0,28*0.D0,
     3  480*0.D0/
     
      DATA EBH/.1149517D0,.1046257D0,.1148757D0,.1186010D0,
     1    .7947389D0,2.281208D0,8.7D0,.4D0,
     1    .1561219D0,17*0.D0,
     1    .1013293D0, .1020966D0,
     1    .1080476D0, 
     1    3*0.001D0,28*0.D0,
     2    3*.1847221D0,
     2    4*0.D0, .4D0, .1150687D0, 17*0.D0,
     2    .09736083D0, .09708310D0,
     2    .1080476D0, 
     2    3*0.001D0,28*0.D0,
     3    480*0.D0/

C---Parameters for GPIXI (I,J), J=1,2 -> piXi isospin=1/2,3/2
C--- s-wave scatt. length, eff. radius in 1/GeV
C--  for isospin 1/2 also p-wave Xi* mass, width, decay mom. in GeV
c      DATA AAPIXI/ .1D-6, .0D0, 1.5318D0, 9.1D-3, .14656D0, 15*0.D0,
c     1             .1D-6, .0D0,     18*0.D0/



C=======< constants >========================

      W=1/.1973D0    ! from fm to 1/GeV
      PI=4*DATAN(1.D0)
      PI2=2*PI
      SPI=DSQRT(PI)
      DR=180.D0/PI   ! from radian to degree
        AC1=1.D10
        AC2=1.D10

      CALL LEDBLDATA


C=======< condition de calculs >=============
      NUNIT=11 ! for IBM or HP
C      NUNIT=4 ! for SUN in Prague
c-ml      CALL readint4(NUNIT,'ITEST     ',ITEST)
c-ml      CALL readint4(NUNIT,'LL        ',LL)        ! Two-particle system
c-ml      CALL readint4(NUNIT,'NS        ',NS)
c      CALL READ_FILE(NUNIT,'ITEST     ',CHAR,ITEST,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'LL        ',CHAR,LL,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'NS        ',CHAR,NS,REAL8,IERR)

C---setting particle masses and charges
      AM1=AM1H(LL)
      AM2=AM2H(LL)
      C1=C1H(LL)
      C2=C2H(LL)
      AM1S=AM1*AM1
      AM2S=AM2*AM2

C-- Switches:
C   ISI=1(0)  the strong interaction between the two particles ON (OFF)
C   IQS=1(0)  the quantum statistics ON (OFF);
C             should be OFF for nonidentical particles
C   I3C=1(0)  the Coulomb interaction with the nucleus ON (OFF)
C   I3S=1(0)  the strong interaction with the nucleus ON (OFF)
C   ICH=1(0)  if C1*C2 is different from 0 (is equal to 0)
C-  to switch off the Coulomb force between the two particles
C   put ICH=0 and substitute the strong amplitude parameters by
C   the ones not affected by Coulomb interaction

c-ml       ICH=0
c-ml       IF(C1*C2.NE.0.D0) ICH=1
c-ml       IQS=0
c-ml       IF(C1+AM1.EQ.C2+AM2) IQS=1
c-ml       I3S=0  ! only this option is available
c-ml       ISI=1
c-ml       I3C=0
c-ml       IF(CN*ICH.NE.0.D0) I3C=1


       write(*,*)'in fsiini AM1 AM2 C1 C2 ICH IQS I3C ISI'
       write(*,*)AM1,AM2,C1,C2,ICH,IQS,I3C,ISI


c-ml      IF(ITEST.EQ.1)THEN
c-ml      CALL readint4(NUNIT,'ICH       ',ICH)
c-ml      CALL readint4(NUNIT,'IQS       ',IQS)
c-ml      CALL readint4(NUNIT,'ISI       ',ISI)
c-ml      CALL readint4(NUNIT,'I3C       ',I3C)
c      CALL READ_FILE(NUNIT,'ICH     ',CHAR,ICH,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'IQS     ',CHAR,IQS,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'ISI     ',CHAR,ISI,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'I3C     ',CHAR,I3C,REAL8,IERR)
c-ml      ENDIF
C==================================================================
C---fm to 1/GeV
      DO 3 J1=1,60
      DO 3 J2=1,10
      FDH(J1,J2)=FDH(J1,J2)*W
      RDH(J1,J2)=RDH(J1,J2)*W
      RBH(J1,J2)=RBH(J1,J2)*W

      write(*,*)FDH(12,J2),RDH(12,j2),RBH(12,j2)

 3    continue
 
 
 
C---calcul. twice the reduced mass (AM), the relative mass difference
C-- (DM) and the Bohr radius (AC)
      AM=2*AM1*AM2/(AM1+AM2)
      AMS=AM*AM
      SM=AM1+AM2
      DM=(AM1-AM2)/SM
      AC=1.D10
      C12=C1*C2
      IF(C12.NE.0.D0)AC=2*137.036D0/(C12*AM) ! Bohr radius AC = a
C                                              a < 0 for opposite charges


C---Setting spin factors
      MSPIN=MSPINH(LL)

c-ml
      write(*,*)'AC=',AC,'MSPIN=',MSPIN

      
      MSP=MSPIN
      DO 91 ISPIN=1,10
  91  RHO(ISPIN)=RHOH(LL,ISPIN)
C---Integration limit AA in the spherical wave approximation
      AA=0.D0
cml
c      IF(NS.EQ.2.OR.NS.EQ.4)AA=.5D0 !!in 1/GeV --> 0.1 fm
c       IF(NS.EQ.2.OR.NS.EQ.4)AA=6.D0 !!in 1/GeV --> 1.2 fm
c      IF(NS.EQ.2.OR.NS.EQ.4)AA=.1D0 !!in 1/GeV --> .02 fm
      IF(NS.EQ.2.OR.NS.EQ.4)AA=5.D0 !!in 1/GeV --> 1.0 fm
C---Setting scatt. length (FD), eff. radius (RD) and, if possible,
C-- also the corresp. square well parameters (EB, RB)
      DO 55 JJ=1,MSP
      ISPIN=JJ
      FD(JJ)=FDH(LL,JJ)
      RD(JJ)=RDH(LL,JJ)
      EB(JJ)=EBH(LL,JJ)
      RB(JJ)=RBH(LL,JJ)
C---Resets FD and RD for a nucleon-deuteron system (LL=8,9)
      IF(LL.EQ.8.OR.LL.EQ.9)THEN
       JH=LL-7+2*JJ-2
       FD(JJ)=AAND(1,JH)
       RD(JJ)=AAND(2,JH)-2*AAND(3,JH)/AAND(1,JH)
      ENDIF
C---Resets FD and RD for a pion-pion system (LL=5,6,7)
       IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)THEN
       IF(LL.EQ.7)FD(JJ)=AAPI(1,2)/AM
       IF(LL.EQ.5)FD(JJ)=(.6667D0*AAPI(1,1)+.3333D0*AAPI(1,2))/AM
       IF(LL.EQ.6)FD(JJ)=(.3333D0*AAPI(1,1)+.6667D0*AAPI(1,2))/AM
       AKS=0.D0
       DAKS=1.D-5
       AKSH=AKS+DAKS
       AKH=DSQRT(AKSH)
       G1H=GPIPI(AKSH,1)
       G2H=GPIPI(AKSH,2)
       H=1/FD(JJ)

       write(*,*)'H', H

       IF(LL.EQ.7)C(JJ)=1/DCMPLX(G2H,-AKH)
       IF(LL.EQ.5)
     + C(JJ)=.6667D0/DCMPLX(G1H,-AKH)+.3333D0/DCMPLX(G2H,-AKH)
       IF(LL.EQ.6)
     + C(JJ)=.3333D0/DCMPLX(G1H,-AKH)+.6667D0/DCMPLX(G2H,-AKH)
       HH=DREAL(1/C(JJ))
       RD(JJ)=2*(HH-H)/DAKS
      ENDIF
C---Resets FD and RD for a pion-nucleon system (LL=12,13)
       IF(LL.EQ.12.OR.LL.EQ.13)THEN
       IF(LL.EQ.12)FD(JJ)=AAPIN(1,2)
       IF(LL.EQ.13)FD(JJ)=(.6667D0*AAPIN(1,1)+.3333D0*AAPIN(1,2))


       AKS=0.D0
       DAKS=1.D-5
       AKSH=AKS+DAKS
       AKH=DSQRT(AKSH)
       G1H=GPIN(AKSH,1)
       G2H=GPIN(AKSH,2)
       H=1/FD(JJ)              
       IF(LL.EQ.12)C(JJ)=1/DCMPLX(G2H,-AKH)
       IF(LL.EQ.13)
     + C(JJ)=.6667D0/DCMPLX(G1H,-AKH)+.3333D0/DCMPLX(G2H,-AKH)
       HH=DREAL(1/C(JJ))
       RD(JJ)=2*(HH-H)/DAKS
      ENDIF
C---Resets FD and RD for a pion-Xi system (LL=31,32)
      IF(LL.EQ.31.OR.LL.EQ.32)THEN
       IF(LL.EQ.32)FD(JJ)=AAPIXI(1,2)
       IF(LL.EQ.31)FD(JJ)=(.6667D0*AAPIXI(1,1)+.3333D0*AAPIXI(1,2))
       AKS=0.D0
       DAKS=1.D-5
       AKSH=AKS+DAKS
       AKH=DSQRT(AKSH)
       G1H=GPIXI(AKSH,1)
       G2H=GPIXI(AKSH,2)
       H=1/FD(JJ)
       IF(LL.EQ.32)C(JJ)=1/DCMPLX(G2H,-AKH)
       IF(LL.EQ.31)
     + C(JJ)=.6667D0/DCMPLX(G1H,-AKH)+.3333D0/DCMPLX(G2H,-AKH)
       HH=DREAL(1/C(JJ))
       RD(JJ)=2*(HH-H)/DAKS
      ENDIF
C---fm to 1/GeV for pp-bar system
      IF(LL.EQ.30)THEN
      DO 4 I3=1,3
      AAPAPR(I3,JJ)=AAPAPR(I3,JJ)*W
 4    AAPAPI(I3,JJ)=AAPAPI(I3,JJ)*W
C---Calculates complex elements M11=M22=C(6), M12=M21=C(7) for I=0
C---   at k*=0                  M11=M22=C(8), M12=M21=C(9) for I=1
        C(7+(JJ-1)*2)=2*DCMPLX(AAPAPR(1,JJ),AAPAPI(1,JJ))*
     *                DCMPLX(AAPAPR(2,JJ),AAPAPI(2,JJ)) ! 2a_0Sa_1S
        C(6+(JJ-1)*2)=DCMPLX(AAPAPR(1,JJ)+AAPAPR(2,JJ),
     ,                     AAPAPI(1,JJ)+AAPAPI(2,JJ))/
     /                                  C(7+(JJ-1)*2)   ! M11=M22
      C(7+(JJ-1)*2)=-DCMPLX(AAPAPR(1,JJ)-AAPAPR(2,JJ),
     ,                      AAPAPI(1,JJ)-AAPAPI(2,JJ))/
     /                                   C(7+(JJ-1)*2)  ! M12=M21
      ENDIF
    
              
c      write(*,*)'C(jj)',C(jj)
c      write(*,*)'C(jjj)',(C(jjj),jjj=1,10)
c      RD(JJ)=0
c       write(*,*)'FD RD',FD(JJ),RD(JJ)
                    
      
C---Calculation continues for any system (any LL)
 55   CONTINUE
      RETURN
      END



      subroutine LEDBLDATA
    
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AADD/AADD(20,3)/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_AAPIXI/AAPIXI(20,2)
      COMMON/FSI_AAKK/AAKK(9)/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
C---Parameters for GPIPI (I,J), J=1,2 -> isospin=0,2
C--- Nagels (79)
c      DATA AAPI/.2600D00, .2500D00, .3480D00,-.0637D00, 16*0.D0,
c     1         -.0280D00,-.0820D00, .2795D00,-.0086D00, 16*0.D0/
C--- Colangelo(01)
      DATA AAPI/.220D0, .268D0,-.139D-1,-.139D-2, 36.77D0, 15*0.D0,
     1         -.444D-1,-.857D-1,-.221D-2,-.129D-3,-21.62D0, 15*0.D0/
C***
C---Parameters for GPIN (I,J), J=1,2 -> piN isospin=1/2,3/2
C-- s-wave scattering length and effective radius in 1/GeV
      DATA AAPIN/ .12265D1, .1563D2,     18*0.D0,
     1            -.750D0,  .7688D2,     18*0.D0/
C***
C---Parameters for GPIXI (I,J), J=1,2 -> piXi isospin=1/2,3/2
C--- s-wave scatt. length, eff. radius in 1/GeV
C--  for isospin 1/2 also p-wave Xi* mass, width, decay mom. in GeV
      DATA AAPIXI/ .1D-6, .0D0, 1.5318D0, 9.1D-3, .14656D0, 15*0.D0,
     1             .1D-6, .0D0,     18*0.D0/

c-non zero scat lengths      DATA AAPIXI/ 0.5, .0D0, 1.5318D0, 9.1D-3, .14656D0, 15*0.D0,
c     1             0.5, .0D0,     18*0.D0/

c-non zero scat lengths      
c        DATA AAPIXI/ 2.5338, .0D0, 1.5318D0, 9.1D-3, .14656D0, 15*0.D0,
c     1             2.5338, .0D0,     18*0.D0/ ! 0.5 fm ~ 2.5338 1/GeV



C***
C---Parameters for GND (I,J), J=1(2) corresp. to nd(pd), S=1/2,
C---                          J=3(4) corresp. to nd(pd), S=3/2
      DATA AAND/-.3295D1,-.8819D3, .4537D4, .28645D5, 16*0.D0,
     1          -.13837D2, .11505D2, .0D0, .10416D2,  16*0.D0,
     2          -.32180D2, .1014D2,  .0D0, .0D0,      16*0.D0,
     3          -.60213D2, .1333D2,  .0D0,-.70309D2,  16*0.D0/
      DATA AADD/ .10617D4, .3194D-2, .56849D3, 17*0.D0,
     1           20*0.D0,
     2          -.1085D4, .1987D5, 18*0.D0/
C***
C--- AAKK= m_K^2, m_pi^2, m_eta^2,
C---       m_S*^2, m_delta^2, gam(S*-->KK-b),
C---       gam(S*-->pipi), gam(delta-->KK-b), gam(delta-->pi eta)
      DATA AAKK/.247677D0,.01947977D0,.2997015D0,
     1 .9565D0,.9487D0,.792D0,.199D0,.333D00,.222D0/ ! Martin'77
c     1 .9467D0,.9698D0,2.763D0,.5283D0,.4038D0,.3711D0/ ! Antonelli'02
c     1 .9920D0,.9841D0,1.305D0,.2684D0,.5555D0,.4401D0/ ! Achasov1'01,03
c     1 .9920D0,1.0060D0,1.305D0,.2684D0,.8365D0,.4580D0/! Achasov2'01,03
C***
C---Parameters for PAP (I,J), j=1,2 -> isospin I=0,2
C---                          i=1-3 -> a_singlet, a_triplet, d [fm]
C---    Im a_IS (I=isospin, S=spin) are fixed by atomic data and
C       n-bar survival time up to one free parameter, e.g. Im a_00
C---    Batty (89), Kerbikov (93):
C--- Ima_10=1.96-Ima_00, Ima_01=0.367-Ima_00/3, Ima_11=0.453+Ima_00/3
C---       In DATA we put Ima_00=0.3.
C---    Re a_IS are fixed by atomic data up to three free parameters
C---    Batty (89):
C---         Rea_aver(pp-bar)=Re[(a_00+a_10)+3(a_01+a_11)]/8=-0.9
C---       In DATA we used Rea_IS from Paris potential Pignone (94)
C---       rescaled by 1.67 to satisfy the atomic constraint.
C---    Effective radius is taken independent of IS from the phase
C---    shift fit by Pirner et al. (91).
      DATA AAPAPR/-0.94D0, -1.98D0,  .1D0,
     1            -1.40D0,  0.37D0,  .1D0/ ! Re
      DATA AAPAPI/ 0.3 D0,  .267D0,-.01D0,
     1             1.66D0,  .553D0,-.01D0/ ! Im
      END
