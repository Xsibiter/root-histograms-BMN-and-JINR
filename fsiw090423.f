c         1         2         3         4         5         6         7         8
c---------|---------|---------|---------|---------|---------|---------|---------|
*-- Author :    R.Lednicky      20/01/95
      SUBROUTINE FSIW(J,WEIF,WEI,WEIN)

C=======================================================================
C    Calculates final state interaction (FSI) weights
C    WEIF = weight due to particle - (effective) nucleus FSI (p-N)
C    WEI  = weight due to p-p-N FSI
C    WEIN = weight due to p-p FSI; note that WEIN=WEI if I3C=0;
C                                  note that if I3C=1 the calculation of
C                                  WEIN can be skipped by putting J=0
C.......................................................................
C    Correlation Functions:
C    CF(p-p-N)   = sum(WEI)/sum(WEIF)
C    CF(p-p)     = sum(WEIN)/sum(1); here the nucleus is completely
C                                    inactive
C    CF(p-p-"N") = sum(WEIN*WEIF')/sum(WEIF'), where WEIN and WEIF'
C                  are not correlated (calculated at different emission
C                  points, e.g., for different events);
C                  thus here the nucleus affects one-particle
C                  spectra but not the correlation
C.......................................................................
C    User must supply data file <fn> on unit NUNIT (e.g. =11) specifying
C    LL   : particle pair
C    NS   : approximation used to calculate Bethe-Salpeter amplitude
C    ITEST: test switch
C           If ITEST=1 then also following parameters are required
C    ICH  : 1(0) Coulomb interaction between the two particles ON (OFF)
C    IQS  : 1(0) quantum statistics for the two particles ON (OFF)
C    ISI  : 1(0) strong interaction between the two particles ON (OFF)
C    I3C  : 1(0) Coulomb interaction with residual nucleus ON (OFF)
C    This data file can contain other information useful for the user.
C    It is read by subroutines READINT4 and READREA8(4) (or READ_FILE).
C----------------------------------------------------------------------
C-   LL       1  2  3  4  5   6   7   8  9 10  11  12  13  14 15 16 17
C-   part. 1: n  p  n  a  pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
C-   part. 2: n  p  p  a  pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p
C   NS=1 y/n: +  +  +  +  +   -   -   -  -  -   -   -   -  -  -  -  -
C----------------------------------------------------------------------
C-   LL       18 19 20 21 22 23  24 25 26 27 28 29 30 31  32
C-   part. 1: d  d  t  t  K0 K0  d  p  p  p  n  /\ p  pi+ pi-
C-   part. 2: d  a  t  a  K0 K0b t  t  a  /\ /\ /\ pb Xi- Xi-
C   NS=1 y/n: -  -  -  -  -  -   -  -  -  +  +  +  -  -   -
C----------------------------------------------------------------------
C   NS=1  Square well potential,
C   NS=3  not used
C   NS=4  scattered wave approximated by the spherical wave,
C   NS=2  same as NS=4 but the approx. of equal emission times in PRF
C         not required (t=0 approx. used in all other cases).
C   Note: if NS=2,4, the B-S amplitude diverges at zero distance r* in
C         the two-particle c.m.s.; user can specify a cutoff AA in
C         SUBROUTINE FSIINI, for example:
C         IF(NS.EQ.2.OR.NS.EQ.4)AA=5.D0 !! in 1/GeV --> AA=1. fm
C---------------------------------------------------------------------
C    ITEST=1 any values of parameters ICH, IQS, ISI, I3C are allowed
C            and should be given in data file <fn>
C    ITEST=0 physical values of these parameters are put automatically
C            in FSIINI (their values are not required in data file)
C=====================================================================
C    At the beginning of calculation user should call FSIINI,
C    which reads LL, NS, ITEST (and eventually ICH, IQS, ISI, I3C)
C    and initializes various parameters.
C    In particular the constants in
C      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
C    may be useful for the user:
C     W=1/.1973D0    ! from fm to 1/GeV
C     PI=4*DATAN(1.D0)
C     PI2=2*PI
C     SPI=DSQRT(PI)
C     DR=180.D0/PI   ! from radian to degree
C      _______________________________________________________
C  !! |Important note: all real quantities are assumed REAL*8 | !!
C      -------------------------------------------------------
C    For each event user should fill in the following information
C    in COMMONs (all COMMONs in FSI calculation start with FSI_):
C    ...................................................................
C     COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
C    Only
C         AMN  = mass of the effective nucleus   [GeV/c**2]
C         CN   = charge of the effective nucleus [elem. charge units]
C    are required
C    ...................................................................
C     COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1, !part. momenta in the rest frame
C    1               P2X,P2Y,P2Z,E2,P2  !of effective nucleus (NRF)
C    Only the components
C                        PiX,PiY,PiZ  [GeV/c]
C    in NRF are required.
C    To make the corresponding Lorentz transformation user can use the
C    subroutines LTRAN and LTRANB
C    ...................................................................
C    COMMON/FSI_COOR/X1,Y1,Z1,T1,R1,     ! 4-coord. of emission
C    1               X2,Y2,Z2,T2,R2      ! points in NRF
C    The componets
C                       Xi,Yi,Zi  [fm]
C    and emission times
C                          Ti   [fm/c]
C    should be given in NRF with the origin assumed at the center
C    of the effective nucleus. If the effect of residual nucleus is
C    not calculated within FSIW, the NRF can be any fixed frame.
C-----------------------------------------------------------------------
C    Before calling FSIW the user must call
C     CALL LTRAN12
C    Besides Lorentz transformation to pair rest frame:
C    (p1-p2)/2 --> k* it also transforms 4-coordinates of
C    emission points from fm to 1/GeV and calculates Ei,Pi and Ri.
C    Note that |k*|=AK in COMMON/FSI_PRF/
C-----------------------------------------------------------------------
C    After making some additional filtering using k* (say k* < k*max)
C    or direction of vector k*,
C    user can finally call FSIW to calculate the FSI weights
C    to be used to construct the correlation function
C=======================================================================

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_JR/JRAT
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  ! particle momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! k*=(p1-p2)/2 and x1-x2
     1               X,Y,Z,T,RP,RPS      ! in pair rest frame (PRF)
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, !4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF12,FF21
C------------------------------------------------------------------
C==> AC1,2 = "relativistic" Bohr radii for particle-nucleus systems
      C1N=C1*CN
      IF(C1N.NE.0.D0)AC1=137.036D0/(C1N*E1) !m1-->E1
      C2N=C2*CN
      IF(C2N.NE.0.D0)AC2=137.036D0/(C2N*E2) !m2-->E2

C-----------------------------------------------------------
      CALL FSIPN(WEIF) !weight due to particle-nucleus FSI
      JRAT=0
      CALL FSIWF(WEI)  !weight due to particle-particle-nucleus FSI
      WEIN=WEI
         IF(I3C*J.NE.0) THEN
      FF12=DCMPLX(1.D0,0.D0)
      FF21=DCMPLX(1.D0,0.D0)
      JRAT=1
      CALL VZ(WEIN) ! weight due to particle-particle FSI
         ENDIF
         RETURN
      END
C=======================================================================

      SUBROUTINE LTRAN(P0,P,PS)
C==>calculating particle 4-momentum PS={PSX,PSY,PSZ,ES}
C   in rest frame of a system 0 with 4-momentum P0={P0X,P0Y,P0Z,E0}
C   from its 4-momentum P={PX,PY,PZ,E}

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P0(4),P(4),PS(4)
C-----------------------------------------------------------------------
      P0S=P0(1)**2+P0(2)**2+P0(3)**2
      AM0=DSQRT(P0(4)**2-P0S)
      EPM=P0(4)+AM0
      PP0=P(1)*P0(1)+P(2)*P0(2)+P(3)*P0(3)
      H=(PP0/EPM-P(4))/AM0
      PS(1)=P(1)+P0(1)*H
      PS(2)=P(2)+P0(2)*H
      PS(3)=P(3)+P0(3)*H
      PS(4)=(P0(4)*P(4)-PP0)/AM0
      RETURN
      END

      SUBROUTINE LTRANB(P0,PS,P)
C==>calculating particle 4-momentum P={PX,PY,PZ,E}
C   from its 4-momentum PS={PSX,PSY,PSZ,ES}
C   in rest frame of a system 0 with 4-momentum P0={P0X,P0Y,P0Z,E0}

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P0(4),P(4),PS(4)
C-----------------------------------------------------------------------
      P0S=P0(1)**2+P0(2)**2+P0(3)**2
      AM0=DSQRT(P0(4)**2-P0S)
      EPM=P0(4)+AM0
      PSP0=PS(1)*P0(1)+PS(2)*P0(2)+PS(3)*P0(3)
      HS=(PSP0/EPM+PS(4))/AM0
      P(1)=PS(1)+P0(1)*HS
      P(2)=PS(2)+P0(2)*HS
      P(3)=PS(3)+P0(3)*HS
      P(4)=(P0(4)*PS(4)+PSP0)/AM0
      RETURN
      END

      SUBROUTINE LTRAN12
C==>calculating particle momentum in PRF {EE,PPX,PPY,PPZ} from
C-  the momentum of the first particle {E1,P1X,P1Y,P1Z) in NRF
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, !4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W

C   fm --> 1/GeV
      X1=X1*W
      Y1=Y1*W
      Z1=Z1*W
      T1=T1*W
      X2=X2*W
      Y2=Y2*W
      Z2=Z2*W
      T2=T2*W
C   calculating Ri, Pi and Ei
      R1=DSQRT(X1*X1+Y1*Y1+Z1*Z1)
      R2=DSQRT(X2*X2+Y2*Y2+Z2*Z2)
      P1S=P1X*P1X+P1Y*P1Y+P1Z*P1Z
      P2S=P2X*P2X+P2Y*P2Y+P2Z*P2Z
      P1=DSQRT(P1S)
      P2=DSQRT(P2S)
      E1=DSQRT(AM1*AM1+P1S)
      E2=DSQRT(AM2*AM2+P2S)
C-----------------------------------------------------------------------
      E12=E1+E2
      P12X=P1X+P2X
      P12Y=P1Y+P2Y
      P12Z=P1Z+P2Z
      P12S=P12X**2+P12Y**2+P12Z**2
      AM12=DSQRT(E12**2-P12S)
      EPM=E12+AM12
      P12=DSQRT(P12S)
      P112=P1X*P12X+P1Y*P12Y+P1Z*P12Z
      H1=(P112/EPM-E1)/AM12
      PPX=P1X+P12X*H1
      PPY=P1Y+P12Y*H1
      PPZ=P1Z+P12Z*H1
      EE=(E12*E1-P112)/AM12
      AKS=EE**2-AM1**2
      AK=DSQRT(AKS)
CW      WRITE(6,38)'AK ',AK,'K ',PPX,PPY,PPZ,EE
38    FORMAT(A7,E11.4,A7,4E11.4)
      RETURN
      END

      SUBROUTINE FSIPN(WEIF)
C  calculating particle-nucleus Coulomb Wave functions FFij
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_ICH1/ICH1
      COMMON/FSI_ETA/ETA
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF1,FF12,FF21
      FF12=DCMPLX(1.D0,0.D0)
      FF21=DCMPLX(1.D0,0.D0)
      ACH=1.D0
      WEIF=1.D0
      IF(I3C.EQ.0)RETURN
      ICH1=IDINT(C1)
      IF(ICH1.EQ.0)GOTO 11
      XH=AC1*P1
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
      RHOS=P1*R1
      HS=X1*P1X+Y1*P1Y+Z1*P1Z
      FF12=FF12*FF1(RHOS,HS)
      IF(IQS.EQ.0)GOTO 11
      RHOS=P1*R2
      HS=X2*P1X+Y2*P1Y+Z2*P1Z
      FF21=FF21*FF1(RHOS,HS)
 11   ICH1=IDINT(C2)
      IF(ICH1.EQ.0)GOTO 10
      XH=AC2*P2
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
      RHOS=P2*R2
      HS=X2*P2X+Y2*P2Y+Z2*P2Z
      FF12=FF12*FF1(RHOS,HS)
CW      WRITE(6,41)'AC2 ',AC2,'ACH ',ACH,'ETA ',ETA,'RHOS ',RHOS,'HS ',HS
41    FORMAT(5(A5,E11.4))
CW      WRITE(6,40)'FF12 ',DREAL(FF12),DIMAG(FF12)
      IF(IQS.EQ.0)GOTO 10
      RHOS=P2*R1
      HS=X1*P2X+Y1*P2Y+Z1*P2Z
      FF21=FF21*FF1(RHOS,HS)
CW      WRITE(6,41)'AC1 ',AC1,'ACH ',ACH,'ETA ',ETA,'RHOS ',RHOS,'HS ',HS
CW      WRITE(6,40)'FF21 ',DREAL(FF21),DIMAG(FF21)
40    FORMAT(A7,2E12.4)
 10   CONTINUE

C  WEIF = the weight due to the Coulomb particle-nucleus interaction
      WEIF=DREAL(FF12)**2+DIMAG(FF12)**2
      IF(IQS.EQ.1)WEIF=0.5D0*(WEIF+DREAL(FF21)**2+DIMAG(FF21)**2)
      RETURN
      END

      FUNCTION GPIPI(X,J)
C--- GPIPI = k*COTG(DELTA), X=k^2
C--  J=1(2) corresponds to isospin=0(2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPI/AAPI(20,2)
      COMMON/FSI_C/HELP(20)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      OMS=X+AMS
      OM=DSQRT(OMS)
      XX=X/AMS
      XXS=XX*XX
c      GOTO 1 ! Nagels (79)
C--- Colangelo(01)
      GPIPI=OM*(4*OMS-AAPI(5,J))/(4*AMS-AAPI(5,J))
      GPIPI=GPIPI/(AAPI(1,J)+AAPI(2,J)*XX+AAPI(3,J)*XXS+
     +             AAPI(4,J)*XXS*XX)
      RETURN
C--- Nagels (79)
 1    continue
      GPIPI=OM/AAPI(1,J)
      GPIPI=GPIPI*(1+(AAPI(3,J)-AAPI(1,J)**2)*XX+AAPI(4,J)*XXS)
      GPIPI=GPIPI/(1+(AAPI(3,J)+AAPI(2,J)/AAPI(1,J))*XX)
      RETURN
      END

      FUNCTION GPIN(X,J)
C--- GPIN = k*COTG(DELTA), X=k^2
C--  J=1(2) corresponds to piN isospin=1/2(3/2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      GPIN=1/AAPIN(1,J)+.5D0*AAPIN(2,J)*X
      RETURN
      END

      FUNCTION GPIXI(X,J)
C--- GPIXI = k*COTG(DELTA), X=k^2
C--  J=1(2) corresponds to piXi isospin=1/2(3/2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPIXI/AAPIXI(20,2)
      GPIXI=1/AAPIXI(1,J)+.5D0*AAPIXI(2,J)*X
      RETURN
      END

      FUNCTION GND(X,J)
C--- GND = k*COTG(DELTA), X=k^2
C--- J=1(2) corresp. to nd(pd), S=1/2,
C--- J=3(4) corresp. to nd(pd), S=3/2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAND/AAND(20,4)
      XX=X
      GND=1/AAND(1,J)+.5D0*AAND(2,J)*X
      DO 1 I=4,4
      XX=XX*X
   1  GND=GND+AAND(I,J)*XX
      GND=GND/(1+AAND(3,J)*X)
      RETURN
      END

      FUNCTION GDD(X,J)
C--- GDD = k*COTG(DELTA), X=k^2
C--- J=1,2,3 corresp. to S=0,1,2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AADD/AADD(20,3)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
      COMPLEX*16 C
      E=X/2/AM
      ER=DSQRT(E)
      IF(J.EQ.1)THEN
       GDD=ER*(AADD(1,1)*DEXP(-E/AADD(2,1))-AADD(3,1))
       GDD=GDD/DR   ! from degree to radian
       TAND=DTAN(GDD)
       IF(TAND.EQ.0.D0)TAND=1.D-10
       GDD=DSQRT(X)/TAND
      END IF
      IF(J.EQ.2)THEN
       GDD=1.D10
      END IF
      IF(J.EQ.3)THEN
       GDD=ER*(AADD(1,3)+AADD(2,3)*E)
       GDD=GDD/DR    ! from degree to radian
       TAND=DTAN(GDD)
       IF(TAND.EQ.0.D0)TAND=1.D-10
       GDD=DSQRT(X)/TAND
      END IF
      RETURN
      END

      BLOCK DATA
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

      SUBROUTINE CPIPI ! calculates pi+pi- s-wave scattering amplitudes
                       ! accounting for pi0pi0->pi+pi- transition
	                 ! C(1)= f_c(pi+pi-->pi+pi-)
                       ! C(2)= f_c(pi+pi-->pi0pi0)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      COMPLEX*16 C
      DATA AM2C/.13498D0/ ! pi0 mass
      AMU2_AMU1=AM2C/AM   ! AM2C=2*mu(pi0pi0), AM=2*mu(pi+pi-)
      AK2S0=AMU2_AMU1*AKS
      AK2S =AK2S0-2*AM2C*(AM2C-AM)
      AK2=DSQRT(AK2S)  !  k2
      G2=GPIPI(AKS,2)
      G1=GPIPI(AKS,1)
      RK11=.6667D0*G1+.3333D0*G2
      RK22=.6667D0*G2+.3333D0*G1
      RK12=-0.47140452D0*(G1-G2)
      C(3)=RK11+DCMPLX(-HCP2,-AAK)      ! (1/f)11
      C(5)=RK22
      C(5)=C(5)-DCMPLX(0.D0,AK2)        ! (1/f)22
      C(7)=RK12                         ! (1/f)12
      C(10)=C(3)*C(5)-C(7)*C(7)
      C(1)=C(5)/C(10)                   ! f11
      C(2)=-C(7)/C(10)                  ! f12
      RETURN
      END

      SUBROUTINE CPIN  ! calculates pi-p s-wave scattering amplitudes
                       ! accounting for pi0n->pi-p transition
	                 ! C(1)= f_c(pi-p->pi-p)
                       ! C(2)= f_c(pi-p->pi0n)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      COMPLEX*16 C
      DATA AM2C1/.13498D0/  ! pi0 mass
	DATA AM2C2/.93956563D0/ ! n mass
      DATA AM2C/.2360487D0/ ! twice the pi0n reduced mass
      AMU2_AMU1=AM2C/AM   ! AM2C=2*mu(pi0n), AM=2*mu(pi-p)
      AK2S0=AMU2_AMU1*AKS
      AK2S =AK2S0-AM2C*(AM2C1+AM2C2-SM) ! here AM2C1+AM2C2 < SM = AM1+AM2
      AK2=DSQRT(AK2S)  !  k2
      G2=GPIN(AKS,2)
      G1=GPIN(AKS,1)
      RK11=.6667D0*G1+.3333D0*G2
      RK22=.6667D0*G2+.3333D0*G1
      RK12=-0.47140452D0*(G1-G2)
      C(3)=RK11+DCMPLX(-HCP2,-AAK)      ! (1/f)11
      C(5)=RK22
      C(5)=C(5)-DCMPLX(0.D0,AK2)        ! (1/f)22
      C(7)=RK12                         ! (1/f)12
      C(10)=C(3)*C(5)-C(7)*C(7)
      C(1)=C(5)/C(10)                   ! f11
      C(2)=-C(7)/C(10)                  ! f12
      RETURN
      END

	SUBROUTINE CPIXI ! calculates pi+Xi- s-wave scattering amplitudes
	                 ! and p-wave Xi*(1532) resonance amplitudes
                       ! accounting for pi0Xi0->pi+Xi- transition
	                 ! C(1)= f_c(pi+Xi-->pi+Xi-)=f(l=0,J=1/2)11/Ac
                       ! C(2)= f_c(pi+Xi-->pi0Xi0)=f(l=0,J=1/2)12/Ac
                       ! C(3)= f(l=1,J=3/2)11 = (2/3) f(l=1,J=3/2,I=1/2)
                       ! C(4)= f(l=1,J=3/2)12 = -(sqrt(2)/3) f(l=1,J=3/2,I=1/2) 
	                 ! C(5)= f(l=1,J=3/2,I=1/2)
	                 ! calculates DW_INNER= inner contribution to the weight
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      COMMON/FSI_AAPIXI/AAPIXI(20,2)
      COMMON/FSI_AA/AA,DW_INNER
      COMPLEX*16 C
      DATA AM1_2/.13498D0/  ! pi0 mass
	DATA AM2_2/1.31483D0/ ! Xi0 mass
      DATA AM1S_2/.0182196D0/  ! pi0 mass squared
	DATA AM2S_2/1.7287779D0/ ! Xi0 mass squared
      DATA AMU2_2/.2448262D0/ ! twice the pi0Xi0 reduced mass
      OM1=DSQRT(AM1S+AKS)
      OM2=DSQRT(AM2S+AKS)
      SQ=OM1+OM2 ! SQ= M= eff. mass of channel 1=pi+Xi-
      S=SQ*SQ
C--      AMU2_AMU1=AMU2_2/AM   ! AMU2_2=2*mu(pi0Xi0), AM=2*mu(pi+Xi-)
C--      AK2S0=AMU2_AMU1*AKS
C--      AK2S =AK2S0-AM2C*(AM2C1+AM2C2-SM) ! here AM2C1+AM2C2 < SM = AM1+AM2
      AMU=OM1*OM2/(OM1+OM2)          ! "relativistic" reduced mass in channel 1
	OM1_2=(S+AM1S_2-AM2S_2)/(2*SQ)
      OM2_2=(S+AM2S_2-AM1S_2)/(2*SQ)
      AMU_2=OM1_2*OM2_2/(OM1_2+OM2_2)! "relativistic" reduced mass in channel 2
      AMU2_AMU1=AMU_2/AMU ! ratio of "relativistic" reduced masses
	AK2S=OM1_2*OM1_2-AM1S_2
      AK2=DSQRT(AK2S)  !  k2 = momentum in channel 2
      G2=GPIXI(AKS,2)
      G1=GPIXI(AKS,1)
      RK11=.6667D0*G1+.3333D0*G2
      RK22=.6667D0*G2+.3333D0*G1
      RK12=-0.47140452D0*(G1-G2)
      C(3)=RK11+DCMPLX(-HCP2,-AAK)      ! (1/f)11
      C(5)=RK22
      C(5)=C(5)-DCMPLX(0.D0,AK2)        ! (1/f)22
      C(7)=RK12                         ! (1/f)12
      C(10)=C(3)*C(5)-C(7)*C(7)
      C(1)=C(5)/C(10)                   ! f11 l=0
      C(2)=-C(7)/C(10)                  ! f12 l=0
	AM0S=AAPIXI(3,1)*AAPIXI(3,1)
      GAM0K=AAPIXI(4,1)*(AKS/AAPIXI(5,1)**3)*AM0S/SQ ! GAMMA*M0/k
C--                                      GAMMA=GAMMA0*(k/k0)^3*(M0/M)
C--                                      GAMMA0= resonace width= AAPIXI(4,1)
C--                                      M0 = resonace mass = AAPIXI(3,1)
C--                                      k = AK = channel 1 momentum
C--                                      k0 = k(M0) = AAPIXI(5,1)
      GAM0=AK*GAM0K                   ! GAMMA*M0     
      C(5)=GAM0K/DCMPLX(AM0S-S,-GAM0) ! f(l=1,J=3/2,I=1/2)  
	C(3)=0.66666667D0*C(5)          ! f(l=1,J=3/2)11= (2/3)f(l=1,J=3/2,I=1/2)
      C(4)=-0.47140452D0*C(5)         ! f(l=1,J=3/2)12= -(sqrt(2)/3) 
c                                                            f(l=1,J=3/2,I=1/2) 
C---  calculating DW_1 = resonance p-wave contribution to the weight for r < AA
	RHOA=AK*AA
      C_0=DCOS(RHOA)
	S_0=DSIN(RHOA)
      C_1=C_0/RHOA+S_0
	S_1=S_0/RHOA-C_0
	AC1S1=RHOA*(C_0*S_0+C_1*S_1)-2*C_0*S_1-C_1*S_0
	AC1C1=RHOA*(C_0*C_0+C_1*C_1)-3*C_0*C_1
      AS1S1=RHOA*(S_0*S_0+S_1*S_1)-3*S_0*S_1
C---  calculating C(6) = df(J=3/2,I=1/2)/dk^2
      DAKS=1.D-5
      AKSH=AKS+DAKS
	AKH=DSQRT(AKSH)
      OM1H=DSQRT(AM1S+AKSH)
      OM2H=DSQRT(AM2S+AKSH)
      SQH=OM1H+OM2H
      SH=SQH*SQH
      GAM0KH=AAPIXI(4,1)*(AKSH/AAPIXI(5,1)**3)*AM0S/SQH 
      GAM0H=AKH*GAM0KH
      C(10)=GAM0KH/DCMPLX(AM0S-SH,-GAM0H)
      C(6)=(C(10)-C(5))/DAKS
	F1R=DREAL(C(5))  ! Ref(J=3/2,I=1/2)
	F1I=DIMAG(C(5))  ! Imf(J=3/2,I=1/2)
	H=2*AK
      DW_INNER= (4/AA**3)*(DREAL(C(6))+(AC1S1*F1R-AS1S1*F1I)/AKS +
     +                    (AC1C1+AS1S1)*(F1R*F1R+F1I*F1I)/H + 
     +                    H*DIMAG(DCONJG(C(5))*C(6)))
      RETURN
      END

      SUBROUTINE CKKB  ! calculates KK-b s-wave scattering amplitude,
                       ! saturated by isospin-0 S*(980) 
					 !          and isospin-1 delta(982) resonances
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_AAKK/AAKK(9)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMPLEX*16 C
      S4=AKS+AAKK(1)
      S=4*S4
      AKPIPI=DSQRT(S4-AAKK(2))
      EETA2=(S+AAKK(3)-AAKK(2))**2/4/S
      AKPIETA=DSQRT(EETA2-AAKK(3))
      C(1)=AAKK(6)/2/DCMPLX(AAKK(4)-S,
     ,-AK*AAKK(6)-AKPIPI*AAKK(7))
      C(1)=C(1)+AAKK(8)/2/DCMPLX(AAKK(5)-S,
     ,-AK*AAKK(8)-AKPIETA*AAKK(9))
      RETURN
      END
	
      SUBROUTINE CPAP  ! calculates pp-bar s-wave scattering amplitudes
                       ! accounting for nn-bar->pp-bar channel 
					 ! C(1)= f_c(pp-bar->pp-bar,S=0)
                       ! C(2)= f_c(pp-bar->pp-bar,S=1)
                       ! C(3)= f_c(pp-bar->nn-bar,S=0)
                       ! C(4)= f_c(pp-bar->nn-bar,S=1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      COMPLEX*16 C
        DATA AM2C/.93956563D0/ ! neutron mass
      AMU2_AMU1=AM2C/AM   ! AM2C=2*mu(nn-bar), AM=2*mu(pp-bar)
      AK2S0=AMU2_AMU1*AKS
        AK2S =AK2S0-2*AM2C*(AM2C-AM)
      IF(AK2S.GE.0.D0)THEN
           AK2=DCMPLX(DSQRT(AK2S),0.D0)  !  k2
       ELSE
           AK2=DCMPLX(0.D0,DSQRT(-AK2S)) !  kappa2
      ENDIF
      C(10)=C(6+(ISPIN-1)*2)+
     + DCMPLX(AAPAPR(3,ISPIN)*AKS/2-0.016D0-HCP2,
     ,        AAPAPI(3,ISPIN)*AKS/2-AAK)         ! (1/f)11
      C(5)=C(6+(ISPIN-1)*2)+
     + DCMPLX(AAPAPR(3,ISPIN)*AK2S0/2,
     ,        AAPAPI(3,ISPIN)*AK2S0/2)
      IF(AK2S.GE.0.D0)THEN
           C(5)=C(5)-DCMPLX(0.D0,AK2)
       ELSE
           C(5)=C(5)+DCMPLX(AK2,0.D0)              ! (1/f)22
        ENDIF
      C(10)=C(10)*C(5)-C(7+(ISPIN-1)*2)*C(7+(ISPIN-1)*2)
      C(ISPIN)=C(5)/C(10)                ! f11
      C(ISPIN+2)=-C(7+(ISPIN-1)*2)/C(10) ! f12
      RETURN
      END

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
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
      COMPLEX*16 C
      COMMON/FSI_AA/AA,DW_INNER
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
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
      DATA FDH/17.0D0,7.8D0,23.7D0,2230.1218D0,.225D0,.081D0,-.063D0,
     1  -.65D0,-2.73D0,
     1  .137D0,-.071D0,-.148D0,.112D0,2*1.D-6,
c     1  -.360D0,   ! Martin'77 (K+p)
     1  -.330D0,   ! Martin'81 (K+p)
     1  2*1.D-6,1.344D0,6*1.D-6,-5.628D0,2.18D0,2.40D0,
     1  2.81D0, ! ND potential /\/\
C     1  0.50D0, ! NSC97e potential /\/\
     1  1*0.001D0,
	1  2*.1D-6,28*0.D0,
cc     2 -10.8D0,2*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2  3*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2  1.93D0,1.84D0,
     2  0.50D0, ! triplet f0 /\/\=singlet f0 ND
                ! not contributing in s-wave FSI approx.
     2  1*0.001D0,
	2  2*0.D0,28*0.D0,
     3  480*0.D0/
c--------|---------|---------|---------|---------|---------|---------|----------
      DATA RDH/2.7D0,2.8D0,2.7D0,1.12139906D0,-44.36D0,64.0D0,784.9D0,
     1  477.9D0, 2.27D0, 9*0.D0,-69.973D0, 6*0.D0,3.529D0,
     1  3.19D0,3.15D0,
     1  2.95D0, ! ND potential /\/\
C     1  10.6D0, ! NSC97e potential /\/\
     1  1*0.D0,
     1  2*0.D0,28*0.D0,
     2  3*1.7D0,4*0.D0,2.0D0,2.63D0, 17*0.D0,3.35D0,3.37D0,
     2  2.95D0, ! triplet d0 /\/\=singlet d0 ND
                ! not contributing in s-wave approx.
     2  1*0.D0,
     2  2*0.D0,28*0.D0,
     3  480*0.D0/
C---Corresponding square well parameters RB (width in fm) and
C-- EB =SQRT(-AM*U) (in GeV/c); U is the well height
      DATA RBH/2.545739D0,   2.779789D0, 2.585795D0, 5.023544D0,
     1 .124673D0, .3925180D0,.09D0, 2.D0, 4.058058D0, 17*0.D0,
     1  2.252623D0, 2.278575D0,
     1  2.234089D0, ! ND potential /\/\
C     1  3.065796D0, ! NSC97e potential /\/\
     1  3*0.001D0,28*0.D0,
     2  3*2.003144D0,
     2  4*0.D0, 2.D0, 4.132163D0, 17*0.D0,
     2  2.272703D0, 2.256355D0,
     2  2.234089D0, ! triplet potential /\/\=singlet ND
                    ! not contributing in s-wave FSI approx.
     2  3*0.001D0,28*0.D0,
     3  480*0.D0/
      DATA EBH/.1149517D0,    .1046257D0,   .1148757D0, .1186010D0,
     1    .7947389D0,2.281208D0,8.7D0,.4D0,.1561219D0,17*0.D0,
     1    .1013293D0, .1020966D0,
     1    .1080476D0, ! ND potential /\/\
C     1    .04115994D0, ! NSC97e potential /\/\
     1    3*0.001D0,28*0.D0,
     2    3*.1847221D0,
     2    4*0.D0, .4D0, .1150687D0, 17*0.D0,
     2    .09736083D0, .09708310D0,
     2    .1080476D0, ! triplet potential /\/\= singlet ND
                      ! not contributing in s-wave FSI approx.
     2    3*0.001D0,28*0.D0,
     3    480*0.D0/
C=======< constants >========================
      W=1/.1973D0    ! from fm to 1/GeV
      PI=4*DATAN(1.D0)
      PI2=2*PI
      SPI=DSQRT(PI)
      DR=180.D0/PI   ! from radian to degree
        AC1=1.D10
        AC2=1.D10
C=======< condition de calculs >=============
      NUNIT=11 ! for IBM or HP
C      NUNIT=4 ! for SUN in Prague
      CALL readint4(NUNIT,'ITEST     ',ITEST)
      CALL readint4(NUNIT,'LL        ',LL)        ! Two-particle system
      CALL readint4(NUNIT,'NS        ',NS)
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

       ICH=0
       IF(C1*C2.NE.0.D0) ICH=1
       IQS=0
       IF(C1+AM1.EQ.C2+AM2) IQS=1
       I3S=0  ! only this option is available
       ISI=1
       I3C=0
       IF(CN*ICH.NE.0.D0) I3C=1

      IF(ITEST.EQ.1)THEN
      CALL readint4(NUNIT,'ICH       ',ICH)
      CALL readint4(NUNIT,'IQS       ',IQS)
      CALL readint4(NUNIT,'ISI       ',ISI)
      CALL readint4(NUNIT,'I3C       ',I3C)
c      CALL READ_FILE(NUNIT,'ICH     ',CHAR,ICH,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'IQS     ',CHAR,IQS,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'ISI     ',CHAR,ISI,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'I3C     ',CHAR,I3C,REAL8,IERR)
      ENDIF
C==================================================================
C---fm to 1/GeV
      DO 3 J1=1,60
      DO 3 J2=1,10
      FDH(J1,J2)=FDH(J1,J2)*W
      RDH(J1,J2)=RDH(J1,J2)*W
 3    RBH(J1,J2)=RBH(J1,J2)*W
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
      MSP=MSPIN
      DO 91 ISPIN=1,10
  91  RHO(ISPIN)=RHOH(LL,ISPIN)
C---Integration limit AA in the spherical wave approximation
      AA=0.D0
cc      IF(NS.EQ.2.OR.NS.EQ.4)AA=.5D0 !!in 1/GeV --> 0.1 fm
ccc      IF(NS.EQ.2.OR.NS.EQ.4)AA=6.D0 !!in 1/GeV --> 1.2 fm
ccc      IF(NS.EQ.2.OR.NS.EQ.4)AA=.1D0 !!in 1/GeV --> .02 fm
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
C---Calculation continues for any system (any LL)
 55   CONTINUE
      RETURN
      END

C=======================================================
C
      SUBROUTINE FSIWF(WEI)
C==> Prepares necessary quantities, e.g. the amplitudes f_c= C(i) 
C    and call VZ(WEI) to calculate the weight due to FSI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_CVK/V,CVK
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_BP/B,P
      COMMON/FSI_ETA/ETA
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_RR/F(10)
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMPLEX*16 C,F
      COMMON/FSI_AA/AA,DW_INNER
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1    ! k* (kappa) for 2-nd channel
C***
C==>calculating relative 4-coordinates of the particles in PRF
C-  {T,X,Y,Z} from the relative coordinates {TS,XS,YS,ZS} in NRF
      XS=X1-X2
      YS=Y1-Y2
      ZS=Z1-Z2
      TS=T1-T2
      RS12=XS*P12X+YS*P12Y+ZS*P12Z
      H1=(RS12/EPM-TS)/AM12
      X=XS+P12X*H1
      Y=YS+P12Y*H1
      Z=ZS+P12Z*H1
      T=(E12*TS-RS12)/AM12
      RPS=X*X+Y*Y+Z*Z
      RP=DSQRT(RPS)
CW      WRITE(6,38)'RP ',RP,'X ',X,Y,Z,T
38    FORMAT(A7,E11.4,A7,4E11.4)

      CVK=(P12X*PPX+P12Y*PPY+P12Z*PPZ)/(P12*AK)
      V=P12/E12

      ACH=1.D0
      IF(ICH.EQ.0)GOTO 21
      XH=AC*AK
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
C---HCP, HPR needed (e.g. in GST) if ICH=1
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
      AAK=ACH*AK    !
      HCP2=2*HCP/AC ! needed to calculate C(JJ) for charged particles
  21  CONTINUE
      MSP=MSPIN
      DO 30 JJ=1,MSP
      ISPIN=JJ
      IF(NS.NE.1)GOTO22
C---Calc. quantities for the square well potential;
C-- for LL=6-26 the square well potential is not possible or available
      IF(LL.EQ.4)GOTO 22
      BK(JJ)=DSQRT(EB(JJ)**2+AKS)
      XRA=2*RB(JJ)/AC
      HRA=BK(JJ)*RB(JJ)
      CALL SEQ(XRA,HRA)
      SBKRB(JJ)=HRA*B
      HRA=AK*RB(JJ)
      CALL GST(XRA,HRA)
      SDK(JJ)=SH
      CDK(JJ)=CHH
      SDKK(JJ)=RB(JJ)
      IF(AK.NE.0.D0)SDKK(JJ)=SH/AK
      IF(ICH.EQ.1)SDK(JJ)=ACH*SDK(JJ)
  22  CONTINUE
C-----------------------------------------------------------------------
C---Calc. the strong s-wave scattering amplitude = C(JJ)
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(NS.NE.1)GOTO 230
      IF(LL.NE.4)GOTO 230 ! SW scat. amplitude used for alfa-alfa only
      GAK=G(AK)
      AKACH=AK
      IF(ICH.EQ.1)AKACH=AK*ACH
      C(JJ)=1/DCMPLX(GAK,-AKACH) ! amplitude for the SW-potential
      GOTO 30
 230  IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)GOTO20    ! pipi
      IF(LL.EQ.12.OR.LL.EQ.13)GOTO20             ! piN
      IF(LL.EQ.31.OR.LL.EQ.32)GOTO20             ! piXi
      IF(LL.EQ.8.OR.LL.EQ.9.OR.LL.EQ.18)GOTO20   ! Nd, dd
      IF(LL.EQ.14.OR.LL.EQ.17.OR.LL.EQ.23)GOTO27 ! K+K-, K-p, K0K0-b
        IF(LL.EQ.30)GOTO 28                        ! pp-bar
      A1=RD(JJ)*FD(JJ)*AKS
      A2=1+.5D0*A1
      IF(ICH.EQ.1)A2=A2-2*HCP*FD(JJ)/AC
      AKF=AK*FD(JJ)
      IF(ICH.EQ.1)AKF=AKF*ACH
      C(JJ)=FD(JJ)/DCMPLX(A2,-AKF)
      GOTO30
 20   CONTINUE
C---Calc. scatt. ampl. C(JJ) for pipi, piN, piXi and Nd, dd
      IF(LL.EQ.5)THEN
       CALL CPIPI  ! pi+pi- C(1)=f_c(pi+pi- -> pi+pi-)
c                           C(2)=f_c(pi+pi- -> pi0pi0)
       GOTO 30
      ENDIF
      IF(LL.EQ.13)THEN
       CALL CPIN  ! pi-p    C(1)=f_c(pi-p -> pi-p)
c                           C(2)=f_c(pi-p -> pi0n)
       GOTO 30
      ENDIF
      IF(LL.EQ.31)THEN
       CALL CPIXI      ! C(1)= f_c(pi+Xi-->pi+Xi-)=f(l=0,J=1/2)11/Ac
                       ! C(2)= f_c(pi+Xi-->pi0Xi0)=f(l=0,J=1/2)12/Ac
                       ! C(3)= f(l=1,J=3/2)11 = (2/3) f(l=1,J=3/2,I=1/2)
                       ! C(4)= f(l=1,J=3/2)12 = -(sqrt(2)/3) f(l=1,J=3/2,I=1/2) 
	                 ! C(5)= f(l=1,J=3/2,I=1/2)
	                 ! DW_INNER
       GOTO 30
      ENDIF
      JH=LL-7+2*JJ-2
      IF(LL.EQ.8.OR.LL.EQ.9)G2=GND(AKS,JH)
      IF(LL.EQ.18)G2=GDD(AKS,JJ)
      IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)G2=GPIPI(AKS,2)
      IF(LL.EQ.12.OR.LL.EQ.13)G2=GPIN(AKS,2)
      IF(LL.EQ.31.OR.LL.EQ.32)G2=GPIXI(AKS,2)
      C(JJ)=1.D0/DCMPLX(G2,-AK) !pi+pi+, nd, pd, pi+p, pi-Xi-, dd
      IF(LL.NE.5.AND.LL.NE.6.AND.LL.NE.13.AND.LL.NE.31)GOTO27
      IF(LL.EQ.5.OR.LL.EQ.6)G1=GPIPI(AKS,1)
      IF(LL.EQ.13)G1=GPIN(AKS,1)
      IF(LL.EQ.31)G1=GPIXI(AKS,1)
      IF(LL.EQ.5.OR.LL.EQ.13.OR.LL.EQ.31)
     c           C(JJ)=.6667D0/DCMPLX(G1,-AK)+.3333D0*C(JJ)!pi+pi-,pi-p,pi+Xi-
      IF(LL.EQ.6)C(JJ)=.3333D0/DCMPLX(G1,-AK)+.6667D0*C(JJ)!pi0pi0
 27   CONTINUE
C---Calc. K+K-, K0K0-b or K-p s-wave scatt. ampl.
      IF(LL.EQ.14.OR.LL.EQ.23)CALL CKKB
c      IF(LL.EQ.17)C(JJ)=DCMPLX(-3.29D0,3.55D0)     ! Martin'76 (K-p)
c      IF(LL.EQ.17)C(JJ)=DCMPLX(-2.585D0,4.156D0)   ! Borasoy'04 (K-p)
      IF(LL.EQ.17)C(JJ)=DCMPLX(-3.371D0,3.244D0)   ! Martin'81 (K-p)
C---Calc. pi+pi-/+, pd, pi+/-p, pi+/-Xi-, K+K- or K-p s-wave scatt. ampl.
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(ICH.EQ.0)GOTO 30
      C(JJ)=1/(1/C(JJ)-HCP2+DCMPLX(0.D0,AK-AAK))
 28   CONTINUE
C---Calc. pp-bar s-wave scatt. ampl.
      CALL CPAP
 30   CONTINUE
C***********************************************************************
      CALL VZ(WEI)
      RETURN
      END

      SUBROUTINE VZ(WEI)
C==> Calculates the weight WEI due to FSI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_JR/JRAT
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ETA/ETA
      COMMON/FSI_AA/AA,DW_INNER
      COMMON/FSI_FFF/F12,F21
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_RR/F(10)
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_COULPH/EIDC
      COMPLEX*16 F,C,G,PSI12,PSI21
      COMPLEX*16 F12,F21
      COMPLEX*16 EIDC
      COMPLEX*8 Z8,CGAMMA
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF12,FF21
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      WEI=0.D0
      IF(JRAT.EQ.1)GOTO 11
      RHOS=AK*RP
      HS=X*PPX+Y*PPY+Z*PPZ
      IF(RHOS.LT.15.D0.AND.RHOS+DABS(HS).LT.20.D0)GOTO 2
C---Calc. EIDC=exp(i*Coul.Ph.);
C-- used in calc. of hypergeom. f-s in SEQA, FAS at k*R > 15, 20
      Z8=CMPLX(1.,SNGL(ETA))
      Z8=CGAMMA(Z8)
      EIDC=Z8/CABS(Z8)
 2    CALL FF(RHOS,HS)
 11   MSP=MSPIN
      IF(ISI.EQ.0)GOTO 4  ! the strong interaction ON (OFF) if ISI=1(0)
      IF(RP.LT.AA)GOTO 4  !
      IF(JRAT.NE.1) CALL FIRT
      IF(IQS.EQ.0)GOTO 5  ! the quantum statistics ON (OFF) if IQS=1(0)
      JSIGN=-1
      DO 1 JJ=1,MSP
      JSIGN=-JSIGN
      G=F(JJ)*C(JJ)
      IF(ICH.EQ.1)G=G*ACHR
      PSI12=FF12*(F12+G)
      PSI21=FF21*(F21+G)
      G=PSI12+JSIGN*PSI21
 1    WEI=WEI+RHO(JJ)*(DREAL(G)**2+DIMAG(G)**2)
      GOTO 8
 5    FF12S=DREAL(FF12)**2+DIMAG(FF12)**2
      DO 6 JJ=1,MSP
      G=F(JJ)*C(JJ)
c---   account for p-wave resonance Xi* in piXi channel
	IF(LL.EQ.31) THEN
	COSTH=1.D0
	IF(RHOS.NE.0.D0) COSTH=HS/RHOS
	SINTHS=1-COSTH*COSTH
	CRHOR=DCOS(RHOS)/RP
	SRHOR=DSIN(RHOS)/RP
c--   F(2)= [C_0(rho)+iS_0(rho)]/r; F(3)= [C_1(rho)+iS_1(rho)]/r
      F(2)=DCMPLX(CRHOR,SRHOR) 
	F(3)=DCMPLX(CRHOR/RHOS+SRHOR,SRHOR/RHOS-CRHOR)
	G=G-2*DCMPLX(0.D0,1.D0)*F(3)*C(3)*COSTH
	ENDIF
      IF(ICH.EQ.1)G=G*ACHR
CW      WRITE(6,38)'JJ ',JJ,'F ',DREAL(F(JJ)),DIMAG(F(JJ))
CW      WRITE(6,38)'JJ ',JJ,'C ',DREAL(C(JJ)),DIMAG(C(JJ))
CW      WRITE(6,38)'JJ ',JJ,'G ',DREAL(G),DIMAG(G)
CW      WRITE(6,38)'JJ ',JJ,'F12+G ',DREAL(F12+G),DIMAG(F12+G)
CW      WRITE(6,38)'JJ ',JJ,'F21+G ',DREAL(F21+G),DIMAG(F21+G)
38    FORMAT(A7,I3,A7,2E11.4)
      PSI12=F12+G
 6    WEI=WEI+RHO(JJ)*FF12S*(DREAL(PSI12)**2+DIMAG(PSI12)**2)
c--- Account for pi0pi0 -> pi+pi- 
c---          or   pi0n -> pi-p
c---          or nn-bar -> pp-bar  
      IF(LL.EQ.5.OR.LL.EQ.13.OR.LL.EQ.30)THEN
      DO 61 JJ=1,MSP
        HH=RHO(JJ)*FF12S*(DREAL(C(JJ+MSP))**2+DIMAG(C(JJ+MSP))**2)*
     *  AMU2_AMU1*ACH/RPS
        IF(AK2S.LT.0)HH=HH*DEXP(-2*RP*AK2)
 61   WEI=WEI+HH
        ENDIF
c--- Account for pi+Xi- -> pi+Xi- spin flip p-wave resonance Xi*
      IF(LL.EQ.31)THEN
	HN1=FF12S*ACH
	G=F(3)*C(3)
      WEI=WEI+HN1*(DREAL(G)**2+DIMAG(G)**2)*SINTHS
c--          and pi0Xi0 -> pi+Xi- s-wave + p-wave resonance Xi*
	F(10)=F(3)*C(4)
	HN2=HN1*AMU2_AMU1
	G=F(2)*C(2)-2*DCMPLX(0.D0,1.D0)*F(10)*COSTH
      WEI=WEI+HN2*(DREAL(G)**2+DIMAG(G)**2)
      WEI=WEI+HN2*(DREAL(F(10))**2+DIMAG(F(10))**2)*SINTHS
	ENDIF
c------------------------------------------------------------------
      RETURN
 4    PSI12=FF12*F12
      IF(IQS.EQ.0)GOTO 50 ! the quantum statistics ON (OFF) if IQS=1(0)
      PSI21=FF21*F21
      JSIGN=-1
      DO 3 JJ=1,MSP
      JSIGN=-JSIGN
      G=PSI12+JSIGN*PSI21
 3    WEI=WEI+RHO(JJ)*(DREAL(G)**2+DIMAG(G)**2)
      GOTO 8
 50   WEI=DREAL(PSI12)**2+DIMAG(PSI12)**2
      FF12S=DREAL(FF12)**2+DIMAG(FF12)**2
c--- Account for the contribution of the inner region r < AA 
C--  to the pi+Xi- weight due to the p-wave resonance Xi*
      IF(LL.EQ.31.AND.ISI.NE.0)THEN
	HN1=FF12S*ACH	
      WEI=WEI+HN1*DW_INNER
	ENDIF
      RETURN
 8    WEI=WEI/2
      RETURN
      END

      SUBROUTINE FIRT
C---CALC. THE F(JJ)
C-- F(JJ)*C(JJ)= DEVIATION OF THE BETHE-SALPETER AMPL. FROM PLANE WAVE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_BP/B,P
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_C/C(10)/FSI_AM/AM,AMS,DM,SM,AM1S,AM2S
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_RR/F(10)
      EQUIVALENCE(RSS,RP),(TSS,T)
      COMPLEX*16 F,C,CH1
      MSP=MSPIN
      DO 10 JJ=1,MSP
      IF(JJ.GT.1)GOTO 3
      XRA=2*RSS/AC
      IF(AK.NE.0.D0)GOTO2
      SHK=1.D0
      SHR=.0D0
      SHHR=SHR
      CHHR=1/RSS
      GOTO3
  2   H=AK*RSS
      CALL GST(XRA,H)
      SHR=SH/RSS
      CHHR=CHH/RSS
      SHHR=SHR
      IF(ICH.EQ.1) SHHR=ACH*SHR
  3   IF(NS.EQ.2)GOTO1
C---F= ASYMPTOTIC FORMULA (T= 0 APPROX.); NS= 4
  6   F(JJ)=DCMPLX(CHHR,SHHR)
      IF(NS.NE.1)GOTO 10
C---F INSIDE THE SQUARE-WELL (T= 0 APPROX.); NS= 1
      IF(RSS.GE.RB(JJ)) GOTO 10
      IF(AK.NE.0.D0.AND.JJ.EQ.1)SHK=SHR/AK
      H=BK(JJ)*RSS
      CALL GST(XRA,H)
      SKR=B*BK(JJ)
      F(JJ)=DCMPLX(CDK(JJ),SDK(JJ))*SKR
      CH1=(SDKK(JJ)*SKR-SHK*SBKRB(JJ))/C(JJ)
      F(JJ)=(F(JJ)+CH1)/SBKRB(JJ)
      GOTO 10
  1   CONTINUE
C---F= ASYMPTOTIC FORMULA (T= 0 NOT REQUIRED); NS= 2
      IF(JJ.GT.1)GOTO 8
      IF(TSS.EQ.0.D0)GOTO6
      TSSA=DABS(TSS)
      IF(DM.NE.0.D0)GOTO 11
      H=AM*.5D0/TSSA
      IF(AK.NE.0.D0)GOTO4
      HM=H*RPS
      IF(HM.GE.3.D15)GOTO6
      FS1=DFRSIN(HM)
      FC1=DFRCOS(HM)
      FC2=FC1
      FS2=FS1
      GOTO5
  4   CONTINUE
      H1=AK*TSSA/AM
      HM=H*(RSS-H1)**2
      HP=H*(RSS+H1)**2
      IF(HP.GE.3.D15)GOTO6
      FS1=DFRSIN(HM)
      FC1=DFRCOS(HM)
      FS2=DFRSIN(HP)
      FC2=DFRCOS(HP)
      GOTO 5
  11  CONTINUE
      FS1=0.D0
      FS2=0.D0
      FC1=0.D0
      FC2=0.D0
      DO 13 I=1,2
      IF(I.EQ.1)TSSH=TSSA*(1+DM)
      IF(I.EQ.2)TSSH=TSSA*(1-DM)
      H=AM*.5D0/TSSH
      IF(AK.NE.0.D0)GOTO 12
      HM=H*RPS
      IF(HM.GE.3.D15)GOTO6
      FS1=FS1+DFRSIN(HM)/2
      FC1=FC1+DFRCOS(HM)/2
      IF(I.EQ.1)GOTO 13
      FC2=FC1
      FS2=FS1
      GOTO 13
  12  CONTINUE
      H1=AK*TSSH/AM
      HM=H*(RSS-H1)**2
      HP=H*(RSS+H1)**2
      IF(HP.GE.3.D15)GOTO6
      FS1=FS1+DFRSIN(HM)/2
      FC1=FC1+DFRCOS(HM)/2
      FS2=FS2+DFRSIN(HP)/2
      FC2=FC2+DFRCOS(HP)/2
  13  CONTINUE
  5   C12=FC1+FS2
      S12=FS1+FC2
      A12=FS1-FC2
      A21=FS2-FC1
      A2=.5D0*(CHHR*(A12+A21)+SHR*(A12-A21))+SHHR
      A1=.5D0*(CHHR*(C12+S12)+SHR*(C12-S12))
      F(JJ)=.3989422D0*DCMPLX(A1,A2)
      GOTO 10
  8   F(JJ)=F(1)
 10   CONTINUE
      RETURN
      END

      FUNCTION EXF(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.LT.-15.D0) GO TO 1
      EXF=DEXP(X)
      RETURN
  1   EXF=.0D0
      RETURN
      END

      SUBROUTINE SEQ(X,H)
C---CALC. FUNCTIONS B, P (EQS. (17) OF G-K-L-L);
C-- NEEDED TO CALC. THE CONFLUENT HYPERGEOMETRIC FUNCTION GST.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_BP/B,P
      DIMENSION BH(3),PH(3)
      DATA ERR/1.D-7/
      BH(1)=1.D0
      PH(1)=1.D0
      PH(2)=.0D0
      BH(2)=.5D0*X
      B=1+BH(2)
      P=1.D0
      HS=H*H
      J=0
  2   J=J+1
      BH(3)=(X*BH(2)-HS*BH(1))/((J+1)*(J+2))
      PH(3)=(X*PH(2)-HS*PH(1)-(2*J+1)*X*BH(2))/(J*(J+1))
      B=B+BH(3)
      P=P+PH(3)
      Z=DABS(BH(2))+DABS(BH(3))+DABS(PH(2))+DABS(PH(3))
      IF(Z.LT.ERR)RETURN
      BH(1)=BH(2)
      BH(2)=BH(3)
      PH(1)=PH(2)
      PH(2)=PH(3)
      GOTO 2
      END

      SUBROUTINE SEQA(X,H)
C---CALC. FUNCTIONS CHH=REAL(GST), SH=IMAG(GST)/ACH, B=SH/H
C-- IN THE ASYMPTOTIC REGION H=K*R >> 1.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_BP/B,P
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_ETA/ETA
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_COULPH/EIDC
      COMPLEX*16 EIDC,GST
      ARG=H-ETA*DLOG(2*H)
      GST=DCMPLX(DCOS(ARG),DSIN(ARG))
      GST=ACHR*EIDC*GST
      CHH=DREAL(GST)
      SH=DIMAG(GST)/ACH
      B=SH/H
      RETURN
      END

      SUBROUTINE FF(RHO,H)
C---Calc. F12, F21;
C-- F12= FF0* plane wave,  FF0=F*ACHR,
C---F is the confluent hypergeometric function,
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMMON/FSI_FFF/F12,F21
      COMPLEX*16 FF0,F12,F21
      C=DCOS(H)
      S=DSIN(H)
      F12=DCMPLX(C,-S)
      F21=DCMPLX(C,S)
      IF(ICH.EQ.0)RETURN
      RHOP=RHO+H
      RHOM=RHO-H
      F12=FF0(RHO,H)*F12
      F21=FF0(RHO,-H)*F21
      RETURN
      END

      FUNCTION FAS(RKS)
C-- FAS=F*ACHR
C---F is the confluent hypergeometric function at k*r >> 1
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 FAS,EIDC,ZZ1
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_ETA/ETA
      COMMON/FSI_COULPH/EIDC
      D1=DLOG(RKS)*ETA
      D2=ETA*ETA/RKS
      ZZ1=DCMPLX(DCOS(D1),DSIN(D1))/EIDC
      FAS=DCMPLX(1.D0,-D2)*ZZ1
      FAS=FAS-DCMPLX(DCOS(RKS),DSIN(RKS))*ETA/RKS/ZZ1
      RETURN
      END

      FUNCTION FF0(RHO,H)
C-- FF0=F*ACHR
C-- F is the confluent hypergeometric function
C-- (Eq. (15) of G-K-L-L), F= 1 at r* << AC
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMPLEX*16 ALF,ALF1,Z,S,A,FF0,FAS
      DATA ERR/1.D-5/
      S=DCMPLX(1.D0,0.D0)
      FF0=S
      RHOP=RHO+H
CC      GOTO 5 ! rejects the approx. calcul. of hyperg. f-ion F
      IF(RHOP.LT.20.D0)GOTO5
      FF0=FAS(RHOP) ! approx. calc.
      RETURN
  5   ALF=DCMPLX(.0D0,-ETA)
      ALF1=ALF-1
      Z=DCMPLX(.0D0,RHOP)
      J=0
  3   J=J+1
      A=(ALF1+J)/(J*J)
      S=S*A*Z
      FF0=FF0+S
      ZR=DABS(DREAL(S))
      ZI=DABS(DIMAG(S))
      IF((ZR+ZI).GT.ERR)GOTO3
      FF0=FF0*ACHR
      RETURN
      END

      FUNCTION HC(XA)
C---HC = h-function of Landau-Lifshitz: h(x)=Re[psi(1-i/x)]+ln(x)
C-- psi(x) is the digamma function (the logarithmic derivative of
C-- the gamma function)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BN(15)
      DATA BN/.8333333333D-1,.8333333333D-2,.396825396825D-2,
     1        .4166666667D-2,.7575757576D-2,.2109279609D-1,
     2        .8333333333D-1,.4432598039D0 ,.305395433D1,
     3        .2645621212D2, .2814601449D3, .3607510546D4,
     4        .5482758333D5, .9749368235D6, .200526958D8/
      X=DABS(XA)
      IF(X.LT..33D0) GOTO 1
CC      IF(X.GE.3.5D0) GO TO 2
      S=.0D0
      N=0
   3  N=N+1
      DS=1.D0/N/((N*X)**2+1)
      S=S+DS
      IF(DS.GT.0.1D-12) GOTO 3
C---Provides 7 digit accuracy
      HC=S-.5772156649D0+DLOG(X)
      RETURN
CC   2  HC=1.2D0/X**2+DLOG(X)-.5772156649 D0
CC      RETURN
   1  X2=X*X
      XP=X2
      HC=0.D0
      IMA=9
      IF(X.LT.0.1D0)IMA=3
      DO 4 I=1,IMA
      HC=HC+XP*BN(I)
   4  XP=X2*XP
      RETURN
      END

      FUNCTION ACP(X)
C--- ACP = COULOMB PENETRATION FACTOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.LT.0.05D0.AND.X.GE.0.D0) GO TO 1
      Y=6.2831853D0/X
      ACP=Y/(EXF(Y)-1)
      RETURN
   1  ACP=1.D-6
      RETURN
      END

      SUBROUTINE GST(X,H)
C---CALC. THE CONFL. HYPERGEOM. F-N: GST = CHH+i*SHH
C-- AND THE COULOMB F-S B, P (CALLS SEQ OR SEQA).
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_BP/B,P
  1   IF(ICH.EQ.1)GOTO 2
  3   SH=DSIN(H)
      CHH=DCOS(H)
      B=1.D0
      IF(H.NE.0.D0)B=SH/H
      P=CHH
      RETURN
  2   CONTINUE
      IF(H.GT.15.D0)GOTO4 ! comment out if you want to reject
                   ! the approximate calculation of hyperg. f-ion G
      CALL SEQ(X,H) ! exact calculation
      SH=H*B        ! = SHH/ACH
      CHH=P+B*X*(DLOG(DABS(X))+HPR)
      RETURN
  4   CALL SEQA(X,H)
      RETURN
      END

      FUNCTION FF1(RHO,H)
C---FF1=FF0; used for particle-nucleus system
C-- FF0=F12*ACHR
C-- F12 is the confluent hypergeometric function
C-- (Eq. (15) of G-K-L-L), F12= 1 at r* << AC
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMMON/FSI_COULPH/EIDC
      COMMON/FSI_ICH1/ICH1
      COMPLEX*16 FF0,FF1
      COMPLEX*16 EIDC
      COMPLEX*8 Z8,CGAMMA
      FF1=DCMPLX(1.D0,0.D0)
      IF(ICH1.EQ.0)GOTO 2
      IF(RHO.LT.15.D0.AND.RHO+H.LT.20.D0)GOTO 2
C---Calc. EIDC=exp(i*Coul.Ph.);
C-- used in calc. of hypergeom. f-s in SEQA, FAS at k*R > 15, 20
      Z8=CMPLX(1.,SNGL(ETA))
      Z8=CGAMMA(Z8)
      EIDC=Z8/CABS(Z8)
 2    FF1=FF0(RHO,H)
      RETURN
      END

      FUNCTION G(AK)
C---Used to calculate SW scattering amplitude for alpa-alpha system
C-- and for sear.f (square well potential search)
C---NOTE THAT SCATT. AMPL.= 1/CMPLX(G(AK),-AK*ACH)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,JJ,MSPIN
      COMMON/FSI_BP/B,P/FSI_DERIV/BPR,PPR/FSI_SHH/SH,CHH
      COMMON/FSI_DAK/DAK,IFUN
      HCP2=.0D0
      ACH=1.D0
      IF(ICH.EQ.0)GOTO 1
      XH=AC*AK
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
      ACH=ACP(XH)
      HCP2=2*HCP/AC
   1  AKS=AK**2
      BK(JJ)=DSQRT(AKS+EB(JJ)**2) ! kappa=kp
      X=2*RB(JJ)/AC
      H=BK(JJ)*RB(JJ)             ! kp*d
      CALL GST(X,H)
      BRHO=B                      ! B(kp,d)
      SBKRB(JJ)=SH                ! kp*d*B(kp,d)
      CALL DERIW(X,H)
      BRHOP=BPR                   ! B'(kp,d)= dB(kp,r)/dln(r) at r=d
      H=AK*RB(JJ)
      CALL GST(X,H)
      CDK(JJ)=CHH                 ! ReG(k,d)
      BRHOS=B                     !  B(k,d)
      PRHOS=P                     !  P(k,d)
      SDK(JJ)=SH
      IF(ICH.EQ.0)GOTO 2
      SDK(JJ)=ACH*SH              ! ImG(k,d)
      IF(AK.EQ.0.D0.AND.AC.LT.0.D0)SDK(JJ)=3.14159*X*B
   2  SDKK(JJ)=RB(JJ)
      IF(AK.NE.0.D0)SDKK(JJ)=SH/AK ! d*B(k,d)
      CALL DERIW(X,H)              ! PPR=P'(k,d)= dP(k,r)/dln(r) at r=d
      ZZ=PPR-PRHOS
      IF(ICH.EQ.1)ZZ=ZZ+X*(BRHOS+BPR*(DLOG(DABS(X))+HPR))
C   ZZ= P'(k,d)-P(k,d)+x*{B(k,d)+B'(k,d)*[ln!x!+2*C-1+h(k*ac)]}
      GG=(BRHOP*CDK(JJ)-BRHO*ZZ)/RB(JJ)
C   GG= [B'(kp,d)*ReG(k,d)-B(kp,d)*ZZ]/d
      G=GG/(BRHO*BPR-BRHOP*BRHOS)
C    G= GG/[B(kp,d)*B'(k,d)-B'(kp,d)*B(k,d)]
      RETURN
      END

      SUBROUTINE DERIW(X,H)
C---CALLED BY F-N G(AK)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_BP/B,P/FSI_DERIV/BPR,PPR
      HH=.1D-3
      CALL GST(X,H-HH)
      Q1=P
      B1=B
      CALL GST(X,H+HH)
      HHH=HH+HH
      BPR=H*(B-B1)/HHH
      PPR=H*(P-Q1)/HHH
      IF(ICH.EQ.0)RETURN
      CALL GST(X-HH,H)
      Q1=P
      B1=B
      CALL GST(X+HH,H)
      BPR=BPR+X*(B-B1)/HHH
      PPR=PPR+X*(P-Q1)/HHH
      RETURN
      END


C=============================================================
      SUBROUTINE READ_FILE(NUNIT,KEY,CH8,INT4,REAL8,IERR)
C     ==========
C
C     Routine to read one parameter of the program in the file
C     DATA NUNIT defined FSIINI
C     NUNIT=11 for IBM or HP, 4 for SUN in Prague
C
C     INPUT  : KEY (CHARACTER*10) :
C     OUTPUT : case of KEY : CH8   : (CHARACTER*8)
C                            INT4  : (INTEGER*4)
C                            REAL8 : (REAL*8)
C                     (only one of them)
C              IERR (INTEGER) : 0 : no error
C                               1 : key not found

      CHARACTER*10 KEY,TEST
      CHARACTER*4  TYPE
      CHARACTER*8  CH8
      INTEGER*4    INT4
      REAL*8       REAL8
      INTEGER*4    IERR

      IERR=0
      REWIND(NUNIT)
1     READ(NUNIT,FMT='(A10,2X,A4,64X)')TEST,TYPE
      IF (TEST.EQ.KEY) THEN
        BACKSPACE(NUNIT)
        IF (TYPE.EQ.'CHAR') READ(NUNIT,FMT='(18X,A8,54X)')CH8
        IF (TYPE.EQ.'INT4') READ(NUNIT,FMT='(18X,I8,54X)')INT4
        IF (TYPE.EQ.'REA8') READ(NUNIT,FMT='(18X,F10.5,52X)')REAL8
      ELSE
        IF (TEST.NE.'* E.O.F. *') THEN
          GOTO 1
        ELSE
          IERR=1
        ENDIF
      ENDIF
c      IF(IERR.EQ.1)STOP
      RETURN
      END
C================================================================

      SUBROUTINE readrea8(lun,key,rea8)
C---> Read a value (rea8) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 text
      REAL*8 rea8
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'REA8') THEN
        READ(lun,FMT='(18X,E10.4,A50)',ERR=13,END=13) rea8,text
        PRINT 2,   key,rea8,text
        WRITE(9,2) key,rea8,text
  2     FORMAT('      REA8|==> ',a10,'= ',E10.4,2x,a50)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readrea4(lun,key,rea4)
C---> Read a value (rea4) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'REA4') THEN
        READ(lun,FMT='(18X,E10.4,A50)',ERR=13,END=13) rea4,text
        PRINT 2,   key,rea4,text
        WRITE(9,2) key,rea4,text
  2     FORMAT('      REA4|==> ',a10,'= ',E10.4,2x,a50)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readint4(lun,key,int4)
C---> Read a value (int4) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'INT4') THEN
        READ(lun,FMT='(18X,I8,A50)',ERR=13,END=13) int4,text
        PRINT 2,   key,int4,text
        WRITE(9,2) key,int4,text
  2     FORMAT('      INT4|==> ',a10,'= ',I9,2x,a50)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readch50(lun,key,ch50)
C---> Read a value (ch50) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 ch50
      CHARACTER*12 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'CH50') THEN
C        READ(lun,FMT='(18X,A50,A12)',ERR=13,END=13) ch50,text
        READ(lun,FMT='(18X,A50,A12)',ERR=13) ch50,text
        PRINT 2,   key,ch50,text
        WRITE(9,2) key,ch50,text
  2     FORMAT('      CH50|==> ',a10,'= ',A50,2x,a12)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readch30(lun,key,ch30)
C---> Read a value (ch30) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*30 ch30
      CHARACTER*32 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'CH30') THEN
C        READ(lun,FMT='(18X,A30,A32)',ERR=13,END=13) ch30,text
        READ(lun,FMT='(18X,A30,A32)',ERR=13) ch30,text
        PRINT 2,   key,ch30,text
        WRITE(9,2) key,ch30,text
  2     FORMAT('      CH30|==> ',a10,'= ',A30,2x,a32)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readtest(mode,lun,key,type,ierr)
C---> Test of reading a variable (key) from unit (lun).
C     =================================================
      CHARACTER*10 key,test
      CHARACTER*4  type
      data istart/1/
      IF(istart.eq.1) THEN
        istart=0
        PRINT 10,    lun
        WRITE (9,10) lun
  10    FORMAT(' READ-PAR.|==> Control values from unit (',i2,')',/,
     &         50('='))
      END IF
      IF(mode.eq.1) THEN
      ierr=0
      REWIND(lun)
  1   READ(lun,FMT='(A10,2X,A4)',ERR=13,END=13) test,type
        IF(test.EQ.key) THEN
          BACKSPACE(lun)
          RETURN
        ELSE
          GOTO 1
        END IF
 13   ierr=1
      RETURN
      ELSE
      PRINT 14, key,test,type
 14   FORMAT(' *** ReadErr: key,test,type=',
     &       2a10,1x,a4,i9,/' *** STOP')
      STOP
      END IF
      END
