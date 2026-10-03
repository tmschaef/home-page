      PROGRAM MAIN
C--------------------------------------------------------------------------C
C     FIT EQUAL TIME NUCLEON CORRELATION FUNCTIONS WITH SIMPLE             C
C     NUCLEON POLE PLUS CONTIUUM MODEL.                                    C
C--------------------------------------------------------------------------C
C     THIS PROGRAM IS SLIGHTLY MODIFIED TO INCORPORATE NEW INPUT FORMAT    C
C     AND DIFFERENT NORMALISATION OF K2-K6.                                C
C--------------------------------------------------------------------------C
C     VERSION            :    3.1                                          C
C     CREATION DATE      :  11-24-92                                       C
C     LAST MODIFICATION  :  07-07-93                                       C
C--------------------------------------------------------------------------C
C     THE CORRELATION FUNCTIONS ARE DEFINED AS FOLLOWS:                    C
C     F1     11-CORRELATOR, SCALAR STRUCTURE                               C
C     F2     11-CORRELATOR, GAMMA.TAU STRUCTURE                            C
C     F3     22-CORRELATOR, SCALAR STRUCTURE                               C
C     F4     22-CORRELATOR, GAMMA.TAU STRUCTURE                            C
C     F5     12-CORRELATOR, SCALAR STRUCTURE                               C
C     F6     12-CORRELATOR, GAMMA.TAU STRUCTURE                            C
C--------------------------------------------------------------------------C
C     N      NUMBER OF POINTS IN THE DATA FILE                             C
C     NP     NUMBER OF POINTS USED IN THE FIT                              C
C     NPAR   NUMBER OF FIT PARAMETERS                                      C
C     IS     FIT STRATEGY                                                  C
C--------------------------------------------------------------------------C
C     INPUT:   FOR015.DAT                                                  C
C     TITLE                                                                C
C     X(INDEX) F1(INDEX) DF1(INDEX)                                        C
C       ...      ...       ...                                             C
C                                                                          C
C--------------------------------------------------------------------------C
C     OUTPUT:  FOR016.DAT                                                  C
C     TAU    FIT(TAU) DATA(TAU) DEL(TAU) OPE(TAU)                          C
C      ...    ...      ...       ...      ...                              C
C                                                                          C
C--------------------------------------------------------------------------C
C     PLOT FILE: NUCP.MGO                                                  C
C--------------------------------------------------------------------------C
C     NEW VERSION USES MINUIT MINIMIZATION AND ERROR ANNALYSIS.            C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      REAL*8 F1(100),F2(100),F3(100),F4(100),F5(100),F6(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100),DF5(100),DF6(100)
      REAL*8 EF1(100),EF2(100),EF3(100),EF4(100),EF5(100),EF6(100)
      REAL*8 X(100)
      REAL*8 C(10),ARGLIS(10)
      INTEGER I,J,K,N,NP,NPAR,NPARI,NPARX,IS
      INTEGER IERFLG,IVAR,ISTAT
      INTEGER IPARAM(7)
      CHARACTER*80 NAME(10)
      EXTERNAL CHISQU
      COMMON /DATA/ X,F1,F2,F3,F4,F5,F6,DF1,DF2,DF3,DF4,DF5,DF6,NP
      COMMON /OPE/  QQ,SS,G2,M0,MS
C--------------------------------------------------------------------------C
C     INPUT                                                                C
C--------------------------------------------------------------------------C
      OPEN(UNIT=15,FILE='NUC.DAT',STATUS='OLD')
      OPEN(UNIT=16,FILE='FITNUC.DAT',STATUS='OLD')
      N = 15
      WRITE(6,*) 'DETERMINE FIT WINDOW'
      WRITE(6,*) 'MINIMUM DISTANCE (LAM^-1))'
      READ(5,*)  X0
      WRITE(6,*) 'MAXIMUM DISTANCE'
      READ(5,*)  X1
      WRITE(6,*) 'DETERMINE FITTING PROCEDURE'
      WRITE(6,*) '  (1) STATISTICAL WEIGHT'
      WRITE(6,*) '  (2) MODIFIED STAT. WEIGHT'
      WRITE(6,*) '  (3) EQUAL WEIGHTS'
      READ(5,*) IS
      IF (IS .EQ. 2) THEN
      WRITE(6,*) 'MINIMUM ERROR'
      READ(5,*)  DEL
      ENDIF
      WRITE(6,*) 'SCALE PARAMETER (GEV)'
      READ(5,*) CONV
C     WRITE(6,*) 'QUARK CONDENSATE (LAM**3)'
C     READ(5,*) QQ
C-------------------------------------------------------------------------C
      NP = 0
      READ(15,*) 
      DO 10 I=1,N
         READ(15,*) XX,FF1,DFF1
         IF (XX .GE. X0 .AND. XX .LE. X1) THEN
            NP = NP + 1
            X(NP) = XX
            F1(NP) = FF1
            EF1(NP)= DFF1
         ENDIF
 10   CONTINUE
C-------------------------------------------------------------------------C
      NP = 0
      READ(15,*) 
      DO 11 I=1,N
         READ(15,*) XX,FF2,DFF2
         IF (XX .GE. X0 .AND. XX .LE. X1) THEN
            NP = NP + 1
            X(NP) = XX
            F2(NP) = FF2
            EF2(NP)= DFF2
         ENDIF
 11   CONTINUE
C-------------------------------------------------------------------------C
      NP = 0
      READ(15,*)
      DO 12 I=1,N
         READ(15,*) XX,FF3,DFF3
         IF (XX .GE. X0 .AND. XX .LE. X1) THEN
            NP = NP + 1
            X(NP) = XX
            F3(NP) = FF3
            EF3(NP)= DFF3
         ENDIF
 12   CONTINUE
C-------------------------------------------------------------------------C
      NP = 0
      READ(15,*)
      DO 13 I=1,N
         READ(15,*) XX,FF4,DFF4
         IF (XX .GE. X0 .AND. XX .LE. X1) THEN
            NP = NP + 1
            X(NP) = XX
            F4(NP) = FF4
            EF4(NP)= DFF4
         ENDIF
 13   CONTINUE
C-------------------------------------------------------------------------C
      NP = 0
      READ(15,*)
      DO 14 I=1,N
         READ(15,*) XX,FF5,DFF5
         IF (XX .GE. X0 .AND. XX .LE. X1) THEN
            NP = NP + 1
            X(NP) = XX
            F5(NP) = FF5
            EF5(NP)= DFF5
         ENDIF
 14   CONTINUE
C--------------------------------------------------------------------------C
      NP = 0
      READ(15,*)
      DO 15 I=1,N
         READ(15,*) XX,FF6,DFF6
         IF (XX .GE. X0 .AND. XX .LE. X1) THEN
            NP = NP + 1
            X(NP) = XX
            F6(NP) = FF6
            EF6(NP)= DFF6
         ENDIF
 15   CONTINUE
C--------------------------------------------------------------------------C
C     DETERMINE STATISTIACL WEIGHT                                         C
C--------------------------------------------------------------------------C
      DO 16 I=1,NP
         IF (IS .EQ. 1) THEN
            DF1(I) = EF1(I)
            DF2(I) = EF2(I)
            DF3(I) = EF3(I)
            DF4(I) = EF4(I)
            DF5(I) = EF5(I)
            DF6(I) = EF6(I)
         ELSE IF (IS .EQ. 2) THEN
            DF1(I) = MAX(EF1(I),DEL)
            DF2(I) = MAX(EF2(I),DEL)
            DF3(I) = MAX(EF3(I),DEL)
            DF4(I) = MAX(EF4(I),DEL)
            DF5(I) = MAX(EF5(I),DEL)
            DF6(I) = MAX(EF6(I),DEL)
         ELSE 
            DF1(I) = 1.D0
            DF2(I) = 1.D0
            DF3(I) = 1.D0
            DF4(I) = 1.D0
            DF5(I) = 1.D0
            DF6(I) = 1.D0
         ENDIF
 16   CONTINUE
C--------------------------------------------------------------------------C
C     VACUUM CONDENSATES                                                   C
C--------------------------------------------------------------------------C
      CONV3= CONV**3
      QQ   =-0.25D0**3/CONV3
      SS   = 0.82D0*QQ
      G2   = 0.36D0**4/CONV**4
      M0   = 0.70D0/CONV
      MS   = 0.15D0/CONV
C--------------------------------------------------------------------------C
C     INITIALIZE MINUIT                                                    C
C--------------------------------------------------------------------------C
      CALL MNINIT(5,6,7)
C--------------------------------------------------------------------------C
C     SYNTAX: NUMBER,NAME,START,ERROR,LBOUND,UBOUND,ERRORFLAG              C
C--------------------------------------------------------------------------C
      CALL MNPARM(1,'MN',4.7D0,1.D0,0.D0,1.D2,IERFLG)
      CALL MNPARM(2,'L1',3.9D0,1.D0,0.D0,1.D4,IERFLG)
      CALL MNPARM(3,'E1',8.5D0,1.D0,0.D0,2.D2,IERFLG)
      CALL MNPARM(4,'L2',9.0D0,1.D0,0.D0,1.D4,IERFLG)
      CALL MNPARM(5,'E2',8.5D0,1.D0,0.D0,2.D2,IERFLG)
      CALL MNPARM(6,'E3',6.0D0,1.D0,0.D0,2.D2,IERFLG)
      CALL MNPARM(7,'E4',6.0D0,1.D0,0.D0,2.D2,IERFLG)
C--------------------------------------------------------------------------C
C     FIX PARAMETERS DURING MINIMIZATION (OPTIONAL)                        C
C--------------------------------------------------------------------------C
C     ARGLIS(1) = 1
C     CALL MNEXCM(CHISQU,'FIX',ARGLIS,1,IERFLG)
C--------------------------------------------------------------------------C
C     MINIMIZATION AND ERROR ANALYSIS                                      C
C--------------------------------------------------------------------------C
      ARGLIS(1) = 0
      CALL MNEXCM(CHISQU,'SET_PRINT',ARGLIS,1,IERFLG)
      CALL MNEXCM(CHISQU,'MIGRAD',ARGLIS,0,IERFLG)
      CALL MNEXCM(CHISQU,'MINOS',ARGLIS,0,IERFLG)
      CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
C--------------------------------------------------------------------------C  
C     RESULTS                                                              C
C--------------------------------------------------------------------------C
      CALL MNPOUT(1,NAME(1),MN,DMN,BND1,BND2,IVAR)
      CALL MNPOUT(2,NAME(2),LAM1,DLAM1,BND1,BND2,IVAR)
      CALL MNPOUT(3,NAME(3),E1,DE1,BND1,BND2,IVAR)
      CALL MNPOUT(4,NAME(4),LAM2,DLAM2,BND1,BND2,IVAR)
      CALL MNPOUT(5,NAME(5),E2,DE2,BND1,BND2,IVAR)
      CALL MNPOUT(6,NAME(6),E3,DE3,BND1,BND2,IVAR)
      CALL MNPOUT(7,NAME(7),E4,DE4,BND1,BND2,IVAR)
      CALL MNERRS(1,MNUP,MNLOW,MNPAR,MNCOR)
      CALL MNERRS(2,L1UP,L1LOW,L1PAR,L1COR)
      CALL MNERRS(3,E1UP,E1LOW,E1PAR,E1COR)
      CALL MNERRS(4,L2UP,L2LOW,L1PAR,L1COR)
      CALL MNERRS(5,E2UP,E2LOW,E2PAR,E2COR)
      CALL MNERRS(6,E3UP,E3LOW,E3PAR,E3COR)
      CALL MNERRS(7,E4UP,E4LOW,E4PAR,E4COR)
C--------------------------------------------------------------------------C
C     NORMALIZE LAMBDA2                                                    C
C--------------------------------------------------------------------------C
      LP   = LAM2/(2.D0*SQRT(3.D0))
      LPUP = L2UP/(2.D0*SQRT(3.D0))
      LPLOW= L2LOW/(2.D0*SQRT(3.D0))
C--------------------------------------------------------------------------C
C     COMPARE WITH DATA                                                    C
C--------------------------------------------------------------------------C
      WRITE(16,*) '11-CORRELATOR, SCALAR PART'
      DO 100 I=1,NP
         TAU = X(I)
         WRITE(16,555) TAU,K1(TAU,MN,LAM1,E3)
     1                  ,F1(I),EF1(I),TH1(TAU)
 100  CONTINUE
      WRITE(16,*) '11-CORRELATOR, GAMMA.TAU PART'
      DO 110 I=1,NP
         TAU = X(I)
         WRITE(16,555) TAU,K2(TAU,MN,LAM1,E1),F2(I),EF2(I),TH2(TAU)
 110  CONTINUE
      WRITE(16,*) '22-CORRELATOR, SCALAR PART'
      DO 120 I=1,NP
         TAU = X(I)
         WRITE(16,555) TAU,K3(TAU,MN,LAM2),F3(I),EF3(I),TH3(TAU)
 120  CONTINUE
      WRITE(16,*) '22-CORRELATOR, GAMMA.TAU PART'
      DO 130 I=1,NP
         TAU = X(I)
         WRITE(16,555) TAU,K4(TAU,MN,LAM2,E2),F4(I),EF4(I),TH4(TAU)
 130  CONTINUE
      WRITE(16,*) '12-CORRELATOR, SCALAR PART'
      DO 140 I=1,NP
         TAU = X(I)
         WRITE(16,555) TAU,K5(TAU,MN,LAM1,LAM2,E4),
     1                                  F5(I),EF5(I),TH5(TAU)
 140  CONTINUE
      WRITE(16,*) '12-CORRELATOR, GAMMA.TAU PART'
      DO 150 I=1,NP
         TAU = X(I)
         WRITE(16,555) TAU,K6(TAU,MN,LAM1,LAM2),F6(I),EF6(I),TH6(TAU)
 150  CONTINUE
C--------------------------------------------------------------------------C
C     OUTPUT                                                               C
C--------------------------------------------------------------------------C
      DO 200 I=6,16,10
         WRITE(I,999)
         WRITE(I,998) NAME(1),MN*CONV,MNUP*CONV,MNLOW*CONV 
         WRITE(I,998) NAME(2),LAM1*CONV3,L1UP*CONV3,L1LOW*CONV3
         WRITE(I,998) NAME(3),E1*CONV,E1UP*CONV,E1LOW*CONV
         WRITE(I,998) NAME(4),LAM2*CONV3,L2UP*CONV3,L2LOW*CONV3
         WRITE(I,998) 'SCAL.',LP*CONV3,LPUP*CONV3,LPLOW*CONV3
         WRITE(I,998) NAME(5),E2*CONV,E2UP*CONV,E2LOW*CONV
         WRITE(I,998) NAME(6),E3*CONV,E3UP*CONV,E3LOW*CONV
         WRITE(I,998) NAME(7),E4*CONV,E4UP*CONV,E4LOW*CONV
         WRITE(I,997) 'CHI2F',FMIN/(6.D0*NP-NPARI)
         WRITE(I,999)
 200  CONTINUE
C--------------------------------------------------------------------------C
 111  FORMAT(1X,F12.5)
 222  FORMAT(1X,F12.5,1X,F12.5)
 333  FORMAT(1X,F12.5,1X,F12.5,1X,F12.5)
 444  FORMAT(1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5)
 555  FORMAT(1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5)
 997  FORMAT(1X,A5,1X,F12.5)
 998  FORMAT(1X,A5,1X,F12.5,1X,F12.5,1X,F12.5)
 999  FORMAT(70('-'))
      END

      SUBROUTINE CHISQU(NPAR,GRAD,FUNC,C,IFLAG)
C--------------------------------------------------------------------------C
C     EVALUATES CHISQU FOR FIT TO MODEL SPECTRUM                           C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER I,NP,NPAR,IS,IFLAG
      REAL*8 C(NPAR),GRAD(NPAR)
      REAL*8 X(100)
      REAL*8 F1(100),F2(100),F3(100),F4(100),F5(100),F6(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100),DF5(100),DF6(100)
      COMMON /DATA/ X,F1,F2,F3,F4,F5,F6,DF1,DF2,DF3,DF4,DF5,DF6,NP
C--------------------------------------------------------------------------C
      MN   = C(1)
      LAM1 = C(2)
      E1   = C(3)
      LAM2 = C(4)
      E2   = C(5)
      E3   = C(6)
      E4   = C(7)
C--------------------------------------------------------------------------C
      FUNC = 0.0D0
      DO 10 I=1,NP
         TAU = X(I)
         FUNC = FUNC + ( F2(I)-K2(TAU,MN,LAM1,E1)     )**2/DF2(I)**2
     1               + ( F1(I)-K1(TAU,MN,LAM1,E3)     )**2/DF1(I)**2
     2               + ( F4(I)-K4(TAU,MN,LAM2,E2)     )**2/DF4(I)**2
     3               + ( F3(I)-K3(TAU,MN,LAM2)        )**2/DF3(I)**2
     4               + ( F6(I)-K6(TAU,MN,LAM1,LAM2)   )**2/DF6(I)**2
     5               + ( F5(I)-K5(TAU,MN,LAM1,LAM2,E4))**2/DF5(I)**2
 10   CONTINUE
      RETURN
      END

      FUNCTION K1(TAU,MN,LAM1,E0)
C-------------------------------------------------------------------------C
C     MODEL FOR FIRST CORRELATION FUNCTION                                C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,GG,M0,MS
      PI = 3.1415926D0
      K1 = LAM1**2 * D(MN,TAU) * MN - QQ/(4.D0*PI**2)*CONTP(TAU,E0)
      K20= 24.D0/(PI**6*TAU**9)
      K1 = K1/K20
      END

      FUNCTION K2(TAU,MN,LAM1,E0)
C--------------------------------------------------------------------------C
C     MODEL FOR SECOND CORRELATION FUNCTION                                C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K2 = LAM1**2 * DP(MN,TAU) + 1.D0/(64.D0*PI**4)*CONT(TAU,E0)
      K20= 24.D0/(PI**6*TAU**9)
      K2 = K2/K20
      END

      FUNCTION K3(TAU,MN,LAM2)
C-------------------------------------------------------------------------C
C     MODEL FOR THIRD CORRELATION FUNCTION                                C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K3 = LAM2**2 * D(MN,TAU) * MN
      K20=  24.D0/(PI**6*TAU**9)
      K40= 144.D0/(PI**6*TAU**9)
      K3 = K3/K40
      END

      FUNCTION K4(TAU,MN,LAM2,E0)
C--------------------------------------------------------------------------C
C     MODEL FOR FORTH  CORRELATION FUNCTION                                C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K4 = LAM2**2 * DP(MN,TAU) + 6.D0/(64.D0*PI**4)*CONT(TAU,E0)
      K20=  24.D0/(PI**6*TAU**9)
      K40= 144.D0/(PI**6*TAU**9)
      K4 = K4/K40
      END

      FUNCTION K5(TAU,MN,LAM1,LAM2,E0)
C-------------------------------------------------------------------------C
C     MODEL FOR FIFTH CORRELATION FUNCTION                                C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,GG,M0,MS
      PI = 3.1415926D0
      K5 = LAM1*LAM2*D(MN,TAU)*MN - 3.D0*QQ/(2.D0*PI**2)*CONTP(TAU,E0)
      K20= 24.D0/(PI**6*TAU**9)
      K5 = K5/K20
      END

      FUNCTION K6(TAU,MN,LAM1,LAM2)
C--------------------------------------------------------------------------C
C     MODEL FOR SIXTH  CORRELATION FUNCTION                                C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K6 = LAM1*LAM2 * DP(MN,TAU) 
      K20= 24.D0/(PI**6*TAU**9)
      K6 = K6/K20
      END

      FUNCTION D(M,TAU)
C-------------------------------------------------------------------------C
C     EQUAL TIME PROPAGATOR FOR A MASSIVE PARTICLE                        C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      X  = M*TAU
      IF (X .GT. 0.D0 .AND. X .LT. 80.D0) THEN
C-------------------------------------------------------------------------C
C     IMSL OR CERNLIB BESSEL FUNCTIONS                                    C
C-------------------------------------------------------------------------C
C        D = M/(4.D0*PI**2*TAU) * DBSK1(M*TAU)
         D = M/(4.D0*PI**2*TAU) * DBESK1(M*TAU)
             ELSE 
         D = 0.0D0
      ENDIF 
      END

      FUNCTION DP(M,TAU)
C------------------------------------------------------------------------C
C     DERIVATIVE OF EQUAL TIME PROPAGATOR                                C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      REAL*8 K(0:2)
      PI= 3.1415926D0
      X = M*TAU
      IF (X .GT. 0.D0 .AND. X .LT. 80.D0) THEN
C-------------------------------------------------------------------------C
C     IMSL OR CERNLIB BESSEL FUNCTIONS                                    C
C-------------------------------------------------------------------------C
C        CALL DBSKS(0.D0,X,3,K)
         K(0) = DBESK0(X)
         K(1) = DBESK1(X)
         K(2) = K(0)+2.0/X*K(1)
         DP = M**2/(4.D0*PI**2*TAU) * K(2)
      ELSE
         DP = 0.0D0
      ENDIF
      END

      FUNCTION CONT(T,E)
C------------------------------------------------------------------------C
C     EVALUATE CONTINUUM INTEGRAL                                        C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER IRULE
      EXTERNAL FF
      COMMON /PAR/ TAU
C------------------------------------------------------------------------C
      TAU= T
      A  = E
      B  = 20.D0*E
      ERRABS = 1.D-5
      ERRREL = 1.D-3
      IRULE  = 2
C------------------------------------------------------------------------C
C     IMSL INTEGRATION                                                   C
C------------------------------------------------------------------------C
C     CALL DQDAG(FF,A,B,ERRABS,ERRREL,IRULE,RESULT,ERREST)
C------------------------------------------------------------------------C
C     CERNLIB INTEGRATION                                                C
C------------------------------------------------------------------------C
      RESULT = DGAUSS(FF,A,B,ERRABS)
C------------------------------------------------------------------------C
C     ADD ASYMPTOTIC PART                                                C
C------------------------------------------------------------------------C
      PI   = 3.1415926D0
      CONT = RESULT + 1.D0/SQRT(8.D0) * PI**(-1.5D0) * TAU**(-2.5D0)
     1           * ABS(B)**6.5D0 * EXP(-ABS(B)*TAU)
      END

      FUNCTION FF(E)
C------------------------------------------------------------------------C
C     INTEGRAND FOR CONTIUUM PART                                        C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /PAR/ TAU
      FF= 2.D0 * E**5 * DP(E,TAU) 
      END

      FUNCTION CONTP(T,E)
C------------------------------------------------------------------------C
C     EVALUATE QQ CONTINUUM                                              C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER IRULE
      EXTERNAL FP
      COMMON /PAR/ TAU
C------------------------------------------------------------------------C
      TAU= T
      A  = E
      B  = 20.D0*E
      ERRABS = 1.D-5
      ERRREL = 1.D-3
      IRULE  = 2
C------------------------------------------------------------------------C
C     IMSL INTEGRATION                                                   C
C------------------------------------------------------------------------C
C     CALL DQDAG(FP,A,B,ERRABS,ERRREL,IRULE,RESULT,ERREST)
C------------------------------------------------------------------------C
C     CERNLIB INTEGRATION                                                C
C------------------------------------------------------------------------C
      RESULT = DGAUSS(FP,A,B,ERRABS)

      CONTP  = RESULT
      END

      FUNCTION FP(E)
C------------------------------------------------------------------------C
C     INTEGRAND FOR QQ CONTINUUM                                         C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /PAR/ TAU
      FP= 2.D0 * E**3 * D(E,TAU) 
      END


      FUNCTION TH1(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR FIRST CORRELATION FUNCTION                      C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,G2,M0,MS
      PI = 3.1415926D0
      TH1=-PI**2/12.D0*QQ * TAU**3 - PI**6*TAU**9/216.D0*QQ**3
      END

      FUNCTION TH2(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR SECOND CORRELATION FUNCTION                     C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,G2,M0,MS
      PI = 3.1415926D0
      C4 = PI**2/192.D0*G2
      C6 = PI**4/72.D0*QQ**2
      C8 =-PI**4/(72.D0*32.D0)*M0**2*QQ**2
      TH2= 1.D0 + C4*TAU**4 + C6*TAU**6 + C8*TAU**8
      END

      FUNCTION TH3(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR THIRD CORRELATION FUNCTION                      C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,G2,M0,MS
      PI = 3.1415926D0
      TH3=-PI**6/216.D0*QQ**3 * TAU**9
      END

      FUNCTION TH4(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR FORTH CORRELATION FUNCTION                      C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,G2,M0,MS
      PI = 3.1415926D0
      TH4= 1.0D0
      END

      FUNCTION TH5(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR FITH CORRELATION FUNCTION                       C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,G2,M0,MS
      PI = 3.1415926D0
      TH5=-PI**2/2.D0*QQ * TAU**3
      END

      FUNCTION TH6(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR SIXTH CORRELATION FUNCTION                      C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,SS,G2,M0,MS
      PI = 3.1415926D0
      TH6= PI**4/12.D0*QQ**2 * TAU**6
      END





