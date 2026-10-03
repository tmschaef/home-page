      PROGRAM MAIN
C--------------------------------------------------------------------------C
C     FIT EQUAL TIME MESON CORRELATION FUNCTIONS WITH SIMPLE POLE PLUS     C
C     CONTIUUM MODELS.                                                     C
C--------------------------------------------------------------------------C
C     NOTE THAT FPI IS DETERMINED IN THE A1 CHANNEL ONLY WHEREAS LPI IS    C
C     FIT INDEPENDENTLY IN THE PI CHANNEL (THE TWO ARE REALTED BY PCAC,    C
C     BUT ONE WOULD HAVE TO KNOW THE VALUE OF THE CONDENSATE IN THIS       C
C     SIMULATION.)                                                         C
C--------------------------------------------------------------------------C
C     VERSION            :    1.0                                          C
C     CREATION DATE      :  07-06-93                                       C
C     LAST MODIFICATION  :  12-06-95                                       C
C--------------------------------------------------------------------------C
C     THE CORRELATION FUNCTIONS ARE DEFINED AS FOLLOWS:                    C
C     F1     SCALAR CHANNEL                                                C
C     F2     PSEUDOSCALAR CHANNEL                                          C
C     F3     AXIAL VECTOR CHANNEL                                          C
C     F4     VECTOR CHANNEL                                                C
C--------------------------------------------------------------------------C
C     N      NUMBER OF POINTS IN THE DATA FILE                             C
C     NP     NUMBER OF POINTS USED IN THE FIT                              C
C     NPAR   NUMBER OF FIT PARAMETERS                                      C
C     IS     FIT STRATEGY                                                  C
C--------------------------------------------------------------------------C
C     INPUT:   MES.DAT                                                     C
C     TITLE                                                                C
C     X(INDEX) F1(INDEX) DF1(INDEX)                                        C
C      ...        ...       ...                                            C
C                                                                          C
C--------------------------------------------------------------------------C
C     OUTPUT:  FITMES.DAT                                                  C
C     TAU    FIT(TAU) DATA(TAU) DEL(TAU) OPE(TAU)                          C
C      ...    ...      ...       ...      ...                              C
C                                                                          C
C--------------------------------------------------------------------------C
C     PLOT FILE: MES.MGO                                                   C
C--------------------------------------------------------------------------C
C     THIS VERSION USES MINUIT MINIMIZATION AND ERROR ANNALYSIS.           C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      REAL*8 F1(100),F2(100),F3(100),F4(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100)
      REAL*8 EF1(100),EF2(100),EF3(100),EF4(100)
      REAL*8 X(100)
      REAL*8 ARGLIS(10)
      INTEGER I,J,K,N,NP,NPARI,NPARX,IS
      INTEGER IERFLG,IVAR,ISTAT
      CHARACTER*80 NAME(20)
      EXTERNAL CHISQU1,CHISQU2,CHISQU3,CHISQU4
      COMMON /DATA/ X,F1,F2,F3,F4,DF1,DF2,DF3,DF4,NP
      COMMON /OPE/  QQ,G2,M0
C--------------------------------------------------------------------------C
C     INPUT                                                                C
C--------------------------------------------------------------------------C
      N = 15
      OPEN(UNIT=15,FILE='mes.dat',STATUS='unknown')
      OPEN(UNIT=16,FILE='fitmes.dat',STATUS='unknown')
      WRITE(6,*) 'DETERMINE FIT WINDOW'
      WRITE(6,*) 'MINIMUM DISTANCE (LAM^-1)'
      READ(5,*)  X0
      WRITE(6,*) 'MAXIMUM DISTANCE'
      READ(5,*)  X1
      WRITE(6,*) 'DETERMINE FITTING PROCEDURE'
      WRITE(6,*) '  (1) STATISTICAL WEIGHT'
      WRITE(6,*) '  (2) MODIFIED STAT. WEIGHT'
      WRITE(6,*) '  (3) EQUAL WEIGHTS'
      READ(5,*) IS
      IF ( IS .EQ. 2 ) THEN
      WRITE(6,*) 'MINIMUM STATISTICAL ERROR'
      READ(5,*)  DEL
      ENDIF
      WRITE(6,*) 'SCALE PARAMETER'
      READ(5,*) CONV
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
C--------------------------------------------------------------------------C
C     DETERMINE STATISTIACL WEIGHT                                         C
C--------------------------------------------------------------------------C
      DO 15 I=1,NP
         IF (IS .EQ. 1) THEN
            DF1(I) = EF1(I)
            DF2(I) = EF2(I)
            DF3(I) = EF3(I)
            DF4(I) = EF4(I)
         ELSE IF (IS .EQ. 2) THEN
            DF1(I) = MAX(EF1(I),DEL)
            DF2(I) = MAX(EF2(I),DEL)
            DF3(I) = MAX(EF3(I),DEL)
            DF4(I) = MAX(EF4(I),DEL)
         ELSE 
            DF1(I) = 1.D0
            DF2(I) = 1.D0
            DF3(I) = 1.D0
            DF4(I) = 1.D0
         ENDIF
 15   CONTINUE
C--------------------------------------------------------------------------C
C     VACUUM CONDENSATES                                                   C
C--------------------------------------------------------------------------C
      CONV2= CONV**2
      CONV3= CONV**3
      QQ   =-0.25D0**3/CONV3
      G2   = 0.36D0**4/CONV**4
      M0   = 0.70D0/CONV
C--------------------------------------------------------------------------C
C     INITIALIZE MINUIT                                                    C
C--------------------------------------------------------------------------C
      CALL MNINIT(5,6,7)
C--------------------------------------------------------------------------C
C     SYNTAX: NUMBER,NAME,START,ERROR,LBOUND,UBOUND,ERRORFLAG              C
C--------------------------------------------------------------------------C
C     START WITH SCALAR DELTA MESON                                        C
C--------------------------------------------------------------------------C
      CALL MNPARM(1,'MS',5.0D0,1.D0,0.D0,1.D2,IERFLG)
      CALL MNPARM(2,'LS',1.0D0,1.D0,0.D0,0.D0,IERFLG)
      CALL MNPARM(3,'E1',7.3D0,1.D0,0.D0,0.D0,IERFLG)
C--------------------------------------------------------------------------C
C     FIX SOME OF THE PARAMETERS:                                          C
C--------------------------------------------------------------------------C
C     ARGLIS(1) = 1
C     CALL MNEXCM(CHISQU,'FIX',ARGLIS,1,IERFLG,0)
C--------------------------------------------------------------------------C
C     MINIMIZATION AND ERROR ANALYSIS                                      C
C--------------------------------------------------------------------------C
      ARGLIS(1) = 0
      CALL MNEXCM(CHISQU1,'SET_PRINT',ARGLIS,1,IERFLG,0)
      CALL MNEXCM(CHISQU1,'MIGRAD',ARGLIS,0,IERFLG,0)
C--------------------------------------------------------------------------C
C     IF NO ERROR ANALYSIS IS NEEDED, DELETE NEXT LINE                     C
C--------------------------------------------------------------------------C
      CALL MNEXCM(CHISQU1,'MINOS',ARGLIS,0,IERFLG,0)
      CALL MNSTAT(FMIN1,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
C--------------------------------------------------------------------------C  
C     RESULTS                                                              C
C--------------------------------------------------------------------------C
      CALL MNPOUT(1,NAME(1),MS,DMS,BND1,BND1,IVAR)
      CALL MNPOUT(2,NAME(2),LS,DLS,BND1,BND2,IVAR)
      CALL MNPOUT(3,NAME(3),E1,DE1,BND1,BND2,IVAR)
C--------------------------------------------------------------------------C
C     ERROR ESTIMATES                                                      C
C--------------------------------------------------------------------------C
      CALL MNERRS(1,MSUP,MSLOW,MSPAR,MSCOR)
      CALL MNERRS(2,LSUP,LSLOW,LSPAR,LSCOR)
      CALL MNERRS(3,E1UP,E1LOW,E1PAR,E1COR)
C--------------------------------------------------------------------------C
C     PION                                                                 C
C--------------------------------------------------------------------------C
      CALL MNINIT(5,6,7)
      CALL MNPARM(1,'MP',1.0D0,1.D0,0.D0,1.D2,IERFLG)
      CALL MNPARM(2,'LP',2.5D0,1.D0,0.D0,0.D0,IERFLG)
      CALL MNPARM(3,'E2',8.5D0,1.D0,0.D0,0.D0,IERFLG)

      ARGLIS(1) = 0
      CALL MNEXCM(CHISQU2,'SET_PRINT',ARGLIS,1,IERFLG,0)
      CALL MNEXCM(CHISQU2,'MIGRAD',ARGLIS,0,IERFLG,0)
      CALL MNEXCM(CHISQU2,'MINOS',ARGLIS,0,IERFLG,0)
      CALL MNSTAT(FMIN2,FEDM,ERRDEF,NPARI,NPARX,ISTAT)

      CALL MNPOUT(1,NAME(4),MP,DMP,BND1,BND2,IVAR)
      CALL MNPOUT(2,NAME(5),LP,DLP,BND1,BND2,IVAR)
      CALL MNPOUT(3,NAME(6),E2,DE2,BND1,BND2,IVAR)

      CALL MNERRS(1,MPUP,MPLOW,MPPAR,MPCOR)
      CALL MNERRS(2,LPUP,LPLOW,LPPAR,LPCOR)
      CALL MNERRS(3,E2UP,E2LOW,E2PAR,E2COR)
C--------------------------------------------------------------------------C
C     A1 MESON                                                             C
C--------------------------------------------------------------------------C
      CALL MNINIT(5,6,7)
      CALL MNPARM(1,'MA',5.0D0,1.D0,2.D0,1.D2,IERFLG)
      CALL MNPARM(2,'GA',9.0D0,1.D0,0.D0,0.D0,IERFLG)
      CALL MNPARM(3,'E3',8.5D0,1.D0,0.D0,0.D0,IERFLG)
      CALL MNPARM(4,'FP',0.5D0,0.1D0,0.D0,0.D0,IERFLG)
      CALL MNPARM(5,'MP2',MP,0.1D0,0.D0,0.D0,IERFLG)
C--------------------------------------------------------------------------C
C     FIX PION MASS FROM ABOVE                                             C
C--------------------------------------------------------------------------C
      ARGLIS(1) = 5
      CALL MNEXCM(CHISQU1,'FIX',ARGLIS,1,IERFLG,0)

      ARGLIS(1) = 0
      CALL MNEXCM(CHISQU3,'SET_PRINT',ARGLIS,1,IERFLG,0)
      CALL MNEXCM(CHISQU3,'MIGRAD',ARGLIS,0,IERFLG,0)
      CALL MNEXCM(CHISQU3,'MINOS',ARGLIS,0,IERFLG,0)
      CALL MNSTAT(FMIN3,FEDM,ERRDEF,NPARI,NPARX,ISTAT)

      CALL MNPOUT(1,NAME(7),MA,DMA,BND1,BND2,IVAR)
      CALL MNPOUT(2,NAME(8),GA,DGA,BND1,BND2,IVAR)
      CALL MNPOUT(3,NAME(9),E3,DE3,BND1,BND2,IVAR)
      CALL MNPOUT(4,NAME(10),FP,DFP,BND1,BND2,IVAR)
      CALL MNPOUT(5,NAME(11),MP2,DMP2,BND1,BND2,IVAR)

      CALL MNERRS(1,MAUP,MALOW,MAPAR,MACOR)
      CALL MNERRS(2,GAUP,GALOW,GAPAR,GACOR)
      CALL MNERRS(3,E3UP,E3LOW,E3PAR,E3COR)
      CALL MNERRS(4,FPUP,FPLOW,FPPAR,FPCOR)
      CALL MNERRS(5,MP2UP,MP2LOW,MP2PAR,FPCOR)
C--------------------------------------------------------------------------C
C     RHO MESON                                                            C
C--------------------------------------------------------------------------C
      CALL MNINIT(5,6,7)
      CALL MNPARM(1,'MV',4.0D0,1.D0,0.D0,1.D2,IERFLG)
      CALL MNPARM(2,'GR',5.0D0,1.D0,0.D0,0.D0,IERFLG)
      CALL MNPARM(3,'E4',8.5D0,1.D0,0.D0,0.D0,IERFLG)

      ARGLIS(1) = 0
      CALL MNEXCM(CHISQU4,'SET_PRINT',ARGLIS,1,IERFLG,0)
      CALL MNEXCM(CHISQU4,'MIGRAD',ARGLIS,0,IERFLG,0)
      CALL MNEXCM(CHISQU4,'MINOS',ARGLIS,0,IERFLG,0)
      CALL MNSTAT(FMIN4,FEDM,ERRDEF,NPARI,NPARX,ISTAT)

      CALL MNPOUT(1,NAME(12),MV,DMV,BND1,BND2,IVAR)
      CALL MNPOUT(2,NAME(13),GV,DGV,BND1,BND2,IVAR)
      CALL MNPOUT(3,NAME(14),E4,DE4,BND1,BND2,IVAR)

      CALL MNERRS(1,MVUP,MVLOW,MVPAR,MVCOR)
      CALL MNERRS(2,GVUP,GVLOW,GVPAR,GVCOR)
      CALL MNERRS(3,E4UP,E4LOW,E4PAR,E4COR)

C--------------------------------------------------------------------------C
C     PHENOMENOLOGY                                                        C
C--------------------------------------------------------------------------C
 
      MSEX  = 0.90/CONV
      LSEX  = (0.0/CONV)**2
      E1EX  = 1.2/CONV
      MPEX  = 0.136/CONV
      LPEX  = (0.550/CONV)**2
      E2EX  = 1.30/CONV
      MAEX  = 1.23/CONV
      GAEX  = 9.1
      E3EX  = 1.40/CONV
      FPEX  = 0.093/CONV
      MVEX  = 0.770/CONV
      GVEX  = 5.44
      E4EX  = 1.40/CONV

C--------------------------------------------------------------------------C
C     COMPARE WITH DATA                                                    C
C--------------------------------------------------------------------------C

      WRITE(16,*) 'SCALAR CHANNEL (TAU,FIT,DATA,ERROR,OPE,PH)' 
      DO 100 I=1,NP   
         TAU = X(I)
         WRITE(16,666) TAU,K1(TAU,MS,LS,E1),F1(I),
     1    EF1(I),TH1(TAU),K1(TAU,MSEX,LSEX,E1EX)
 100  CONTINUE
      WRITE(16,*) 'PSEUDOSCALAR CHANNEL '
      DO 110 I=1,NP
         TAU = X(I)
         WRITE(16,666) TAU,K2(TAU,MP,LP,E2),F2(I),
     1    EF2(I),TH2(TAU),K2(TAU,MPEX,LPEX,E2EX)
 110  CONTINUE
      WRITE(16,*) 'AXIAL VECTOR CHANNEL'
      DO 120 I=1,NP
         TAU = X(I)
         WRITE(16,666) TAU,K3(TAU,MA,MP,GA,FP,E3),F3(I),
     1    EF3(I),TH3(TAU),K3(TAU,MAEX,MPEX,GAEX,FPEX,E3EX)
 120  CONTINUE
      WRITE(16,*) 'VECTOR CHANNEL'
      DO 130 I=1,NP
         TAU = X(I)
         WRITE(16,666) TAU,K4(TAU,MV,GV,E4),F4(I),
     1    EF4(I),TH4(TAU),K4(TAU,MVEX,GVEX,E4EX)
 130  CONTINUE

C--------------------------------------------------------------------------C
C     OUTPUT                                                               C
C--------------------------------------------------------------------------C

      FMIN = FMIN1+FMIN2+FMIN3+FMIN4
      DO 200 I=6,16,10
         WRITE(I,999)
         WRITE(I,998) NAME(1),MS*CONV,MSUP*CONV,MSLOW*CONV
         WRITE(I,998) NAME(2),LS*CONV2,LSUP*CONV2,LSLOW*CONV2
         WRITE(I,998) NAME(3),E1*CONV,E1UP*CONV,E1LOW*CONV
         WRITE(I,998) NAME(4),MP*CONV,MPUP*CONV,MPLOW*CONV
         WRITE(I,998) NAME(5),LP*CONV2,LPUP*CONV2,LPLOW*CONV2
         WRITE(I,998) NAME(6),E2*CONV,E2UP*CONV,E2LOW*CONV
         WRITE(I,998) NAME(7),MA*CONV,MAUP*CONV,MALOW*CONV 
         WRITE(I,998) NAME(8),GA,GAUP,GALOW
         WRITE(I,998) NAME(9),E3*CONV,E3UP*CONV,E3LOW*CONV
         WRITE(I,998) NAME(10),FP*CONV,FPUP*CONV,FPLOW*CONV
         WRITE(I,998) NAME(11),MP2*CONV,MP2UP*CONV,MP2LOW*CONV
         WRITE(I,998) NAME(12),MV*CONV,MVUP*CONV,MVLOW*CONV
         WRITE(I,998) NAME(13),GV,GVUP,GVLOW
         WRITE(I,998) NAME(14),E4*CONV,E4UP*CONV,E4LOW*CONV
         WRITE(I,999)
         WRITE(I,*)   'CHI2/NDF'
         WRITE(I,111) FMIN/(4.D0*NP-NPARI)
         WRITE(I,999)
 200  CONTINUE
C--------------------------------------------------------------------------C
 111  FORMAT(1X,F12.5)
 222  FORMAT(1X,F12.5,1X,F12.5)
 333  FORMAT(1X,F12.5,1X,F12.5,1X,F12.5)
 444  FORMAT(1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5)
 555  FORMAT(1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5)
 666  FORMAT(6(1X,F12.5))
 998  FORMAT(1X,A5,1X,F12.5,1X,F12.5,1X,F12.5)
 999  FORMAT(70('-'))
      END

      SUBROUTINE CHISQU1(NPAR,GRAD,FUNC,C,IFLAG)
C--------------------------------------------------------------------------C
C     EVALUATES CHISQU FOR FIT TO MODEL SPECTRUM                           C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER I,NP,NPAR,IS,IFLAG
      REAL*8 C(NPAR),GRAD(NPAR)
      REAL*8 X(100)
      REAL*8 F1(100),F2(100),F3(100),F4(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100)
      COMMON /DATA/ X,F1,F2,F3,F4,DF1,DF2,DF3,DF4,NP
C--------------------------------------------------------------------------C
      MS = C(1)
      LS = C(2)
      E1 = C(3)
C--------------------------------------------------------------------------C
      FUNC = 0.0D0
      DO 10 I=1,NP
         TAU = X(I)
         FUNC = FUNC + ( F1(I)-K1(TAU,MS,LS,E1) )**2/DF1(I)**2
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CHISQU2(NPAR,GRAD,FUNC,C,IFLAG)
C--------------------------------------------------------------------------C
C     EVALUATES CHISQU FOR FIT TO MODEL SPECTRUM                           C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER I,NP,NPAR,IS,IFLAG
      REAL*8 C(NPAR),GRAD(NPAR)
      REAL*8 X(100)
      REAL*8 F1(100),F2(100),F3(100),F4(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100)
      COMMON /DATA/ X,F1,F2,F3,F4,DF1,DF2,DF3,DF4,NP
C--------------------------------------------------------------------------C
      MP = C(1)
      LP = C(2)
      E2 = C(3)
C--------------------------------------------------------------------------C
      FUNC = 0.0D0
      DO 10 I=1,NP
         TAU = X(I)
         FUNC = FUNC + ( F2(I)-K2(TAU,MP,LP,E2) )**2/DF2(I)**2
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CHISQU3(NPAR,GRAD,FUNC,C,IFLAG)
C--------------------------------------------------------------------------C
C     EVALUATES CHISQU FOR FIT TO MODEL SPECTRUM                           C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER I,NP,NPAR,IS,IFLAG
      REAL*8 C(NPAR),GRAD(NPAR)
      REAL*8 X(100)
      REAL*8 F1(100),F2(100),F3(100),F4(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100)
      COMMON /DATA/ X,F1,F2,F3,F4,DF1,DF2,DF3,DF4,NP
C--------------------------------------------------------------------------C
      MA = C(1)
      GA = C(2)
      E3 = C(3)
      FP = C(4)
      MP = C(5)
C--------------------------------------------------------------------------C
      FUNC = 0.0D0
      DO 10 I=1,NP
         TAU = X(I)
         FUNC = FUNC + ( F3(I)-K3(TAU,MA,MP,GA,FP,E3))**2/DF3(I)**2
 10   CONTINUE
      RETURN
      END

      SUBROUTINE CHISQU4(NPAR,GRAD,FUNC,C,IFLAG)
C--------------------------------------------------------------------------C
C     EVALUATES CHISQU FOR FIT TO MODEL SPECTRUM                           C
C--------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER I,NP,NPAR,IS,IFLAG
      REAL*8 C(NPAR),GRAD(NPAR)
      REAL*8 X(100)
      REAL*8 F1(100),F2(100),F3(100),F4(100)
      REAL*8 DF1(100),DF2(100),DF3(100),DF4(100)
      COMMON /DATA/ X,F1,F2,F3,F4,DF1,DF2,DF3,DF4,NP
C--------------------------------------------------------------------------C
      MV = C(1)
      GV = C(2)
      E4 = C(3)
C--------------------------------------------------------------------------C
      FUNC = 0.0D0
      DO 10 I=1,NP
         TAU = X(I)
         FUNC = FUNC + ( F4(I)-K4(TAU,MV,GV,E4) )**2/DF4(I)**2
 10   CONTINUE
      RETURN
      END

      FUNCTION K1(TAU,MS,LS,E1)
C-------------------------------------------------------------------------C
C     MODEL FOR SCALAR (ISOVECTOR) MESON CORRELATION FUNCTION             C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K1 = 0.5*LS**2 * D(MS,TAU)  + 3.D0/(8.D0*PI**2) * CONT(TAU,E1)
      K10= 3.D0/(PI**4*TAU**6)
      K1 = K1/K10
      END

      FUNCTION K2(TAU,MP,LP,E2)
C-------------------------------------------------------------------------C
C     MODEL FOR PSEUDOSCALAR MESON CORRELATION FUNCTION                   C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K2 = 0.5*LP**2 * D(MP,TAU)  + 3.D0/(8.D0*PI**2) * CONT(TAU,E2)
      K20= 3.D0/(PI**4*TAU**6)
      K2 = K2/K20
      END

      FUNCTION K3(TAU,MA,MP,GA,FP,E3)
C-------------------------------------------------------------------------C
C     MODEL FOR AXIAL VECTOR MESON CORRELATION FUNCTION                   C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K3 = 6.D0*MA**4/GA**2 * D(MA,TAU) - 2.0*FP**2*MP**2 * D(MP,TAU)  
     1    + 6.D0/(8.D0*PI**2) * CONT(TAU,E3)
      K30= 6.D0/(PI**4*TAU**6)
      K3 = K3/K30
      END

      FUNCTION K4(TAU,MV,GV,E4)
C-------------------------------------------------------------------------C
C     MODEL FOR VECTOR MESON CORRELATION FUNCTION                         C
C-------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      PI = 3.1415926D0
      K4 = 6*MV**4/GV**2 * D(MV,TAU) + 6.D0/(8.D0*PI**2) * CONT(TAU,E4)
      K40= 6.D0/(PI**4*TAU**6)
      K4 = K4/K40
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

      FUNCTION CONT(T,E)
C------------------------------------------------------------------------C
C     EVALUATE CONTINUUM                                                 C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      INTEGER IRULE
      EXTERNAL FF
      COMMON /PAR/ TAU
C------------------------------------------------------------------------C
      TAU= T
      A  = E
      B  = 200.d0
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
      CONT   = RESULT
      END

      FUNCTION FF(E)
C------------------------------------------------------------------------C
C     INTEGRAND FOR CONTINUUM INTEGRAL                                   C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /PAR/ TAU
      FF = 2.D0 * E**3 * D(E,TAU) 
      END

      FUNCTION TH1(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR SCALAR MESON CORRELATOR                         C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,G2,M0
      PI = 3.1415926D0
      TH1= 1.D0 - PI**4*TAU**6/36.D0*QQ**2
      END

      FUNCTION TH2(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR PSEUDOSCALAR MESON CORRELATOR                   C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,G2,M0
      PI = 3.1415926D0
      TH2= 1.D0 + PI**4*TAU**6/36.D0*QQ**2
      END

      FUNCTION TH3(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR AXIAL VECTOR CORRELATOR                         C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,G2,M0
      PI = 3.1415926D0
      TH3= 1.D0 - PI**4*TAU**6/18.D0*QQ**2
      END

      FUNCTION TH4(TAU)
C------------------------------------------------------------------------C
C     OPE PREDICTION FOR VECTOR CORRELATOR                               C
C------------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-Z)
      COMMON /OPE/ QQ,G2,M0
      PI = 3.1415926D0
      TH4= 1.D0 + PI**4*TAU**6/18.D0*QQ**2
      END
