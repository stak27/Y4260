      implicit real*8(a-h,o-z)

      call BCODIN(52)

! spin of quark
      sq = 0.5d0
! 3P0
      aL12 = 0
      as12 = 1
      aj12 = 1
! 3S1
      aL34 = 1
      as34 = 1
      aj34 = 0

      totj = 1

! orbital max
      almax = 1


      sum0 = 0
      suma = 0


      do 10 alr = 0d0,almax,1d0
! relative orbital ang mom relative (13)-(24)

      do 10 al13 = 0d0,almax,1d0
      do 10 al24 = 0d0,almax,1d0

      do 10 as13 = 0d0,1d0,1d0
      do 10 as24 = 0d0,1d0,1d0

      do 10 aj13 = abs(al13-as13),al13+as13,1d0
      do 10 aj24 = abs(al24-as24),al24+as24,1d0

      do 10 ajh = abs(totj-alr),totj+alr,1d0
! total spin of final hadron sector

      sum = 0
      do 20 tots = 0d0,2d0,1d0
! total spin of quarks
      do 20 totl = abs(al12-al34),al12+al34,1d0
! total ang mom of initial hadrons
      do 20 alh = abs(al13-al24),al13+al24,1d0
! total ang mom of final hadrons w/o relative

        if(isankaku(alh,alr,totl).eq.0)then
!        write(6,*)alh,alr,totl
        go to 20
        endif

      fac1 = f9j(al12,as12,aj12,al34,as34,aj34,totl,tots,totj)
      fac2 = f9j(sq,sq,as12, sq,sq,as34, as13,as24,tots)
      fac3 = f9j(alh,alr,totl,tots,0d0,tots, ajh,alr,totj)

      fac4 = f9j(al13,al24,alh,as13,as24,tots,aj13,aj24,ajh)

      fac  = fac1*fac2*fac3*fac4
      sum = sum + fac
      if(abs(fac).gt.1d-8)
     & write(6,601)al13,as13,aj13,al24,as24,aj24,ajh,alr,totl,tots,alh,
     & fac,fac**2
601   format(3f5.1,3x,3f5.1,3x,2f5.1,3x,3f5.1,1p2e20.10)
20    continue
      if(abs(sum).gt.1d-8) then
      write(6,600)al13,as13,aj13,al24,as24,aj24,ajh,alr,sum,sum**2
      write(8,600)al13,as13,aj13,al24,as24,aj24,ajh,alr,sum,sum**2
600   format(3f5.1,3x,3f5.1,3x,2f5.1,18x,1p2e20.10)
      endif
      sum0 = sum0 + sum
      suma = suma + sum*sum

10    continue
      write(6,*)sum0,suma
!      write(6,*)nint(0.8d0),nint(-0.9d0)
      call check
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine check
      implicit real*8(a-h,j,l,o-z)
! do again

      hlf = 0.5d0

      l12 = 1
      s12 = 1
      j12 = 0
      l34 = 0
      s34 = 1
      j34 = 1

      l12 = 0
      s12 = 1
      j12 = 1
      l34 = 1
      s34 = 1
      j34 = 0

      l12 = 1
      s12 = 1
      j12 = 1
      l34 = 0
      s34 = 1
      j34 = 1

      l12 = 0
      s12 = 1
      j12 = 1
      l34 = 1
      s34 = 1
      j34 = 1

      l12 = 1
      s12 = 1
      j12 = 2
      l34 = 0
      s34 = 1
      j34 = 1

      l12 = 0
      s12 = 1
      j12 = 1
      l34 = 1
      s34 = 1
      j34 = 2

      l12 = 0
      s12 = 0
      j12 = 0
      l34 = 1
      s34 = 1
      j34 = 1

      l12 = 1
      s12 = 1
      j12 = 1
      l34 = 0
      s34 = 0
      j34 = 0

      l12 = 0
      s12 = 1
      j12 = 1
      l34 = 1
      s34 = 0
      j34 = 1

      l12 = 1
      s12 = 0
      j12 = 1
      l34 = 0
      s34 = 1
      j34 = 1

      l12 = 0
      s12 = 0
      j12 = 0
      l34 = 1
      s34 = 0
      j34 = 1

      l12 = 1
      s12 = 0
      j12 = 1
      l34 = 0
      s34 = 0
      j34 = 0

      jtot = 1
      
      do 30 lr =0d0,1d0,1d0

      do 30 l14 =0d0,1d0,1d0
      do 30 l32 =0d0,1d0,1d0
      if(abs(lr+l14+l32-1).lt.1d-3) then

      sum0 = 0
      do 10 s14=0d0,1d0,1d0
      do 10 s32=0d0,1d0,1d0

      do 10 j14=0d0,2d0,1d0
      do 10 j32=0d0,2d0,1d0

      do 10 jf = 0d0,3d0,1d0
      sum = 0
      
      do 20 li=0d0,2d0,1d0
      do 20 lf=0d0,2d0,1d0
      
      do 20 si = 0d0,2d0,1d0
      
      fac1 = fugo(s34+s32)
      fac2 = f9j(l12,s12,j12,l34,s34,j34,li,si,jtot)
      fac3 = f9j(hlf,hlf,s12,hlf,hlf,s34,s14,s32,si)
      fac4 = f9j(lf,lr,li,si,0d0,si,jf,lr,jtot)
      fac5 = f9j(l14,s14,j14,l32,s32,j32,lf,si,jf)
      
      sum = sum + fac1*fac2*fac3*fac4*fac5
20    continue
      sum0 = sum0+sum**2

      if(abs(sum).gt.1d-8) then
      write(6,600)l14,s14,j14,l32,s32,j32,jf,lr,sum,sum**2
      write(8,600)l14,s14,j14,l32,s32,j32,jf,lr,sum,sum**2
600   format(3f5.1,3x,3f5.1,3x,2f5.1,18x,1p2e20.10)
      endif
      
10    continue
      write(6,*)sum0
      end if

30    continue
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function f9j(a1,a2,a3,a4,a5,a6,a7,a8,a9)
      implicit real*8(a-h,o-z)

      j1 = nint(a1*2)
      j2 = nint(a2*2)
      j3 = nint(a3*2)
      j4 = nint(a4*2)
      j5 = nint(a5*2)
      j6 = nint(a6*2)
      j7 = nint(a7*2)
      j8 = nint(a8*2)
      j9 = nint(a9*2)
      
      fac = (j3+1)*(j6+1)*(j7+1)*(j8+1)
      d9 = d9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
      f9j = sqrt(fac) * d9
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function isankaku(a,b,c)
      implicit real*8(a-h,o-z)

      isankaku = 1
      if(nint(a+b).lt.nint(c)) then
        isankaku = 0
        return
      endif
      if(nint(b+c).lt.nint(a)) then
        isankaku = 0
        return
      endif
      if(nint(c+a).lt.nint(b)) then
        isankaku = 0
        return
      endif
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function fugo(a)
      implicit real*8(a-h,o-z)

      fugo = 1
      ia = nint(a)
      if(abs(mod(ia,2)).eq.1) fugo = -1

      end
*    D3J.FOR FOR KOHCHAN LIBARARY
*    OCT 19 1989
*    MINOR REVISION FOR IMPROVING PORTABILITY BY S ICHII
      FUNCTION D9J(J11A,J12A,J13A,J21A,J22A,J23A,J31A,J32A,J33A)
      IMPLICIT REAL*8 (A-H,O-Z)
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
C** ** FAST?  9J COEF.   WRITTEN  BY  HIDETOSHI NISHIMURA
C** ** REVISED BY M. OKA FOR THE VOS3 SYSTEM (AUGUST 1980)
C** ** ARGUMENT  INTEGER*4 (TWICE)
C---- ---- ---- ---- ---- ---- --- REVISED BY S.TAKEUCHI (APR.1988) -- -
      PARAMETER (IERUNI=6)
      DIMENSION RF(102)
      REAL*4 A,B,C,D,E,F,G,H,T
*      REAL*8 A,B,C,D,E,F,G,H,T
C FACOM
C     DATA IEVEN/Z80000001/
C HITAC
C     DATA IEVEN/'80000001'Z/
C VAX
      DATA IEVEN/'80000001'X/
C UNIX
c      DATA IEVEN/X'80000001'/
c#define IOR OR
c#define IAND AND

      COMMON/BCODCM/PP(52,52)

C STATEMENT FUNCTIONS
      Q(K,KK)=PP(KK+1,K+1)
      P(K,KK)=PP(K+1,KK+1)
      IE(K)=IAND(K,IEVEN)
      R(K)=RF(K+1)
      LPM(L1,L2,L3)=IOR(IOR(L1+L2-L3,L1-L2+L3),-L1+L2+L3)
      DELT(L1,L2,L3)=P(L1,(L1+L2+L3)/2-L3)*P((L1+L2+L3)/2,L1)*R((L1+L2+
     $L3)/2)
C     IFUGO IS NOW A STATEMENT FUNCTION
      IFUGO(KDUMMY)=1-2*MOD(ABS(KDUMMY),2)

      J11 = J11A
      J12 = J12A
      J13 = J13A
      J21 = J21A
      J22 = J22A
      J23 = J23A
      J31 = J31A
      J32 = J32A
      J33 = J33A
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
 10   CONTINUE
      IF(ABS(RF(1)-1D0).GT.1D-8) THEN
      X=0.D0
      DO 20 I=1,102
      X=X+1.D0
      RF(I)=1.D0/SQRT(X)
 20   CONTINUE
      IF(ABS(Q(12,2)-66.D0).GT.1D-8) CALL BCODIN(52)
      END IF

      IEE=IE(IOR(IOR(IOR(IOR(IOR(LPM(J11,J12,J13),
c      IEE=IE(OR(OR(OR(OR(OR(LPM(J11,J12,J13),
     $                           LPM(J21,J22,J23)),
     $                           LPM(J31,J32,J33)),
     $                           LPM(J11,J21,J31)),
     $                           LPM(J12,J22,J32)),
     $                           LPM(J13,J23,J33)))
      IF(IEE.NE.0) THEN
      D9J=0D0
      RETURN
      END IF

      MAXJ=J11+J21+J31
      MAXJ=MAX(MAXJ,J12+J22+J32)
      MAXJ=MAX(MAXJ,J13+J23+J33)
      MAXJ=MAX(MAXJ,J11+J12+J13)
      MAXJ=MAX(MAXJ,J21+J22+J23)
      MAXJ=MAX(MAXJ,J31+J32+J33)
      IF(MAXJ.GT.100)THEN
      WRITE(IERUNI,601)J11,J12,J13,J21,J22,J23,J31,J32,J33
  601 FORMAT(' ARG. RANGE OVER IN  (D9J) (',2(I5,'/2,'),I5,'/2:',
     & 2(I5,'/2,'),I5,'/2:',2(I5,'/2,'),I5,'/2);;')
      D9J=0.D0
      RETURN
      END IF

   30 CONTINUE
      JS=MAX(ABS(J21-J32),ABS(J11-J33),ABS(J12-J23))
      JE=MIN(J21+J32,J11+J33,J12+J23)
      JJN=(JE-JS)/2
      JT1=(J11+J21+J32+J33)/2+1
      JT2=(J32+J12+J23+J21)/2+1
      JT3=(J23+J33+J11+J12)/2+1
      JU1=(J11+J33+JS)/2
      JU2=(J32+J21+JS)/2
      JU3=(J23+J12+JS)/2
      JV1=JU1-JS
      JV2=JU2-JS
      JV3=JU3-JS
      JW1=JU1-J11
      JW2=JU2-J32
      JW3=JU3-J23
      JX1=JU1-J33
      JX2=JU2-J21
      JX3=JU3-J12
      JY1=(-J11-J32+J31+JS)/2
      JY2=(-J32-J23+J22+JS)/2
      JY3=(-J23-J11+J13+JS)/2
      JZ1=(-J21-J33+J31+JS)/2
      JZ2=(-J12-J21+J22+JS)/2
      JZ3=(-J33-J12+J13+JS)/2

   40 CONTINUE
      Y=0.D0
      AJS=JS+1
      DO  50 JJ=0,JJN
      K1S=MAX(0,-JY1-JJ,-JZ1-JJ)
      K2S=MAX(0,-JY2-JJ,-JZ2-JJ)
      K3S=MAX(0,-JY3-JJ,-JZ3-JJ)
      K1E=MIN(JT1-JU1-JJ-1,JW1-JY1,JX1-JZ1,JV1-JJ)
      K2E=MIN(JT2-JU2-JJ-1,JW2-JY2,JX2-JZ2,JV2-JJ)
      K3E=MIN(JT3-JU3-JJ-1,JW3-JY3,JX3-JZ3,JV3-JJ)
      X1=0.D0
      DO 60 K1=K1S,K1E
      X1=-X1+Q(JT1-K1,JU1+JJ+1)*Q(JV1-JJ,K1)*Q(JW1+JJ,JY1+JJ+K1)*
     $Q(JX1+JJ,JZ1+JJ+K1)
 60   CONTINUE
      X2=0.D0
      DO 70 K2=K2S,K2E
      X2=-X2+Q(JT2-K2,JU2+JJ+1)*Q(JV2-JJ,K2)*Q(JW2+JJ,JY2+JJ+K2)*
     $Q(JX2+JJ,JZ2+JJ+K2)
 70   CONTINUE
      X3=0.D0
      DO 80 K3=K3S,K3E
      X3=-X3+Q(JT3-K3,JU3+JJ+1)*Q(JV3-JJ,K3)*Q(JW3+JJ,JY3+JJ+K3)*
     $Q(JX3+JJ,JZ3+JJ+K3)
 80   CONTINUE
      Y=Y+X1*X2*X3*AJS*IFUGO(K1E+K2E+K3E)
      AJS=AJS+2D0
 50   CONTINUE

      D9J=Y*
     $DELT(J11,J12,J13)*
     $DELT(J21,J22,J23)*
     $DELT(J31,J32,J33)*
     $DELT(J11,J21,J31)*
     $DELT(J12,J22,J32)*
     $DELT(J13,J23,J33)

      RETURN
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      ENTRY DNINJ(A,B,C,D,E,F,G,H,T)
      J11=NINT(2.*A)
      J12=NINT(2.*B)
      J13=NINT(2.*C)
      J21=NINT(2.*D)
      J22=NINT(2.*E)
      J23=NINT(2.*F)
      J31=NINT(2.*G)
      J32=NINT(2.*H)
      J33=NINT(2.*T)
      GO TO 10

      END
*    BCODIN.FOR FOR KOHCHAN LIBRARY
*    OCT 19 1989
*    MINOR REVISION FOR IMPROVING PORTABILITY BY S ICHII
      SUBROUTINE BCODI2(NM)
*      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IERUNI=6)
*      COMMON /BCODCM/Q(1)
*    ORIGINAL DEF CAUSES COMPILER WARNING ON SOME IMPLIMENTATION (SI)
      COMMON /BCODCM/Q(52*52)
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
C   ** MAIN ROUTINE DE HITSUYO NA COMMON   CODED BY M.OKA
C   **  COMMON /BCODCM/Q(NM,NM)
C   ** Q(M,I) (I.GE.M) BINOMIAL COEF. C(I-1,M-1)
C   ** Q(I,M) (I.GE.M) Q(M,I)**(-1/2)
C   ** NM>52 DEWA Q GA KETAOCHI SURU KANOUSEI GA ARIMASU.
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      DATA NMAX/0/
      NA(I,J)=I+(J-1)*NM

      IF(NM.LT.3) THEN
      WRITE(IERUNI,100)NM
  100 FORMAT('NMAX.LT.3 IN (BCODIN)  NMAX=',I4,';')
      RETURN
      END IF

      IF(NMAX.EQ.NM) RETURN

      IF(NMAX.NE.0) THEN
      WRITE(IERUNI,101)NMAX,NM
  101 FORMAT('/BCODCM/ LENGTH MISMATCHED  OLD NMAX=',I4,', NEW NMAX=',
     +I4,';')
      RETURN
      END IF

      NMAX=NM
      Q(1)=1.D0
      Q(2)=1.D0
      Q(NM+1)=1.D0
      Q(NM+2)=1.D0
      DO 10 I=3,NM
      Q(NA(I,1))=1.D0
      Q(NA(I,I))=1.D0
      Q(NA(1,I))=1.D0
      DO 10 J=2,I-1
      Q(NA(J,I))=Q(NA(J-1,I-1))+Q(NA(J,I-1))
   10 CONTINUE

      DO 20 J=3,NM
      DO 20 I=1,J-1
      Q(NA(J,I))=1.D0/SQRT(Q(NA(I,J)))
   20 CONTINUE

      END
C==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      SUBROUTINE BCODIN(NM)
*      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      COMMON /BCODCM/ Q(52,52)
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      IF(NM.NE.52) THEN
   99 CALL BCODI2(NM)
      RETURN
      END IF
C---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      Q(1,1)=1.0D0
      Q(2,1)=1.0D0
      Q(1,2)=1.0D0
      Q(2,2)=1.0D0
      DO 10 I=3,52
      Q(I,1)=1.0D0
      Q(1,I)=1.0D0
      Q(I,I)=1.0D0
      DO 10 J=2,I-1
      Q(J,I)=Q(J-1,I-1)+Q(J,I-1)
   10 CONTINUE

      DO 20 J=3,52
      DO 20 I=1,J-1
      Q(J,I)=1.0D0/SQRT(Q(I,J))
   20 CONTINUE

      END
