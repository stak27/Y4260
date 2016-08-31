      implicit real*8(a-h,o-z)

      call BCODIN(52)


      call check
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine check
      implicit real*8(a-h,j,l,o-z)
      character*100 af
      character*23 a1
      character*23 aa(1000,18)
! do again
      dimension iqni(9),ihit(1000)
      hlf = 0.5d0

      ihit = 0
! total j = 1
      do 100 ic1 = 1,2
      
      do 100 ic = 1,18
      if(ic.eq.1) then
! 1s0 1p1  l_rela=0
      iqni=(/0,0,0, 1,0,1, 1,0,1/)
      elseif(ic.eq.2) then
! 1p1 1s0  l_rela=0
      iqni=(/1,0,1, 0,0,0, 1,0,1/)
      elseif(ic.eq.3) then
! 1s0 3p1  l_rela=0
      iqni=(/0,0,0, 1,1,1, 1,0,1/)
      elseif(ic.eq.4) then
! 3p1 1s0  l_rela=0
      iqni=(/1,1,1, 0,0,0, 1,0,1/)

      elseif(ic.eq.5) then
! 3s1 1p1  l_rela=0
      iqni=(/0,1,1, 1,0,1, 1,0,1/)
      elseif(ic.eq.6) then
! 1p1 3s1  l_rela=0
      iqni=(/1,0,1, 0,1,1, 1,0,1/)
      elseif(ic.eq.7) then
! 3s1 3p0  l_rela=0
      iqni=(/0,1,1, 1,1,0, 1,0,1/)
      elseif(ic.eq.8) then
! 3p0 3s1  l_rela=0
      iqni=(/1,1,0, 0,1,1, 1,0,1/)
      elseif(ic.eq.9) then
! 3s1 3p1  l_rela=0
      iqni=(/0,1,1, 1,1,1, 1,0,1/)
      elseif(ic.eq.10) then
! 3p1 3s1  l_rela=0
      iqni=(/1,1,1, 0,1,1, 1,0,1/)
      elseif(ic.eq.11) then
! 3s1 3p2  l_rela=0
      iqni=(/0,1,1, 1,1,2, 1,0,1/)
      elseif(ic.eq.12) then
! 3p2 3s1  l_rela=0
      iqni=(/1,1,2, 0,1,1, 1,0,1/)

      elseif(ic.eq.13) then
! 1s0 1s0  l_rela=1
      iqni=(/0,0,0, 0,0,0, 0,1,1/)
      elseif(ic.eq.14) then
! 1s0 3s1  l_rela=1
      iqni=(/0,0,0, 0,1,1, 1,1,1/)
      elseif(ic.eq.15) then
! 3s1 1s0  l_rela=1
      iqni=(/0,1,1, 0,0,0, 1,1,1/)
      elseif(ic.eq.16) then
! 3s1 3s1 ji=0 l_rela=1
      iqni=(/0,1,1, 0,1,1, 0,1,1/)
      elseif(ic.eq.17) then
! 3s1 3s1 ji=1 l_rela=1
      iqni=(/0,1,1, 0,1,1, 1,1,1/)
      elseif(ic.eq.18) then
! 3s1 3s1 ji=2 l_rela=1
      iqni=(/0,1,1, 0,1,1, 2,1,1/)
      endif

      l12 = iqni(1)
      s12 = iqni(2)
      j12 = iqni(3)
      l34 = iqni(4)
      s34 = iqni(5)
      j34 = iqni(6)
      ji  = iqni(7)
      lri = iqni(8)
      jtot = iqni(9)

      write(6,'(3i3,3x,3i3,3x,3i3)')(iqni(i),i=1,9)
      write(8,'(3i3,3x,3i3,3x,3i3)')(iqni(i),i=1,9)
      write(18,'(3x,3i2,1x,3i2,1x,3i2)')(iqni(i),i=1,9)


      ii=0
      do 30 lrf =0d0,1d0,1d0

      sum0 = 0
      do 15 l14 =0d0,1d0,1d0
      
      do 15 s14=0d0,1d0,1d0
      do 15 j14=0d0,2d0,1d0

      do 15 l32 =0d0,1d0,1d0
      if(abs(lrf+l14+l32-1).lt.1d-3) then
      do 10 s32=0d0,1d0,1d0
      do 10 j32=0d0,2d0,1d0

      do 10 jf = 0d0,3d0,1d0
      sum = 0
      
      do 20 li=0d0,2d0,1d0
      do 20 lf=0d0,2d0,1d0
      do 20 ll=0d0,2d0,1d0
      
      do 20 si = 0d0,2d0,1d0
      
      fac1 = fugo(s34+s32)
      fac2 = f9j(l12,s12,j12,l34,s34,j34,li,si,ji)
      fac3 = f9j(hlf,hlf,s12,hlf,hlf,s34,s14,s32,si)
      fac4 = f9j(l14,s14,j14,l32,s32,j32,lf,si,jf)
      fac5 = f9j(li,lri,ll,si,0d0,si,ji,lri,jtot)
      fac6 = f9j(lf,lrf,ll,si,0d0,si,jf,lrf,jtot)
      
      sum = sum + fac1*fac2*fac3*fac4*fac5*fac6
20    continue
      sum0 = sum0+sum**2
      ii=ii+1

      




      if(ic1.eq.1 .and. abs(sum).gt.1d-8) then
      ihit(ii) = ihit(ii)+1
      end if
      if(ic1.eq.2 .and. ihit(ii).gt.0) then
        call aform(sum,af)

      write(6,600)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & trim(af),sum,sum**2
      write(8,600)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & trim(af),sum,sum**2
600   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a20,1p2e20.10)
      write(18,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & '$',trim(af),'$&',sum,sum**2
610   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a,a20,a,1pe18.8,e12.4)
      write(a1,'(a,a20,a)')'$',trim(af),'$&'
      aa(ii,ic) = a1
      endif
      
10    continue
      end if
15    continue

      write(6,*)sum0
30    continue

100   continue

      do ii=1,1000
      if(ihit(ii).gt.0) then
      write(28,'(18a)')(aa(ii,ic),ic=1,18)
      endif
      end do

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
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine aform(a,af)
      implicit real*8(a-h,o-z)
      character*1 ap
      character*100 af,ad
      eps = 1d-8
      jmx = 1024

      iflatex = 1
      
      do 20 j=1,jmx+1


      if(abs(a).lt.eps) then
!        write(6,'(4x,i6,a,i5,1p2e18.10)')0,'/',1,a,0d0
        write(af,'(i6)') 0
        go to 21
      endif


        ia = nint(a)
      
        if(abs((ia-a)).lt.eps) then
!        write(6,'(4x,i6,1p2e18.10)')ia,a
        write(af,'(i6)') ia
        go to 21
        endif

      if(a.gt.0) then
      ap = ' '
      else
      ap = '-'
      endif

      do ix = 2,1000

        ia = nint(a*ix)
      
        if(abs((ia-a*ix)/1).lt.eps) then
!        write(6,'(4x,i6,a,i5,1p2e18.10)')ia,'/',ix,a,ia*1d0/ix
        write(af,'(i6,a,i5,a)') ia,'/',ix
      if(iflatex.eq.1) then
        write(af,'(a,i3,a,i3,a)') ap//'{',abs(ia),'\over',ix,'}'
      endif
        go to 21
        endif

      end do

      do ix = 1,1000

        ip = a/abs(a)
        ia = nint(a**2*ix)
      
        if(abs((ia-a**2*ix)/ia).lt.eps) then
!        write(6,'(a,i6,a,i5,1p2e18.10)')ap//' rt',ia,'/',ix,a,
!     & ip*sqrt(ia*1d0/ix)
        write(af,'(a,i6,a,i5,a)') ap//'Sqrt[',abs(ia),'/',ix,']'
      if(iflatex.eq.1) then
        write(af,'(a,i3,a,i3,a)') ap//'\sqrt{',ia,'\over',ix,'}'
      endif
        go to 21
        endif

      end do

      write(6,*)' not found for',a

21    continue
20    continue
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
