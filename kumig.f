      implicit real*8(a-h,o-z)
      dimension iqni0(9,18)
!      dimension c13(12,12,12,12,0:4,0:1)
!      dimension c14(12,12,12,12,0:4,0:1)
      dimension c13(18,12,12,0:4,0:1)
      dimension c14(18,12,12,0:4,0:1)
      dimension amat13(18,18),amat14(18,18)
      dimension bmatmm13(18,18),bmatmm14(18,18)
      dimension bmatmp13(18,18),bmatmp14(18,18)
! q qbar c cbar
      
      amu = 300
      amc = 300
      amuc = amc/(amu+amc)
      amuu = amu/(amu+amc)
      write(6,*)' mu mc',amu,amc
      write(48,*)' mu mc',amu,amc

      call BCODIN(52)

! create index-list of the 12 states 
      do ij=0,2
      do is=0,1
      do il=0,1
      write(6,*)il,is,ij,iix(1d0*il,1d0*is,1d0*ij)
      end do
      end do
      end do

      call cha(iqni0)

!      call check(iqni0,0,c14,amuu,amuc)
!      write(28,*)'=========='
!      call check(iqni0,1,c13,amuu,amuc)
!      write(6,*)'====='
!      write(28,*)'====='
      do amuu = 0d0,1d0,1d0
      amuc = 1-amuu
      write(28,*)'==========',amuu,amuc
      call checkfx(iqni0,0,c14,amuu,amuc,amat14)
      write(28,*)'=========='
      call cparity(amat14,0,bmatmm14,bmatmp14)

      write(28,*)'==========',amuu,amuc
      call checkfx(iqni0,1,c13,amuu,amuc,amat13)
      write(28,*)'=========='
      call cparity(amat13,1,bmatmm13,bmatmp13)

      call matchk(amat14,18)
      call matchk(amat13,18)

      end do

      write(6,*)'====='



      call calcmi(iqni0,c14,c13,amuu,amuc)

      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine cparity(amat,if13,bmatmm,bmatmp)
      implicit real*8(a-h,j,l,o-z)
      character*100 af
      character*23 a1
      character*23 aa(18,18)
      dimension amat(18,18)
!      dimension iqni0(9,18)
      dimension bmatmm(18,18),bmatmp(18,18)
      dimension iflgmm(10),iflgmp(8),cmm(10,18),cmp(10,18)

      dimension facmu(0:1,0:1),facnu(0:1,0:1)
      character*1 amu(0:1,0:1),anu(0:1,0:1)

      iflatex=1

      iflgmm = (/2,8,10,12,1,7,9,11,14,15/)
      iflgmp = (/3,5,4,6,13,16,17,18/)
      nmxmm = 10
      nmxmp = 8
      nmx = 18
      
      facmu=0
      facnu=0
      amu=' '
      anu=' '

      ip = 1
      if(if13.eq.1) ip = -1
      
      rt2i = 1/sqrt(2d0)
      cmm = 0
      cmm(1, 1) =  rt2i
      cmm(1, 2) =  rt2i*ip
      cmm(2, 3) =  rt2i
      cmm(2, 4) = -rt2i
      cmm(3, 5) =  rt2i
      cmm(3, 6) =  rt2i
      cmm(4, 7) =  rt2i
      cmm(4, 8) =  rt2i*ip
      cmm(5, 9) =  rt2i
      cmm(5,10) = -rt2i*ip
      cmm(6,11) =  rt2i
      cmm(6,12) =  rt2i*ip
      cmm(7,13) =  1
      cmm(8,14) =  rt2i
      cmm(8,15) = -rt2i*ip
      cmm(9,16) =  1
      cmm(10,18)=  1

      cmp = 0
      cmp(1, 1) =  rt2i
      cmp(1, 2) = -rt2i*ip
      cmp(2, 3) =  rt2i
      cmp(2, 4) =  rt2i
      cmp(3, 5) =  rt2i
      cmp(3, 6) = -rt2i
      cmp(4, 7) =  rt2i
      cmp(4, 8) = -rt2i*ip
      cmp(5, 9) =  rt2i
      cmp(5,10) =  rt2i*ip
      cmp(6,11) =  rt2i
      cmp(6,12) = -rt2i*ip
      cmp(7,14) =  rt2i
      cmp(7,15) =  rt2i*ip
      cmp(8,17) =  1


      write(28,*)amat(1,1),amat(18,18)

      bmatmm = 0
      do i2=1,10
      ix = iflgmm(i2)
      
      do i1=1,10
      do i3=1,nmx
      bmatmm(i1,i2) = bmatmm(i1,i2)+ cmm(i1,i3)*amat(i3,ix)
      end do
      end do

      end do

      do 10 iflatex=0,1

      do i1=1,10
      do i2=1,10
        call aform(bmatmm(i1,i2),af,iflatex)

!      write(8,600)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
!     & trim(af),sum,sum**2
!600   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a20,1p2e20.10)
!      write(18,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
!     & '$',trim(af),'$&',sum,sum**2
!
!
!610   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a,a20,a,1pe18.8,e12.4)
      if(iflatex.eq.1) then
        write(a1,'(a,a20,a)')'$',trim(af),'$&'
      else
        write(a1,'(a,a20,a)')'',trim(af),','
      endif
      aa(i1,i2) = a1

      end do
      end do

      if(iflatex.eq.1) then
        do ii=1,10
          write(28,'(18a)')(aa(ii,ici),ici=1,10)
        end do
      else
        do ii=1,10
          write(29,'(18a)')(aa(ii,ici),ici=1,10)
        end do
        write(29,*)' '
      endif
      
10    continue

      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine matchk(amat,nmx)
      implicit real*8(a-h,j,l,o-z)
      character*100 af
      dimension amat(18,18)
      
      do i1=1,nmx
      do i2=1,nmx
      sum = 0
      do i3=1,nmx
      sum = sum + amat(i1,i3)*amat(i2,i3)
      end do
      if(abs(sum).gt.1d-9) write(9,*)i1,i2,sum
      end do
      end do

      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine cha(iqni0)
      implicit real*8(a-h,j,l,o-z)
!      character*100 af
!      character*23 a1
!      character*23 aa(1000,18)
! do again
      dimension iqni(9),iqni0(9,18)
      hlf = 0.5d0

! total j = 1
! iqni = l1 s1 j1 l2 s2 j2 j L J
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

!        ix1 = ixx(iqni)
      do i=1,9
        iqni0(i,ic) = iqni(i)
      end do

!      write(6,'(5(3i5,2x))')iqni,ixx(iqni)

100   continue
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine calcmi(iqni0,c14,c13,amuu,amuc)
      implicit real*8(a-h,j,l,o-z)
      dimension iqnia(9),iqnib(9),iqni0(9,18)
      dimension c13(18,12,12,0:4,0:1)
      dimension c14(18,12,12,0:4,0:1)

      hlf = 0.5d0
      cc12 = -16d0/3
      cc34 = -16d0/3
      cc13 = -16d0/3
      cc24 = -16d0/3
      cc14 = -16d0/3
      cc32 = -16d0/3
      ihit = 0
      write(6,*)'imhere'

! total j = 1
      do 100 ica = 1,18
      do i=1,9
         iqnia(i) = iqni0(i,ica)
      end do
        l12a = iqnia(1)
        s12a = iqnia(2)
        j12a = iqnia(3)
        l34a = iqnia(4)
        s34a = iqnia(5)
        j34a = iqnia(6)
        jia  = iqnia(7)
        lria = iqnia(8)
        jtot = iqnia(9)

      do 100 icb = 1,18
      do i=1,9
         iqnib(i) = iqni0(i,icb)
      end do

      ss12 = 0
      ss34 = 0
      ss14 = 0
      ss23 = 0
      ss13 = 0
      ss24 = 0
      cmi12 = 0
      cmi34 = 0
      cmi14 = 0
      cmi23 = 0
      cmi13 = 0
      cmi24 = 0





        sum12 = 0
        sum13 = 0
        sum14 = 0

        if(abs(jtot-iqnib(9)).gt.1d-5) go to 101
        l12b = iqnib(1)
        s12b = iqnib(2)
        j12b = iqnib(3)
        l34b = iqnib(4)
        s34b = iqnib(5)
        j34b = iqnib(6)
        jib  = iqnib(7)
        lrib = iqnib(8)

! o12+o34
      sum = 0
        do i=1,9
          sum = sum+abs(iqnia(i)-iqnib(i))
        end do
      if(sum.lt.1d-8)then
        sum12 = sum12+1

        ss12 = sigsig(s12a,s12b)
        ss34 = sigsig(s34a,s34b)
        cmi12 = ss12 * cc12
        cmi34 = ss34 * cc34
      endif

! o14+o32
      sum = 0
      sum1 = 0
      do 51 lrf =0d0,1d0,1d0

      do 52 l14 =0d0,1d0,1d0
      do 52 s14=0d0,1d0,1d0
      do 52 j14=0d0,2d0,1d0

      do 52 l32 =0d0,1d0,1d0
      if(abs(lrf+l14+l32-1).lt.1d-3) then
      do 53 s32=0d0,1d0,1d0
      do 53 j32=0d0,2d0,1d0

      ixb1 = iix(l14,s14,j14)
      ixb2 = iix(l32,s32,j32)

      do 53 jf = 0d0,3d0,1d0
!      sum14 =sum14 + c14(ica,ixb1,ixb2,jf,lrf)*c14(icb,ixb1,ixb2,jf,lrf)
      ww = c14(ica,ixb1,ixb2,jf,lrf)*c14(icb,ixb1,ixb2,jf,lrf)
     & *f14(l12a,l34a,lria,l14,l32,lrf,amuu,amuc)
     & *f14(l12b,l34b,lrib,l14,l32,lrf,amuu,amuc)
      sum14 = sum14 + ww
      ss14 = ss14 + sigsig(s14,s14)*ww
      ss23 = ss23 + sigsig(s32,s32)*ww
53    continue
      endif
52    continue
51    continue
! o13+o24
      sum = 0
      sum1 = 0
      do 61 lrf =0d0,1d0,1d0

      do 62 l13 =0d0,1d0,1d0
      
      do 62 s13=0d0,1d0,1d0
      do 62 j13=0d0,2d0,1d0

      do 62 l24 =0d0,1d0,1d0
      if(abs(lrf+l13+l24-1).lt.1d-3) then
      do 63 s24=0d0,1d0,1d0
      do 63 j24=0d0,2d0,1d0

      ixb1 = iix(l13,s13,j13)
      ixb2 = iix(l24,s24,j24)

      do 63 jf = 0d0,3d0,1d0
!      sum13 =sum13 + c13(ica,ixb1,ixb2,jf,lrf)*c13(icb,ixb1,ixb2,jf,lrf)
       ww = c13(ica,ixb1,ixb2,jf,lrf)*c13(icb,ixb1,ixb2,jf,lrf)
     & *f13(l12a,l34a,lria,l13,l24,lrf,amuu,amuc)
     & *f13(l12b,l34b,lrib,l13,l24,lrf,amuu,amuc)
      sum13 = sum13 + ww
      ss13 = ss13 + sigsig(s13,s13)*ww
      ss24 = ss24 + sigsig(s24,s24)*ww
63    continue
      endif
62    continue
61    continue




101   continue
      if(abs(ss12)+abs(ss13)+abs(ss14)+abs(ss23)+abs(ss24)+abs(ss34)
     & .gt.1d-8) then 
!      write(6,*)'ica icb',ica,icb
       write(6,'(a,2i3,1p9e15.6)')'ss12 13 14 23 24 34:',
     & ica,icb,ss12,ss13,ss14,ss23,ss24,ss34
      endif
      if(abs(sum12)+abs(sum13)+abs(sum14).gt.1d-8) 
     & write(6,'(a,2i3,1p9e15.6)')'n12 13 14          :',
     & ica,icb,sum12,sum13,sum14
100   continue
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function sigsig(s1,s2)
      implicit real*8(a-h,o-z)
      sigsig = 0
      if(abs(s1-s2).lt.1d-9) then
       if(abs(s1-1).lt.1d-8) then
       sigsig = 1
       elseif(abs(s1).lt.1d-8) then
       sigsig = -3
       endif
      endif
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function iix(t1,t2,t3)
      implicit real*8(a-h,o-z)
      iix = nint(t1+t2*2+t3*4 + 1)
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function f14(t12,t34,tri,t14,t32,trf,amuu,amuc)
      implicit real*8(a-h,o-z)
      dimension a14(0:1,0:1,0:1,0:1,0:1,0:1)
      a14 = 0
!     a14(l12,l34,lri,l14,l32,lrf)
      a14(1,0,0,1,0,0) = sqrt(amuc/2)
      a14(1,0,0,0,1,0) = sqrt(amuc/2)
      a14(1,0,0,0,0,1) = sqrt(amuu)
      a14(0,1,0,1,0,0) = sqrt(amuu/2)
      a14(0,1,0,0,1,0) = sqrt(amuu/2)
      a14(0,1,0,0,0,1) =-sqrt(amuc)
      a14(0,0,1,1,0,0) = sqrt(1d0/2)
      a14(0,0,1,0,1,0) =-sqrt(1d0/2)
      a14(0,0,1,0,0,1) = 0

      l12 = nint(t12)
      l34 = nint(t34)
      lri = nint(tri)
      l14 = nint(t14)
      l32 = nint(t32)
      lrf = nint(trf)
      
      f14 = a14(l12,l34,lri,l14,l32,lrf)
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      function f13(t12,t34,tri,t13,t24,trf,amuu,amuc)
      implicit real*8(a-h,o-z)
      dimension a13(0:1,0:1,0:1,0:1,0:1,0:1)
      a13 = 0
!     a13(l12,l34,lri,l13,l24,lrf)
      a13(1,0,0,1,0,0) = sqrt(amuc/2)
      a13(1,0,0,0,1,0) =-sqrt(amuc/2)
      a13(1,0,0,0,0,1) = sqrt(amuu)
      a13(0,1,0,1,0,0) =-sqrt(amuu/2)
      a13(0,1,0,0,1,0) = sqrt(amuu/2)
      a13(0,1,0,0,0,1) = sqrt(amuc)
      a13(0,0,1,1,0,0) = sqrt(1d0/2)
      a13(0,0,1,0,1,0) = sqrt(1d0/2)
      a13(0,0,1,0,0,1) = 0

      l12 = nint(t12)
      l34 = nint(t34)
      lri = nint(tri)
      l13 = nint(t13)
      l24 = nint(t24)
      lrf = nint(trf)
      
      f13 = a13(l12,l34,lri,l13,l24,lrf)
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine checkfx(iqni0,if13,c14,amuu,amuc,amat)
      implicit real*8(a-h,j,l,o-z)
      character*100 af
      character*23 a1
      character*23 aa(1000,18)
      character*23 ab(1000,18)
      dimension iqni(9),iqnf(9),iqni0(9,18)
      dimension c14(18,12,12,0:4,0:1),amat(18,18)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      hlf = 0.5d0

      iflatex = 1
      c14 = 0
      amat = 0
!      ihit = 0
! total j = 1
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      do 100 ici = 1,18
      do i=1,9
         iqni(i) = iqni0(i,ici)
      end do

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
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      sum1 = 0

      do 100 icf = 1,18
      do i=1,9
         iqnf(i) = iqni0(i,icf)
      end do

      l14 = iqnf(1)
      s14 = iqnf(2)
      j14 = iqnf(3)
      l32 = iqnf(4)
      s32 = iqnf(5)
      j32 = iqnf(6)
      jf  = iqnf(7)
      lrf = iqnf(8)
!      jtot = iqnf(9)

      write(6,'(3i3,3x,3i3,3x,3i3)')(iqnf(i),i=1,9)
      write(8,'(3i3,3x,3i3,3x,3i3)')(iqnf(i),i=1,9)
      write(18,'(3x,3i2,1x,3i2,1x,3i2)')(iqnf(i),i=1,9)

!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      ixb1 = iix(l14,s14,j14)
      ixb2 = iix(l32,s32,j32)


      sum = 0

      do 20 li=0d0,2d0,1d0
      do 20 lf=0d0,2d0,1d0
      do 20 ll=0d0,2d0,1d0
      
      do 20 si = 0d0,2d0,1d0
      
      if(if13.eq.1) then
        fac1=1
      else
        fac1 = fugo(s34+s32)
      endif
      fac2 = f9j(l12,s12,j12,l34,s34,j34,li,si,ji)
      fac3 = f9j(hlf,hlf,s12,hlf,hlf,s34,s14,s32,si)
      fac4 = f9j(l14,s14,j14,l32,s32,j32,lf,si,jf)
      fac5 = f9j(li,lri,ll,si,0d0,si,ji,lri,jtot)
      fac6 = f9j(lf,lrf,ll,si,0d0,si,jf,lrf,jtot)
      
      sum = sum + fac1*fac2*fac3*fac4*fac5*fac6

20    continue

      sum0 = sum0+sum**2
      ii=icf
      c14(ici,ixb1,ixb2,jf,lrf) = sum

      if(if13.eq.1) then
      sum1 = sum1 + (sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc))**2
      else
      sum1 = sum1 + (sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc))**2
      endif

      if(if13.eq.1) then
        sum = sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc)
      else
        sum = sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc)
      endif
      amat(icf,ici) = sum
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! for tex files
        call aform(sum,af,iflatex)

      write(8,600)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & trim(af),sum,sum**2
600   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a20,1p2e20.10)
      write(18,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & '$',trim(af),'$&',sum,sum**2


610   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a,a20,a,1pe18.8,e12.4)
      write(a1,'(a,a20,a)')'$',trim(af),'$&'
      aa(icf,ici) = a1

!      if(if13.eq.1) then
!       call aform(sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc),af,iflatex)
!        write(38,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
!     & '$',trim(af),'$&',sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc)
!      write(a1,'(a,a20,a)')'$',trim(af),'$&'
!      ab(ii,ici) = a1
!      else
!       call aform(sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc),af,iflatex)
!        write(38,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
!     & '$',trim(af),'$&',sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc)
!      write(a1,'(a,a20,a)')'$',trim(af),'$&'
!      ab(ii,ici) = a1
!      endif


      write(6,*)sum0,sum1
30    continue
      write(6,*)' --- sum1',ici,sum1

!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
100   continue
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -

      do ii=1,18
!      if(ihit(ii).gt.0) then
      write(28,'(18a)')(aa(ii,ici),ici=1,18)
      write(48,'(18a)')(ab(ii,ici),ici=1,18)
!      endif
      end do

!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      end
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine check(iqni0,if13,c14,amuu,amuc)
      implicit real*8(a-h,j,l,o-z)
      character*100 af
      character*23 a1
      character*23 aa(1000,18)
      character*23 ab(1000,18)
! do again
      dimension iqni(9),iqni0(9,18)
      dimension ihit(1000)
!      dimension c13(12,12,12,12,0:4,0:1)
      dimension c14(18,12,12,0:4,0:1)
      hlf = 0.5d0

      c14 = 0
      ihit = 0
! total j = 1
      do 100 ic1 = 1,2
      
      do 100 ic = 1,18

      do i=1,9
         iqni(i) = iqni0(i,ic)
      end do

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

      sum1 = 0
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

      ixb1 = iix(l14,s14,j14)
      ixb2 = iix(l32,s32,j32)

      do 10 jf = 0d0,3d0,1d0
      sum = 0

      do 20 li=0d0,2d0,1d0
      do 20 lf=0d0,2d0,1d0
      do 20 ll=0d0,2d0,1d0
      
      do 20 si = 0d0,2d0,1d0
      
      if(if13.eq.1) then
        fac1=1
      else
        fac1 = fugo(s34+s32)
      endif
      fac2 = f9j(l12,s12,j12,l34,s34,j34,li,si,ji)
      fac3 = f9j(hlf,hlf,s12,hlf,hlf,s34,s14,s32,si)
      fac4 = f9j(l14,s14,j14,l32,s32,j32,lf,si,jf)
      fac5 = f9j(li,lri,ll,si,0d0,si,ji,lri,jtot)
      fac6 = f9j(lf,lrf,ll,si,0d0,si,jf,lrf,jtot)
      
      sum = sum + fac1*fac2*fac3*fac4*fac5*fac6

20    continue
      sum0 = sum0+sum**2
      ii=ii+1
      c14(ic,ixb1,ixb2,jf,lrf) = sum
      if(if13.eq.1) then
      sum1 = sum1 + (sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc))**2
      else
      sum1 = sum1 + (sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc))**2
      endif

      if(ic1.eq.1 .and. abs(sum).gt.1d-8) then
      ihit(ii) = ihit(ii)+1
      end if

      if(ic1.eq.2 .and. ihit(ii).gt.0) then
        call aform(sum,af,iflatex)

!      write(6,600)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
!     & trim(af),sum,sum**2
      write(8,600)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & trim(af),sum,sum**2
600   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a20,1p2e20.10)
      write(18,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & '$',trim(af),'$&',sum,sum**2


610   format(3f5.1,3x,3f5.1,3x,3f5.1,2x,a,a20,a,1pe18.8,e12.4)
      write(a1,'(a,a20,a)')'$',trim(af),'$&'
      aa(ii,ic) = a1

      if(if13.eq.1) then
       call aform(sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc),af,iflatex)
        write(38,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & '$',trim(af),'$&',sum*f13(l12,l34,lri,l14,l32,lrf,amuu,amuc)
      write(a1,'(a,a20,a)')'$',trim(af),'$&'
      ab(ii,ic) = a1
      else
       call aform(sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc),af,iflatex)
        write(38,610)l14,s14,j14,l32,s32,j32,jf,lrf,jtot,
     & '$',trim(af),'$&',sum*f14(l12,l34,lri,l14,l32,lrf,amuu,amuc)
      write(a1,'(a,a20,a)')'$',trim(af),'$&'
      ab(ii,ic) = a1
      endif



      endif
      
10    continue
      end if
15    continue

      write(6,*)sum0,sum1
30    continue
      write(6,*)' --- sum1',ic,sum1

100   continue

      do ii=1,1000
      if(ihit(ii).gt.0) then
      write(28,'(18a)')(aa(ii,ic),ic=1,18)
      write(48,'(18a)')(ab(ii,ic),ic=1,18)
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
      subroutine aform(a,af,iflatex)
      implicit real*8(a-h,o-z)
      character*1 ap
      character*100 af,ad
      eps = 1d-8
      jmx = 1024

!      iflatex = 1
      
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
