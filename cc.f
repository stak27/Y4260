c==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      module module_phys_params
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      implicit none
      integer,parameter :: wp = selected_real_kind(p=12)
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! pi
      real(wp),parameter::pi=4*atan(1d0)
      real(wp),parameter::pi4=4*pi,pi2cube=(2*pi)**3,rtpi=sqrt(pi)
! i
      complex(wp),parameter::zi=(0d0,1d0)
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! matrices of su(2),su(3)
      complex(wp),dimension(2,2,0:3)::ztau,ztaub,zsig
      complex(wp),dimension(3,3,0:8)::zlam,zlamb
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! hbarc
!      real(wp),parameter::hbarc=197.327053d0
      real(wp),parameter::hbarc   = 197.3269788d0
      real(wp),parameter::alphaem = 1/137.035999139d0
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! constituent masses of the quarks in MeV
!      real(wp),parameter,dimension(0:6)::
!     & qmass0=(/0d0,313d0,313d0,500d0,1500d0,4200d0,170000d0/)
! in fm-1
!      real(wp),dimension(0:6)::qmass=qmass0/hbarc
! name of the quarks 
      character(len=1),parameter,dimension(0:6)::
     & aqname=(/'q','u','d','s','c','b','t'/)
! name of the antiquarks
      character(len=4),dimension(0:6)::aqbarname=aqname//'bar'
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! date
      character(len=8)::adate
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      contains
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine prepare_phys_params(io,iwrite)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      integer :: io,iwrite
      integer :: i
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      call date_and_time(adate)
      write(io,'(a)')
     & ' -- calculate on yyyymmdd: '//adate//' -- by Sachiko'
      write(io,*)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      write(io,'(a)')' -- fundamental parameters'
      write(io,'(8x,a,f20.12)') 'pi     =',pi
!      write(io,'(8x,a,f20.12)') '4pi    =',pi4
!      write(io,'(8x,a,f20.12)') '(2pi)^3=',pi2cube
!      write(io,'(8x,a,f20.12)') 'rtpi   =',rtpi
      write(io,'(8x,a,f20.12)') 'hbarc  =',hbarc
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
!      if(iwrite.ge.1) then
!        write(io,'(a)')' -- quark names and masses'
!        do i=0,6
!          write(io,'(8x,a,2f15.5)') 
!     &    aqname(i)//', '//aqbarname(i),qmass0(i),qmass(i)
!        end do
!      write(io,*)
!      end if
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      call make_lamsig(io,iwrite)
      write(io,'(8x,a)') 'lambda and sigma matrices are calculated.'
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      end subroutine
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      subroutine make_lamsig(io,iwrite)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      integer :: io,iwrite
      integer :: i,j,k
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
! Gell-Mann matrices lambda
! zlam(n1,n2,m)=   <n1| lambda_m|n2> for m=1-8
!              =  identity matrix for m=0
!
! zlamb lambda for antiparticle
! zlamb = -^t zlam
!
! tau taub sig, isospin spin matrices
!
!(n1,n2=1,2,3 or r g b)
      zlam = 0
      zlam(1,1,0) = 1
      zlam(2,2,0) = 1
      zlam(3,3,0) = 1
      zlam(1,2,1) = 1
      zlam(2,1,1) = 1
      zlam(1,2,2) =-zi
      zlam(2,1,2) = zi
      zlam(1,1,3) = 1
      zlam(2,2,3) =-1
      zlam(1,3,4) = 1
      zlam(3,1,4) = 1
      zlam(1,3,5) =-zi
      zlam(3,1,5) = zi
      zlam(2,3,6) = 1
      zlam(3,2,6) = 1
      zlam(2,3,7) =-zi
      zlam(3,2,7) = zi
      zlam(1,1,8) = 1/sqrt(3d0)
      zlam(2,2,8) = 1/sqrt(3d0)
      zlam(3,3,8) =-2/sqrt(3d0)

! for antiquarks
      do 10 j=1,3
      do 10 k=1,3
      zlamb(j,k,0) =  zlam(j,k,0)
      do 10 i=1,8
      zlamb(j,k,i) = -zlam(k,j,i)
10    continue

!Pauli matrices
      do 20 i=0,3
      do 20 j=1,2
      do 20 k=1,2
      ztau(j,k,i) = zlam(j,k,i)
      ztaub(j,k,i) = zlamb(j,k,i)
20    continue

! spin matrices
      zsig = ztau

! Pauli matrices for ubar and dbar isospin
!      ztaub = 0
!      do 21 i=0,3
!      do 21 j=1,2
!      do 21 k=1,2
!      do 21 j1=1,2
!      do 21 k1=1,2
!      ztaub(j,k,i) = ztaub(j,k,i)
!     & + ztau(j,j1,2)*ztau(j1,k1,i)*ztau(k1,k,2)
!21    continue
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      if(iwrite.ge.3)then
      write(io,*)'-- Pauli spin op: tau, taub'
      do k=1,3
      do i=1,2
      write(io,'(2(2(2f5.1,x),a))')
     & (ztau(i,j,k),j=1,2),',',(ztaub(i,j,k),j=1,2)
      end do
      write(io,*)''
      end do

      write(io,'(a,25x,a)')'-- Gell-Mann color op for q,','qbar'
      do k=1,8
      do i=1,3
      write(io,'(2(3(2f8.3,x),a))') 
     & (zlam(i,j,k),j=1,3),',',(zlamb(i,j,k),j=1,3)
      end do
      write(io,*)''
      end do
      endif
c---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      end subroutine
c==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      end module
!==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== =
      use module_phys_params
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      dimension zcc(3,3,3,3)
      dimension z121(3,3,3,3),z128(3,3,3,3)
      dimension z141(3,3,3,3),z148(3,3,3,3)
      dimension z133(3,3,3,3),z136(3,3,3,3)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      call prepare_phys_params(6,3)

      write(6,*)' color for q(1)qbar(2)q(3)qbar(4) systems'

      write(6,*)' 12 34 is color singlet'
      zcc = 0
      do i1=1,3
      do i2=1,3
      zcc(i1,i1,i2,i2) = 1d0/3
      end do
      end do
      call alamlam(zcc,zcc)
      z121 = zcc

      write(6,*)' 12 34 is color octet'

      zcc = 0
      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
      do ia = 1,8
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i2,ia)*zlam(i3,i4,ia)/sqrt(32d0)
      end do
      end do
      end do
      end do
      end do

      call alamlam(zcc,zcc)
      z128 = zcc

      write(6,*)' 14 32 is color singlet'

      zcc = 0
      do i1=1,3
      do i2=1,3
      zcc(i1,i2,i2,i1) = 1d0/3
      end do
      end do

      call alamlam(zcc,zcc)
      z141 = zcc


      write(6,*)' 14 32 is color octet'

      zcc = 0
      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
      do ia = 1,8
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i4,ia)*zlam(i3,i2,ia)/sqrt(32d0)
      end do
      end do
      end do
      end do
      end do

      call alamlam(zcc,zcc)
      z148 = zcc

      write(6,*)' 13 24 is color 3bar'

      zcc = 0
      f = 1/sqrt(12d0)
      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
      ia = 2
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      ia = 5
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      ia = 7
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      end do
      end do
      end do
      end do

      call alamlam(zcc,zcc)
      z133 = zcc

      write(6,*)' 13 24 is color 6'

      zcc = 0
      f = -1/sqrt(24d0)
!      f = -1/sqrt(16d0)
      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
      ia = 1
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      ia = 3
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      ia = 4
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      ia = 6
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      ia = 8
      zcc(i1,i2,i3,i4) = zcc(i1,i2,i3,i4) 
     & +zlam(i1,i3,ia)*zlamb(i2,i4,ia)*f
      end do
      end do
      end do
      end do
      do i1=1,3
      do i2=1,3
      zcc(i1,i2,i1,i2) = zcc(i1,i2,i1,i2) -f*2/3
      end do
      end do
      call alamlam(zcc,zcc)
      z136 = zcc

      zcc = sqrt(2/3d0)*z121+sqrt(1/3d0)*z128
      write(6,*)z136(1,1,1,1),z136(1,1,2,2),z136(1,3,1,3)
      write(6,*)zcc(1,1,1,1),zcc(1,1,2,2),zcc(1,3,1,3)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      write(6,*)' 121-'
      zcc = z121
      call anorm(zcc,z121,zsum1)
      write(6,'(a,9f15.8)')' 121',zsum1,zsum1**2

      call anorm(zcc,z128,zsum1)
      write(6,'(a,9f15.8)')' 128',zsum1,zsum1**2

      call anorm(zcc,z141,zsum1)
      write(6,'(a,9f15.8)')' 141',zsum1,zsum1**2

      call anorm(zcc,z148,zsum1)
      write(6,'(a,9f15.8)')' 148',zsum1,zsum1**2

      call anorm(zcc,z133,zsum1)
      write(6,'(a,9f15.8)')' 133',zsum1,zsum1**2

      call anorm(zcc,z136,zsum1)
      write(6,'(a,9f15.8)')' 136',zsum1,zsum1**2

      write(6,*)' 128-'
      zcc = z128
      call anorm(zcc,z121,zsum1)
      write(6,'(a,9f15.8)')' 121',zsum1,zsum1**2

      call anorm(zcc,z128,zsum1)
      write(6,'(a,9f15.8)')' 128',zsum1,zsum1**2

      call anorm(zcc,z141,zsum1)
      write(6,'(a,9f15.8)')' 141',zsum1,zsum1**2

      call anorm(zcc,z148,zsum1)
      write(6,'(a,9f15.8)')' 148',zsum1,zsum1**2

      call anorm(zcc,z133,zsum1)
      write(6,'(a,9f15.8)')' 133',zsum1,zsum1**2

      call anorm(zcc,z136,zsum1)
      write(6,'(a,9f15.8)')' 136',zsum1,zsum1**2

      write(6,*)' 141-'
      zcc = z141
      call anorm(zcc,z121,zsum1)
      write(6,'(a,9f15.8)')' 121',zsum1,zsum1**2

      call anorm(zcc,z128,zsum1)
      write(6,'(a,9f15.8)')' 128',zsum1,zsum1**2

      call anorm(zcc,z141,zsum1)
      write(6,'(a,9f15.8)')' 141',zsum1,zsum1**2

      call anorm(zcc,z148,zsum1)
      write(6,'(a,9f15.8)')' 148',zsum1,zsum1**2

      call anorm(zcc,z133,zsum1)
      write(6,'(a,9f15.8)')' 133',zsum1,zsum1**2

      call anorm(zcc,z136,zsum1)
      write(6,'(a,9f15.8)')' 136',zsum1,zsum1**2

      write(6,*)' 148-'
      zcc = z148
      call anorm(zcc,z121,zsum1)
      write(6,'(a,9f15.8)')' 121',zsum1,zsum1**2

      call anorm(zcc,z128,zsum1)
      write(6,'(a,9f15.8)')' 128',zsum1,zsum1**2

      call anorm(zcc,z141,zsum1)
      write(6,'(a,9f15.8)')' 141',zsum1,zsum1**2

      call anorm(zcc,z148,zsum1)
      write(6,'(a,9f15.8)')' 148',zsum1,zsum1**2

      call anorm(zcc,z133,zsum1)
      write(6,'(a,9f15.8)')' 133',zsum1,zsum1**2

      call anorm(zcc,z136,zsum1)
      write(6,'(a,9f15.8)')' 136',zsum1,zsum1**2

      write(6,*)' 133-'
      zcc = z133
      call anorm(zcc,z121,zsum1)
      write(6,'(a,9f15.8)')' 121',zsum1,zsum1**2

      call anorm(zcc,z128,zsum1)
      write(6,'(a,9f15.8)')' 128',zsum1,zsum1**2

      call anorm(zcc,z141,zsum1)
      write(6,'(a,9f15.8)')' 141',zsum1,zsum1**2

      call anorm(zcc,z148,zsum1)
      write(6,'(a,9f15.8)')' 148',zsum1,zsum1**2

      call anorm(zcc,z133,zsum1)
      write(6,'(a,9f15.8)')' 133',zsum1,zsum1**2

      call anorm(zcc,z136,zsum1)
      write(6,'(a,9f15.8)')' 136',zsum1,zsum1**2

      write(6,*)' 136-'
      zcc = z136
      call anorm(zcc,z121,zsum1)
      write(6,'(a,9f15.8)')' 121',zsum1,zsum1**2

      call anorm(zcc,z128,zsum1)
      write(6,'(a,9f15.8)')' 128',zsum1,zsum1**2

      call anorm(zcc,z141,zsum1)
      write(6,'(a,9f15.8)')' 141',zsum1,zsum1**2

      call anorm(zcc,z148,zsum1)
      write(6,'(a,9f15.8)')' 148',zsum1,zsum1**2

      call anorm(zcc,z133,zsum1)
      write(6,'(a,9f15.8)')' 133',zsum1,zsum1**2

      call anorm(zcc,z136,zsum1)
      write(6,'(a,9f15.8)')' 136',zsum1,zsum1**2

!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
      write(6,*)' lamlam z121'
      zcc=z121
      call alamlam(zcc,z121)
      call alamlam(zcc,z128)
      call alamlam(zcc,z141)
      call alamlam(zcc,z148)
      call alamlam(zcc,z133)
      call alamlam(zcc,z136)
      write(6,*)' lamlam z128'
      zcc=z128
      call alamlam(zcc,z121)
      call alamlam(zcc,z128)
      call alamlam(zcc,z141)
      call alamlam(zcc,z148)
      call alamlam(zcc,z133)
      call alamlam(zcc,z136)
      write(6,*)' lamlam z141'
      zcc=z141
      call alamlam(zcc,z121)
      call alamlam(zcc,z128)
      call alamlam(zcc,z141)
      call alamlam(zcc,z148)
      call alamlam(zcc,z133)
      call alamlam(zcc,z136)
      write(6,*)' lamlam z148'
      zcc=z148
      call alamlam(zcc,z121)
      call alamlam(zcc,z128)
      call alamlam(zcc,z141)
      call alamlam(zcc,z148)
      call alamlam(zcc,z133)
      call alamlam(zcc,z136)
      write(6,*)' lamlam z133'
      zcc=z133
      call alamlam(zcc,z121)
      call alamlam(zcc,z128)
      call alamlam(zcc,z141)
      call alamlam(zcc,z148)
      call alamlam(zcc,z133)
      call alamlam(zcc,z136)
      write(6,*)' lamlam z136'
      zcc=z136
      call alamlam(zcc,z121)
      call alamlam(zcc,z128)
      call alamlam(zcc,z141)
      call alamlam(zcc,z148)
      call alamlam(zcc,z133)
      call alamlam(zcc,z136)
!---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -

      end

      subroutine anorm(zcc1,zcc2,zsum1)
      use module_phys_params
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      dimension zcc1(3,3,3,3),zcc2(3,3,3,3)
!1
      zsum1 = 0
      do 10 i1 = 1,3
      do 10 i2 = 1,3
      do 10 i3 = 1,3
      do 10 i4 = 1,3
      zsum1 = zsum1 + conjg(zcc1(i1,i2,i3,i4))*zcc2(i1,i2,i3,i4)
10    continue

      zc = 0
      ia = 1
      do 20 i1 = 1,3
      do 20 i2 = 1,3
      do 20 i3 = 1,3
      do 20 i4 = 1,3
      do 20 j1 = 1,3
      do 20 j2 = 1,3
      do 20 j3 = 1,3
      do 20 j4 = 1,3
      zc = zc + conjg(zcc1(i1,i2,i3,i4))*zcc1(j1,j2,j3,j4)*(
     & zlam(i1,j1,ia)*zlamb(i2,j2, 0)*zlam(i3,j3, 0)*zlamb(i4,j4, 0)
     &+zlam(i1,j1, 0)*zlamb(i2,j2,ia)*zlam(i3,j3, 0)*zlamb(i4,j4, 0)
     &+zlam(i1,j1, 0)*zlamb(i2,j2, 0)*zlam(i3,j3,ia)*zlamb(i4,j4, 0)
     &+zlam(i1,j1, 0)*zlamb(i2,j2, 0)*zlam(i3,j3, 0)*zlamb(i4,j4,ia))
20    continue
!      write(6,*)' ia=3',zc


      end



      subroutine alamlam(zcc1,zcc2)
      use module_phys_params
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      dimension zcc1(3,3,3,3),zcc2(3,3,3,3)
      dimension zsum(4,4)
!1
      zsum1 = 0
      do 10 i1 = 1,3
      do 10 i2 = 1,3
      do 10 i3 = 1,3
      do 10 i4 = 1,3
      zsum1 = zsum1 + conjg(zcc1(i1,i2,i3,i4))*zcc2(i1,i2,i3,i4)
10    continue
      write(6,*)' 1',zsum1

!1
      zsum = 0
      do 20 i1 = 1,3
      do 20 i2 = 1,3
      do 20 i3 = 1,3
      do 20 i4 = 1,3
      do 20 j1 = 1,3
      do 20 j2 = 1,3
      do 20 j3 = 1,3
      do 20 j4 = 1,3
      ww = conjg(zcc1(i1,i2,i3,i4))*zcc2(j1,j2,j3,j4)
      zsum(1,1) = zsum(1,1) + ww*
     & zlam(i1,j1, 0)*zlamb(i2,j2, 0)*zlam(i3,j3, 0)*zlamb(i4,j4, 0)
      do ia = 1,8
      zsum(1,2) = zsum(1,2) + ww*
     & zlam(i1,j1,ia)*zlamb(i2,j2,ia)*zlam(i3,j3, 0)*zlamb(i4,j4, 0)
      zsum(1,3) = zsum(1,3) + ww*
     & zlam(i1,j1,ia)*zlamb(i2,j2, 0)*zlam(i3,j3,ia)*zlamb(i4,j4, 0)
      zsum(1,4) = zsum(1,4) + ww*
     & zlam(i1,j1,ia)*zlamb(i2,j2, 0)*zlam(i3,j3, 0)*zlamb(i4,j4,ia)
      zsum(2,3) = zsum(2,3) + ww*
     & zlam(i1,j1, 0)*zlamb(i2,j2,ia)*zlam(i3,j3,ia)*zlamb(i4,j4, 0)
      zsum(2,4) = zsum(2,4) + ww*
     & zlam(i1,j1, 0)*zlamb(i2,j2,ia)*zlam(i3,j3, 0)*zlamb(i4,j4,ia)
      zsum(3,4) = zsum(3,4) + ww*
     & zlam(i1,j1, 0)*zlamb(i2,j2, 0)*zlam(i3,j3,ia)*zlamb(i4,j4,ia)
      end do
20    continue
      do i=1,4
      write(6,'(a,1p9e15.7)')' lamlam',(zsum(i,j),j=1,4)
      end do

      end
