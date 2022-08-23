      module non_advection_m
      use common_hh
      use common_geom
      use common_hyd
      use ss_nu_m
      implicit none
      contains
!*******************************
      subroutine sor(sorerr,lsor,soralpha,lmax)
!*******************************
      real(8)::sorerr,soralpha
      integer::lsor,lmax
      real(8)::div,err,ww,errx,dxh,dyh
      real(8),dimension(0:im,0:jm,0:km)::work
      integer::l

      div=0.
      work=0.;a_n=0.;a_s=0.;a_w=0.;a_e=0.;a_u=0.;a_d=0.;a_p=0.;a_f=0.
      dxh=0.;dyh=0.

      do i=1,nx
       do j=1,ny
        dxh=2.*(dx(i,j)+dx(i,j-1))
        dyh=2.*(dy(i,j)+dy(i-1,j))
        do k=1,nz
         if(obst3d(i,j,k)==0) then

          work(i,j,k)=ypn(i,j,k)

          if(i==nx.and.j_east==3) then !East Periodic
           a_n(i,j,k)=a_n0(i,j)*(hs(i,j)+hs(1,j)) 
          else if(((i==nx).and.(j_east<=2)).or.obst3d(i+1,j,k)==1) then
           a_n(i,j,k)=0.
          else
           a_n(i,j,k)=a_n0(i,j)*(hs(i,j)+hs(i+1,j))
          end if

          if(i==1.and.j_west==3) then !West Periodic
           a_s(i,j,k)=a_s0(i,j)*(hs(i,j)+hs(nx,j)) 
          else if(((i==1).and.(j_west<=2)).or.obst3d(i-1,j,k)==1) then
           a_s(i,j,k)=0.
          else
           a_s(i,j,k)=a_s0(i,j)*(hs(i,j)+hs(i-1,j))
          end if

          if(j==ny.and.j_north==3) then !North Periodic
           a_w(i,j,k)=a_w0(i,j)*(hs(i,j)+hs(i,1))
          else if(((j==ny).and.(j_north<=2)).or.obst3d(i,j+1,k)==1) then
           a_w(i,j,k)=0.
          else
           a_w(i,j,k)=a_w0(i,j)*(hs(i,j)+hs(i,j+1))
          end if

          if(j==1.and.j_south==3) then !South Periodic
           a_e(i,j,k)=a_e0(i,j)*(hs(i,j)+hs(i,ny))
          else if(((j==1).and.(j_south<=2)).or.obst3d(i,j-1,k)==1) then
           a_e(i,j,k)=0.
          else
           a_e(i,j,k)=a_e0(i,j)*(hs(i,j)+hs(i,j-1))
          end if
!
          if(k==nz .or. obst3d(i,j,k+1)==1) then
           a_u(i,j,k)=0.
          else
           a_u(i,j,k)=a_u0(k)/hs(i,j)
          end if
          if(k==1 .or.obst3d(i,j,k-1)==1) then
           a_d(i,j,k)=0.
          else
           a_d(i,j,k)=a_d0(k)/hs(i,j)
          end if

          a_p(i,j,k)=a_n(i,j,k)+a_s(i,j,k)+a_w(i,j,k)+a_e(i,j,k) &
                     +a_u(i,j,k)+a_d(i,j,k)
          if(k==nz) a_p(i,j,nz)=a_p(i,j,nz)+2./(hs(i,j)*dz(nz)**2)
          div=((hs(i+1,j)+hs(i,j))*yu(i,j,k)*dn(i,j)  &
              -(hs(i,j)+hs(i-1,j))*yu(i-1,j,k)*dn(i-1,j))*2./ &
              ((dxi(i,j)+dxi(i-1,j))*(dn(i,j)+dn(i-1,j)))  &
             +((hs(i,j+1)+hs(i,j))*yv(i,j,k)- &
              (hs(i,j)+hs(i,j-1))*yv(i,j-1,k))/(dyj(i,j)+dyj(i,j-1)) &
             +(w1(i,j,k)-w1(i,j,k-1))/dz(k)

!         a_f(i,j,k)=-div/dt*rho
          a_f(i,j,k)=-div/dt*rho &
           -(omega(i,j,k)*(ypn(i+1,j,k+1)+ypn(i,j,k+1)-ypn(i+1,j,k-1)-ypn(i,j,k-1)) &
           +omega(i-1,j,k)*(ypn(i,j,k+1)+ypn(i-1,j,k+1)-ypn(i,j,k-1)-ypn(i-1,j,k-1))) &
             /(dxh*dz(k)) &
           -(psi(i,j,k)*(ypn(i,j+1,k+1)+ypn(i,j,k+1)-ypn(i,j+1,k-1)-ypn(i,j,k-1)) &
           +psi(i,j-1,k)*(ypn(i,j,k+1)+ypn(i,j-1,k+1)-ypn(i,j,k-1)-ypn(i,j-1,k-1))) &
             /(dyh*dz(k)) 
          if(k==1) a_f(i,j,k)=a_f(i,j,k)-rho*g*(-(yc(i,j,1)-c0))/dz(1)
         end if
         if(j==nym) then
!         write(44,'(2i3,7e12.3)')&
!            i,k,a_n(i,j,k),a_s(i,j,k),a_w(i,j,k),a_e(i,j,k), &
!                 a_u(i,j,k),a_d(i,j,k),a_f(i,j,k)
         end if
        end do
       end do
      end do
!
      do l=1,lsor
       err=0.0
       do i=1,nx
        do j=1,ny
         do k=1,nz
          if(obst3d(i,j,k)==1) then
           work(i,j,k)=0.
          else
           ww=(a_n(i,j,k)*ypn(i+1,j,k)+a_s(i,j,k)*ypn(i-1,j,k) &
            +a_w(i,j,k)*ypn(i,j+1,k)+a_e(i,j,k)*ypn(i,j-1,k) &
           +a_u(i,j,k)*ypn(i,j,k+1)+a_d(i,j,k)*ypn(i,j,k-1)+a_f(i,j,k)) &
           /a_p(i,j,k)
           work(i,j,k)=(1.0-soralpha)*ypn(i,j,k)+soralpha*ww
           errx=abs(work(i,j,k)-ypn(i,j,k))
           err=err+errx
           ypn(i,j,k)=work(i,j,k)
!          write(44,'(3i4,8e12.3)') &
!            i,j,k,a_n(i,j,k),a_s(i,j,k),a_w(i,j,k), &
!            a_e(i,j,k),a_u(i,j,k),a_d(i,j,k),a_f(i,j,k),a_p(i,j,k)
!          write(44,'(8e12.3)') ypn(i+1,j,k),ypn(i-1,j,k), &
!                  ypn(i,j+1,k),ypn(i,j-1,k), &
!                  ypn(i,j,k+1),ypn(i,j,k-1),ww,work(i,j,k)
          end if
         end do
        end do
       end do
       if (err.lt.sorerr) goto 40
!底面と水面
       do i=1,nx
        do j=1,ny
         ypn(i,j,0)=ypn(i,j,1)-rho*g*hs(i,j)*(-(yc(i,j,1)-c0))*dz(1)
         if(j_surf==0 .or. (j_surf==1 .and. time<stime_surf)) &
          ypn(i,j,nz+1)=-ypn(i,j,nz)
        end do
       end do
! Sounth and North (i=0,nx+1),(j=0,ny+1),(k=0,nz+1)
       do i=1,nx
        do k=0,nz+1
         if(j_south==1) then !south closed
          ypn(i,0,k)=ypn(i,1,k)
         else if(j_south==2) then !south open
          ypn(i,0,k)=-ypn(i,1,k)
         else if(j_south==3) then
          ypn(i,0,k)=ypn(i,ny,k)
          ypn(i,ny+1,k)=ypn(i,1,k)
         end if
         if(j_north==1) then !north closed
          ypn(i,ny+1,k)=ypn(i,ny,k)
         else if(j_north==2) then !north open 
          ypn(i,ny+1,k)=-ypn(i,ny,k)
         end if
        end do
       end do
!West and East  (i=0 and nx+1),(j=1,ny),(k=0,nz+1)
       do j=0,ny+1
        do k=0,nz+1
         if(j_west==1) then ! west closed
          ypn(0,j,k)=ypn(1,j,k)
         else if(j_west==2) then !west open 
          ypn(0,j,k)=-ypn(1,j,k) !0圧力条件 
         else if(j_west==3) then
          ypn(0,j,k)=ypn(nx,j,k)
          ypn(nx+1,j,k)=ypn(1,j,k)
         end if
         if(j_east==1) then !east closed
          ypn(nx+1,j,k)=ypn(nx,j,k)
         else if(j_east==2) then !east open
          ypn(nx+1,j,k)=-ypn(nx,j,k) !0圧力条件
         end if
!        write(44,'(a2,2i4,2e12.4)') 'EW',j,k,ypn(0,j,k),ypn(nx+1,j,k)
        end do
       end do
! Four corners
!
!      do i=0,nx+1,nx+1
!       j=ny
!       do k=0,nz+1
!        ypn(i,j+1,k)=ypn(i,j,k)
!        write(44,'(3i4,f12.4)') i,j,k,ypn(i,j,k)
!       end do
!       j=1
!       do k=0,nz+1
!        ypn(i,j-1,k)=ypn(i,j,k)
!        write(44,'(3i4,f12.4)') i,j,k,ypn(i,j,k)
!       end do
!      end do

! Obstacles
       do i=1,nx
        do j=1,ny
         do k=1,nz
          if(obst3d(i,j,k)==0) then
           if(i<nx.and.obst3d(i+1,j,k)==1) ypn(i+1,j,k)=ypn(i,j,k)
           if(i>1 .and.obst3d(i-1,j,k)==1) ypn(i-1,j,k)=ypn(i,j,k)
           if(j<ny.and.obst3d(i,j+1,k)==1) ypn(i,j+1,k)=ypn(i,j,k)
           if(j>1 .and.obst3d(i,j-1,k)==1) ypn(i,j-1,k)=ypn(i,j,k)
           if(k>1.and.obst3d(i,j,k-1)==1) &
            ypn(i,j,k-1)=ypn(i,j,k)-rho*g*hs(i,j)*(-(yc(i,j,k)-c0))*dz(k)
           if(k<nz.and.obst3d(i,j,k+1)==1) &
            ypn(i,j,k+1)=ypn(i,j,k)+rho*g*hs(i,j)*(-(yc(i,j,k)-c0))*dz(k)
          end if
         end do
        end do
       end do
      end do

 40   continue
      lmax=l
!
!     write(44,*) 'sor',time
!     j=nym
!     do i=1,nx
!      write(44,'(i3,20e12.3)') i,(ypn(i,j,k),k=1,nz)
!     end do

      end subroutine sor

!*******************************
      subroutine rhs
!*******************************
!
!     real(8)::deds,dhsds
      real(8)::dhds,hs_upp,dpds,dpdg,dd,cwp,ff
!     real(8)::dedn,dhsdn
      real(8)::dhdn,hs_vp,dpdn,crf
      real(8)::yvup,yuvp,crv_up,crv_vp

!     omega=0.; psi=0.
      dpdg=0.

      do i=0,nx
       do j=1,ny
!       deds=(eta(i+1,j)-eta(i,j))/dxi(i,j)
!       dhsds=(hs(i+1,j)-hs(i,j))/dxi(i,j)
        dhds=(h(i+1,j)-h(i,j))/dxi(i,j)
        hs_upp=(hs(i+1,j)+hs(i,j))*.5
        crv_up=(curv(i,j)+curv(i,j-1))*.5
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i+1,j,k)==1) then
          yun(i,j,k)=0.
         else
          yvup=(yv(i,j,k)+yv(i,j-1,k)+yv(i+1,j,k)+yv(i+1,j-1,k))*.25
!         omega(i,j,k)=deds+xxi(k)*dhsds
          dpds=(ypn(i+1,j,k)-ypn(i,j,k))/dxi(i,j)
          dpdg=(ypn(i+1,j,k+1)+ypn(i,j,k+1)  &
           -ypn(i+1,j,k-1)-ypn(i,j,k-1))/(4.*dz(k))
          dd=dpds-dpdg*omega(i,j,k)/hs_upp
          crf=-yu(i,j,k)*yvup*crv_up
          yun(i,j,k)=yu(i,j,k)-(crf+dd/rho+g*dhds)*dt
         end if
        end do
       end do
      end do
!
      do i=1,nx
       do j=0,ny
!       dedn=(eta(i,j+1)-eta(i,j))/dyj(i,j)
!       dhsdn=(hs(i,j+1)-hs(i,j))/dyj(i,j)
        dhdn=(h(i,j+1)-h(i,j))/dyj(i,j)
        hs_vp=(hs(i,j+1)+hs(i,j))*.5
        crv_vp=(curv(i,j)+curv(i-1,j))*.5
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i,j+1,k)==1) then
          yvn(i,j,k)=0.
         else
          yuvp=(yu(i,j,k)+yu(i-1,j,k)+yu(i,j+1,k)+yu(i-1,j+1,k))*.25
!         psi(i,j,k)=dedn+xxi(k)*dhsdn
          dpdn=(ypn(i,j+1,k)-ypn(i,j,k))/dyj(i,j)
          dpdg=(ypn(i,j+1,k+1)+ypn(i,j,k+1) &
           -ypn(i,j+1,k-1)-ypn(i,j,k-1))/(4.*dz(k))
          dd=dpdn-dpdg*psi(i,j,k)/hs_vp
          crf=yuvp**2*crv_vp
          yvn(i,j,k)=yv(i,j,k)-(crf+dd/rho+g*dhdn)*dt
         end if
        end do
       end do
      end do
!
      do i=1,nx
       do j=1,ny
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i,j,k+1)==1) then
          ywn(i,j,k)=0.
         else
          cwp=(yc(i,j,k)+yc(i,j,k+1))*.5
          ff=cwp-c0
          ywn(i,j,k)=yw(i,j,k) &
            -(ypn(i,j,k+1)-ypn(i,j,k))*dt/(rho*dzk(k)*hs(i,j)) &
            -g*ff*dt
         end if
        end do
       end do
      end do
!
      end subroutine rhs


!*******************************
      subroutine omgpsi
!*******************************
!
      real(8)::deds,dhsds
      real(8)::dedn,dhsdn

      omega=0.; psi=0.
      deds=0.; dhsds=0.; dedn=0.; dhsdn=0.

      do i=0,nx
       do j=1,ny
        deds=(eta(i+1,j)-eta(i,j))/dxi(i,j)
        dhsds=(hs(i+1,j)-hs(i,j))/dxi(i,j)
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i+1,j,k)==1) then
          omega(i,j,k)=0.
         else
          omega(i,j,k)=deds+xxi(k)*dhsds
         end if
        end do
       end do
      end do
!
      do i=1,nx
       do j=0,ny
        dedn=(eta(i,j+1)-eta(i,j))/dyj(i,j)
        dhsdn=(hs(i,j+1)-hs(i,j))/dyj(i,j)
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i,j+1,k)==1) then
          psi(i,j,k)=0.
         else
          psi(i,j,k)=dedn+xxi(k)*dhsdn
         end if
        end do
       end do
      end do
!
      end subroutine omgpsi

!*******************************
      subroutine qcal
!*******************************
      real(8)::hs_upp
      hs_upp=0.

      do i=0,nx
       q(i)=0.
       do j=1,ny
        hs_upp=(hs(i+1,j)+hs(i,j))*.5
        qu(i,j)=0.
        do k=1,nz
         if(obst3d(i,j,k)==0.and.obst3d(i+1,j,k)==0) &
          qu(i,j)=qu(i,j)+yu(i,j,k)*hs_upp*dz(k)
        end do
        qu(i,j)=qu(i,j)*dy(i,j)
        q(i)=q(i)+qu(i,j)
       end do
      end do
      end subroutine qcal

!*******************************
      subroutine uupcal
!*******************************
      real(8)::hs_u

      q(0)=0.
      do j=1,ny
       hs_u=(hs(1,j)+hs(0,j))*.5
       qu(0,j)=0.
       do k=1,nz
        if(obst3d(0,j,k)==0.and.obst3d(1,j,k)==0) then
         uup0(j,k)=u00(k)
        else
         uup0(j,k)=0.
        end if
        qu(0,j)=qu(0,j)+uup0(j,k)*hs_u*dz(k)
       end do
       qu(0,j)=qu(0,j)*dy(0,j)
       q(0)=q(0)+qu(0,j)
      end do
      if(qp>0.) then
       qadjust=1.+(q(0)-qp0)/qp
      else
       qadjust=0.
      end if

      do j=1,ny
       do k=1,nz
        uup0(j,k)=uup0(j,k)*qadjust
       end do
      end do
      end subroutine uupcal

!*******************************
      subroutine ustacal
!*******************************
      real(8)::ubx,vbx,uv2,hsx,usx

      do i=1,nx
       do j=1,ny
        ubx=(yu(i,j,1)+yu(i-1,j,1))*.5
        vbx=(yv(i,j,1)+yu(i,j-1,1))*.5
        uv2=ubx**2+vbx**2
        usta(i,j)=sqrt(uv2)*cd0
        do k=1,nz
         snu_t(i,j,k)=ss_nu(usta(i,j),hs(i,j),xxi(k)) !セル中心
        end do
       end do
      end do

      do k=1,nz
       do i=1,nx
        do j=0,ny
         if(j==0) then
          hsx=hs(i,j+1)
          usx=usta(i,j+1)
         else if(j==ny) then
          hsx=hs(i,j)
          usx=usta(i,j)
         else
          hsx=(hs(i,j)+hs(i,j+1))*.5
          usx=(usta(i,j)+usta(i,j+1))*.5
         end if
         snu_xi_eg(i,j,k)=ss_nu(usx,hsx,xi(k)) !xi 辺上
        end do
       end do
       do i=0,nx
        do j=1,ny
         if(i==0) then
          hsx=hs(i+1,j)
          usx=usta(i+1,j)
         else if(i==nx) then
          hsx=hs(i,j)
          usx=usta(i,j)
         else
          hsx=(hs(i,j)+hs(i+1,j))*.5
          usx=(usta(i,j)+usta(i+1,j))*.5
         end if
         snu_et_eg(i,j,k)=ss_nu(usx,hsx,xi(k)) !eta 辺上
        end do
       end do
      end do

! Boundary Condition
!
!West and East

       do j=1,ny
        do k=0,nz+1
         if(j_west<=2) then
          snu_t(0,j,k)=snu_t(1,j,k)
         else if(j_west==3) then
          snu_t(0,j,k)=snu_t(nx,j,k)
          snu_t(nx+1,j,k)=snu_t(1,j,k)
         end if
         if(j_east<=2) then
          snu_t(nx+1,j,k)=snu_t(nx,j,k)
         end if
        end do
       end do
! Sounth and North
       do i=1,nx
        do k=0,nz+1
         if(j_south<=2) then
          snu_t(i,0,k)=snu_t(i,1,k)
         else if(j_south==3) then
          snu_t(i,0,k)=snu_t(i,ny,k)
          snu_t(i,ny+1,k)=snu_t(i,1,k)
         end if
         if(j_north<=2) then
          snu_t(i,ny+1,k)=snu_t(i,ny,k)
         end if
        end do
       end do

      end subroutine ustacal

!*******************************
      subroutine hcal
!*******************************
      real(8)::dhdt
!     real(8)::h_dw_t,delh
!     integer::j_t

!
      do i=1,nx
       do j=1,ny
        if(obst3d(i,j,nz)==0) then
         dhdt=w1(i,j,nz)
         hn(i,j)=h(i,j)+dhdt*dt*alpha_surf
         hs(i,j)=hn(i,j)-eta(i,j)
        end if
       end do
      end do

!
! Boundary Condition of Water Surface Elevation
!
!     write(*,*) j_west,j_hup,h_up,eta(0,1)
!     write(*,*) j_east,j_hdw,h_dw,eta(nx+1,1)
! West 
      do j=1,ny
       if(j_west==1) then !上流閉鎖
        if(j_hup==1) then !一定値
         write(*,*) 'Impossible(1)'!閉鎖で与えるのは矛盾
         stop
        end if
        hn(0,j)=hn(1,j)
        hs(0,j)=hs(1,j)
       else if(j_west==2) then !上流開放
        if(j_hup==1) then !一定値
         hn(0,j)=h_up
         hs(0,j)=h_up-eta(0,j)
        else if(j_hup==2) then !水平
         hn(0,j)=hn(1,j)
         hs(0,j)=hs(1,j)
        else if(j_hup==3) then
         hn(0,j)=h_up
         hs(0,j)=h_up-eta(0,j)
        else if(j_hup==4) then
         hn(0,j)=h_up+up_slope*ds0_center
         hs(0,j)=h_up-eta(0,j)
        end if
       else if(j_west==3) then !周期 
        hn(0,j)=hn(nx,j)
        hs(0,j)=hs(nx,j)
        hn(nx+1,j)=hn(1,j)
        hs(nx+1,j)=hs(1,j)
       end if
      end do
      dh_up0=hn(0,1)
! East
      do j=1,ny
       if(j_east==1) then !閉鎖
        if(j_hdw==1 .or. j_hdw>=3) then !一定値
         write(*,*) 'Impossible(2)'!閉鎖で与えるのは矛盾
         stop
        end if
        hn(nx+1,j)=hn(nx,j)
        hs(nx+1,j)=hs(nx,j)
       else if(j_east==2) then !下流開放
        if(j_hdw==1.or.j_hdw>=3) then !数値を与える
!        hn(nx+1,j)=hn(nx,j)
         hn(nx+1,j)=h_dw
         hs(nx+1,j)=h_dw-eta(nx+1,j)
        else if(j_hdw==2) then !水平
         hn(nx+1,j)=hn(nx,j)
         hs(nx+1,j)=hs(nx,j)
        end if
       end if
      end do

!     if(j_hdw==1.or.j_hdw>=3) then
!      j_t=0
!      h_dw_t=0.
!      do j=1,ny
!       if(obst(nx,j)==0) then
!        j_t=j_t+1
!        h_dw_t=h_dw_t+hn(nx,j)
!       end if
!      end do
!      h_dw_t=h_dw_t/dble(j_t)
!      delh=h_dw_t-h_dw
!      do j=1,ny
!       hn(nx+1,j)=hn(nx+1,j)-delh
!      end do
!     end if

!     j=2
!     write(44,'(5e12.3)') time,hn(nx-2,j),hn(nx-1,j),hn(nx,j),hn(nx+1,j)

! South
      do i=0,nx+1
       if(j_south<=2) then
        hn(i,0)=hn(i,1)
        hs(i,0)=hs(i,1)
       else if(j_south==3) then
        hn(i,0)=hn(i,ny)
        hs(i,0)=hs(i,ny)
        hn(i,ny+1)=hn(i,1)
        hs(i,ny+1)=hs(i,1)
       end if
      end do
! North
      do i=0,nx+1
       if(j_north<=2) then
        hn(i,ny+1)=hn(i,ny)
        hs(i,ny+1)=hs(i,ny)
       end if
      end do

!Obstacle
      do i=1,nx
       do j=1,ny
        if(obst3d(i,j,nz)==0) then
         if(obst3d(i+1,j,nz)==1) then
          hn(i+1,j)=hn(i,j)
          hs(i+1,j)=hs(i,j)
         end if
         if(obst3d(i-1,j,nz)==1) then
          hn(i-1,j)=hn(i,j)
          hs(i-1,j)=hs(i,j)
         end if
         if(obst3d(i,j+1,nz)==1) then
          hn(i,j+1)=hn(i,j)
          hs(i,j+1)=hs(i,j)
         end if
         if(obst3d(i,j-1,nz)==1) then
          hn(i,j-1)=hn(i,j)
          hs(i,j-1)=hs(i,j)
         end if
        end if
       end do
      end do

      end subroutine hcal

!*******************************
      subroutine hshift
!*******************************
      real(8)::cell_center_height
      
      do i=1,nx
       do j=1,ny
        h(i,j)=hn(i,j)
        hs(i,j)=hn(i,j)-eta(i,j)
        if(obst(i,j)==1) then
         do k=1,nz
          cell_center_height=eta(i,j)+hs(i,j)*xxi(k)
          if(cell_center_height>hobst(i,j)) then
           obst3d(i,j,k)=0
          else
           obst3d(i,j,k)=1
           if(k==1) obst3d(i,j,0)=1
          end if
         end do
        end if
       end do
      end do

      do i=1,nx
       j=0
       h(i,j)=hn(i,j)
       hs(i,j)=hn(i,j)-eta(i,j)
       j=ny+1
       h(i,j)=hn(i,j)
       hs(i,j)=hn(i,j)-eta(i,j)
      end do

      do j=0,ny+1
       i=0
       h(i,j)=hn(i,j)
       hs(i,j)=hn(i,j)-eta(i,j)
       i=nx+1
       h(i,j)=hn(i,j)
       hs(i,j)=hn(i,j)-eta(i,j)
      end do

!Obstacle
      do i=1,nx
       do j=1,ny
        if(obst3d(i,j,nz)==0) then
         if(obst3d(i+1,j,nz)==1) then
          hs(i+1,j)=hs(i,j)
          h(i+1,j)=h(i,j)
         end if
         if(obst3d(i-1,j,nz)==1) then
          hs(i-1,j)=hs(i,j)
          h(i-1,j)=h(i,j)
         end if
         if(obst3d(i,j+1,nz)==1) then
          hs(i,j+1)=hs(i,j)
          h(i,j+1)=h(i,j)
         end if
         if(obst3d(i,j-1,nz)==1) then
          hs(i,j-1)=hs(i,j)
          h(i,j-1)=h(i,j)
         end if
        end if
       end do
      end do
      end subroutine hshift
      end module non_advection_m
