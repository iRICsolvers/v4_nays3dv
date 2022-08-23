      module bound_m
      use common_hh
      use common_geom
      use common_hyd
      implicit none
      contains
!*******************************
      subroutine bound_u(u,w)
!*******************************
!
      real(8),dimension(-1:im,-1:jm,0:km)::u,w
      real(8)::dwds,dhds,hs_up,dwdxi,duds,dd
      real(8)::dxii
!
      do i=0,nx
       do j=1,ny
!       if(i==0.or.i==nx.or.j_surf==0.or.(j_surf==1.and.time<stime_surf)) then
        if(j_surf==0.or.(j_surf==1.and.time<stime_surf)) then
         u(i,j,nz+1)=u(i,j,nz)
        else
         dwds=(w(i+1,j,nz)-w(i,j,nz))/dxi(i,j)
         dhds=(hn(i+1,j)-hn(i,j))/dxi(i,j)
         hs_up=(hs(i+1,j)+hs(i,j))*.5
         dwdxi=(w(i+1,j,nz)+w(i,j,nz)-w(i+1,j,nz-1)-w(i,j,nz-1)) &
            /(dz(nz)*2.)
         if(i==0) then
          dxii=(dx(i+1,j)+dx(i+1,j-1))*.5
          duds=(u(i+1,j,nz)-u(i,j,nz))/dxii
         else if(i==nx) then
          dxii=(dx(i,j)+dx(i,j-1))*.5
          duds=(u(i,j,nz)-u(i-1,j,nz))/dxii
         else
          duds=(u(i+1,j,nz)-u(i-1,j,nz))/(2.*dxi(i,j))
         end if
         dd=-dwds+dhds*(duds-dwdxi/hs_up)
         u(i,j,nz+1)=u(i,j,nz)+dd*dz(nz)*hs_up
        end if
        u(i,j,0)=-u(i,j,1)  !Bottom Nonskip
       end do
      end do
!
! North and South Boundary
!
      do i=0,nx
       do k=0,nz+1
        if(j_south<=2) then  !South
         u(i,0,k)=u(i,1,k) !Slip
        else if(j_south==3) then !South and North Periodic
         u(i,0,k)=u(i,ny,k)
         u(i,ny+1,k)=u(i,1,k)
        end if
        if(j_north<=2) then !North Closed or Open
         u(i,ny+1,k)=u(i,ny,k) !Slip
        end if
       end do
      end do
!
! West and East Boundary
!
      do j=0,ny+1
       do k=0,nz+1
        if(j_west==1) then !West Closed
         u(0,j,k)=0.
         u(-1,j,k)=0.
        else if(j_west==2) then !West Open è„ó¨äJï˙
         if(j_qin==0) then !è„ó¨ÇÕÉtÉäÅ[
          u(-1,j,k)=u(0,j,k)
         else if(j_qin>=1) then !ó¨ó Çó^Ç¶ÇÈ-->ó¨ë¨Çó^Ç¶ÇÈÅD
          if(obst(1,j)==1) then
           u(-1,j,k)=0.
           u(0,j,k)=0.
          else
           u(0,j,k)=uup0(j,k)
           u(-1,j,k)=u(0,j,k)
          end if
         end if
        else if(j_west==3) then !West and East Periodic
         u(-1,j,k)=u(nx-1,j,k)
         u(nx+1,j,k)=u(1,j,k)
        end if
        if(j_east==1) then !East Closed
         u(nx,j,k)=0.
         u(nx+1,j,k)=0.
        else if(j_east==2) then
!        if(u(nx,j,k)<0) u(nx,j,k)=0. !â∫ó¨Ç≈u<0ÇÕîFÇﬂÇ»Ç¢
         u(nx+1,j,k)=u(nx,j,k)
        end if
       end do
      end do
!
! Obstacles 
      do i=0,nx
       do j=1,ny
        do k=1,nz+1
         if(obst3d(i,j,k)==1 .or. obst3d(i+1,j,k)==1) then
          if(j==1) then
           u(i,j,k)=0.
           u(i,j-1,k)=0.
          else if(j==ny) then
           u(i,j,k)=0.
           u(i,j+1,k)=0.
          else
           u(i,j,k)=0.
          end if
         end if
        end do
       end do
      end do

      end subroutine bound_u

!*******************************
      subroutine bound_v(v,w)
!*******************************
!
      real(8),dimension(-1:im,-1:jm,0:km)::v,w
      real(8)::dwdn,dhdn,hs_vp,dwdxi,dvdn,dd
      real(8)::dyjj
!
      do i=1,nx
       do j=0,ny
!       if(j==0.or.j==ny.or.j_surf==0.or.(j_surf==1.and.time<stime_surf)) then
        if(j_surf==0.or.(j_surf==1.and.time<stime_surf)) then
         v(i,j,nz+1)=v(i,j,nz)
        else
         dwdn=(w(i,j+1,nz)-w(i,j,nz))/dyj(i,j)
         dhdn=(hn(i,j+1)-hn(i,j))/dyj(i,j)
         hs_vp=(hs(i,j+1)+hs(i,j))*.5
         dwdxi=(w(i,j+1,nz)+w(i,j,nz)-w(i,j+1,nz-1)-w(i,j,nz-1)) &
           /(dz(nz)*2.)
         if(j==0) then
          dyjj=(dy(i,j+1)+dy(i-1,j+1))*.5
          dvdn=(v(i,j+1,nz)-v(i,j,nz))/dyjj
         else if(j==ny) then
          dyjj=(dy(i,j)+dy(i-1,j))*.5
          dvdn=(v(i,j,nz)-v(i,j-1,nz))/dyjj
         else
          dvdn=(v(i,j+1,nz)-v(i,j-1,nz))/(dyj(i,j)*2.)
         end if
         dd=-dwdn+dhdn*(dvdn-dwdxi/hs_vp)
         v(i,j,nz+1)=v(i,j,nz)+dd*dz(nz)*hs_vp
        end if
        v(i,j,0)=-v(i,j,1) !Bottom Nonslip
       end do
      end do
!
! West and East Boundary
!
      do j=0,ny
       do k=0,nz+1
        if(j_west==1) then !West Closed or Open
         v(0,j,k)=v(1,j,k) !Slip
        else if(j_west==2) then 
!        v(0,j,k)=v(1,j,k)
         v(0,j,k)=0.
        else if(j_west==3) then !West and East Periodic
         v(0,j,k)=v(nx,j,k)
         v(nx+1,j,k)=v(1,j,k)
        end if
        if(j_east==1) then !East Closed or Open
         v(nx+1,j,k)=-v(nx,j,k)
        else if(j_east==2) then 
!        v(nx+1,j,k)=v(nx,j,k)
         v(nx+1,j,k)=0.
        end if
       end do
      end do
!
! South and North
!
      do i=1,nx
       do k=0,nz+1
        if(j_south==1) then !South Closed
         v(i,0,k)=0.
         v(i,-1,k)=0.
        else if(j_south==2) then !South Open
         v(i,-1,k)=v(i,0,k)
        else if(j_south==3) then !South and North Periodic
         v(i,0,k)=v(i,ny,k)
         v(i,ny+1,k)=v(i,1,k)
        end if
        if(j_north==1) then !North Closed
         v(i,ny,k)=0.
         v(i,ny+1,k)=0.
        else if(j_north==2) then !North Open
         v(i,ny+1,k)=v(i,ny,k)
        end if
       end do
      end do
!
! Obstacles 
!
      do i=1,nx
       do j=1,ny-1
        do k=1,nz+1
         if(obst3d(i,j,k)==1 .or. obst3d(i,j+1,k)==1) v(i,j,k)=0.
        end do
       end do
      end do

      end subroutine bound_v

!*******************************
      subroutine bound_w(u,v,w)
!*******************************
!
      real(8),dimension(-1:im,-1:jm,0:km)::u,v,w
      real(8)::dhsds,dhsdn,deds,dedn
      real(8)::u_wp,v_wp,ome_wp,psi_wp
!
! Water Surface Condition
!
      do i=1,nx
       do j=1,ny
        if(j_surf==1 .and. time>stime_surf) then
         w(i,j,nz+1)=3.*w(i,j,nz)-3.*w(i,j,nz-1)+w(i,j,nz-2)
        else
         w(i,j,nz+1)=0.
         w(i,j,nz)=0.
        end if
        dhsds=(hs(i+1,j)-hs(i-1,j))/(dxi(i,j)+dxi(i-1,j))
        dhsdn=(hs(i,j+1)-hs(i,j-1))/(dyj(i,j)+dyj(i,j-1))
        deds=(eta(i+1,j)-eta(i-1,j))/(dxi(i,j)+dxi(i-1,j))
        dedn=(eta(i,j+1)-eta(i,j-1))/(dyj(i,j)+dyj(i,j-1))
        k=0
        u_wp=(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1))*.25
        v_wp=(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1))*.25
        ome_wp=xi(k)*dhsds+deds
        psi_wp=xi(k)*dhsdn+dedn
        w(i,j,k)=ome_wp*u_wp+psi_wp*v_wp
       end do
      end do
!
! West and East Condition
!
      do j=1,ny
       do k=0,nz+1
        if(j_west<=2) then ! West Closed or Open
         w(0,j,k)=w(1,j,k)
        else if(j_west==3) then !Periodic Condition
         w(0,j,k)=w(nx,j,k)
         w(nx+1,j,k)=w(1,j,k)
        end if
        if(j_east==1) then !East Closed
         w(nx+1,j,k)=w(nx,j,k)
        else if(j_east==2) then
         if(j_hdw==1) then
          w(nx+1,j,k)=0.
         else if(j_hdw>=2) then
           w(nx+1,j,k)=w(nx,j,k)
         else
          if(k==nz+1) then
           w(nx+1,j,k)=-w(nx,j,k)+dhdt_dw
          else
           w(nx+1,j,k)=-w(nx,j,k)+dhdt_dw*xi(k)
          end if
         end if
        end if
       end do
      end do
!
! North and South Condition
!
      do i=1,nx
       do k=0,nz+1
        if(j_south<=2) then
         w(i,0,k)=w(i,1,k)
        else if(j_south==3) then
         w(i,0,k)=w(i,ny,k)
         w(i,ny+1,k)=w(i,1,k)
        end if
        if(j_north<=2) then
         w(i,ny+1,k)=w(i,ny,k)
        end if
       end do
      end do
!
! Obstacles 
!
      do i=1,nx
       do j=1,ny
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i,j,k+1)==1) w(i,j,k)=0.
        end do
       end do
      end do
      end subroutine bound_w

!*******************************
      subroutine bound_c(c,u,v,w)
!*******************************
!
      real(8),dimension(-1:im,-1:jm,0:km)::u,v,w
      real(8),dimension(-1:im,-1:jm,0:km)::c
      real(8)::elh
!
! North and South Boundary
!
      do i=1,nx
       do k=1,nz
        if(j_south==1) then
         c(i,0,k)=c(i,1,k)
        else if(j_south==2) then
         c(i,0,k)=c0
        else if(j_south==3) then
         c(i,ny+1,k)=c(i,1,k)
         c(i,0,k)=c(i,ny,k)
        end if
        if(j_north==1) then
         c(i,ny+1,k)=c(i,ny,k)
        else if(j_north==2) then
         c(i,ny+1,k)=c0
        end if
       end do
      end do
!
! East and West Boundary
! 
      do j=1,ny
       do k=1,nz
        if(j_west==1) then !West closed
         c(0,j,k)=c(1,j,k)
        else if(j_west==2) then !West open
         c(0,j,k)=c0
        else if(j_west==3) then !periodic
         c(0,j,k)=c(nx,j,k)
         c(nx+1,j,k)=c(1,j,k)
        end if
        if(j_east==1) then !East closed
         c(nx+1,j,k)=c(nx,j,k)
        else if(j_east==2) then !East open
         c(nx+1,j,k)=c0
        end if
       end do
      end do
!
! Surface and Bottom Boundary
!
      do i=1,nx
       do j=1,ny
        c(i,j,nz+1)=c(i,j,nz)
        c(i,j,0)=c(i,j,1)
       end do
      end do
!
! Obstacles
!
      do i=1,nx
       do j=1,ny
        do k=1,nz
         if(obst3d(i,j,k)==1) then
          if(i<nx .and. obst3d(i+1,j,k)==0) c(i,j,k)=c(i+1,j,k)
          if(i>1  .and. obst3d(i-1,j,k)==0) c(i,j,k)=c(i-1,j,k)
          if(j<ny .and. obst3d(i,j+1,k)==0) c(i,j,k)=c(i,j+1,k)
          if(j>1  .and. obst3d(i,j-1,k)==0) c(i,j,k)=c(i,j-1,k)
          if(k<nz .and. obst3d(i,j,k+1)==0) c(i,j,k)=c(i,j,k+1)
          if(k>1  .and. obst3d(i,j,k-1)==0) c(i,j,k)=c(i,j,k-1)
         end if
        end do
       end do
      end do
!
! ñßìxã´äEèåè
!
      if(j_dens==1.and.j_bc_dens==1) then

      if(jc_in_west>0) then
       do m=1,jc_in_west
        do j=jc_west_s(m),jc_west_e(m)
         do k=1,nz
!         if(k<=k_bc_west(m)) then
!          c(0,j,k)=c_bound_west(m)
!          if(yu(0,j,k)>0) c(1,j,k)=c(0,j,k)
!         else
!          c(0,j,k)=c0
!         end if
          elh=eta(1,j)+hs(1,j)*xxi(k)
          if(elh<=c_up_west(m).and.obst3d(1,j,k)==0) then
           c(0,j,k)=c_bound_west(m)
           if(j_west==1) c(1,j,k)=c_bound_west(m)
          else
           c(0,j,k)=c0
          end if
         end do
        end do
       end do
      end if  

      if(jc_in_east>0) then
       do m=1,jc_in_east
        do j=jc_east_s(m),jc_east_e(m)
         do k=1,nz
!         if(k<=k_bc_east(m)) then
!          c(nx+1,j,k)=c_bound_east(m)
!          if(yu(nx,j,k)<0) c(nx,j,k)=c(nx+1,j,k)
!         else
!          c(nx+1,j,k)=c0
!         end if
          elh=eta(nx,j)+hs(nx,j)*xxi(k)
          if(elh<=c_up_east(m).and.obst3d(nx,j,k)==0) then
           c(nx+1,j,k)=c_bound_east(m)
           if(j_east<=2) c(nx,j,k)=c_bound_east(m)
          else
           c(nx+1,j,k)=c0
          end if
!         write(47,'(2i5,3f10.4)') j,k,c(nx+1,j,k),c(nx,j,k),c(nx-1,j,k)
         end do
        end do
       end do
      end if  
!
      if(jc_in_south>0.and.j_south==2) then
       do m=1,jc_in_south
        do i=ic_south_s(m),ic_south_e(m)
         do k=1,nz
!         if(k<=k_bc_south(m)) then
!          c(i,0,k)=c_bound_south(m)
!          if(yv(i,0,k)>0) c(i,1,k)=c(i,0,k)
!         else
!          c(i,0,k)=c0
!         end if
          elh=eta(i,1)+hs(i,1)*xxi(k)
          if(elh<=c_up_south(m).and.obst3d(i,1,k)==0) then
           c(i,0,k)=c_bound_south(m)
           if(j_south==1) c(i,1,k)=c_bound_south(m)
          else
           c(i,0,k)=c0
          end if
         end do
        end do
       end do
      end if  

      if(jc_in_north>0.and.j_north==2) then
       do m=1,jc_in_north
        do i=ic_north_s(m),ic_north_e(m)
         do k=1,nz
!         if(k<=k_bc_north(m)) then
!          c(i,ny+1,k)=c_bound_north(m)
!          if(yv(i,ny,k)<0) c(i,ny,k)=c(i,ny+1,k)
!         else
!          c(i,ny+1,k)=c0
!         end if
          elh=eta(i,ny)+hs(i,ny)*xxi(k)
          if(elh<=c_up_north(m).and.obst3d(i,ny,k)==0) then
           c(i,ny+1,k)=c_bound_north(m)
           if(j_north==1) c(i,ny,k)=c_bound_north(m)
          else
           c(i,ny+1,k)=c0
          end if
         end do
        end do
       end do
      end if  

      end if

!     do j=1,ny
!      do k=1,nk
!         write(47,'(2i5,3f10.4)') j,k,c(nx+1,j,k),c(nx,j,k),c(nx-1,j,k)
!      end do
!     end do

      end subroutine bound_c

!*******************************
      subroutine bound_p(smg_g)
!*******************************
      real(8)::smg_g,dx1,dy1,ypi,ypj
   
      do i=1,nx
       do j=1,ny
        if(obst3d(i,j,nz)==1) then
         ypsurf(i,j)=0.
         ypn(i,j,nz+1)=0.
        else
         dx1=(dxi(i,j)+dxi(i-1,j))*.5
         ypi=(hn(i-1,j)-2.*hn(i,j)+hn(i+1,j))/dx1**2
         dy1=(dyj(i,j)+dyj(i,j-1))*.5
         ypj=(hn(i,j-1)-2.*hn(i,j)+hn(i,j+1))/dx1**2
         ypsurf(i,j)=smg_g*(ypi+ypj)
         ypn(i,j,nz+1)=2.*ypsurf(i,j)-ypn(i,j,nz)
        end if
       end do
      end do

! Obstacle

      do i=1,nx
       do j=1,ny
        if(obst3d(i,j,nz)==1) then
         if(i<nx .and. obst3d(i+1,j,nz)==0) ypn(i,j,nz+1)=ypn(i+1,j,nz)
         if(i>1  .and. obst3d(i-1,j,nz)==0) ypn(i,j,nz+1)=ypn(i-1,j,nz)
         if(j<nx .and. obst3d(1,j+1,nz)==0) ypn(i,j,nz+1)=ypn(i,j+1,nz)
         if(j>1  .and. obst3d(1,j-1,nz)==0) ypn(i,j,nz+1)=ypn(i,j-1,nz)
        end if
       end do
      end do
      end subroutine bound_p
!*******************************
      subroutine shift(yn,y)
!*******************************
      real(8),dimension(-1:im,-1:jm,0:km)::yn,y
      do j=0,ny+1
       do i=0,nx+1
        do k=0,nz+1
         y(i,j,k)=yn(i,j,k)
        end do
       end do
      end do
      end subroutine shift
      end module bound_m
