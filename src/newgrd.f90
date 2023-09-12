module newgrd_m
	use common_hh
	use common_geom
	implicit none
contains
	!*******************************
	subroutine newgrd_u(yun,yu,gux_n,guy_n,guz_n,gux,guy,guz)
	!*******************************
		real(8),dimension(-1:im,-1:jm,0:km)::yun,yu,gux,guy,guz,gux_n,guy_n,guz_n

		do i=0,nx
			do j=1,ny
				do k=1,nz
					gux_n(i,j,k)=gux(i,j,k) &
								+(yun(i+1,j,k)-yun(i-1,j,k)-yu(i+1,j,k)+yu(i-1,j,k)) &
								*0.5/dxi(i,j)
					guy_n(i,j,k)=guy(i,j,k) &
								+(yun(i,j+1,k)-yun(i,j-1,k)-yu(i,j+1,k)+yu(i,j-1,k)) &
								*0.5/dy(i,j)
					guz_n(i,j,k)=guz(i,j,k) &
								+(yun(i,j,k+1)-yun(i,j,k-1)-yu(i,j,k+1)+yu(i,j,k-1))/ &
								(dzk(k-1)+dzk(k))
				end do
			end do
		end do

		!
		! North and South Boundary
		!
		do i=0,nx
			do k=0,nz+1
				if(j_south<=2) then  !South Closed or Open
					gux_n(i,0,k)=0.
					guy_n(i,0,k)=0.
					guz_n(i,0,k)=0.
				else if(j_south==3) then !South and North Periodic
					gux_n(i,0,k)=gux_n(i,ny,k)
					guy_n(i,0,k)=guy_n(i,ny,k)
					guz_n(i,0,k)=guz_n(i,ny,k)
					gux_n(i,ny+1,k)=gux_n(i,1,k)
					guy_n(i,ny+1,k)=guy_n(i,1,k)
					guz_n(i,ny+1,k)=guz_n(i,1,k)
				end if
				if(j_north<=2) then !North Closed or Open
					gux_n(i,ny+1,k)=0.
					guy_n(i,ny+1,k)=0.
					guz_n(i,ny+1,k)=0.
				end if
			end do
		end do

		! West and East Boundary
		do j=1,ny
			do k=0,nz+1
				if(j_west==1) then ! West Closed
					gux_n(0,j,k)=0.
					guy_n(0,j,k)=0.
					guz_n(0,j,k)=0.
				else if(j_west==2) then !West Open
					gux_n(-1,j,k)=0.
					guy_n(-1,j,k)=0.
					guz_n(-1,j,k)=0.
				else if(j_west==3) then !West and East Periodic
					gux_n(0,j,k)   =gux_n(nx,j,k)
					gux_n(nx+1,j,k)=gux_n(1,j,k)
					guy_n(0,j,k)   =guy_n(nx,j,k)
					guy_n(nx+1,j,k)=guy_n(1,j,k)
					guz_n(0,j,k)   =guz_n(nx,j,k)
					guz_n(nx+1,j,k)=guz_n(1,j,k)
				end if
				if(j_east==1) then !East Closed
					gux_n(nx,j,k)=0.
					guy_n(nx,j,k)=0.
					guz_n(nx,j,k)=0.
				else if(j_east==2) then !East Open
					gux_n(nx+1,j,k)=0.
					guy_n(nx+1,j,k)=0.
					guz_n(nx+1,j,k)=0.
				end if
			end do
		end do
		! Obstacles
		do i=0,nx
			do j=1,ny
				do k=0,nz+1
					if(obst3d(i,j,k)==1.or.obst3d(i+1,j,k)==1) then
						gux_n(i,j,k)=0.; guy_n(i,j,k)=0.; guz_n(i,j,k)=0.
					end if
				end do
			end do
		end do
	end subroutine newgrd_u

	!*******************************
	subroutine newgrd_v(yvn,yv,gvx_n,gvy_n,gvz_n,gvx,gvy,gvz)
	!*******************************
		real(8),dimension(-1:im,-1:jm,0:km)::yvn,yv,gvx,gvy,gvz,gvx_n,gvy_n,gvz_n

		do i=1,nx
			do j=0,ny
				do k=1,nz
					gvx_n(i,j,k)=gvx(i,j,k) &
								+(yvn(i+1,j,k)-yvn(i-1,j,k)-yv(i+1,j,k)+yv(i-1,j,k)) &
								*.5/dx(i,j)
					gvy_n(i,j,k)=gvy(i,j,k) &
								+(yvn(i,j+1,k)-yvn(i,j-1,k)-yv(i,j+1,k)+yv(i,j-1,k)) &
								*.5/dyj(i,j)
					gvz_n(i,j,k)=gvz(i,j,k) &
								+(yvn(i,j,k+1)-yvn(i,j,k-1)-yv(i,j,k+1)+yv(i,j,k-1)) &
								/(dzk(k-1)+dzk(k))
				end do
			end do
		end do

		!
		! West and East Boundary
		!
		do j=0,ny
			do k=0,nz+1
				if(j_west<=2) then !West Closed or Open
					gvx_n(0,j,k)=0.
					gvy_n(0,j,k)=0.
					gvz_n(0,j,k)=0.
				else if(j_west==3) then !West and East Periodic
					gvx_n(0,j,k)   =gvx_n(nx,j,k)
					gvx_n(nx+1,j,k)=gvx_n(1,j,k)
					gvy_n(0,j,k)   =gvy_n(nx,j,k)
					gvy_n(nx+1,j,k)=gvy_n(1,j,k)
					gvz_n(0,j,k)   =gvz_n(nx,j,k)
					gvz_n(nx+1,j,k)=gvz_n(1,j,k)
				end if
				if(j_east<=2) then !East Closed or Open
					gvx_n(nx+1,j,k)=0.
					gvy_n(nx+1,j,k)=0.
					gvz_n(nx+1,j,k)=0.
				end if
			end do
		end do
		!
		! South and North
		!
		do i=1,nx
			do k=0,nz+1
				if(j_south==1) then !South Closed
					gvx_n(i,0,k)=0.
					gvy_n(i,0,k)=0.
					gvz_n(i,0,k)=0.
				else if(j_south==2) then !South Open
					gvx_n(i,-1,k)=0.
					gvy_n(i,-1,k)=0.
					gvz_n(i,-1,k)=0.
				else if(j_south==3) then !South and North Periodic
					gvx_n(i,0,k)   =gvx_n(i,ny,k)
					gvx_n(i,ny+1,k)=gvx_n(i,1,k)
					gvy_n(i,0,k)   =gvy_n(i,ny,k)
					gvy_n(i,ny+1,k)=gvy_n(i,1,k)
					gvz_n(i,0,k)   =gvz_n(i,ny,k)
					gvz_n(i,ny+1,k)=gvz_n(i,1,k)
				end if
				if(j_north==1) then ! North Closed
					gvx_n(i,ny,k)=0.
					gvy_n(i,ny,k)=0.
					gvz_n(i,ny,k)=0.
				else if(j_north==2) then ! North Open
					gvx_n(i,ny+1,k)=0.
					gvy_n(i,ny+1,k)=0.
					gvz_n(i,ny+1,k)=0.
				end if
			end do
		end do
		! Obstacles
		do i=1,nx
			do j=0,ny
				do k=0,nz+1
					if(obst3d(i,j,k)==1.or.obst3d(i,j+1,k)==1) then
						gvx_n(i,j,k)=0.; gvy_n(i,j,k)=0.; gvz_n(i,j,k)=0.
					end if
				end do
			end do
		end do

	end subroutine newgrd_v

	!*******************************
	subroutine newgrd_w(ywn,yw,gwx_n,gwy_n,gwz_n,gwx,gwy,gwz)
	!*******************************
		real(8),dimension(-1:im,-1:jm,0:km)::ywn,yw,gwx,gwy,gwz,gwx_n,gwy_n,gwz_n

		do i=1,nx
			do j=1,ny
				do k=1,nz-1
					gwx_n(i,j,k)=gwx(i,j,k) &
								+(ywn(i+1,j,k)-ywn(i-1,j,k)-yw(i+1,j,k)+yw(i-1,j,k)) &
								/(dxi(i,j)+dxi(i-1,j))
					gwy_n(i,j,k)=gwy(i,j,k) &
								+(ywn(i,j+1,k)-ywn(i,j-1,k)-yw(i,j+1,k)+yw(i,j-1,k)) &
								/(dyj(i,j)+dyj(i,j-1))
					gwz_n(i,j,k)=gwz(i,j,k) &
								+(ywn(i,j,k+1)-ywn(i,j,k-1)-yw(i,j,k+1)+yw(i,j,k-1)) &
								/(dz(k)+dz(k+1))
				end do
				gwx_n(i,j,0)=0.; gwx_n(i,j,nz)=0.
				gwy_n(i,j,0)=0.; gwy_n(i,j,nz)=0.
				gwz_n(i,j,0)=0.; gwz_n(i,j,nz)=0.
			end do
		end do
		!
		! West and East Condition
		!
		do j=1,ny
			do k=0,nz+1
				if(j_west<=2) then
					gwx_n(0,j,k)=0.
					gwy_n(0,j,k)=0.
					gwz_n(0,j,k)=0.
				else if(j_west==3) then
					gwx_n(0,j,k)   =gwx_n(nx,j,k)
					gwx_n(nx+1,j,k)=gwx_n(1,j,k)
					gwy_n(0,j,k)   =gwy_n(nx,j,k)
					gwy_n(nx+1,j,k)=gwy_n(1,j,k)
					gwz_n(0,j,k)   =gwz_n(nx,j,k)
					gwz_n(nx+1,j,k)=gwz_n(1,j,k)
				end if
				if(j_east<=2) then
					gwx_n(nx+1,j,k)=0.
					gwy_n(nx+1,j,k)=0.
					gwy_n(nx+1,j,k)=0.
				end if
			end do
		end do
		!
		! North and South Condition
		!
		do i=1,nx
			do k=0,nz+1
				if(j_south<=2) then
					gwx_n(i,0,k)=0.
					gwy_n(i,0,k)=0.
					gwz_n(i,0,k)=0.
				else if(j_south==3) then
					gwx_n(i,0,k)   =gwx_n(i,ny,k)
					gwx_n(i,ny+1,k)=gwx_n(i,1,k)
					gwy_n(i,0,k)   =gwy_n(i,ny,k)
					gwy_n(i,ny+1,k)=gwy_n(i,1,k)
					gwz_n(i,0,k)   =gwz_n(i,ny,k)
					gwz_n(i,ny+1,k)=gwz_n(i,1,k)
				end if
				if(j_north<=2) then
					gwx_n(i,ny+1,k)=0.
					gwy_n(i,ny+1,k)=0.
					gwz_n(i,ny+1,k)=0.
				end if
			end do
		end do
		! Obstacles
		do i=1,nx
			do j=1,ny
				do k=0,nz
					if(obst3d(i,j,k)==1.or.obst3d(i,j,k+1)==1) then
						gwx_n(i,j,k)=0.; gwy_n(i,j,k)=0.; gwz_n(i,j,k)=0.
					end if
				end do
			end do
		end do
	end subroutine newgrd_w

	!*******************************
	subroutine newgrd_c(ycn,yc,gcx_n,gcy_n,gcz_n,gcx,gcy,gcz)
	!*******************************
		real(8),dimension(-1:im,-1:jm,0:km)::ycn,yc,gcx,gcy,gcz,gcx_n,gcy_n,gcz_n

		do i=1,nx
			do j=1,ny
				do k=1,nz
					gcx_n(i,j,k)=gcx(i,j,k) &
								+(ycn(i+1,j,k)-ycn(i-1,j,k)-yc(i+1,j,k)+yc(i-1,j,k)) &
								/(dxi(i,j)+dxi(i-1,j))
					gcy_n(i,j,k)=gcy(i,j,k) &
								+(ycn(i,j+1,k)-ycn(i,j-1,k)-yc(i,j+1,k)+yc(i,j-1,k)) &
								/(dyj(i,j)+dyj(i,j-1))
					gcz_n(i,j,k)=gcz(i,j,k) &
								+(ycn(i,j,k+1)-ycn(i,j,k-1)-yc(i,j,k+1)+yc(i,j,k-1)) &
								/(dzk(k)+dzk(k-1))
				end do
			end do
		end do
		!
		! North and South Boundary
		!
		do i=1,nx
			do k=1,nz
				if(j_south==3) then
					gcx_n(i,ny+1,k)=gcx_n(i,1,k)
					gcx_n(i,0,k)=gcx_n(i,ny,k)
					gcy_n(i,ny+1,k)=gcy_n(i,1,k)
					gcy_n(i,0,k)=gcy_n(i,ny,k)
					gcz_n(i,ny+1,k)=gcz_n(i,1,k)
					gcz_n(i,0,k)=gcz_n(i,ny,k)
				else
					gcx_n(i,ny+1,k)=0.
					gcx_n(i,0,k)=0.
					gcy_n(i,ny+1,k)=0.
					gcy_n(i,0,k)=0.
					gcz_n(i,ny+1,k)=0.
					gcz_n(i,0,k)=0.
				end if
			end do
		end do
		!
		! East and West Boundary
		!
		do j=1,ny
			do k=1,nz
				if(j_west==3) then
					gcx_n(0,j,k)=gcx_n(nx,j,k)
					gcx_n(nx+1,j,k)=gcx_n(1,j,k)
					gcy_n(0,j,k)=gcy_n(nx,j,k)
					gcy_n(nx+1,j,k)=gcy_n(1,j,k)
					gcz_n(0,j,k)=gcz_n(nx,j,k)
					gcz_n(nx+1,j,k)=gcz_n(1,j,k)
				else
					gcx_n(0,j,k)=0.
					gcx_n(nx+1,j,k)=0.
					gcy_n(0,j,k)=0.
					gcy_n(nx+1,j,k)=0.
					gcz_n(0,j,k)=0.
					gcz_n(nx+1,j,k)=0.
				end if
			end do
		end do
		!
		! Surface and Bottom Boundary
		!
		do i=1,nx
			do j=1,ny
				gcx_n(i,j,nz+1)=0.
				gcx_n(i,j,0)=0.
				gcy_n(i,j,nz+1)=0.
				gcy_n(i,j,0)=0.
				gcz_n(i,j,nz+1)=0.
				gcz_n(i,j,0)=0.
			end do
		end do

	end subroutine newgrd_c
end module newgrd_m
