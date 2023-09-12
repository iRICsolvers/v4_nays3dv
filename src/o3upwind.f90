module o3upwind_m

	use common_hh
	use common_geom
	use common_hyd
	implicit none
	real(8)::dx1,dy1,dz1
	real(8)::hs_upp,hs_vp
	real(8)::u_dfdx,v_dfdy,w_dfdz

contains

	!***************************************************
	subroutine o3upwind_u(f,gx,gy,gz)
	!***************************************************

		real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
		real(8),dimension(-1:im,-1:jm,0:km)::fn,u,v,w

		fn=0.;  u=0.; v=0.; w=0.

		do i=0,nx
			do j=1,ny
				hs_upp=(hs(i,j)+hs(i+1,j))*.5
				do k=1,nz
					u(i,j,k)=yu(i,j,k)
					v(i,j,k)=(yv(i,j,k)+yv(i+1,j,k)+yv(i,j-1,k)+yv(i+1,j-1,k))*.25
					w(i,j,k)=(w2(i,j,k)+w2(i+1,j,k)+w2(i,j,k-1)+w2(i+1,j,k-1))*.25/hs_upp
				end do
			end do
		end do

		do i=0,nx
			do j=1,ny
				dx1=dxi(i,j)
				dy1=dy(i,j)
				do k=1,nz
					if(i==0) then
						u_dfdx=(u(i,j,k)-abs(u(i,j,k)))*(f(i+1,j,k)-f(i,j,k))/dx1
					else if(i==nx) then
						u_dfdx=(u(i,j,k)+abs(u(i,j,k)))*(f(i,j,k)-f(i-1,j,k))/dx1
					else if(i==1.or.i==nx-1) then
						u_dfdx=((u(i,j,k)+abs(u(i,j,k)))*(f(i,j,k)-f(i-1,j,k)) &
							  +(u(i,j,k)-abs(u(i,j,k)))*(f(i+1,j,k)-f(i,j,k)))/(2.*dx1)
					else
						u_dfdx=((u(i,j,k)+abs(u(i,j,k))) &
							  *(2.*f(i+1,j,k)+3.*f(i,j,k)-6.*f(i-1,j,k)+f(i-2,j,k)) &
							  +(u(i,j,k)-abs(u(i,j,k)))  &
							  *(2.*f(i+2,j,k)+3.*f(i+1,j,k)-6.*f(i,j,k)+f(i-1,j,k)))/(12.*dx1)
					end if
					if(j==1.or.j==ny) then
						v_dfdy=((v(i,j,k)+abs(v(i,j,k)))*(f(i,j,k)-f(i,j-1,k)) &
							  +(v(i,j,k)-abs(v(i,j,k)))*(f(i,j+1,k)-f(i,j,k)))/(2.*dy1)
					else
						v_dfdy=((v(i,j,k)+abs(v(i,j,k))) &
							  *(2.*f(i,j+1,k)+3.*f(i,j,k)-6.*f(i,j-1,k)+f(i,j-2,k)) &
							  +(v(i,j,k)-abs(v(i,j,k))) &
							  *(2.*f(i,j+2,k)+3.*f(i,j+1,k)-6.*f(i,j,k)+f(i,j-1,k)))/(12.*dy1)
					end if
					if(k==1.or.k==nz) then
						w_dfdz=((w(i,j,k)+abs(w(i,j,k)))*(f(i,j,k)-f(i,j,k-1)) &
							  +(w(i,j,k)-abs(w(i,j,k)))*(f(i,j,k+1)-f(i,j,k)))/(2.*dz(k))
					else
						w_dfdz=((w(i,j,k)+abs(w(i,j,k))) &
							  *(2.*f(i,j,k+1)+3.*f(i,j,k)-6.*f(i,j,k-1)+f(i,j,k-2)) &
							  +(w(i,j,k)-abs(w(i,j,k))) &
							  *(2.*f(i,j,k+2)+3.*f(i,j,k+1)-6.*f(i,j,k)+f(i,j,k-1)))/(2.*dz(k))
					end if
					fn(i,j,k)=f(i,j,k)-(u_dfdx+v_dfdy+w_dfdz)*dt
				end do
			end do
		end do

		do j=1,ny
			do i=0,nx
				dx1=dxi(i,j)
				dy1=dy(i,j)
				do k=1,nz
					f(i,j,k)=fn(i,j,k)
					if(i==0) then
						gx(i,j,k)=(fn(i+1,j,k)-fn(i,j,k))/dx1
					else if(i==nx) then
						gx(i,j,k)=(fn(i,j,k)-fn(i-1,j,k))/dx1
					end if
					gx(i,j,k)=(fn(i+1,j,k)-fn(i-1,j,k))/(2.*dx1)
					gy(i,j,k)=(fn(i,j+1,k)-fn(i,j-1,k))/(2.*dy1)
					gz(i,j,k)=(fn(i,j,k+1)-fn(i,j,k-1))/(2.*dz(k))
				end do
			end do
		end do

	end subroutine o3upwind_u

	!***************************************************
	subroutine o3upwind_v(f,gx,gy,gz)
	!***************************************************

		real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
		real(8),dimension(-1:im,-1:jm,0:km)::fn,u,v,w
		real(8)::hs_vp

		fn=0.; u=0.; v=0.; w=0.

		do i=1,nx
			do j=0,ny
				hs_vp=(hs(i,j)+hs(i,j+1))*.5
				do k=1,nz
					u(i,j,k)=(yu(i,j,k)+yu(i-1,j,k)+yu(i,j+1,k)+yu(i-1,j+1,k))*.25
					v(i,j,k)=yv(i,j,k)
					w(i,j,k)=(w2(i,j,k)+w2(i,j+1,k)+w2(i,j,k-1)+w2(i,j+1,k-1))*.25/hs_vp
				end do
			end do
		end do

		do i=1,nx
			do j=0,ny
				dx1=dx(i,j)
				dy1=dyj(i,j)
				do k=1,nz
					if(i==1.or.i==nx) then
						u_dfdx=((u(i,j,k)+abs(u(i,j,k)))*(f(i,j,k)-f(i-1,j,k)) &
							  +(u(i,j,k)-abs(u(i,j,k)))*(f(i+1,j,k)-f(i,j,k)))/(2.*dx1)
					else
						u_dfdx=((u(i,j,k)+abs(u(i,j,k)))  &
							  *(2.*f(i+1,j,k)+3.*f(i,j,k)-6.*f(i-1,j,k)+f(i-2,j,k)) &
							  +(u(i,j,k)-abs(u(i,j,k)))  &
							  *(2.*f(i+2,j,k)+3.*f(i+1,j,k)-6.*f(i,j,k)+f(i-1,j,k)))/(12.*dx1)
					end if
					if(j==0) then
						v_dfdy=(v(i,j,k)-abs(v(i,j,k)))*(f(i,j+1,k)-f(i,j,k))/dy1
					else if(j==ny) then
						v_dfdy=(v(i,j,k)+abs(v(i,j,k)))*(f(i,j,k)-f(i,j-1,k))/dy1
					else if(j==1.or.j==ny-1) then
						v_dfdy=((v(i,j,k)+abs(v(i,j,k)))*(f(i,j,k)-f(i,j-1,k)) &
							  +(v(i,j,k)-abs(v(i,j,k)))*(f(i,j+1,k)-f(i,j,k)))/(2.*dy1)
					else
						v_dfdy=((v(i,j,k)+abs(v(i,j,k))) &
							  *(2.*f(i,j+1,k)+3.*f(i,j,k)-6.*f(i,j-1,k)+f(i,j-2,k)) &
							  +(v(i,j,k)-abs(v(i,j,k))) &
							  *(2.*f(i,j+2,k)+3.*f(i,j+1,k)-6.*f(i,j,k)+f(i,j-1,k)))/(12.*dy1)
					end if
					if(k==1.or.k==nz) then
						w_dfdz=((w(i,j,k)+abs(w(i,j,k)))*(f(i,j,k)-f(i,j,k-1)) &
							  +(w(i,j,k)-abs(w(i,j,k)))*(f(i,j,k+1)-f(i,j,k)))/(2.*dz(k))
					else
						w_dfdz=((w(i,j,k)+abs(w(i,j,k))) &
							  *(2.*f(i,j,k+1)+3.*f(i,j,k)-6.*f(i,j,k-1)+f(i,j,k-2)) &
							  +(w(i,j,k)-abs(w(i,j,k))) &
							  *(2.*f(i,j,k+2)+3.*f(i,j,k+1)-6.*f(i,j,k)+f(i,j,k-1)))/(2.*dz(k))
					end if
					fn(i,j,k)=f(i,j,k)-(u_dfdx+v_dfdy+w_dfdz)*dt
				end do
			end do
		end do

		do j=0,ny
			do i=1,nx
				dx1=dx(i,j)
				dy1=dyj(i,j)
				do k=1,nz
					f(i,j,k)=fn(i,j,k)
					gx(i,j,k)=(fn(i+1,j,k)-fn(i-1,j,k))/(2.*dx1)
					if(j==0) then
						gy(i,j,k)=(fn(i,j+1,k)-fn(i,j,k))/dy1
					else if(j==ny) then
						gy(i,j,k)=(fn(i,j,k)-fn(i,j-1,k))/dy1
					end if
					gy(i,j,k)=(fn(i,j+1,k)-fn(i,j-1,k))/(2.*dy1)
					gz(i,j,k)=(fn(i,j,k+1)-fn(i,j,k-1))/(2.*dz(k))
				end do
			end do
		end do

	end subroutine o3upwind_v

	!***************************************************
	subroutine o3upwind_w(f,gx,gy,gz)
	!***************************************************

		real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
		real(8),dimension(-1:im,-1:jm,0:km)::fn,u,v,w

		fn=0.; u=0.; v=0.; w=0.

		do i=1,nx
			do j=1,ny
				do k=1,nz
					u(i,j,k)=(yu(i,j,k)+yu(i-1,j,k)+yu(i,j,k+1)+yu(i-1,j,k+1))*.25
					v(i,j,k)=(yv(i,j,k)+yv(i,j-1,k)+yv(i,j,k+1)+yv(i,j-1,k+1))*.25
					w(i,j,k)=w2(i,j,k)/hs(i,j)
				end do
			end do
		end do

		do i=1,nx
			do j=1,ny
				dx1=(dxi(i,j)+dxi(i-1,j))*.5
				dy1=(dyj(i,j)+dyj(i,j-1))*.5
				do k=1,nz
					dz1=(dz(k)+dz(k+1))*.5
					if(i==1.or.i==nx) then
						u_dfdx=((u(i,j,k)+abs(u(i,j,k)))*(f(i,j,k)-f(i-1,j,k)) &
							  +(u(i,j,k)-abs(u(i,j,k)))*(f(i+1,j,k)-f(i,j,k)))/(2.*dx1)
					else
						u_dfdx=((u(i,j,k)+abs(u(i,j,k)))  &
							  *(2.*f(i+1,j,k)+3.*f(i,j,k)-6.*f(i-1,j,k)+f(i-2,j,k)) &
							  +(u(i,j,k)-abs(u(i,j,k)))  &
							  *(2.*f(i+2,j,k)+3.*f(i+1,j,k)-6.*f(i,j,k)+f(i-1,j,k)))/(12.*dx1)
					end if
					if(j==1.or.j==ny) then
						v_dfdy=((v(i,j,k)+abs(v(i,j,k)))*(f(i,j,k)-f(i,j-1,k)) &
							  +(v(i,j,k)-abs(v(i,j,k)))*(f(i,j+1,k)-f(i,j,k)))/(2.*dy1)
					else
						v_dfdy=((v(i,j,k)+abs(v(i,j,k))) &
							  *(2.*f(i,j+1,k)+3.*f(i,j,k)-6.*f(i,j-1,k)+f(i,j-2,k)) &
							  +(v(i,j,k)-abs(v(i,j,k))) &
							  *(2.*f(i,j+2,k)+3.*f(i,j+1,k)-6.*f(i,j,k)+f(i,j-1,k)))/(12.*dy1)
					end if
					if(k==1.or.k==nz) then
						w_dfdz=((w(i,j,k)+abs(w(i,j,k)))*(f(i,j,k)-f(i,j,k-1)) &
							  +(w(i,j,k)-abs(w(i,j,k)))*(f(i,j,k+1)-f(i,j,k)))/(2.*dz1)
					else
						w_dfdz=((w(i,j,k)+abs(w(i,j,k))) &
							  *(2.*f(i,j,k+1)+3.*f(i,j,k)-6.*f(i,j,k-1)+f(i,j,k-2)) &
							  +(w(i,j,k)-abs(w(i,j,k))) &
							  *(2.*f(i,j,k+2)+3.*f(i,j,k+1)-6.*f(i,j,k)+f(i,j,k-1)))/(2.*dz(k))
					end if
					fn(i,j,k)=f(i,j,k)-(u_dfdx+v_dfdy+w_dfdz)*dt
				end do
			end do
		end do

		do j=1,ny
			do i=1,nx
				dx1=(dxi(i,j)+dxi(i-1,j))*.5
				dy1=(dyj(i,j)+dyj(i,j-1))*.5
				do k=1,nz
					dz1=(dz(k)+dz(k+1))*.5
					f(i,j,k)=fn(i,j,k)
					gx(i,j,k)=(fn(i+1,j,k)-fn(i-1,j,k))/(2.*dx1)
					gy(i,j,k)=(fn(i,j+1,k)-fn(i,j-1,k))/(2.*dy1)
					gz(i,j,k)=(fn(i,j,k+1)-fn(i,j,k-1))/(2.*dz1)
				end do
			end do
		end do

	end subroutine o3upwind_w

	!***************************************************
	subroutine o3upwind_c(f,gx,gy,gz)
	!***************************************************

		real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
		real(8),dimension(-1:im,-1:jm,0:km)::fn,u,v,w

		fn=0.; u=0.; v=0.; w=0.

		do i=1,nx
			do j=1,ny
				do k=1,nz
					u(i,j,k)=(yu(i,j,k)+yu(i-1,j,k))*.5
					v(i,j,k)=(yv(i,j,k)+yv(i,j-1,k))*.5
					w(i,j,k)=(w2(i,j,k)+w2(i,j,k-1))*.5/hs(i,j)
				end do
			end do
		end do

		do i=1,nx
			do j=1,ny
				dx1=(dxi(i,j)+dxi(i-1,j))*.5
				dy1=(dyj(i,j)+dyj(i,j-1))*.5
				do k=1,nz
					if(obst3d(i,j,k)==0) then
						dz1=(dzk(k)+dzk(k-1))*.5
						if(i==1.or.i==nx) then
							u_dfdx=((u(i,j,k)+abs(u(i,j,k)))*(f(i,j,k)-f(i-1,j,k)) &
								  +(u(i,j,k)-abs(u(i,j,k)))*(f(i+1,j,k)-f(i,j,k)))/(2.*dx1)
						else
							u_dfdx=((u(i,j,k)+abs(u(i,j,k)))  &
								  *(2.*f(i+1,j,k)+3.*f(i,j,k)-6.*f(i-1,j,k)+f(i-2,j,k)) &
								  +(u(i,j,k)-abs(u(i,j,k)))  &
								  *(2.*f(i+2,j,k)+3.*f(i+1,j,k)-6.*f(i,j,k)+f(i-1,j,k)))/(12.*dx1)
						end if
						if(j==1.or.j==ny) then
							v_dfdy=((v(i,j,k)+abs(v(i,j,k)))*(f(i,j,k)-f(i,j-1,k)) &
								  +(v(i,j,k)-abs(v(i,j,k)))*(f(i,j+1,k)-f(i,j,k)))/(2.*dy1)
						else
							v_dfdy=((v(i,j,k)+abs(v(i,j,k))) &
								  *(2.*f(i,j+1,k)+3.*f(i,j,k)-6.*f(i,j-1,k)+f(i,j-2,k)) &
								  +(v(i,j,k)-abs(v(i,j,k))) &
								  *(2.*f(i,j+2,k)+3.*f(i,j+1,k)-6.*f(i,j,k)+f(i,j-1,k)))/(12.*dy1)
						end if
						if(k==1.or.k==nz) then
							w_dfdz=((w(i,j,k)+abs(w(i,j,k)))*(f(i,j,k)-f(i,j,k-1)) &
								  +(w(i,j,k)-abs(w(i,j,k)))*(f(i,j,k+1)-f(i,j,k)))/(2.*dz1)
						else
							w_dfdz=((w(i,j,k)+abs(w(i,j,k))) &
								  *(2.*f(i,j,k+1)+3.*f(i,j,k)-6.*f(i,j,k-1)+f(i,j,k-2)) &
								  +(w(i,j,k)-abs(w(i,j,k))) &
								  *(2.*f(i,j,k+2)+3.*f(i,j,k+1)-6.*f(i,j,k)+f(i,j,k-1)))/(2.*dz(k))
						end if
						fn(i,j,k)=f(i,j,k)-(u_dfdx+v_dfdy+w_dfdz)*dt
					else
						fn(i,j,k)=0.
					end if
				end do
			end do
		end do

		do j=1,ny
			do i=1,nx
				dx1=(dxi(i,j)+dxi(i-1,j))*.5
				dy1=(dyj(i,j)+dyj(i,j-1))*.5
				do k=1,nz
					if(obst3d(i,j,k)==0) then
						dz1=(dz(k)+dz(k+1))*.5
						f(i,j,k)=fn(i,j,k)
						gx(i,j,k)=(fn(i+1,j,k)-fn(i-1,j,k))/(2.*dx1)
						gy(i,j,k)=(fn(i,j+1,k)-fn(i,j-1,k))/(2.*dy1)
						gz(i,j,k)=(fn(i,j,k+1)-fn(i,j,k-1))/(2.*dz1)
					end if
				end do
			end do
		end do

	end subroutine o3upwind_c

end module o3upwind_m
