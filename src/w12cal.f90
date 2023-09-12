module w12cal_m
	use common_hh
	use common_geom
	use common_hyd
	use common_grad
	implicit none
contains
	!*******************************
	subroutine w1cal(u,v,w)
	!*******************************
		real(8),dimension(-1:im,-1:jm,0:km)::u,v,w
		real(8)::hs_up,hs_vp,u_wp,v_wp,ome_wp,psi_wp
		real(8)::dhsds,dhsdn,deds,dedn

		do j=1,ny
			do i=0,nx
				hs_up=(hs(i,j)+hs(i+1,j))*.5
				u(i,j,nz+1)=u(i,j,nz)-(w(i+1,j,nz)-w(i,j,nz))*dz(nz)*hs_up/dxi(i,j)
			end do
		end do
		do i=1,nx
			do j=0,ny
				hs_vp=(hs(i,j)+hs(i,j+1))*.5
				v(i,j,nz+1)=v(i,j,nz)-(w(i,j+1,nz)-w(i,j,nz))*dz(nz)*hs_vp/dyj(i,j)
			end do
		end do

		do i=1,nx
			do j=1,ny
				dhsds=(hs(i+1,j)-hs(i-1,j))/(dxi(i,j)+dxi(i-1,j))
				dhsdn=(hs(i,j+1)-hs(i,j-1))/(dyj(i,j)+dyj(i,j-1))
				deds=(eta(i+1,j)-eta(i-1,j))/(dxi(i,j)+dxi(i-1,j))
				dedn=(eta(i,j+1)-eta(i,j-1))/(dyj(i,j)+dyj(i,j-1))
				do k=0,nz
					u_wp=(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1))*.25
					v_wp=(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1))*.25
					ome_wp=xi(k)*dhsds+deds
					psi_wp=xi(k)*dhsdn+dedn
					if(k==0) then
						w1(i,j,k)=0.
						w(i,j,k)=ome_wp*u_wp+psi_wp*v_wp
					else
						w1(i,j,k)=w(i,j,k)-ome_wp*u_wp-psi_wp*v_wp
					end if
				end do
			end do
		end do

	end subroutine w1cal

	!*******************************
	subroutine w2cal
	!*******************************

		do i=1,nx
			do j=1,ny
				do k=0,nz
					w2(i,j,k)=w1(i,j,k)-xi(k)*w1(i,j,nz)
				end do
			end do
		end do

	end subroutine w2cal

end module w12cal_m
