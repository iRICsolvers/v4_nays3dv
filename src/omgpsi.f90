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
