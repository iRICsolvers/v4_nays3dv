      module diffusion_m
      use common_hh
      use common_geom
      use common_hyd
      implicit none
      contains
!*******************************
      subroutine diffusion
!*******************************
      real(8),dimension(0:im,0:jm,0:km)::uvis,vvis,wvis,cvis
      real(8)::hs_up,hs_vp,d1,d2,d3,d4,d5,dx1,dy1
      real(8)::d11,d12,d21,d22,d31,d32,d3x,d3x1,d3x2
      real(8)::snu_up,snu_vp,snu_wp,vup,uvp

      d11=0.;d12=0.;d21=0.;d22=0.;d31=0.;d32=0.
      d1=0.;d2=0.;d3=0.;d4=0.;d5=0.; d3x=0.
      snu_up=0.; snu_vp=0.;snu_wp=0. ;vup=0.

!--------------------------- U--------------------------
      do i=0,nx
       do j=1,ny
        hs_up=(hs(i,j)+hs(i+1,j))*.5
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i+1,j,k)==1) then
          yun(i,j,k)=0.
         else
          snu_up=(snu_t(i,j,k)+snu_t(i+1,j,k))*.5
          if(i==0) then
           d1=(yun(i+1,j,k)-yun(i,j,k))/dxi(i,j)**2
           d4=0.
           d5=0.
          else if(i==nx) then
           d1=(-yun(i,j,k)+yun(i-1,j,k))/dxi(i,j)**2
           d4=0.
           d5=0.
          else
           d1=(yun(i+1,j,k)-2.*yun(i,j,k)+yun(i-1,j,k))/dxi(i,j)**2
           if(k==1) then
            d4=-omega(i,j,k)/(hs_up*dz(k)*dxi(i,j))*( &
               yun(i+1,j,k+1)-yun(i-1,j,k+1)-yun(i+1,j,k)+yun(i-1,j,k))
           else
            d4=-omega(i,j,k)/(hs_up*2.*dz(k)*dxi(i,j))*( &
               yun(i+1,j,k+1)-yun(i-1,j,k+1)-yun(i+1,j,k-1)+yun(i-1,j,k-1))
           end if
           d5=-1./(hs_up*4.*dz(k)*dxi(i,j))* &
              (omega(i+1,j,k)-omega(i-1,j,k))*(yun(i,j,k+1)-yun(i,j,k-1))
          end if
          if(j==1.or.obst3d(i,j-1,k)==1) then
           d2=(yun(i,j+1,k)-yun(i,j,k))/dy(i,j)**2
          else if(j==ny.or.obst3d(i,j+1,k)==1) then
           d2=(-yun(i,j,k)+yun(i,j-1,k))/dy(i,j)**2
          else
           d2=(yun(i,j+1,k)-2.*yun(i,j,k)+yun(i,j-1,k))/dy(i,j)**2
          end if
          if(k==1) then
           if(j==1) then
            vup=(yvn(i,j,k)+yvn(i+1,j,k))*.5
           else if(j==ny) then
            vup=(yvn(i,j-1,k)+yvn(i+1,j-1,k))*.5
           else
            vup=(yvn(i,j,k)+yvn(i+1,j,k)+yvn(i,j-1,k)+yvn(i+1,j-1,k))*.25
           end if
           d3x1=snu_et_eg(i,j,k)/(dz(k)**2*hs_up**2)*(yun(i,j,k+1)-yun(i,j,k))
           d3x2=-cd2*yun(i,j,k)*sqrt(yun(i,j,k)**2+vup**2)/(dz(k)*hs_up)
           d3x=d3x1+d3x2
           uvis(i,j,k)=snu_up*(d1+d2+d4+d5)+d3x
!          if(time.gt.10000..and.j==1) then
!          write(44,'(i5,5e12.3)') i,yun(i,j,k),d3x1,d3x2,d3x,uvis(i,j,k)
!          nd if
          else
           d3=(snu_et_eg(i,j,k)*(yun(i,j,k+1)-yun(i,j,k)) &
               -snu_et_eg(i,j,k-1)*(yun(i,j,k)-yun(i,j,k-1))) &
                /(dz(k)**2*hs_up**2)
           uvis(i,j,k)=snu_up*(d1+d2+d4+d5)+d3
          end if
          yun(i,j,k)=yun(i,j,k)+uvis(i,j,k)*dt
         end if
        end do
       end do
      end do
!------------------------------V------------------------------
      do i=1,nx
       do j=0,ny
        hs_vp=(hs(i,j)+hs(i,j+1))*.5
        do k=1,nz
         if(obst3d(i,j,k)==1 .or. obst3d(i,j+1,k)==1) then
          yvn(i,j,k)=0.
         else
          snu_vp=(snu_t(i,j,k)+snu_t(i,j+1,k))*.5
          if(i==1.or.obst3d(i-1,j,k)==1) then
           d1=(yvn(i+1,j,k)-yvn(i,j,k))/dx(i,j)**2
          else if(i==nx.or.obst3d(i+1,j,k)==1) then
           d1=(-yvn(i,j,k)+yvn(i-1,j,k))/dx(i,j)**2
          else
           d1=(yvn(i+1,j,k)-2.*yvn(i,j,k)+yvn(i-1,j,k))/dx(i,j)**2
          end if
          if(j==0) then
           d2=(yvn(i,j+1,k)-yvn(i,j,k))/dyj(i,j)**2
           d4=0.
           d5=0.
          else if(j==ny) then
           d2=(-yvn(i,j,k)+yvn(i,j-1,k))/dyj(i,j)**2
           d4=0.
           d5=0.
          else
           d2=(yvn(i,j+1,k)-2.*yvn(i,j,k)+yvn(i,j-1,k))/dyj(i,j)**2
           if(k==1) then
            d4=-psi(i,j,k)/(hs_vp*dz(k)*dyj(i,j))*( &
             yvn(i,j+1,k+1)-yvn(i,j-1,k+1)-yvn(i,j+1,k)+yvn(i,j-1,k))
           else
            d4=-psi(i,j,k)/(hs_vp*2.*dz(k)*dyj(i,j))*( &
             yvn(i,j+1,k+1)-yvn(i,j-1,k+1)-yvn(i,j+1,k-1)+yvn(i,j-1,k-1))
           end if
           d5=-1./(hs_vp*4.*dz(k)*dyj(i,j))* &
           (psi(i,j+1,k)-psi(i,j-1,k))*(yvn(i,j,k+1)-yvn(i,j,k-1))
          end if
          if(k==1) then
           if(j==0) then
            uvp=(yun(i,j+1,k)+yun(i-1,j+1,k))*.5
           else if(j==ny) then
            uvp=(yun(i,j,k)+yun(i-1,j,k))*.5
           else
            uvp=(yun(i,j,k)+yun(i,j+1,k)+yun(i-1,j,k)+yun(i-1,j+1,k))*.25
           end if
           d3x1=snu_xi_eg(i,j,k)/(dz(k)**2*hs_vp**2)*(yvn(i,j,k+1)-yvn(i,j,k))
           d3x2=-cd2*yvn(i,j,k)*sqrt(yvn(i,j,k)**2+uvp**2)/(dz(k)*hs_vp)
           d3x=d3x1+d3x2
           vvis(i,j,k)=snu_vp*(d1+d2+d4+d5)+d3x
          else
           d3=(snu_xi_eg(i,j,k)*(yvn(i,j,k+1)-yvn(i,j,k)) &
               -snu_xi_eg(i,j,k-1)*(yvn(i,j,k)-yvn(i,j,k-1))) &
               /(dz(k)**2*hs_vp**2)
           vvis(i,j,k)=snu_vp*(d1+d2+d4+d5)+d3
          end if
          yvn(i,j,k)=yvn(i,j,k)+vvis(i,j,k)*dt
         end if
        end do
       end do
      end do
!------------------------------W------------------------------
      do i=1,nx
       do j=1,ny
        do k=1,nz-1
         if(obst3d(i,j,k)==1 .or. obst3d(i,j,k+1)==1) then
          ywn(i,j,k)=0.
         else
          snu_wp=(snu_t(i,j,k)+snu_t(i,j,k+1))*.5
          wvis(i,j,k)=snu_wp*( &
           (ywn(i+1,j,k)-2.*ywn(i,j,k)+ywn(i-1,j,k)) &
            *(2./(dxi(i,j)+dxi(i-1,j)))**2 &
          +(ywn(i,j+1,k)-2.*ywn(i,j,k)+ywn(i,j-1,k)) &
            *(2./(dyj(i,j)+dyj(i,j-1)))**2 &
          +(ywn(i,j,k+1)-2.*ywn(i,j,k)+ywn(i,j,k-1)) &
            /(dzk(k)**2*hs(i,j)**2))
          ywn(i,j,k)=ywn(i,j,k)+wvis(i,j,k)*dt
         end if
        end do
       end do
      end do
!------------------------------C------------------------------
      if(j_dens==1) then
      do i=1,nx
       do j=1,ny
        dx1=(dxi(i,j)+dxi(i-1,j))*.5
        dy1=(dyj(i,j)+dyj(i,j-1))*.5
        do k=1,nz
         if(obst3d(i,j,k)==1) then
          cvis(i,j,k)=0.
          ycn(i,j,k)=yc(i,j,k)
         else
          skc=snu_t(i,j,k)
          if(i==nx) then !下流端境界
           if(j_east==1) then !下流閉鎖
            d11=0.
           else  !下流開放
            if(obst3d(i,j,k)==1) then
             d11=0.
            else
             d11=(yc(i+1,j,k)-yc(i,j,k))/dxi(i,j)
            end if
           end if
          else !下流端以外
           if(obst3d(i,j,k)==1.or.obst3d(i+1,j,k)==1) then
            d11=0.
           else
            d11=(yc(i+1,j,k)-yc(i,j,k))/dxi(i,j)
           end if
          end if


          if(i==1) then !上流端境界
           if(j_west==1) then !上流閉鎖
            d12=0.
           else !上流開放
            if(obst3d(i,j,k)==1) then
             d12=0.
            else
             d12=(yc(i,j,k)-yc(i-1,j,k))/dxi(i-1,j)
            end if
           end if
          else !上流端以外
           if(obst3d(i,j,k)==1.or.obst3d(i-1,j,k)==1) then
            d12=0.
           else
            d12=(yc(i,j,k)-yc(i-1,j,k))/dxi(i-1,j)
           end if        
          end if
          d1=(d11-d12)/dx1

!         if(j==nym.and.k==2.and.i>=nx-3) &
!              write(44,'(a3,i5,5f10.5)')'d1',i,d11,d12,d1,yc(i,j,k),yc(i+1,j,k)

          if(j==ny) then
           if(j_north==1) then
            d21=0.
           else
            if(obst3d(i,j,k)==1) then
             d21=0.
            else
             d21=(yc(i,j+1,k)-yc(i,j,k))/dyj(i,j)
            end if
           end if
          else
           if(obst3d(i,j,k)==1.or.obst3d(i,j+1,k)==1) then
            d21=0.
           else
            d21=(yc(i,j+1,k)-yc(i,j,k))/dyj(i,j)
           end if
          end if

          if(j==1) then
           if(j_south==1) then
            d22=0.
           else
            if(obst3d(i,j,k)==1) then
             d22=0.
            else
             d22=(yc(i,j,k)-yc(i,j-1,k))/dyj(i,j)
            end if
           end if
          else
           if(obst3d(i,j,k)==1.or.obst3d(i,j-1,k)==1) then
            d22=0.
           else
            d22=(yc(i,j,k)-yc(i,j-1,k))/dyj(i,j)
           end if
          end if
          d2=(d21-d22)/dy1

          if(k==nz) then
           d31=0.
          else
           d31=(yc(i,j,k+1)-yc(i,j,k))/(dzk(k)*hs(i,j))
          end if
          if(k==1) then
           d32=0.
          else
           d32=(yc(i,j,k)-yc(i,j,k-1))/(dzk(k-1)*hs(i,j))
          end if
          d3=(d31-d32)/(dz(k)*hs(i,j))

          cvis(i,j,k)=skc*(d1+d2+d3)
          ycn(i,j,k)=yc(i,j,k)+cvis(i,j,k)*dt
         end if
        end do
       end do
      end do

!     k=2;j=nym
!     write(44,'(a3,5e12.3)') 'yc ',(yc(i,j,k),i=nx-3,nx+1)
!     write(44,'(a3,5e12.3)') 'ycn',(ycn(i,j,k),i=nx-3,nx+1)
!
      end if

      end subroutine diffusion
      end module diffusion_m
