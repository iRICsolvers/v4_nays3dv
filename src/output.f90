      module output_m
      use common_hh
      use common_geom
      use common_hyd
      implicit none
      contains
! **************************************
      subroutine output
! **************************************
       integer::ier
       call CG_IRIC_WRITE_SOL_TIME_F(time, ier)
       call cg_iric_write_sol_baseiterative_real_f("Input Discharge",qp0,IER)
       call cg_iric_write_sol_baseiterative_real_f("Downstream Stage",h_dw,IER)
       call cg_iric_write_sol_baseiterative_real_f("Upstream Stage",h_up,IER)
       call cg_iric_write_sol_baseiterative_real_f("Upstream Depth",hs_up,IER)
       call cg_iric_write_sol_baseiterative_real_f("Upstream Velocity",uave_up,IER)
       call cg_iric_write_sol_baseiterative_real_f("Upstream Dsicharge",q(0),IER)
       call cg_iric_write_sol_baseiterative_real_f("Qadjust",qadjust,IER)
!      write(45,'(2f10.3)') time,h_dw
       call CG_IRIC_WRITE_SOL_GRIDCOORD3D_F(x38,y38,z38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("3dVelocityX",u38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("3dVelocityY",v38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("3dVelocityZ",w38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Concentration",c38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Pressure",p38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Eddy Viscosity",snu_t38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Sigma",Sigma38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Position",z38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("3dObstacle",ob38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Discharge",q38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("X-Velocity",u38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Y-Velocity",v38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Z-Velocity",w38,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("X-Vorticity",vol_x,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Y-Vorticity",vol_y,ier)
       call CG_IRIC_WRITE_SOL_REAL_F("Z-Vorticity",vol_z,ier)
!      call cg_iric_write_sol_cell_real_f("Concentration_cell",c38_cell,ier)
!      call cg_iric_write_sol_cell_real_f("Pressure_cell",p38_cell,ier)

!      if(time>00.) then
!      write(44,*) 'time=',time
!      write(44,*) 'c(nx+1)'
!      do k=nz+1,0,-1
!       write(44,'(20f8.4)')(ycn(nx+1,j,k),j=0,ny+1)
!      end do
!      write(44,*) 'c(nx)'
!      do k=nz+1,0,-1
!       write(44,'(20f8.4)')(ycn(nx,j,k),j=0,ny+1)
!      end do
!      write(44,*) 'u(nx+1)'
!      do k=nz+1,0,-1
!       write(44,'(20f8.4)')(yun(nx+1,j,k),j=0,ny+1)
!      end do
!      write(44,*) 'u(nx)'
!      do k=nz+1,0,-1
!       write(44,'(20f8.4)')(yun(nx,j,k),j=0,ny+1)
!      end do
!      write(44,*) 'p(nx)'
!      do k=nz+1,0,-1
!       write(44,'(20f8.4)')(ypn(nx+1,j,k),j=0,ny+1)
!      end do
!      write(44,*) 'p(nx+1)'
!      do k=nz+1,0,-1
!       write(44,'(20f8.4)')(ypn(nx,j,k),j=0,ny+1)
!      end do
!      write(44,*) 'h(nx+1)'
!      write(44,'(20f8.5)')(h(nx+1,j),j=0,ny+1)
!      write(44,*) 'h(nx)'
!      write(44,'(20f8.5)')(h(nx,j),j=0,ny+1)
!      end if

       end subroutine output
      end module output_m
     

