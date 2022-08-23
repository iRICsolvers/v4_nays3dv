      module output_m
      use common_hh
      use common_geom
      use common_hyd
      use iric
      implicit none
      contains
! **************************************
      subroutine output(fid)
! **************************************
       integer::fid, ier
       call cg_iric_write_sol_start(fid, ier)
       call CG_IRIC_WRITE_SOL_TIME(time, ier)
       call cg_iric_write_sol_baseiterative_real(fid, "Input Discharge",qp0,IER)
       call cg_iric_write_sol_baseiterative_real(fid, "Downstream Stage",h_dw,IER)
       call cg_iric_write_sol_baseiterative_real(fid, "Upstream Stage",h_up,IER)
       call cg_iric_write_sol_baseiterative_real(fid, "Upstream Depth",hs_up,IER)
       call cg_iric_write_sol_baseiterative_real(fid, "Upstream Velocity",uave_up,IER)
       call cg_iric_write_sol_baseiterative_real(fid, "Upstream Dsicharge",q(0),IER)
       call cg_iric_write_sol_baseiterative_real(fid, "Qadjust",qadjust,IER)
!      write(45,'(2f10.3)') time,h_dw
       call cg_iRIC_Write_Sol_Grid3d_Coords(fid, x38,y38,z38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "3dVelocityX",u38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "3dVelocityY",v38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "3dVelocityZ",w38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Concentration",c38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Pressure",p38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Eddy Viscosity",snu_t38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Sigma",Sigma38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Position",z38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "3dObstacle",ob38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Discharge",q38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "X-Velocity",u38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Y-Velocity",v38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Z-Velocity",w38,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "X-Vorticity",vol_x,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Y-Vorticity",vol_y,ier)
       call cg_iRIC_Write_Sol_Node_Real(fid, "Z-Vorticity",vol_z,ier)
!      call cg_iric_write_sol_cell_real(fid, "Concentration_cell",c38_cell,ier)
!      call cg_iric_write_sol_cell_real(fid, "Pressure_cell",p38_cell,ier)

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
       call cg_iric_write_sol_end(fid,ier)
       end subroutine output
      end module output_m
     

