!*************************************************
!
! Main 
!
!*************************************************
      Program nays3dv
      use common_hh
      use common_geom
      use common_hyd
      use common_grad
      use alloc_var_m
      use initial_m
      use output_m
      use cell2node_m
      use w12cal_m
      use non_advection_m
      use bound_m
      use diffusion_m
      use newgrd_m
      use cip3d_m
      use upwind3d_m
      use o3upwind_m
      use iric

      implicit none
      
      real(8)::sorerr,soralpha
      real(8)::total_1,total_2,diff_ratio,ss,elh
      integer::lsor,lmax,n

      integer, parameter :: strMax = 1000
      character(len = strMax) :: condFile

!------ Values for iRIC ------
      integer::ircount,ier,fid,istatus
!------------------------------------------
!     open(44,file='c:\tmp\tmp3.d',status='unknown')
!     open(45,file='c:\tmp\stage.d',status='unknown')
!     open(47,file='c:\tmp\hd_prof.d',status='unknown')
!------------------------------------------

      ircount=0; ier=0; fid=0; ni=0; nj=0; nk=0
      tuk=0.;etime=0.;dt=0.; st_dens=0.
      sorerr=0.; lsor=0; soralpha=0.;lmax=0
      surf_tension=0.
      hloop=0;hloop_err=0.
!=========================================
      g=9.8
      kappa=0.4

!==========================================================
! 格子生成データファイルを開く
!     call cg_open("./nays3dv/iRICZone/gridcreate.cgn", CG_MODE_READ, fin, ier)
!     if (ier /=0) stop "*** Open error of CGNS file ***"



!--- Read from grid generator file -----
!  格子生成条件データファイルからの読み出し
!     call cg_iric_read_real(fin,'xlen',xlen,ier)
!     call cg_iric_read_real(fin,'ylen',ylen,ier)
!     call cg_iric_read_integer(fin,'nx',nx, ier)
!     call cg_iric_read_integer(fin,'ny',ny, ier)    機能してなさそうだったのでコメントアウト（22.08.23星野)
!     write(44,*) 'nx, xlen=',nx,xlen
!     write(44,*) 'ny, ylen=',ny,ylen
!     nxm=nx/2; nym=ny/2
!     call cg_iric_close(fin, ier)
!=========================================================

      ircount = nargs()
      if(ircount == 2) then
       call getarg(1,condFile,ier)
      else
       write(*,*) "Input File not specified."
       stop
      end if

! 計算データファイルを開く
      call cg_iric_open(condFile,IRIC_MODE_MODIFY, fid, ier)
      if (ier /=0) then
       write(*,*) "*** Open error of CGNS file ***"
      end if
 
!
!guiにcgnsファイルを読込みであることを知らせるファイルを生成
      call iric_initoption(IRIC_OPTION_CANCEL, ier)

!
! 2次元格子データーのサイズを読み込む ni,nj
!
      call cg_iRIC_Read_Grid2d_Str_Size(fid, ni, nj, ier)
!     write(44,*) 'ni,nj=',ni,nj
      nx=ni-1; ny=nj-1
      nxm=nx/2; nym=ny/2

      if (ier /=0) then
       write(*,*) "*** error:cg_iRIC_Read_Grid2d_Str_Size ***"
       stop
      end if
!     
      allocate(x8(ni,nj),y8(ni,nj),z8(ni,nj),h8(ni,nj))
      allocate(hobst38(ni-1,nj-1),obst4(ni-1,nj-1))
      allocate(hobst_node8(ni,nj))
      x8=0.;y8=0.;z8=0.;h8=0.;obst4=0; hobst38=0.; hobst_node8=0.

  

! 2次元格子データーを読み込む x8(ni,nj),y8(ni,nj)
      call cg_iRIC_Read_Grid2d_Coords(fid, x8,y8,ier)
!
      if (ier /=0) then
       write(*,*) "*** error:cg_iRIC_Read_Grid2d_Coords ***"
       stop
      end if

!read bed elevation data
! 河床データの読み込み
!     call cg_iric_read_grid_real_node("BedElevation_n", z8, ier)
      call cg_iric_read_grid_real_node(fid, "Elevation", z8, ier)
      if (ier /=0) then
         write(*,*) "error:cg_iric_read_grid_real_node-BedElevation"
         read(*,*)
         call calc_stop(fid)
      end if

!初期水位データの読み込み

      call cg_iric_read_integer(fid, 'j_hinit',j_hinit,ier)  !初期水位条件
!     write(*,*) 'j_hinit=',j_hinit
!   j_hinit==1: 水平  j_hinit==2:地形ファイルからの読み込み
      if(j_hinit==1) &
       call cg_iric_read_real(fid, 'h_horizontal',h_horizontal,ier) !水平水位

      if(j_hinit==2) then
!      call cg_iric_read_grid_real_node(fid, "WaterSurface_n", h8, ier)
       call cg_iric_read_grid_real_node(fid, "WaterSurface", h8, ier)
       if (ier /=0) then
          write(*,*) "error:cg_iric_read_grid_real_node_f-WaterSurface"
          read(*,*)
          call calc_stop(fid)
       end if
      else if(j_hinit==1) then
       h8=h_horizontal
      end if

!最小水深の指定

      call cg_iric_read_integer(fid, 'j_mindep',j_mindep,ier)  !最少水深条件
      if(j_mindep==1) then
       call cg_iric_read_real(fid, 'dep_min',dep_min,ier) !最少水深
      end if

!     write(44,*) 'j_hinit=',j_hinit
!     do i=1,ni
!      do j=1,nj
!       write(44,*) i,j,h8(i,j)
!      end do
!     end do
!     stop

!障害物データの読み込み
      call cg_iric_read_grid_integer_cell(fid, "Obstacle", obst4, ier)
      if (ier /=0) then
         write(*,*) "error:cg_iric_read_grid_integer_cell_f-Obstacle"
         read(*,*)
         call calc_stop(fid)
      end if

      call cg_iric_read_grid_real_cell(fid, "ObstacleTop", hobst38, ier)
      if (ier /=0) then
         write(*,*) "error:cg_iric_read_grid_real_cell_f-ObstacleTop"
         read(*,*)
         call calc_stop(fid)
      end if
!
      call cg_iric_read_integer(fid, 'j_uadvec',j_uadvec, ier)
      call cg_iric_read_integer(fid, 'j_cadvec',j_cadvec, ier)
      call cg_iric_read_integer(fid, 'nz',nz, ier)
      call cg_iric_read_integer(fid, 'j_zgrid',j_zgrid, ier)
      if(j_zgrid==2) call cg_iric_read_real(fid, 'dz1',dz10, ier)
      call cg_iric_read_real(fid, 'diam',diam, ier)
      ks=2.*diam; z0=ks/30.; snm=0.0486*diam**(1./6.)

      call cg_iric_read_integer(fid, 'j_snu',j_snu, ier)
      call cg_iric_read_real(fid, 'al_ep',al_ep, ier)
!     call cg_iric_read_real(fid, 'qp',qp, ier)
!
! 平面上の境界条件
! j_west(i=0), j_east(i=nx), j_south(i=0), j_north(j=ny)
! 1==Closeed  2==Open  3=Periodic
      call cg_iric_read_integer(fid, 'j_west',j_west,ier)
!     call cg_iric_read_integer(fid, 'jp_west',jp_west,ier)
      if(j_west==3) then
       j_east=3
      else
       call cg_iric_read_integer(fid, 'j_east',j_east,ier)
!      call cg_iric_read_integer(fid, 'jp_east',jp_east,ier)
      end if

      call cg_iric_read_integer(fid, 'j_south',j_south,ier)
      if(j_south==3) then
       j_north=3
      else
       call cg_iric_read_integer(fid, 'j_north',j_north,ier)
      end if

!     write(*,'(4i5)') j_west,j_east,j_south,j_north
!     write(44,'(4i5)') j_west,j_east,j_south,j_north
!
! Hydraulic Boundary Condition
! 上流流量
      if(j_west<=2) then
       call cg_iric_read_integer(fid, 'j_qin',j_qin,ier)  !流量条件
       if(j_qin==0) then
        q_up_const=0.
       else if(j_qin==1) then
        call cg_iric_read_real(fid, 'q_up_const',q_up_const, ier) !一定流量
       else if(j_qin==2) then
        call cg_iric_read_functionalsize(fid, "q_hyd",n_qsize,ier) !流量の個数
        allocate(time_q(n_qsize),qt_up(n_qsize))
        call cg_iric_read_functionalwithname(fid, "q_hyd","time_q",time_q,ier) !時間
        call cg_iric_read_functionalwithname(fid, "q_hyd","qt_up",qt_up,ier) !流量時系列
        etime_q=time_q(n_qsize)
       end if
       if(j_qin>=1) then
!       call cg_iric_read_real(fid, 'dh_ref',dh_ref, ier)
!       call cg_iric_read_real(fid, 'dh_alpha',dh_alpha, ier)
!流量補正に用いる水位変動単位に対する倍率
        call cg_iric_read_real(fid, 'q_stt',q_stt, ier) !流量補正の開始時間 
        call cg_iric_read_real(fid, 'q_trn',q_trn, ier) !流量補正の開始時間 
       end if
      end if
!上流水位
      if(j_west==2) then
       call cg_iric_read_integer(fid, 'j_hup',j_hup,ier)  !水位条件
! j_hup==1: 一定  j_hup==2: 水平  j_hup==3: 勾配を与える  j_hup=4:等流計算
       if(j_hup==1) then
        call cg_iric_read_real(fid, 'h_up_const',h_up_const, ier) !一定水位
       end if
       if(j_hup==3) then
        call cg_iric_read_real(fid, 'up_wslope',up_wslope, ier) !上流端の水面勾配
       end if
       if(j_hup==4) then
        call cg_iric_read_real(fid, 'up_slope',up_slope, ier)
!上流の水深流速を等流計算するための勾配
       end if
      end if
!下流水位
      if(j_east==2) then
       call cg_iric_read_integer(fid, 'j_hdw',j_hdw,ier)  !水位条件
       if(j_hdw==1) then
        call cg_iric_read_real(fid, 'h_dw_const',h_dw_const, ier) !一定水位
       else if(j_hdw==3) then !sine Oscillation
        call cg_iric_read_real(fid, 'hd_amp',hd_amp,ier) !波高
        call cg_iric_read_real(fid, 'hd_wl',hd_wl,ier) !周期
        call cg_iric_read_real(fid, 'hd_st',hd_st,ier) !開始時刻
        call cg_iric_read_real(fid, 'hd_ap',hd_ap,ier) !発達時間
       else if(j_hdw==4) then ! 下流端水位を時系列で与える
        call cg_iric_read_functionalsize(fid, "h_hyd",n_hsize,ier) !流量の個数
        allocate(time_h(n_qsize),ht_up(n_hsize))
        call cg_iric_read_functionalwithname(fid, "h_hyd","time_h",time_h,ier) !時間
        call cg_iric_read_functionalwithname(fid, "h_hyd","ht_up",ht_up,ier) !流量時系列
        etime_h=time_h(n_hsize)
       end if
      end if


      if(j_west==1) then !上流閉鎖
       qp=0.
      else !上流自由またはPeriodicの場合
       if(j_qin==0) then !フリー
        qp=0.
       else if(j_qin==1) then !一定流量
        qp=q_up_const
       else if(j_qin==2) then !流量時系列
        qp=qt_up(1)
       end if
      end if
!   
! アウトプット用3次元変数のアロケート
!
      nk=nz+1
      allocate(x38(ni,nj,nk),y38(ni,nj,nk),z38(ni,nj,nk))
      allocate(u38(ni,nj,nk),v38(ni,nj,nk),w38(ni,nj,nk))
      allocate(p38(ni,nj,nk),c38(ni,nj,nk),sigma38(ni,nj,nk))
      allocate(snu_t38(ni,nj,nk))
      allocate(ob38(ni,nj,nk),q38(ni,nj,nk))
      allocate(vol_x(ni,nj,nk),vol_y(ni,nj,nk),vol_z(ni,nj,nk))
!     allocate(c38_cell(ni-1,nj-1,nk-1),p38_cell(ni-1,nj-1,nk-1))

      x38=0.; y38=0.; z38=0.; u38=0.; v38=0.; w38=0.; p38=0.; c38=0.
      vol_x=0.; vol_y=0.; vol_z=0.
      sigma38=0.; ob38=0. ;q38=0.; snu_t38=0.
!
! Read computational parameters related in time 
!
      call cg_iric_read_real(fid, 'tuk',tuk,ier)
      call cg_iric_read_real(fid, 'etime',etime,ier)
      if(j_west<=2 .and. j_qin==2) then
       etime=min(etime0,etime_q)
       etime0=etime
      else if(j_east==2 .and. j_hdw==4) then
       etime=min(etime0,etime_h)
      end if
      call cg_iric_read_real(fid, 'dt',dt,ier)
      call cg_iric_read_real(fid, 'st_dens',st_dens,ier)
!     write(44,'(a13,3f10.5)') 'tuk,etime,dt=',tuk,etime,dt
!
!Read parameters related with iteration
      call cg_iric_read_real(fid, 'sorerr',sorerr,ier)
      call cg_iric_read_integer(fid, 'lsor',lsor,ier)
! lsor=sorの反復回数(5回程度) ----------
      call cg_iric_read_real(fid, 'soralpha',soralpha,ier)
! soralpha = sor法の緩和回数(=1.5で加速緩和) ------
!     write(44,'(a21,f10.5,i5,f10.5)') 'sorerr,lsor,soralpha=',sorerr,lsor,soralpha

      call cg_iric_read_integer(fid, 'j_surf',j_surf,ier)
! 自由水面計算 j_surf=1 (Yes)  j_suef=0 (no)
! j_surf=0--->固定水面
      call cg_iric_read_real(fid, 'alpha_surf',alpha_surf,ier)
      call cg_iric_read_real(fid, 'stime_surf',stime_surf,ier)

!     write(44,*) 'j_surf,al_suef,st_surf=',j_surf,alpha_surf,stime_surf

      call cg_iric_read_integer(fid, 'hloop',hloop,ier)
! 自由水面計算 繰り返し回数
      if(j_surf==0) hloop=1
      call cg_iric_read_real(fid, 'hloop_err',hloop_err,ier)
!     write(44,*) 'hloop,hloop_err=',hloop,hloop_err

! 自由水面計算 繰り返し計算打ち切り誤差
      call cg_iric_read_real(fid, 'snu',snu,ier)
!     call cg_iric_read_real(fid, 'skt',skt,ier)
      call cg_iric_read_real(fid, 'skc',skc,ier)
!     call cg_iric_read_real(fid, 'beta_t',beta_t,ier)
      call cg_iric_read_real(fid, 'rho',rho,ier)
      call cg_iric_read_real(fid, 'surf_tension',surf_tension,ier)
      smg_g=-surf_tension/(rho*g)
!初期の濃度条件の読み込み
      call cg_iric_read_integer(fid, 'j_dens',j_dens0,ier)
      j_dens=j_dens0
      if(j_dens0==0) then
       j_ini_dens=0
       j_bc_dens=0
      else
       call cg_iric_read_integer(fid, 'j_ini_dens',j_ini_dens,ier)
       call cg_iric_read_integer(fid, 'j_bc_dens',j_bc_dens,ier)
      end if
      if(j_dens0==1) then
       call cg_iric_read_real(fid, 'c0',c0,ier)
      else
       c0=0.
      end if
      if(j_dens0==0.or.j_ini_dens==0) then
       ic1=1;ic2=1;jc1=1;jc2=1;kc1=1;kc2=1
       c1=0.
      else
       call cg_iric_read_real(fid, 'c1',c1,ier)
       call cg_iric_read_integer(fid, 'ic1',ic1,ier)
       call cg_iric_read_integer(fid, 'ic2',ic2,ier)
       call cg_iric_read_integer(fid, 'jc1',jc1,ier)
       call cg_iric_read_integer(fid, 'jc2',jc2,ier)
       call cg_iric_read_integer(fid, 'kc1',kc1,ier)
       call cg_iric_read_integer(fid, 'kc2',kc2,ier)
       if(ic1<1 .or. ic1>nx.or.ic2<1.or.ic2>nx.or.ic1>ic2.or. & 
        jc1<1 .or. jc1>ny.or.jc2<1.or.jc2>ny.or.jc1>jc2.or. & 
        kc1<1 .or. kc1>nz.or.kc2<1.or.kc2>nz.or.kc1>kc2) then
        write(*,'(a9,3i5)') 'nx,ny,nz=',nx,ny,nz
        write(*,*) 'ic1,ic2,jc1,jc2,kc1,kc2='
        write(*,'(6i5)') ic1,ic2,jc1,jc2,kc1,kc2
        write(*,*) 'initial condition domain error'
        stop
       end if
      end if
!     write(44,'(2f10.5)') c0,c1
!     write(44,'(6i5)')ic1,ic2,jc1,jc2,kc1,kc2

!***********************************************************
!
!  変数のアロケーション
!***********************************************************

      im=nx+3; jm=ny+3; km=nz+3
      call alloc_var(im,jm,km) 

      call initial

      call bound_u(yu,yw)
      yun=yu

!     do i=1,nx
!      write(44,'(20i1)') (obst(i,j),j=1,ny)
!     end do
!     write(44,*)
!     do i=0,nx
!      write(44,'(20f8.3)')(z(i,j),j=0,ny)
!     end do
!     write(44,*)
!     do i=1,nx
!      write(44,'(20f8.3)')(eta(i,j),j=1,ny)
!     end do

!     write(44,*) ni,nj,nk
!     do i=1,ni
!      do j=1,nj
!       do k=1,nk
!        write(44,'(3i5,3f10.5)') i,j,k,x38(i,j,k),y38(i,j,k),z38(i,j,k)
!       end do
!      end do
!     end do

!     call cg_iRIC_Write_Grid3d_Coords(ni,nj,nk,x38,y38,z38,ier)
!------------------------------------------
! 濃度 j_in_c 境界条件の個数
      call cg_iric_read_bc_count(fid, 'b_con',jc_in)
!     write(44,*) 'jc_in=',jc_in
! -----------濃度境界条件の読み込み -----

      allocate(jc_inlen(jc_in))  ! ede number*2
      allocate(boundary_con_value(jc_in),boundary_con_up(jc_in))  !Boundary Concentration Value 

      indexmax_c = 0 ! max of jc_inlen(i)

      do i = 1, jc_in
       call cg_iric_read_bc_indicessize(fid, 'b_con',i,jc_inlen(i),ier)
!      write(44,*) 'i,jc_inlen(i)=',i,jc_inlen(i)
       if(indexmax_c.lt. jc_inlen(i)) indexmax_c= jc_inlen(i)
!      call cg_iric_read_bc_string('b_con',i,'_caption',flowname_c(i),ier)
!      write(44,*) flowname_c(i)(1:30)
      end do

!     write(44,*) 'indexmax_c=',indexmax_c
      allocate(indices_c(jc_in, 2, indexmax_c))

      do i = 1, jc_in
       call cg_iric_read_bc_indices(fid, 'b_con',i, indices_c(i:i,:,:), ier)
       call cg_iric_read_bc_real(fid, 'b_con', i, 'b_con_val', &
       boundary_con_value(i), ier)
!      write(44,*) 'concentration=',boundary_con_value(i)

       call cg_iric_read_bc_real(fid, 'b_con', i, 'e_up_con_val', &
       boundary_con_up(i), ier)
!      write(44,*) 'upper limit=',boundary_con_up(i)
      end do
!
!-------境界条件の場所を探す-----
! Assign Boudary Condition
!

      do i=1,jc_in
       do j=1,jc_inlen(i)-1
! Check i=1 Boundary (West Boundary)
        if(indices_c(i,1,1)==1 .and. indices_c(i,1,jc_inlen(i))==1) then
         if(j==1) then
          jc_in_west=jc_in_west+1
          c_bound_west(jc_in_west)=boundary_con_value(i)
          c_up_west(jc_in_west)=boundary_con_up(i)
!         write(44,*) '[con]Bound West',jc_in_west,c_bound_west(jc_in_west)
          jc_west_s(jc_in_west)=indices_c(i,2,1)
          jc_west_e(jc_in_west)=indices_c(i,2,jc_inlen(i)-1)
!         write(44,*) '[con]Bound West start end=',jc_west_s(jc_in_west),jc_west_e(jc_in_west)
         end if
! Check i=ni Boundary (East Boundary)
        else if(indices_c(i,1,1)==ni .and. indices_c(i,1,jc_inlen(i))==ni) then
         if(j==1) then
          jc_in_east=jc_in_east+1
          c_bound_east(jc_in_east)=boundary_con_value(i)
          c_up_east(jc_in_east)=boundary_con_up(i)
!         write(44,*) '[con]Bound East',jc_in_east,c_bound_east(jc_in_east)
          jc_east_s(jc_in_east)=indices_c(i,2,1)
          jc_east_e(jc_in_east)=indices_c(i,2,jc_inlen(i)-1)
!         write(44,*) '[con]Bound East start end=',jc_east_s(jc_in_east),jc_east_e(jc_in_east)
         end if
! Check j=1 Boundary (South Boundary)
        else if(indices_c(i,2,j)==1) then
         if(j==1) then
          jc_in_south=jc_in_south+1
          c_bound_south(jc_in_south)=boundary_con_value(i)
          c_up_south(jc_in_south)=boundary_con_up(i)
!         write(44,*) '[con]Bound South',jc_in_south,c_bound_south(jc_in_south)
          ic_south_s(jc_in_south)=indices_c(i,1,1)
          ic_south_e(jc_in_south)=indices_c(i,1,jc_inlen(i)-1)
!         write(44,*) '[con]Bound South start end=',ic_south_s(jc_in_south),ic_south_e(jc_in_south)
         end if
! Check j=nj Boundary (NorthBoundary)
        else if(indices_c(i,2,j)==nj) then
         if(j==1) then
          jc_in_north=jc_in_north+1
          c_bound_north(jc_in_north)=boundary_con_value(i)
          c_up_north(jc_in_north)=boundary_con_up(i)
!         write(44,*) '[con]Bound North',jc_in_north,c_bound_north(jc_in_north)
          ic_north_s(jc_in_north)=indices_c(i,1,1)
          ic_north_e(jc_in_north)=indices_c(i,1,jc_inlen(i)-1)
!         write(44,*) '[con]Bound North start end=',ic_north_s(jc_in_north),ic_north_e(jc_in_north)
         end if
        end if
       end do
      end do
!============================================================
!  初期状態での境界濃度セルの指定
!============================================================
      if(j_dens0==1.and.j_bc_dens==1) then
!      write(47,*) jc_in_west,jc_in_east,jc_in_south,jc_in_north

      if(jc_in_west>0) then
       do m=1,jc_in_west
        do j=jc_west_s(m),jc_west_e(m)
         do k=1,nz
          elh=eta(1,j)+hs(1,j)*xxi(k)
          if(elh<=c_up_west(m)) k_bc_west(m)=k
         end do
        end do
       end do
      end if

      if(jc_in_east>0) then
       do m=1,jc_in_east
        do j=jc_east_s(m),jc_east_e(m)
         do k=1,nz
          elh=eta(nx,j)+hs(nx,j)*xxi(k)
          if(elh<=c_up_east(m)) k_bc_east(m)=k
         end do
        end do
       end do
      end if  
!
      if(jc_in_south>0) then
       do m=1,jc_in_south
        do i=ic_south_s(m),ic_south_e(m)
         do k=1,nz
          elh=eta(i,1)+hs(i,1)*xxi(k)
          if(elh<=c_up_south(m)) k_bc_south(m)=k
         end do
        end do
       end do
      end if  

       if(jc_in_north>0) then
       do m=1,jc_in_north
        do i=ic_north_s(m),ic_north_e(m)
         do k=1,nz
          elh=eta(i,ny)+hs(i,ny)*xxi(k)
          if(elh<=c_up_north(m)) k_bc_north(m)=k
         end do
        end do
       end do
      end if  

      end if

!=============================================================
      call cg_iRIC_Write_Grid3d_Coords(fid, ni,nj,nk,x38,y38,z38,ier)
!=============================================================
      itout=int(tuk/dt)
      time=0.
      itt=itout
      icount=0
      pi=4.*atan(1.0d0)
      m=0

      h_dw_old=h(nx,nym)
      h_dw=h_dw_old

!
! Starting Point of Main Roop 
!
      do while(time<etime) 
      icount=icount+1

      if(j_dens0==1 .and. time>st_dens) then
       j_dens=1
      else
       j_dens=0
      end if
!
! ---- Nan Check ------
!
      do i=1,nx
       do j=1,ny
        do k=1,nz,nz-1
         if(isnan(yu(i,j,k)).or.isnan(ypn(i,j,k))) then 
          write(*,*) 'Calculation is Falure (Nan found)',time
          call cg_iric_close(fid,ier); stop
         end if
        end do
        if(isnan(hs(i,j))) then
         write(*,*) 'Calculation is Falure (Nan found)',time
         call cg_iric_close(fid,ier); stop
        end if
       end do
      end do

      call iric_check_cancel(istatus)
      if(istatus==1) then
       write(*,*) &
        "Solver is stopped because the STOP button was clicked."
         call cg_iric_close(fid,ier); stop
      end if
!
! utsaの計算
!
      call ustacal
!
! 流量の計算

      call qcal
!
! Set input discharge qp
! 上流から与える流量の計算
!

!     write(*,*) 'j_west,j_qin,qp=',j_west,j_qin,qp
!
! j_qin==0:Free,  j_qin==1:Constant,  j_qin==2: Hydrograph
! j_hup==1:Constant,  j_hup==2:Horizontal,  j_hup==3:Give a Slope,
! j_hup=4:Uniform Flow
! 

      if(j_west==2) then !上流自由のとき
!
! 上流端の有効幅
!
       width_in=0.
       do j=1,ny
        if(obst(i,1)==0) then
         width_in=width_in+dy(0,j)
        end if
       end do
!
! 上流端に与える流量
!
       if(j_qin==2) then !流量時系列を与える場合
        do n=2,n_qsize
         if(time>=time_q(n-1).and.time<=time_q(n)) then
          ss=(time-time_q(n-1))/(time_q(n)-time_q(n-1))
          qp=qt_up(n-1)+ss*(qt_up(n)-qt_up(n-1))
         end if
        end do
       else if(j_qin==1) then
        if(time<=q_stt) then
         qp0=0.
        else if(time<=(q_stt+q_trn)) then
         qp0=qp*0.5*(1-cos(pi/q_trn*(time-q_stt)))
        else if(time>(q_stt+q_trn)) then
         qp0=qp
        else
         qp0=0.
        end if
       else
        qp0=0.
       end if
! 
! 上流端の水位
!
       if(j_hup==1) then
        h_up=h_up_const
        hs_up=h_up-eave_up
       else if(j_hup==2 .or. j_hup==3) then
        h_up=0.
        jj=0
        do j=1,ny
         if(obst(i,1)==0) then
          jj=jj+1
          h_up=h_up+h(1,j)
         end if
        end do
        if(jj>0) then
         h_up=h_up/float(jj) !上流平均水位
        else
         h_up=have_dw_ini
        end if
        hs_up=h_up-eave_up
       else if(j_hup==4) then
        hs_up=(snm*qp0/(width_in*sqrt(up_slope)))**0.6 !等流水深
        h_up=eave_up+hs_up
        if(h_up.lt.have_dw_ini) h_up=have_dw_ini
        hs_up=h_up-eave_up
       end if
! 
! 上流端の流速
!
       if(j_qin>=1) then
        uave_up=qp0/(width_in*hs_up)
        sb=log(hs_up/z0)-1.
        do k=1,nz
         zz=xxi(k)*hs_up
         u00(k)=uave_up*log(zz/z0)/sb
        end do
        u00(0)=-u00(1)
        u00(nz+1)=u00(nz)
       end if
      end if 
! 下流端の水位
!
      if(j_east==2) then
!      write(47,*) time
       if(j_hdw==1) then
        h_dw=h_dw_const
       else if(j_hdw==3) then
        if(time<hd_st) then
         h_dw=have_dw_ini
         hd_amp0=0.
        else
         if(time<(hd_st+hd_ap)) then
          amp_alpha=0.5*(1.-cos(pi/hd_ap*(time-hd_st)))
         else
          amp_alpha=1.0
         end if
         hd_amp0=hd_amp*sin(2.*pi*(time-hd_st)/hd_wl)
         h_dw=have_dw_ini+hd_amp0*amp_alpha
        end if
!       write(47,'(5f12.3)') time,have_dw_ini,amp_alpha,hd_amp0,h_dw
       else if(j_hdw==4) then
        write(*,*) 'not ready yet'
        stop
       end if
      end if

      if(j_hdw==1) then
       dhdt_dw=0.
      else
       dhdt_dw=(h_dw-h_dw_old)/dt
       h_dw_old=h_dw
      end if
!     if(time.gt.1000) then
!      write(47,'(21f10.5,e12.3)') time,h_dw,dhdt_dw
!     end if


!     write(44,*)'time=',time
!----------------------------
!     i=10;  k=nz
!     write(44,'(a2,4f8.4)') 'u ', yu(0,j,nz),yu(1,j,nz),yu(nx,j,nz),yu(nx+1,j,nz)
!     write(44,'(a2,4f8.0)') 'p ', ypn(0,j,nz),ypn(1,j,nz),ypn(nx,j,nz),ypn(nx+1,j,nz)
!     write(44,'(a2,4f8.0)') 'p ', ypn(0,j,nz-1),ypn(1,j,nz-1),ypn(nx,j,nz-1),ypn(nx+1,j,nz-1)
!     write(44,'(a2,4f8.0)') 'p ', ypn(0,j,nz-2),ypn(1,j,nz-2),ypn(nx,j,nz-2),ypn(nx+1,j,nz-2)
!     write(44,'(a2,4f8.4)') 'h ', h(0,j),h(1,j),h(nx,j),h(nx+1,j)
!     write(44,'(a2,4f12.3)') 'p ', ypn(i,0,nz),ypn(i,1,nz),ypn(i,ny,nz),ypn(i,ny+1,nz)
!     write(44,'(a2,4f12.3)') 'p ', ypn(i,0,nz-1),ypn(i,1,nz-1),ypn(i,ny,nz-1),ypn(i,ny+1,nz-1)
!     write(44,'(a2,4f12.3)') 'p ', ypn(i,0,nz-2),ypn(i,1,nz-2),ypn(i,ny,nz-2),ypn(i,ny+1,nz-2)
!     write(44,'(a2,4f8.4)') 'h ', h(i,0),h(i,1),h(i,ny),h(i,ny+1)

      if(itt==itout) then
       itt=0
!      write(*,'(a5,f8.4,f8.5,2i5)') 'time=',time,qp,lmax,m
!      write(44,'(a5,f8.4,f8.5,2i5)') 'time=',time,qp,lmax,m
       write(*,'(a5,f10.3,2f13.3,f10.3)') 'time=',time,qp,q(0),h_dw
!      write(44,'(a5,f10.3,3f8.5)') 'time=',time,qp,q(0),h_dw
!      write(45,'(7f10.4)') time,h_dw,hs(0,nym),hs(5,nym),q(5),yu(5,nym,nz),u0_in(nym,nz)
       call cell2node_2d(h,h_node)
       call cell2node_2d(hs,hs_node)
       call cell2node_3dout
! h,hs,eta ---> cell
! h_node,hs_node,z --->node (0--nx)
! 
!      call cg_iric_write_sol_baseiterative_real_f &
!       (fid, "Input Discharge", qp, IER)
!      call cg_iric_write_sol_baseiterative_real_f &
!       (fid, "Downstream Stage", h_dw, IER)
       call output(fid)
!      write(44,'(20f10.5)')(h(i,2),i=1,8)
!      write(44,'(20f10.5)')(h_node(i,2),i=0,8)
!      write(44,'(20i10)')(int(ypn(i,3,nz+1)),i=0,5),&
!                         (int(ypn(i,3,nz+1)),i=nx-5,nx+1)
!      write(44,'(20i10)')(int(ypn(i,3,nz)),i=0,5),&
!                         (int(ypn(i,3,nz)),i=nx-5,nx+1)
!      write(44,'(20i10)')(int(ypn(i,3,nz-1)),i=0,5),&
!                         (int(ypn(i,3,nz-1)),i=nx-5,nx+1)
!      write(44,'(20f10.5)')(yun(i,3,nz),i=0,5),(yun(i,3,nz),i=nx-5,nx+1)
      end if
      itt=itt+1
!-----------------------------------------

      if(j_surf==0 .or. (j_surf==1 .and. time<stime_surf)) then
       hloop0=1
      else
       hloop0=hloop
      end if
!     write(*,*) 'hloop0=',hloop0;stop
      yun=yu; yvn=yv
!     call w1cal(yu,yv,yw)
      call omgpsi
      do m=1,hloop0
       call w1cal(yun,yvn,ywn) 
       if(j_surf==1 .and. time>stime_surf) then
        call hcal              ! h-->hn,hs
        call bound_p(smg_g)    ! hn,ypn--->ypn(i,ny+1)
       end if
       if(j_west==2.and.j_qin>=1) call uupcal !上流端の流速断面分布計算
       call omgpsi
       call sor(sorerr,lsor,soralpha,lmax)  !ypn--->ypn
       call rhs ! yu,yv,yw -->yun,yvn,ywn
       call bound_u(yun,ywn)
       call bound_v(yvn,ywn)
       call bound_w(yun,yvn,ywn)
!-----------------------------------------
       call diffusion  ! yun,yvn,ywn ---> yun,yvn,ywn
                       ! yc ----> ycn
!-----------------------------------------
       call bound_u(yun,ywn)
       call bound_v(yvn,ywn)
       call bound_w(yun,yvn,ywn)
       call newgrd_u(yun,yu,gux_n,guy_n,guz_n,gux,guy,guz)
       call newgrd_v(yvn,yv,gvx_n,gvy_n,gvz_n,gvx,gvy,gvz)
       call newgrd_w(ywn,yw,gwx_n,gwy_n,gwz_n,gwx,gwy,gwz)
!
       if(j_dens==1) then
        call bound_c(ycn,yun,yvn,ywn)
        call newgrd_c(ycn,yc,gcx_n,gcy_n,gcz_n,gcx,gcy,gcz)
       end if
!-----------------------------------------
       call w1cal(yun,yvn,ywn)
       call w2cal
!-----------------------------------------
       if(j_uadvec==1) then
        call dcip3d_u(yun,gux_n,guy_n,guz_n)
       else if(j_uadvec==2) then
        call upwind3d_u(yun,gux_n,guy_n,guz_n)
       else
        call o3upwind_u(yun,gux_n,guy_n,guz_n)
       end if
       call bound_u(yun,ywn)

       if(j_uadvec==1) then
        call dcip3d_v(yvn,gvx_n,gvy_n,gvz_n)
       else if(j_uadvec==2) then
        call upwind3d_v(yvn,gvx_n,gvy_n,gvz_n)
       else
        call o3upwind_v(yvn,gvx_n,gvy_n,gvz_n)
       end if
       call bound_v(yvn,ywn)

       if(j_uadvec==1) then
        call dcip3d_w(ywn,gwx_n,gwy_n,gwz_n)
       else if(j_uadvec==2) then
        call upwind3d_w(ywn,gwx_n,gwy_n,gwz_n)
       else
        call o3upwind_w(ywn,gwx_n,gwy_n,gwz_n)
       end if
       call bound_w(yun,yvn,ywn)

       if(j_dens==1) then
        if(j_cadvec==1) then
         call dcip3d_c(ycn,gcx_n,gcy_n,gcz_n)
        else if(j_uadvec==2) then
         call upwind3d_c(ycn,gcx_n,gcy_n,gcz_n)
        else
         call o3upwind_c(ycn,gcx_n,gcy_n,gcz_n)
        end if
        call bound_c(ycn,yun,yvn,ywn)
       end if
!-----------------------------------------
       total_1=0.; total_2=0.
       do i=1,nx
        do j=1,ny
         do k=1,nz+1
          total_1=total_1+abs(yun(i,j,k)-yuo(i,j,k))
          total_2=total_2+abs(yun(i,j,k))
          yuo(i,j,k)=yun(i,j,k)
         end do
        end do
       end do
       if(total_2==0.) then
        diff_ratio=0.
       else
        diff_ratio=total_1/total_2
       end if
       if(diff_ratio < hloop_err) exit
!---------------------------------------
      end do !--end of hloop0
!---------------------------------------
      call shift(yun,yu)
      call shift(yvn,yv)
      call shift(ywn,yw)
!     call shift(ypn,yp)
      if(j_dens==1) call shift(ycn,yc)
      call shift(gux_n,gux)
      call shift(guy_n,guy)
      call shift(guz_n,guz)
      call shift(gvx_n,gvx)
      call shift(gvy_n,gvy)
      call shift(gvz_n,gvz)
      call shift(gwx_n,gwx)
      call shift(gwy_n,gwy)
      call shift(gwz_n,gwz)
      if(j_dens==1) then
       call shift(gcx_n,gcx)
       call shift(gcy_n,gcy)
       call shift(gcz_n,gcz)
      end if
!--------------------------------------------
      if(j_surf==1 .and. time>stime_surf) then
       call hshift
      end if
!--------------------------------------------
      time=time+dt
      end do ! roop end of Main

      call cg_iric_close(fid,ier)
!     close(44)
!     close(45)
!     close(47)
      stop
      end program nays3dv
