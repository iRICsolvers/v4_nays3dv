module common_hh
    implicit none

    integer :: nx,ny,nz,i,j,k,ni,nj,nk,im,jm,km,nxm,nym,jj
    real(8)::g,snu,time,rho,skt,beta_t,skc,stime_surf,alpha_surf
    real(8)::dt,dz0,dz10,dz20
    real(8)::z0,sb,hs_ave,diam,ks,snm,kappa,xi1,al_ep
    real(8)::tmp0,con0
    real(8)::c0,c1,q_diff
    real(8)::h_dw_old,dhdt_dw,dh_up,dh_up0,dh_ref,q_stt,q_trn
    real(8)::dh_alpha,up_slope,up_wslope,ds0_center
    real(8)::eave_up,eave_dw,have_up_ini,have_dw_ini,hsave_up_ini,hsave_dw_ini
    real(8)::h_horizontal,dep_min
    integer::j_zgrid
    integer::j_dens,j_dens0,ic1,ic2,jc1,jc2,kc1,kc2,j_ini_dens,j_bc_dens
    integer::j_surf,j_snu,j_uadvec,j_cadvec
    integer::j_west,j_east,j_south,j_north
    integer::j_hinit,j_mindep

    real(8)::tuk,etime0,surf_tension,etime_q,etime_h,etime,st_dens
    real(8)::smg_g,pi
    real(8)::hloop_err
    real(8)::qp,qp0,width_in,u_ave,cd0,cd2,h_up,h_dw,amp_alpha,qadjust
    real(8)::hs_up,uave_up
    integer::icount,itout,itt
    integer::hloop,hloop0,m

    integer::j_qin,j_hup,j_hdw,n_qsize,n_hsize
    real(8)::q_up_const,h_up_const,h_dw_const,hd_amp,hd_wl,hd_st,hd_ap,hd_amp0

	! ======= Values for Bounday Concentration(CGNS I/O) ======
    integer::jc_in,indexmax_c
    integer,dimension(:),allocatable :: jc_inlen
    integer, dimension(:,:,:), allocatable:: indices_c
    double precision,dimension(:),allocatable :: boundary_con_value,boundary_con_up
    character(250),dimension(:),allocatable :: flowname_c
	!-----------------------------------------
    integer::jc_in_west,jc_in_east,jc_in_south,jc_in_north
    integer::jc_west_s(20),jc_west_e(20),jc_east_s(20),jc_east_e(20), &
             ic_south_s(20),ic_south_e(20),ic_north_s(20),ic_north_e(20)
    real(8)::c_bound_west(20),c_bound_east(20),c_bound_south(20), &
             c_bound_north(20)
    real(8)::c_up_west(20),c_up_east(20),c_up_south(20),c_up_north(20)
    integer::k_bc_west(20),k_bc_east(20),k_bc_south(20),k_bc_north(20)
end module common_hh

module common_geom
    implicit none
    real(8),dimension(:,:),allocatable::eta,eta_x
    real(8),dimension(:,:),allocatable::x8,y8,z8,h8,x,y,z,h_node,hs_node,dx,dy,dxi,dyj
    real(8),dimension(:,:),allocatable::ds,dn,xds2,yds2,coss,sins
    real(8),dimension(:,:),allocatable::curv,radi
    integer,dimension(:,:),allocatable::obst,obst4
    integer,dimension(:,:,:),allocatable::obst3d
    real(8),dimension(:,:),allocatable::hobst,hobst38,hobst_node8
    real(8),dimension(:),allocatable::dz,xi,dzk,xxi
    real(8),dimension(:,:,:),allocatable::x38,y38,z38
    real(8),dimension(:,:,:),allocatable::u38,v38,w38,p38,c38,sigma38,ob38,q38,snu_t38
    real(8),dimension(:,:,:),allocatable::vol_x,vol_y,vol_z
    ! real(8),dimension(:,:,:),allocatable::p38_cell,c38_cell
    real(8),dimension(:),allocatable::time_q,qt_up,time_h,ht_up
end module common_geom

module common_hyd
    implicit none
    real(8),dimension(:,:),allocatable::h,hn,hs,ypsurf,usta,uup0
    real(8),dimension(:,:,:),allocatable::yu,yun,yuo,yv,yvn,yw,ywn,yp,ypn
    real(8),dimension(:,:),allocatable::qu,qv
    real(8),dimension(:,:,:),allocatable::yc,ycn,w1,w2,a_p
    real(8),dimension(:,:,:),allocatable::omega,psi,snu_t,snu_xi_eg,snu_et_eg
    real(8),dimension(:,:),allocatable::a_n0,a_s0,a_w0,a_e0
    real(8),dimension(:,:,:),allocatable::a_s,a_n,a_e,a_w
    real(8),dimension(:,:,:),allocatable::work,ap,a_u,a_d,a_f
    real(8),dimension(:),allocatable::a_u0,a_d0
    real(8),dimension(:),allocatable::u00,q
end module common_hyd

module common_grad
    implicit none
    real(8),dimension(:,:,:),allocatable::gux,guy,guz,gvx,gvy,gvz,gwx,gwy,gwz
    real(8),dimension(:,:,:),allocatable::gcx,gcy,gcz
    real(8),dimension(:,:,:),allocatable::gux_n,guy_n,guz_n
    real(8),dimension(:,:,:),allocatable::gvx_n,gvy_n,gvz_n
    real(8),dimension(:,:,:),allocatable::gwx_n,gwy_n,gwz_n
    real(8),dimension(:,:,:),allocatable::gcx_n,gcy_n,gcz_n
end module common_grad
!================================================

module alloc_var_m
    use common_geom
    use common_hyd
    use common_grad
    implicit none
contains
!------------------------------------------------
    subroutine alloc_var(im,jm,km)
        integer::im,jm,km
        integer::i,j,k
        i=im; j=jm; k=km
        allocate(x(0:i,0:j),y(0:i,0:j),z(0:i,0:j),dx(0:i,0:j),dy(0:i,0:j) &
                ,dxi(0:i,0:j),dyj(0:i,0:j),h_node(0:i,0:j),hs_node(0:i,0:j))
        allocate(ds(0:i,0:j),dn(0:i,0:j),xds2(0:i,0:j),yds2(0:i,0:j))
        allocate(coss(0:i,0:j),sins(0:i,0:j),radi(0:i,0:j),curv(0:i,0:j))
        allocate(h(0:i,0:j),hn(0:i,0:j),hs(0:i,0:j),ypsurf(0:i,0:j),usta(0:i,0:j))
        allocate(eta(0:i,0:j),eta_x(0:i,0:j),obst(0:i,0:j),obst3d(0:i,0:j,0:k),hobst(0:i,0:j))
        allocate(yu(-1:i,-1:j,0:k),yun(-1:i,-1:j,0:k),yuo(0:i,0:j,0:k))
        allocate(yv(-1:i,-1:j,0:k),yvn(-1:i,-1:j,0:k))
        allocate(yw(-1:i,-1:j,0:k),ywn(-1:i,-1:j,0:k))
        allocate(w1(0:i,0:j,0:k),w2(0:i,0:j,0:k))

        allocate(qu(0:i,0:j),qv(0:i,0:j),q(0:i))
        allocate(omega(0:i,0:j,0:k),psi(0:i,0:j,0:k),snu_t(0:i,0:j,0:k),snu_xi_eg(0:i,0:j,0:k),snu_et_eg(0:i,0:j,0:k))
        allocate(yp(0:i,0:j,0:k),ypn(0:i,0:j,0:k))
        allocate(yc(-1:i,-1:j,0:k),ycn(-1:i,-1:j,0:k))

        allocate(gux(-1:i,-1:j,0:k),guy(-1:i,-1:j,0:k),guz(-1:i,-1:j,0:k))
        allocate(gvx(-1:i,-1:j,0:k),gvy(-1:i,-1:j,0:k),gvz(-1:i,-1:j,0:k))
        allocate(gwx(-1:i,-1:j,0:k),gwy(-1:i,-1:j,0:k),gwz(-1:i,-1:j,0:k))
        allocate(gcx(-1:i,-1:j,0:k),gcy(-1:i,-1:j,0:k),gcz(-1:i,-1:j,0:k))

        allocate(gux_n(-1:i,-1:j,0:k),guy_n(-1:i,-1:j,0:k),guz_n(-1:i,-1:j,0:k))
        allocate(gvx_n(-1:i,-1:j,0:k),gvy_n(-1:i,-1:j,0:k),gvz_n(-1:i,-1:j,0:k))
        allocate(gwx_n(-1:i,-1:j,0:k),gwy_n(-1:i,-1:j,0:k),gwz_n(-1:i,-1:j,0:k))
        allocate(gcx_n(-1:i,-1:j,0:k),gcy_n(-1:i,-1:j,0:k),gcz_n(-1:i,-1:j,0:k))

        allocate(dz(0:k),dzk(0:k),xi(0:k),xxi(0:k),u00(0:k))
        allocate(uup0(0:j,0:k))
        allocate(a_n0(0:i,0:j),a_s0(0:i,0:j),a_w0(0:i,0:j) &
                ,a_e0(0:i,0:j),a_u0(0:k),a_d0(0:k))
        allocate(work(0:i,0:j,0:k),a_p(0:i,0:j,0:k) &
                ,a_e(0:i,0:j,0:k),a_w(0:i,0:j,0:k),a_s(0:i,0:j,0:k),a_n(0:i,0:j,0:k) &
                ,a_u(0:i,0:j,0:k),a_d(0:i,0:j,0:k),a_f(0:i,0:j,0:k))

    end subroutine alloc_var

end module alloc_var_m

!================================================
module initial_m
    use common_hh
    use common_geom
    use common_hyd
    use common_grad
    implicit none
contains
!------------------------------------------------
    subroutine initial
        real(8)::cell_center_height,zz
        real(8)::ds0,dx0,dy0
        integer::jj,inum

        h=0.; hn=0.; hs=0.;obst=0; obst3d=0; hobst=0.; ypsurf=0.;z=0.;h_node=0.; hs_node=0.
        usta=0.
        yu=0.;yun=0.;yv=0.;yvn=0.;yw=0.;ywn=0.; yp=0.; ypn=0.
        yc=0.; ycn=0.
        yuo=0.
        omega=0.; psi=0.; snu_t=0.; snu_xi_eg=0.; snu_et_eg=0.
        gux=0.;guy=0.;guz=0.;gvx=0.;gvy=0.;gvz=0.;gwx=0.;gwy=0.;gwz=0.
        gcx=0.;gcy=0.;gcz=0.

        gux_n=0.;guy_n=0.;guz_n=0.
        gvx_n=0.;gvy_n=0.;gvz_n=0.
        gwx_n=0.;gwy_n=0.;gwz_n=0.
        gcx_n=0.;gcy_n=0.;gcz_n=0.

        inum=0; hs_ave=0.
        jc_in_west=0; jc_in_east=0; jc_in_south=0; jc_in_north=0

        dh_up=0.
        eave_up=0.;eave_dw=0.;have_up_ini=0.;have_dw_ini=0.;hsave_up_ini=0.;hsave_dw_ini=0.
        !
        !Šiqƒf[ƒ^[‚ğŒvZ•Ï”‚Ö‚ÌŠ„‚è•t‚¯
        !
        do i=0,nx
            do j=0,ny
                x(i,j)=x8(i+1,j+1)
                y(i,j)=y8(i+1,j+1)
                z(i,j)=z8(i+1,j+1)
                h_node(i,j)=h8(i+1,j+1)
                if(j_mindep==1.and.(z(i,j).gt.h_node(i,j)-dep_min)) then
                    z(i,j)=h_node(i,j)-dep_min
                    z8(i+1,j+1)=z(i,j)
                end if
            end do
        end do
        ! z,h_node ---> node
        ! eta,h ---> cell


        do i=1,nx
            do j=1,ny
                obst(i,j)=obst4(i,j)
                hobst(i,j)=hobst38(i,j)
                eta(i,j)=(z(i,j)+z(i-1,j-1)+z(i,j-1)+z(i-1,j))*.25
                h(i,j)=(h_node(i,j)+h_node(i-1,j-1)+h_node(i,j-1)+h_node(i-1,j))*.25
                if(obst(i,j)==0) then
                    hobst(i,j)=eta(i,j)
                end if
            end do
        end do

        eave_up=0.;eave_dw=0.
        have_up_ini=0.;have_dw_ini=0.
        do j=1,ny
            if(obst(1,j)==0) then
                eave_up=eave_up+eta(1,j)
                have_up_ini=have_up_ini+h(1,j)
                eave_dw=eave_dw+eta(nx,j)
                have_dw_ini=have_dw_ini+h(nx,j)
            end if
        end do
        eave_up=eave_up/float(ny) !‰Šúã—¬•½‹Ï‰Í°
        eave_dw=eave_dw/float(ny) !‰Šú‰º—¬•½‹Ï‰Í°
        have_up_ini=have_up_ini/float(ny) !‰Šúã—¬•½‹Ï…ˆÊ
        have_dw_ini=have_dw_ini/float(ny) !‰Šú‰º—¬•½‹Ï…ˆÊ
        hsave_up_ini=have_up_ini-eave_up  !‰Šúã—¬•½‹Ï…[
        hsave_dw_ini=have_dw_ini-eave_dw  !‰Šú‰º—¬•½‹Ï…[

        do i=1,nx
            eta(i,0)=eta(i,1)
            h(i,0)=h(i,1)
            eta(i,ny+1)=eta(i,ny)
            h(i,ny+1)=h(i,ny)
        end do

        do j=0,ny+1
            eta(0,j)=eta(1,j)
            h(0,j)=h(1,j)
            eta(nx+1,j)=eta(nx,j)
            h(nx+1,j)=h(nx,j)
        end do

        do i=0,nx+1
            do j=0,ny+1
                hs(i,j)=h(i,j)-eta(i,j)
            end do
        end do

        do i=1,nx+1
            do j=1,ny+1
                if(i==1)then
                    if(j==1) then
                        hobst_node8(i,j)=hobst(i,j)
                    else if(j==ny+1) then
                        hobst_node8(i,j)=hobst(i,j-1)
                    else
                        hobst_node8(i,j)=(hobst(i,j)+hobst(i,j-1))*.5
                    end if

                else if(i==nx+1) then
                    if(j==1) then
                        hobst_node8(i,j)=hobst(i-1,j)
                    else if(j==ny+1) then
                        hobst_node8(i,j)=hobst(i-1,j-1)
                    else
                        hobst_node8(i,j)=(hobst(i-1,j)+hobst(i-1,j-1))*.5
                    end if

                else if(j==1) then
                    hobst_node8(i,j)=(hobst(i,j)+hobst(i-1,j))*.5
                else if(j==ny+1) then
                    hobst_node8(i,j)=(hobst(i,j-1)+hobst(i-1,j-1))*.5

                else if(i>=2.and.i<=nx.and.j>=2.and.j<=ny) then
                    hobst_node8(i,j)=(hobst(i,j)+hobst(i-1,j)&
                            +hobst(i,j-1)+hobst(i-1,j-1))*.25
                end if
            end do
        end do

        if(j_zgrid==1) then
            do k=0,nz+1
                dz(k)=1./float(nz)
            end do
        else
            dz20=(1.-dz10)/float(nz-1)
            ! write(*,*) dz10,dz20
            do k=0,1
                dz(k)=dz10
            end do
            do k=2,nz+1
                dz(k)=dz20
            end do
        end if

        xi(0)=0.
        do k=1,nz+1
            xi(k)=xi(k-1)+dz(k)
        end do

        do k=0,nz
            dzk(k)=(dz(k)+dz(k+1))*.5
        end do

        dzk(nz+1)=dz(nz)
        do k=1,nz+1
            xxi(k)=(xi(k)+xi(k-1))*.5
        end do
        xxi(0)=xxi(1)

        xi1=dz(1)*.5

        do i=1,nx
            do j=1,ny
                if(obst(i,j)==1) then
                    do k=1,nz
                        cell_center_height=eta(i,j)+hs(i,j)*xxi(k)
                        if(cell_center_height>hobst(i,j)) then
                            obst3d(i,j,k)=0
                        else
                            obst3d(i,j,k)=1
                            if(k==nz) then
                                obst3d(i,j,k+1)=1
                                if(i==1) obst3d(0,j,k+1)=1
                                if(i==nx) obst3d(nx+1,j,k+1)=1
                            end if
                            if(k==1) then
                                obst3d(i,j,0)=1
                                if(i==1)  obst3d(0,j,0)=1
                                if(i==nx) obst3d(nx+1,j,0)=1
                                if(j==1)  obst3d(i,0,0)=0
                                if(j==ny) obst3d(i,ny+1,0)=1
                            end if
                            if(i==1) obst3d(0,j,k)=1
                            if(i==nx) obst3d(nx+1,j,k)=1
                            if(j==1) obst3d(i,0,k)=1
                            if(j==ny) obst3d(i,ny+1,k)=1
                            if(i==1.and.j==1) obst3d(0,0,k)=1;obst3d(0,0,0)=1
                            if(i==1.and.j==ny) obst3d(0,ny+1,k)=1;obst3d(0,ny+1,0)=1
                            if(i==nx.and.j==1) obst3d(nx+1,0,k)=1;obst3d(nx+1,0,0)=1
                            if(i==nx.and.j==ny) obst3d(nx+1,ny+1,k)=1;obst3d(nx+1,ny+1,0)=1
                        end if
                    end do
                end if
            end do
        end do

        !
        ! ds,dn‚ÌŒvZ
        !
        do i=0,nx
            do j=0,ny
                if(i>0) ds(i,j)=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)
                if(j>0) dn(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
            end do
        end do
        ds0_center=ds(1,nym)
        !
        ! coss,sins‚ÌŒvZ
        !
        do i=0,nx
            do j=0,ny
                if(i==0) then
                    coss(i,j)=(x(i+1,j)-x(i,j))/ds(i+1,j)
                    sins(i,j)=(y(i+1,j)-y(i,j))/ds(i+1,j)
                else if(i==nx) then
                    coss(i,j)=(x(i,j)-x(i-1,j))/ds(i,j)
                    sins(i,j)=(y(i,j)-y(i-1,j))/ds(i,j)
                else
                    coss(i,j)=((x(i,j)-x(i-1,j))/ds(i,j)+(x(i+1,j)-x(i,j))/ds(i+1,j))*.5
                    sins(i,j)=((y(i,j)-y(i-1,j))/ds(i,j)+(y(i+1,j)-y(i,j))/ds(i+1,j))*.5
                end if
            end do
        end do

        !
        ! xds2,yds2‚ÌŒvZ
        !
        do j=0,ny
            do i=1,nx-1
                xds2(i,j)=((x(i+1,j)-x(i,j))/ds(i+1,j) &
                         -(x(i,j)-x(i-1,j))/ds(i,j))*2./(ds(i+1,j)+ds(i,j))
                yds2(i,j)=((y(i+1,j)-y(i,j))/ds(i+1,j) &
                         -(y(i,j)-y(i-1,j))/ds(i,j))*2./(ds(i+1,j)+ds(i,j))
            end do
        end do

        !
        ! ‹È—¦‚ÌŒvZ@curv  ‹È—¦”¼Œa radi
        !
        do j=0,ny
            do i=1,nx-1
                ds0=ds(i+1,j)+ds(i,j)
                dx0=(x(i+1,j)-x(i-1,j))
                dy0=(y(i+1,j)-y(i-1,j))
                if (abs(dy0)>1e-7) then
                    curv(i,j)=-xds2(i,j)*ds0/dy0
                else if(abs(dx0)>1e-7) then
                    curv(i,j)=yds2(i,j)*ds0/dx0
                else
                    curv(i,j)=0.
                end if
                if (abs(curv(i,j))>1e-7) then
                    radi(i,j)=1./curv(i,j)
                else
                    radi(i,j)=0.
                end if
            end do
        end do

        !
        ! x,y,dx==ds,dy==dn‚È‚ÇÀ•W‚Ìİ’è
        !
        do i=0,nx
            do j=0,ny
                if(i>0) then
                    dx(i,j)=ds(i,j)
                    if(i==1) dx(0,j)=dx(1,j)
                    if(i==nx) dx(nx+1,j)=dx(nx,j)
                end if
                if(j>0) then
                    dy(i,j)=dn(i,j)
                    if(j==1) dy(i,0)=dy(i,1)
                    if(j==ny) dy(i,ny+1)=dy(i,ny)
                end if
            end do
        end do
        !
        ! dxi ----> u‚ÌŒvZ“_‚Ìds
        !
        do i=0,nx
            do j=1,ny
                dxi(i,j)=(dx(i,j)+dx(i+1,j)+dx(i,j-1)+dx(i+1,j-1))*.25
            end do
            dxi(i,0)=dxi(i,1)
            dxi(i,ny+1)=dxi(i,ny)
        end do
        !
        ! dyj ----> v‚ÌŒvZ“_‚Ìdn
        do i=1,nx
            do j=0,ny
                dyj(i,j)=(dy(i,j)+dy(i,j+1)+dy(i-1,j)+dy(i-1,j+1))*.25
            end do
        end do
        do j=0,ny
            dyj(0,j)=dyj(1,j)
            dyj(nx+1,j)=dyj(nx,j)
        end do
        !
        ! —¬“ü•”‚Ì…——Ê‚Ì‰Šúİ’è‚¨‚æ‚Ñ•½‹Ï‘Î”‘¥•ª•z‚Ìİ’è
        !
        width_in=0.; hs_up=0.
        jj=0
        do j=1,ny
            if(obst(i,1)==0) then
                jj=jj+1
                width_in=width_in+dy(0,j)
                hs_up=hs_up+hs(1,j)
            end if
        end do
        if(jj>0) then
            hs_up=hs_up/float(jj)
        else
            write(*,*) 'No enterance found'
            stop
        end if

        ks=2.*diam; z0=ks/30.; snm=0.0486*diam**(1./6.)
        sb=log(hs_up/z0)-1.
        u_ave=qp/(hs_up*width_in)
        ! do k=1,nz
        !     zz=xxi(k)*hs_up
        !     u00(k)=u_ave*log(zz/z0)/sb
        ! end do
        ! u00(0)=-u00(1); u00(nz+1)=u00(nz)

        !
        ! cd0 –€CŒW”‚ÌŒˆ’è
        !
        inum=0
        hs_ave=0.
        do i=1,nx
            do j=1,ny
                if(obst(i,j)==0) then
                    hs_ave=hs_ave+hs(i,j)
                    inum=inum+1
                end if
            end do
        end do
        hs_ave=hs_ave/dble(inum)
        ! cd0=1./(1./kappa*log(xi1*hs_ave/z0))**2
        cd0=1./(1./kappa*log(xi1*hs_ave/z0))**2
        cd2=cd0**2

        write(*,*) diam,z0,cd0,cd2

        yu=0.
        qu=0.
        q=0.

        do i=1,nx
            do j=1,ny
                if(i.eq.nx) then
                    a_n0(i,j)=0.
                else
                    a_n0(i,j)=1./(dxi(i,j)*(dxi(i,j)+dxi(i-1,j)))
                end if
                if(i.eq.1) then
                    a_s0(i,j)=0.
                else
                    a_s0(i,j)=1./(dxi(i-1,j)*(dxi(i,j)+dxi(i-1,j)))
                end if
                if(j.eq.ny) then
                    a_w0(i,j)=0.
                else
                    a_w0(i,j)=1./(dyj(i,j)*(dyj(i,j)+dyj(i,j-1)))
                end if
                if(j.eq.1) then
                    a_e0(i,j)=0.
                else
                    a_e0(i,j)=1./(dyj(i,j-1)*(dyj(i,j)+dyj(i,j-1)))
                end if
            end do
        end do

        do k=1,nz
            if(k.eq.1) then
                a_d0(k)=0.
            else
                a_d0(k)=1./dz(k)**2
            end if
            if(k.eq.nz) then
                a_u0(k)=0.
            else
                a_u0(k)=1./dz(k)**2
            end if
        end do

        do i=1,nx
            do j=1,ny
                do k=1,nz
                    if(obst3d(i,j,k)==1) then
                        a_n(i,j,k)=0.; a_s(i,j,k)=0.; a_w(i,j,k)=0.; a_e(i,j,k)=0.
                    else
                        if(i==nx .or. obst3d(i+1,j,k)==1) then
                            a_n(i,j,k)=0.
                        else
                            a_n(i,j,k)=a_n0(i,j)*(hs(i,j)+hs(i+1,j))
                        end if
                        if(i==1 .or. obst3d(i-1,j,k)==1) then
                            a_s(i,j,k)=0.
                        else
                            a_s(i,j,k)=a_s0(i,j)*(hs(i,j)+hs(i-1,j))
                        end if
                        if(j==ny .or. obst3d(i,j+1,k)==1) then
                            a_w(i,j,k)=0.
                        else
                            a_w(i,j,k)=a_w0(i,j)*(hs(i,j)+hs(i,j+1))
                        end if
                        if(j==1 .or.obst3d(i,j-1,k)==1) then
                            a_e(i,j,k)=0.
                        else
                            a_e(i,j,k)=a_e0(i,j)*(hs(i,j)+hs(i,j-1))
                        end if
                        if(k==nz .or. obst3d(i,j,k+1)==1) then
                            a_u(i,j,k)=0.
                        else
                            a_u(i,j,k)=a_u0(k)/hs(i,j)
                        end if
                        if(k==1 .or. obst3d(i,j,k-1)==1) then
                            a_d(i,j,k)=a_d0(k)/hs(i,j)
                        end if
                            a_p(i,j,k)=a_n(i,j,k)+a_s(i,j,k)+a_w(i,j,k)+a_e(i,j,k) &
                                      +a_u(i,j,k)+a_d(i,j,k)
                        if(k==nz) a_p(i,j,nz)=a_p(i,j,nz)+2./(hs(i,j)*dz(nz)**2)
                    end if
                end do
            end do
        end do

        if(j_dens==1) then
            do i=0,nx+1
                do j=0,ny+1
                    do k=0,nz+1
                        yc(i,j,k)=c0
                    end do
                end do
            end do
        end if

        if(j_dens==1.and.j_ini_dens==1) then
            do i=1,nx
                do j=1,ny
                    do k=1,nz
                        if(i>=ic1 .and. i<=ic2 .and. j>=jc1 .and. j<=jc2 .and. k>=kc1 .and. k<=kc2) then
                            if(obst3d(i,j,k)==0) then
                                yc(i,j,k)=c1
                            else
                                yc(i,j,k)=c0
                            end if
                        else
                            yc(i,j,k)=c0
                        end if
                        ycn(i,j,k)=yc(i,j,k)
                    end do
                end do
            end do
            !
            ! Boundary
            !
            ! West and East
            do j=0,ny+1
                do k=0,nz+1
                    yc(0,j,k)=c0;ycn(0,j,k)=c0
                    yc(nx+1,j,k)=c0;ycn(nx+1,j,k)=c0
                end do
            end do
            ! South and North
            do i=0,nx+1
                do k=0,nz+1
                    yc(i,0,k)=c0;ycn(i,0,k)=c0
                    yc(i,ny+1,k)=c0;ycn(i,ny+1,k)=c0
                end do
            end do
            ! Top and Bottom
            do i=0,nx+1
                do j=0,ny+1
                    yc(i,j,0)=c0;ycn(i,j,0)=c0
                    yc(i,j,nz+1)=c0;ycn(i,j,nz+1)=c0
                end do
            end do
            !
            ! Gradients
            !
            do i=1,nx
                do j=1,ny
                    do k=1,nz
                        gcx(i,j,k)=(yc(i+1,j,k)-yc(i-1,j,k))/(dxi(i,j)+dxi(i-1,j))
                        gcy(i,j,k)=(yc(i,j+1,k)-yc(i,j-1,k))/(dyj(i,j)+dyj(i,j-1))
                        gcz(i,j,k)=(yc(i,j,k+1)-yc(i,j,k-1))/(dzk(k)+dzk(k-1))
                        gcx_n(i,j,k)=gcx(i,j,k)
                        gcy_n(i,j,k)=gcy(i,j,k)
                        gcz_n(i,j,k)=gcz(i,j,k)
                    end do
                end do
            end do

        end if
        !
        ! Values for 3doutput
        !
        x38=0.;y38=0.;z38=0.
        u38=0.;v38=0.;w38=0.;c38=0.;p38=0.;snu_t38=0.
        ! c38_cell=0.; p38_cell=0.


        do i=1,ni
            do j=1,nj
                do k=1,nk
                    x38(i,j,k)=x8(i,j)
                    y38(i,j,k)=y8(i,j)
                    sigma38(i,j,k)=xi(k-1)
                    z38(i,j,k)=z8(i,j)+(h8(i,j)-z8(i,j))*xi(k-1)
                end do
            end do
        end do

    end subroutine initial
end module initial_m


    !--------------------------------------------------
    !Šù‘¶‚Ì"stop"‚ğ‚±‚ÌŠÖ”‚É’u‚«Š·‚¦‚Ä‚ ‚°‚é
    !--------------------------------------------------
    subroutine calc_stop(fid)
        use common_hh
        use iric
        integer :: ier,fid
        call cg_iric_close(fid,ier)
        stop
    end subroutine
