module cip3d_m
    use common_hh
    use common_geom
    use common_hyd
    implicit none
    real(8)::dx1,dx2,dx3,dtdx,dy1,dy2,dy3,dtdy,dz1,dz2,dz3,dtdz
    real(8)::cx,cy,cz,xx,yy,zz,x1,y1,z1
    real(8)::a01,a11,a02,a12,a03,a13,b1,b2,b3,a05,a04,a14,a09,a08,a15
    real(8)::a06,a07,a16,a10
    real(8)::gxo,gyo,gzo
    integer::im1,jm1,km1
contains

    !***************************************************
    subroutine dcip3d_u(f,gx,gy,gz)
    !***************************************************

        real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
        real(8),dimension(-1:im,-1:jm,0:km)::fn,gxn,gyn,gzn,u,v,w
        real(8)::hs_up

        fn=0.; gxn=0.; gyn=0.; u=0.; v=0.; w=0.

        do i=0,nx
            do j=1,ny
                hs_up=(hs(i,j)+hs(i+1,j))*.5
                do k=1,nz
                    u(i,j,k)=yu(i,j,k)
                    v(i,j,k)=(yv(i,j,k)+yv(i+1,j,k)+yv(i,j-1,k)+yv(i+1,j-1,k))*.25
                    w(i,j,k)=(w2(i,j,k)+w2(i+1,j,k)+w2(i,j,k-1)+w2(i+1,j,k-1))*.25/hs_up
                end do
            end do
        end do

        do i=0,nx
            do j=1,ny
                dx1=dxi(i,j)
                dx2 =dx1*dx1
                dx3 =dx2*dx1
                dtdx=dt/dx1
                dy1=dy(i,j)
                dy2 =dy1*dy1
                dy3 =dy2*dy1
                dtdy=dt/dy1
                do k=1,nz
                    dz1=dz(k)
                    dtdz=dt/dz1
                    dz2 =dz1*dz1
                    dz3 =dz2*dz1

                    cx=u(i,j,k)
                    cy=v(i,j,k)
                    cz=w(i,j,k)
                    xx = - cx*dt
                    yy = - cy*dt
                    zz = - cz*dt

                    x1 = -sign(1.0,cx) !u>0‚Å•‰, u<0‚Å³
                    y1 = -sign(1.0,cy)
                    z1 = -sign(1.0,cz)
                    im1=i+int(x1) !u<0‚Åi+1, u>0‚Åi-1
                    jm1=j+int(y1)
                    ! if(j_north==1.and.jm1>ny) jm1=ny
                    ! if(j_south==1.and.jm1<1)  jm1=1
                    km1=k+int(z1)

                    a01 = ( ( gx(i,j,k)+gx(im1,j,k) )*dx1*x1 &
                        +2.0*( f(i,j,k)-f(im1,j,k) ) )/(dx3*x1)
                    a11 = ( 3.0*( f(im1,j,k)-f(i,j,k) ) &
                        - ( gx(im1,j,k)+2.*gx(i,j,k) )*dx1*x1 )/dx2
                    a02 = ( ( gy(i,j,k)+gy(i,jm1,k) )*dy1*y1 &
                        +2.0*( f(i,j,k)-f(i,jm1,k) ) )/(dy3*y1)
                    a12 = ( 3.0*( f(i,jm1,k)-f(i,j,k) ) &
                        - ( gy(i,jm1,k)+2.*gy(i,j,k) )*dy1*y1 )/dy2
                    a03 = ( ( gz(i,j,k)+gz(i,j,km1) )*dz1*z1 &
                        +2.0*( f(i,j,k)-f(i,j,km1) ) )/(dz3*z1)
                    a13 = ( 3.0*( f(i,j,km1)-f(i,j,k) ) &
                        - ( gz(i,j,km1)+2.*gz(i,j,k) )*dz1*z1 )/dz2
                    b1  = f(i,j,k)-f(im1,j,k)-f(i,jm1,k)+f(im1,jm1,k)
                    b2  = f(i,j,k)-f(i,jm1,k)-f(i,j,km1)+f(i,jm1,km1)
                    b3  = f(i,j,k)-f(im1,j,k)-f(i,j,km1)+f(im1,j,km1)
                    a05 = ( b1-(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy2)
                    a04 = ( b1-(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1)/(dx2*dy1*y1)
                    a14 = (-b1+(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1 &
                        +(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy1*y1)
                    a09 = ( b2-(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz2)
                    a08 = ( b2-(-gy(i,j,k)+gy(i,j,km1))*dy1*y1)/(dy2*dz1*z1)
                    a15 = (-b2+(-gy(i,j,k)+gy(i,j,km1))*dy1*y1+(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz1*z1)
                    a06 = ( b3-(-gz(i,j,k)+gz(im1,j,k))*dz1*z1)/(dx1*x1*dz2)
                    a07 = ( b3-(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx2*dz1*z1)
                    a16 = (-b3+(-gz(i,j,k)+gz(im1,j,k))*dz1*z1+(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx1*x1*dz1*z1)
                    a10 = ( -f(i,j,k) + (f(im1,j,k)+f(i,jm1,k)+f(i,j,km1)) &
                        - (f(im1,jm1,k)+f(i,jm1,km1)+f(im1,j,km1)) &
                        + f(im1,jm1,km1) ) / (dx1*x1*dy1*y1*dz1*z1)

                    ! if((i==0.and.im1==-1).or.(i==nx.and.im1==nx+1)) then
                    !     fn(i,j,k) = ((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gy(i,j,k))*yy &
                    !               +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gz(i,j,k))*zz &
                    !               +a10*xx*yy*zz+f(i,j,k)
                    ! else

                        fn(i,j,k) = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gx(i,j,k))*xx &
                                  +((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gy(i,j,k))*yy &
                                  +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gz(i,j,k))*zz &
                                  +a10*xx*yy*zz+f(i,j,k)
                    ! end if

                    gxn(i,j,k)= (3.*a01*xx+2.*(a04*yy+a07*zz+a11))*xx+(a05*yy+a10*zz+a14)*yy+(a06*zz+a16)*zz+gx(i,j,k)
                    gyn(i,j,k)= (3.*a02*yy+2.*(a05*xx+a08*zz+a12))*yy+(a09*zz+a10*xx+a15)*zz+(a04*xx+a14)*xx+gy(i,j,k)
                    gzn(i,j,k)= (3.*a03*zz+2.*(a06*xx+a09*yy+a13))*zz+(a07*xx+a10*yy+a16)*xx+(a08*yy+a15)*yy+gz(i,j,k)
                end do
            end do
        end do

        do j=1,ny
            do i=0,nx
                dtdx=dt/dxi(i,j)
                dtdy=dt/dy(i,j)
                do k=1,nz
                    dtdz=dt/dz(k)
                    f(i,j,k)  = fn(i,j,k)
                    if(i==0) then
                        gx(i,j,k) = 0.
                    else if(i==nx) then
                        gx(i,j,k) = 0.
                    else
                        gxo=(f(i+1,j,k)-f(i-1,j,k))/(2.*dx1)
                        gx(i,j,k) = gxn(i,j,k) &
                                  -(gxo*(u(i+1,j,k)-u(i-1,j,k))+gyo*(v(i+1,j,k)-v(i-1,j,k)) &
                                  +gzn(i,j,k)*(w(i+1,j,k)-w(i-1,j,k)))*0.5*dtdx
                    end if
                    if(j==1.or.j==ny) then
                        gyo=0.
                    else
                        gyo=(f(i,j+1,k)-f(i,j-1,k))/(2.*dy1)
                    end if

                    gy(i,j,k) = gyn(i,j,k) &
                              -(gxo*(u(i,j+1,k)-u(i,j-1,k))+gyo*(v(i,j+1,k)-v(i,j-1,k)) &
                              +gzo*(w(i,j+1,k)-w(i,j-1,k)))*0.5*dtdy
                    gzo=(f(i,j,k+1)-f(i,j,k-1))/(2.*dz1)
                    gz(i,j,k) = gzn(i,j,k) &
                              -(gxo*(u(i,j,k+1)-u(i,j,k-1))+gyo*(v(i,j,k+1)-v(i,j,k-1)) &
                              +gzo*(w(i,j,k+1)-w(i,j,k-1)))*0.5*dtdz
                end do
            end do
        end do
    end subroutine dcip3d_u

    !***************************************************
    subroutine dcip3d_v(f,gx,gy,gz)
    !***************************************************

        real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
        real(8),dimension(-1:im,-1:jm,0:km)::fn,gxn,gyn,gzn,u,v,w
        real(8)::hs_vp

        fn=0.; gxn=0.; gyn=0.; gzn=0.; u=0.; v=0.; w=0.

        do i=1,nx
            do j=0,ny
                hs_vp=(hs(i,j)+hs(i,j+1))*.5
                do k=1,nz
                    u(i,j,k)=(yu(i,j,k)+yu(i-1,j,k)+yu(i,j+1,k)+yu(i-1,j+1,k))*.25
                    v(i,j,k)=yv(i,j,k)
                    w(i,j,k)=(w2(i,j,k)+w2(i,j+1,k)+w2(i,j,k-1)+w2(i,j+1,k-1)) &
                            *.25/hs_vp
                end do
            end do
        end do

        do i=1,nx
            do j=0,ny
                dx1=dx(i,j)
                dx2 =dx1*dx1
                dx3 =dx2*dx1
                dtdx=dt/dx1
                dy1=dyj(i,j)
                dy2 =dy1*dy1
                dy3 =dy2*dy1
                dtdy=dt/dy1
                do k=1,nz
                    dz1=dz(k)
                    dtdz=dt/dz1
                    dz2 =dz1*dz1
                    dz3 =dz2*dz1

                    cx=u(i,j,k)
                    cy=v(i,j,k)
                    cz=w(i,j,k)
                    xx = - cx*dt
                    yy = - cy*dt
                    zz = - cz*dt

                    x1 = -sign(1.0,cx)
                    y1 = -sign(1.0,cy)
                    z1 = -sign(1.0,cz)
                    im1=i+int(x1)
                    jm1=j+int(y1)
                    km1=k+int(z1)

                    a01 = ( ( gx(i,j,k)+gx(im1,j,k) )*dx1*x1+2.0*( f(i,j,k)-f(im1,j,k) ) )/(dx3*x1)
                    a11 = ( 3.0*( f(im1,j,k)-f(i,j,k))-( gx(im1,j,k)+2.*gx(i,j,k) )*dx1*x1 )/dx2
                    a02 = ( ( gy(i,j,k)+gy(i,jm1,k) )*dy1*y1+2.0*( f(i,j,k)-f(i,jm1,k) ) )/(dy3*y1)
                    a12 = ( 3.0*( f(i,jm1,k)-f(i,j,k) )-( gy(i,jm1,k)+2.*gy(i,j,k) )*dy1*y1 )/dy2
                    a03 = ( ( gz(i,j,k)+gz(i,j,km1) )*dz1*z1+2.0*( f(i,j,k)-f(i,j,km1) ) )/(dz3*z1)
                    a13 = ( 3.0*( f(i,j,km1)-f(i,j,k) )-( gz(i,j,km1)+2.*gz(i,j,k) )*dz1*z1 )/dz2
                    b1  = f(i,j,k)-f(im1,j,k)-f(i,jm1,k)+f(im1,jm1,k)
                    b2  = f(i,j,k)-f(i,jm1,k)-f(i,j,km1)+f(i,jm1,km1)
                    b3  = f(i,j,k)-f(im1,j,k)-f(i,j,km1)+f(im1,j,km1)
                    a05 = ( b1-(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy2)
                    a04 = ( b1-(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1)/(dx2*dy1*y1)
                    a14 = (-b1+(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1+(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy1*y1)
                    a09 = ( b2-(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz2)
                    a08 = ( b2-(-gy(i,j,k)+gy(i,j,km1))*dy1*y1)/(dy2*dz1*z1)
                    a15 = (-b2+(-gy(i,j,k)+gy(i,j,km1))*dy1*y1+(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz1*z1)
                    a06 = ( b3-(-gz(i,j,k)+gz(im1,j,k))*dz1*z1)/(dx1*x1*dz2)
                    a07 = ( b3-(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx2*dz1*z1)
                    a16 = (-b3+(-gz(i,j,k)+gz(im1,j,k))*dz1*z1+(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx1*x1*dz1*z1)
                    a10 = ( -f(i,j,k) + (f(im1,j,k)+f(i,jm1,k)+f(i,j,km1)) &
                        - (f(im1,jm1,k)+f(i,jm1,km1)+f(im1,j,km1)) &
                        + f(im1,jm1,km1) ) / (dx1*x1*dy1*y1*dz1*z1)

                    ! if((j==0.and.jm1==-1).or.(j==ny.and.jm1==ny+1)) then
                    !     fn(i,j,k) = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gx(i,j,k))*xx &
                    !               +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gz(i,j,k))*zz &
                    !               +a10*xx*yy*zz+f(i,j,k)
                    ! else
                        fn(i,j,k) = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gx(i,j,k))*xx &
                                  +((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gy(i,j,k))*yy &
                                  +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gz(i,j,k))*zz &
                                  +a10*xx*yy*zz+f(i,j,k)
                    ! end if
                    gxn(i,j,k)= (3.*a01*xx+2.*(a04*yy+a07*zz+a11))*xx &
                              +(a05*yy+a10*zz+a14)*yy+(a06*zz+a16)*zz+gx(i,j,k)
                    gyn(i,j,k)= (3.*a02*yy+2.*(a05*xx+a08*zz+a12))*yy &
                              +(a09*zz+a10*xx+a15)*zz+(a04*xx+a14)*xx+gy(i,j,k)
                    gzn(i,j,k)= (3.*a03*zz+2.*(a06*xx+a09*yy+a13))*zz &
                              +(a07*xx+a10*yy+a16)*xx+(a08*yy+a15)*yy+gz(i,j,k)
                end do
            end do
        end do

        do j=0,ny
            do i=1,nx
                dtdx=dt/dx(i,j)
                dtdy=dt/dyj(i,j)
                do k=1,nz
                    dtdz=dt/dz(k)
                    f(i,j,k)  = fn(i,j,k)
                    gxo=(f(i+1,j,k)-f(i-1,j,k))/(2.*dx1)
                    gx(i,j,k) = gxn(i,j,k) &
                              -(gxo*(u(i+1,j,k)-u(i-1,j,k))+gyo*(v(i+1,j,k)-v(i-1,j,k)) &
                              +gzo*(w(i+1,j,k)-w(i-1,j,k)))*0.5*dtdx
                    if(j==0) then
                        gy(i,j,k) = 0.
                    else if(j==ny) then
                        gy(i,j,k) = 0.
                    else
                        gyo=(f(i,j+1,k)-f(i,j-1,k))/(2.*dy1)
                        gy(i,j,k) = gyn(i,j,k) &
                                  -(gxo*(u(i,j+1,k)-u(i,j-1,k))+gyo*(v(i,j+1,k)-v(i,j-1,k)) &
                                  +gzo*(w(i,j+1,k)-w(i,j-1,k)))*0.5*dtdy
                    end if
                        gzo=(f(i,j,k+1)-f(i,j,k-1))/(2.*dz1)
                        gz(i,j,k) = gzn(i,j,k) &
                                  -(gxo*(u(i,j,k+1)-u(i,j,k-1))+gyo*(v(i,j,k+1)-v(i,j,k-1)) &
                                  +gzo*(w(i,j,k+1)-w(i,j,k-1)))*0.5*dtdz
                end do
            end do
        end do

    end subroutine dcip3d_v

    !***************************************************
    subroutine dcip3d_w(f,gx,gy,gz)
    !***************************************************

        real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
        real(8),dimension(-1:im,-1:jm,0:km)::fn,gxn,gyn,gzn,u,v,w

        fn=0.; gxn=0.; gyn=0.; gzn=0.; u=0.; v=0.; w=0.

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
                dx2 =dx1*dx1
                dx3 =dx2*dx1
                dtdx=dt/dx1
                dy1=(dyj(i,j)+dyj(i,j-1))*.5
                dy2 =dy1*dy1
                dy3 =dy2*dy1
                dtdy=dt/dy1
                ! do k=1,nz-1
                do k=1,nz
                    dz1=(dz(k)+dz(k+1))*.5
                    dtdz=dt/dz1
                    dz2 =dz1*dz1
                    dz3 =dz2*dz1

                    cx=u(i,j,k)
                    cy=v(i,j,k)
                    cz=w(i,j,k)
                    xx = - cx*dt
                    yy = - cy*dt
                    zz = - cz*dt

                    x1 = -sign(1.0,cx)
                    y1 = -sign(1.0,cy)
                    z1 = -sign(1.0,cz)
                    im1=i+int(x1)
                    jm1=j+int(y1)
                    km1=k+int(z1)

                    a01 = ( ( gx(i,j,k)+gx(im1,j,k) )*dx1*x1+2.0*( f(i,j,k)-f(im1,j,k) ) )/(dx3*x1)
                    a11 = ( 3.0*( f(im1,j,k)-f(i,j,k) )-( gx(im1,j,k)+2.*gx(i,j,k) )*dx1*x1 )/dx2
                    a02 = ( ( gy(i,j,k)+gy(i,jm1,k) )*dy1*y1+2.0*( f(i,j,k)-f(i,jm1,k) ) )/(dy3*y1)
                    a12 = ( 3.0*( f(i,jm1,k)-f(i,j,k) )-( gy(i,jm1,k)+2.*gy(i,j,k) )*dy1*y1 )/dy2
                    a03 = ( ( gz(i,j,k)+gz(i,j,km1) )*dz1*z1+2.0*( f(i,j,k)-f(i,j,km1) ) )/(dz3*z1)
                    a13 = ( 3.0*( f(i,j,km1)-f(i,j,k) )-( gz(i,j,km1)+2.*gz(i,j,k) )*dz1*z1 )/dz2
                    b1  = f(i,j,k)-f(im1,j,k)-f(i,jm1,k)+f(im1,jm1,k)
                    b2  = f(i,j,k)-f(i,jm1,k)-f(i,j,km1)+f(i,jm1,km1)
                    b3  = f(i,j,k)-f(im1,j,k)-f(i,j,km1)+f(im1,j,km1)
                    a05 = ( b1-(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy2)
                    a04 = ( b1-(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1)/(dx2*dy1*y1)
                    a14 = (-b1+(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1 &
                        +(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy1*y1)
                    a09 = ( b2-(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz2)
                    a08 = ( b2-(-gy(i,j,k)+gy(i,j,km1))*dy1*y1)/(dy2*dz1*z1)
                    a15 = (-b2+(-gy(i,j,k)+gy(i,j,km1))*dy1*y1 &
                        +(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz1*z1)
                    a06 = ( b3-(-gz(i,j,k)+gz(im1,j,k))*dz1*z1)/(dx1*x1*dz2)
                    a07 = ( b3-(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx2*dz1*z1)
                    a16 = (-b3+(-gz(i,j,k)+gz(im1,j,k))*dz1*z1 &
                        +(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx1*x1*dz1*z1)
                    a10 = ( -f(i,j,k) + (f(im1,j,k)+f(i,jm1,k)+f(i,j,km1)) &
                        - (f(im1,jm1,k)+f(i,jm1,km1)+f(im1,j,km1)) &
                        + f(im1,jm1,km1) ) / (dx1*x1*dy1*y1*dz1*z1)

                    ! if(k==nz.and.km1==nz+1) then
                    !     fn(i,j,k) = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gx(i,j,k))*xx &
                    !               +((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gy(i,j,k))*yy &
                    !               +a10*xx*yy*zz+f(i,j,k)
                    ! else
                        fn(i,j,k) = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gx(i,j,k))*xx &
                                  +((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gy(i,j,k))*yy &
                                  +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gz(i,j,k))*zz &
                                  +a10*xx*yy*zz+f(i,j,k)
                    ! end if

                    gxn(i,j,k)= (3.*a01*xx+2.*(a04*yy+a07*zz+a11))*xx &
                              +(a05*yy+a10*zz+a14)*yy+(a06*zz+a16)*zz+gx(i,j,k)
                    gyn(i,j,k)= (3.*a02*yy+2.*(a05*xx+a08*zz+a12))*yy &
                              +(a09*zz+a10*xx+a15)*zz+(a04*xx+a14)*xx+gy(i,j,k)
                    gzn(i,j,k)= (3.*a03*zz+2.*(a06*xx+a09*yy+a13))*zz &
                              +(a07*xx+a10*yy+a16)*xx+(a08*yy+a15)*yy+gz(i,j,k)
                end do
            end do
        end do

        do j=1,ny
            do i=1,nx
                dtdx=dt/(dxi(i,j)+dxi(i-1,j))*2.
                dtdy=dt/(dyj(i,j)+dyj(i,j-1))*2.
                do k=1,nz
                    dtdz=dt/(dz(k)+dz(k+1))*2.
                    gxo=(f(i+1,j,k)-f(i-1,j,k))/(2.*dx1)
                    gyo=(f(i,j+1,k)-f(i,j-1,k))/(2.*dy1)
                    gzo=(f(i,j,k+1)-f(i,j,k-1))/(2.*dz1)
                    f(i,j,k)  = fn(i,j,k)
                    gx(i,j,k) = gxn(i,j,k) &
                              -(gxo*(u(i+1,j,k)-u(i-1,j,k))+gyo*(v(i+1,j,k)-v(i-1,j,k)) &
                              +gzo*(w(i+1,j,k)-w(i-1,j,k)))*0.5*dtdx
                    gy(i,j,k) = gyn(i,j,k) &
                              -(gxo*(u(i,j+1,k)-u(i,j-1,k))+gyo*(v(i,j+1,k)-v(i,j-1,k)) &
                              +gzo*(w(i,j+1,k)-w(i,j-1,k)))*0.5*dtdy
                    if(k==nz) then
                        gz(i,j,k)=0.
                    else
                        gz(i,j,k) = gzn(i,j,k) &
                                  -(gxo*(u(i,j,k+1)-u(i,j,k-1))+gyo*(v(i,j,k+1)-v(i,j,k-1)) &
                                  +gzo*(w(i,j,k+1)-w(i,j,k-1)))*0.5*dtdz
                    end if
                end do
            end do
        end do

    end subroutine dcip3d_w

    !***************************************************
    subroutine dcip3d_c(f,gx,gy,gz)
    !***************************************************

        real(8),dimension(-1:im,-1:jm,0:km)::f,gx,gy,gz
        real(8),dimension(-1:im,-1:jm,0:km)::fn,gxn,gyn,gzn,u,v,w

        fn=0.; gxn=0.; gyn=0.; gzn=0.; u=0.; v=0.; w=0.

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
                dx2 =dx1*dx1
                dx3 =dx2*dx1
                dtdx=dt/dx1
                dy1=(dyj(i,j)+dyj(i,j-1))*.5
                dy2 =dy1*dy1
                dy3 =dy2*dy1
                dtdy=dt/dy1
                do k=1,nz
                    if(obst3d(i,j,k)==0) then
                        dz1=(dzk(k)+dzk(k-1))*.5
                        dtdz=dt/dz1
                        dz2 =dz1*dz1
                        dz3 =dz2*dz1

                        cx=u(i,j,k)
                        cy=v(i,j,k)
                        cz=w(i,j,k)
                        xx = - cx*dt
                        yy = - cy*dt
                        zz = - cz*dt

                        x1 = -sign(1.0,cx)
                        y1 = -sign(1.0,cy)
                        z1 = -sign(1.0,cz)
                        im1=i+int(x1)
                        if(j_west==1.and.i==1.and.im1==0) im1=1
                        if(j_east==1.and.i==nx.and.im1==nx+1) im1=nx
                        jm1=j+int(y1)
                        if(j_south==1.and.j==1.and.jm1==0) jm1=1
                        if(j_north==1.and.j==ny.and.jm1==ny+1) jm1=ny
                        km1=k+int(z1)
                        if(k==1.and.km1==0) km1=1
                        if(k==nz.and.km1==nz+1) km1=nz

                        a01 = ( ( gx(i,j,k)+gx(im1,j,k) )*dx1*x1+2.0*( f(i,j,k)-f(im1,j,k) ) )/(dx3*x1)
                        a11 = ( 3.0*( f(im1,j,k)-f(i,j,k) )-( gx(im1,j,k)+2.*gx(i,j,k) )*dx1*x1 )/dx2
                        a02 = ( ( gy(i,j,k)+gy(i,jm1,k) )*dy1*y1+2.0*( f(i,j,k)-f(i,jm1,k) ) )/(dy3*y1)
                        a12 = ( 3.0*( f(i,jm1,k)-f(i,j,k) )-( gy(i,jm1,k)+2.*gy(i,j,k) )*dy1*y1 )/dy2
                        a03 = ( ( gz(i,j,k)+gz(i,j,km1) )*dz1*z1+2.0*( f(i,j,k)-f(i,j,km1) ) )/(dz3*z1)
                        a13 = ( 3.0*( f(i,j,km1)-f(i,j,k) )-( gz(i,j,km1)+2.*gz(i,j,k) )*dz1*z1 )/dz2
                        b1  = f(i,j,k)-f(im1,j,k)-f(i,jm1,k)+f(im1,jm1,k)
                        b2  = f(i,j,k)-f(i,jm1,k)-f(i,j,km1)+f(i,jm1,km1)
                        b3  = f(i,j,k)-f(im1,j,k)-f(i,j,km1)+f(im1,j,km1)
                        a05 = ( b1-(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy2)
                        a04 = ( b1-(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1)/(dx2*dy1*y1)
                        a14 = (-b1+(-gx(i,j,k)+gx(i,jm1,k))*dx1*x1+(-gy(i,j,k)+gy(im1,j,k))*dy1*y1)/(dx1*x1*dy1*y1)
                        a09 = ( b2-(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz2)
                        a08 = ( b2-(-gy(i,j,k)+gy(i,j,km1))*dy1*y1)/(dy2*dz1*z1)
                        a15 = (-b2+(-gy(i,j,k)+gy(i,j,km1))*dy1*y1+(-gz(i,j,k)+gz(i,jm1,k))*dz1*z1)/(dy1*y1*dz1*z1)
                        a06 = ( b3-(-gz(i,j,k)+gz(im1,j,k))*dz1*z1)/(dx1*x1*dz2)
                        a07 = ( b3-(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx2*dz1*z1)
                        a16 = (-b3+(-gz(i,j,k)+gz(im1,j,k))*dz1*z1+(-gx(i,j,k)+gx(i,j,km1))*dx1*x1)/(dx1*x1*dz1*z1)
                        a10 = ( -f(i,j,k) + (f(im1,j,k)+f(i,jm1,k)+f(i,j,km1)) &
                            - (f(im1,jm1,k)+f(i,jm1,km1)+f(im1,j,km1)) &
                            + f(im1,jm1,km1) ) / (dx1*x1*dy1*y1*dz1*z1)

                        fn(i,j,k) = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gx(i,j,k))*xx &
                                  +((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gy(i,j,k))*yy &
                                  +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gz(i,j,k))*zz &
                                  +a10*xx*yy*zz+f(i,j,k)
                        gxn(i,j,k)= (3.*a01*xx+2.*(a04*yy+a07*zz+a11))*xx &
                                  +(a05*yy+a10*zz+a14)*yy+(a06*zz+a16)*zz+gx(i,j,k)
                        gyn(i,j,k)= (3.*a02*yy+2.*(a05*xx+a08*zz+a12))*yy &
                                  +(a09*zz+a10*xx+a15)*zz+(a04*xx+a14)*xx+gy(i,j,k)
                        gzn(i,j,k)= (3.*a03*zz+2.*(a06*xx+a09*yy+a13))*zz &
                                  +(a07*xx+a10*yy+a16)*xx+(a08*yy+a15)*yy+gz(i,j,k)
                    end if
                end do
            end do
        end do

        do j=1,ny
            do i=1,nx
                dtdx=dt/(dxi(i,j)+dxi(i-1,j))*2.
                dtdy=dt/(dyj(i,j)+dyj(i,j-1))*2.
                do k=1,nz
                    if(obst3d(i,j,k)==0) then
                        dtdz=dt/(dzk(k)+dzk(k-1))*2.
                        gxo=(f(i+1,j,k)-f(i-1,j,k))/(2.*dx1)
                        gyo=(f(i,j+1,k)-f(i,j-1,k))/(2.*dy1)
                        gzo=(f(i,j,k+1)-f(i,j,k-1))/(2.*dz1)
                        f(i,j,k)  = fn(i,j,k)
                        gx(i,j,k) = gxn(i,j,k) &
                                  -(gxo*(u(i+1,j,k)-u(i-1,j,k))+gyo*(v(i+1,j,k)-v(i-1,j,k)) &
                                  +gzo*(w(i+1,j,k)-w(i-1,j,k)))*0.5*dtdx
                        gy(i,j,k) = gyn(i,j,k) &
                                  -(gxo*(u(i,j+1,k)-u(i,j-1,k))+gyo*(v(i,j+1,k)-v(i,j-1,k)) &
                                  +gzo*(w(i,j+1,k)-w(i,j-1,k)))*0.5*dtdy
                        gz(i,j,k) = gzn(i,j,k) &
                                  -(gxo*(u(i,j,k+1)-u(i,j,k-1))+gyo*(v(i,j,k+1)-v(i,j,k-1)) &
                                  +gzo*(w(i,j,k+1)-w(i,j,k-1)))*0.5*dtdz
                    end if
                end do
            end do
        end do

    end subroutine dcip3d_c

end module cip3d_m
