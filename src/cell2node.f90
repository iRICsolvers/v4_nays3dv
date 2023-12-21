module cell2node_m
  use common_hh
  use common_geom
  use common_hyd
  implicit none
contains
  !----------------------------------------------
  subroutine cell2node_3dout
    !----------------------------------------------
    real(8)::up, vp

    do i = 1, ni
      do j = 1, nj
        do k = 1, nk
          up = (yu(i - 1, j, k) + yu(i - 1, j - 1, k) &
                + yu(i - 1, j, k - 1) + yu(i - 1, j - 1, k - 1))*.25
          vp = (yv(i, j - 1, k) + yv(i - 1, j - 1, k) &
                + yv(i, j - 1, k - 1) + yv(i - 1, j - 1, k - 1))*.25
          u38(i, j, k) = up*coss(i - 1, j - 1) - vp*sins(i - 1, j - 1)
          v38(i, j, k) = up*sins(i - 1, j - 1) + vp*coss(i - 1, j - 1)
          w38(i, j, k) = (yw(i, j, k - 1) + yw(i - 1, j, k - 1) + yw(i, j - 1, k - 1) + yw(i - 1, j - 1, k - 1))*.25
          if (k == 1 .and. j == 1) q38(i, j, k) = q(i - 1)
        end do
      end do
    end do

    do i = 0, ni - 1
      do j = 0, nj - 1
        do k = 0, nk - 1
          if (i == ni - 1 .and. j == nj - 1) then
            c38(i + 1, j + 1, k + 1) = (yc(i, j, k) + yc(i, j, k + 1))*.5
          else if (i == ni - 1 .and. j == 0) then
            c38(i + 1, j + 1, k + 1) = (yc(i, j + 1, k) + yc(i, j + 1, k + 1))*.5
          else if (i == 0 .and. j == nj - 1) then
            c38(i + 1, j + 1, k + 1) = (yc(i + 1, j, k) + yc(i + 1, j, k + 1))*.5
          else if (i == 0 .and. j == 0) then
            c38(i + 1, j + 1, k + 1) = (yc(i + 1, j + 1, k) + yc(i + 1, j + 1, k + 1))*.5
          else if (i == ni - 1) then
            c38(i + 1, j + 1, k + 1) = (yc(i, j, k) + yc(i, j + 1, k) + yc(i, j, k + 1) &
                                        + yc(i, j + 1, k + 1))*.25
          else if (i == 0) then
            c38(i + 1, j + 1, k + 1) = (yc(i + 1, j, k) + yc(i + 1, j + 1, k) + yc(i + 1, j, k + 1) &
                                        + yc(i + 1, j + 1, k + 1))*.25
          else if (j == nj - 1) then
            c38(i + 1, j + 1, k + 1) = (yc(i, j, k) + yc(i + 1, j, k) + yc(i, j, k + 1) &
                                        + yc(i + 1, j, k + 1))*.25
          else if (j == 0) then
            c38(i + 1, j + 1, k + 1) = (yc(i, j + 1, k) + yc(i + 1, j + 1, k) + yc(i, j + 1, k + 1) &
                                        + yc(i + 1, j + 1, k + 1))*.25
          else
            c38(i + 1, j + 1, k + 1) = (yc(i, j, k) + yc(i + 1, j, k) + yc(i, j + 1, k) + yc(i, j, k + 1) &
                                        + yc(i + 1, j + 1, k) + yc(i, j + 1, k + 1) + yc(i + 1, j, k + 1) + yc(i + 1, j + 1, k + 1))*.125
          end if
          p38(i + 1, j + 1, k + 1) = (ypn(i, j, k) + ypn(i + 1, j, k) + ypn(i, j + 1, k) + ypn(i, j, k + 1) &
                                      + ypn(i + 1, j + 1, k) + ypn(i, j + 1, k + 1) + ypn(i + 1, j, k + 1) + ypn(i + 1, j + 1, k + 1))*.125
          snu_t38(i + 1, j + 1, k + 1) = (snu_t(i, j, k) + snu_t(i + 1, j, k) &
                                          + snu_t(i, j + 1, k) + snu_t(i, j, k + 1) + snu_t(i + 1, j + 1, k) &
                                          + snu_t(i, j + 1, k + 1) + snu_t(i + 1, j, k + 1) + snu_t(i + 1, j + 1, k + 1))*.125
          if (k == nk - 1) then
            ob38(i + 1, j + 1, k + 1) = dble(obst3d(i, j, k) + obst3d(i + 1, j, k) &
                                             + obst3d(i, j + 1, k) + obst3d(i + 1, j + 1, k))*.25
          else
            ob38(i + 1, j + 1, k + 1) = dble(obst3d(i, j, k) + obst3d(i + 1, j, k) &
                                             + obst3d(i, j + 1, k) + obst3d(i, j, k + 1) + obst3d(i + 1, j + 1, k) &
                                             + obst3d(i, j + 1, k + 1) + obst3d(i + 1, j, k + 1) + obst3d(i + 1, j + 1, k + 1))*.125
          end if
        end do
      end do
    end do

    ! do i=1,ni
    !     do j=1,nj
    !         do k=1,nk
    !             c38_cell(i,j,k)=yc(i,j,k)
    !             p38_cell(i,j,k)=yp(i,j,k)
    !         end do
    !     end do
    ! end do

    do i = 1, ni
      do j = 1, nj
        h8(i, j) = h_node(i - 1, j - 1)
        do k = 1, nk
          z38(i, j, k) = z8(i, j) + (h_node(i - 1, j - 1) - z8(i, j))*xi(k - 1)
          ! if((k==nk).and.(hobst_node8(i,j)>=z38(i,j,k))) then
          !     z38(i,j,k)=hobst_node8(i,j)
          ! end if
        end do
      end do
    end do

    ! Vorticity
    do i = 0, nx
      do j = 0, ny
        do k = 0, nz
          vol_x(i + 1, j + 1, k + 1) = (yw(i + 1, j + 1, k) + yw(i, j + 1, k) - yw(i + 1, j, k) - yw(i, j, k)) &
                                       /(dy(i, j + 1) + dy(i, j)) &
                                       - (yv(i + 1, j, k + 1) + yv(i, j, k + 1) - yv(i + 1, j, k) - yv(i, j, k)) &
                                       /(dzk(k)*hs_node(i, j))
          vol_y(i + 1, j + 1, k + 1) = (yu(i, j + 1, k + 1) + yu(i, j, k + 1) - yu(i, j + 1, k) - yu(i, j, k)) &
                                       /(dzk(k)*hs_node(i, j)) &
                                       - (yw(i + 1, j + 1, k) + yw(i + 1, j, k) - yw(i, j + 1, k) - yw(i, j, k)) &
                                       /(dx(i + 1, j) + dx(i, j))
          vol_z(i + 1, j + 1, k + 1) = (yv(i + 1, j, k + 1) + yv(i + 1, j, k) - yv(i, j, k + 1) - yv(i, j, k)) &
                                       /(dx(i + 1, j) + dx(i, j)) &
                                       - (yu(i, j + 1, k + 1) + yu(i, j + 1, k) - yu(i, j, k + 1) - yu(i, j, k)) &
                                       /(dy(i, j + 1) + dy(i, j))
        end do
      end do
    end do

  end subroutine cell2node_3dout
  !----------------------------------------------
  subroutine cell2node_2d(y_cell, y_node)
    real(8), dimension(0:im, 0:jm)::y_cell
    real(8), dimension(0:im, 0:jm)::y_node
    do i = 0, nx
      do j = 0, ny
        if (j == 0) then
          if (i == 0) then
            y_node(i, j) = y_cell(i + 1, j + 1)
          else if (i == nx) then
            y_node(i, j) = y_cell(i, j + 1)
          else
            y_node(i, j) = (y_cell(i, j + 1) + y_cell(i + 1, j + 1))*0.5d0
          end if
        else if (j == ny) then
          if (i == 0) then
            y_node(i, j) = y_cell(i + 1, j)
          else if (i == nx) then
            y_node(i, j) = y_cell(i, j)
          else
            y_node(i, j) = (y_cell(i, j) + y_cell(i + 1, j))*0.5d0
          end if
        else
          if (i == 0) then
            y_node(i, j) = (y_cell(i + 1, j) + y_cell(i + 1, j + 1))*0.5d0
          else if (i == nx) then
            y_node(i, j) = (y_cell(i, j) + y_cell(i, j + 1))*0.5d0
          else
            y_node(i, j) = (y_cell(i, j) + y_cell(i + 1, j) &
                            + y_cell(i, j + 1) + y_cell(i + 1, j + 1))*0.25d0
          end if
        end if
      end do
    end do
  end subroutine cell2node_2d
  !----------------------------------------------

end module cell2node_m
