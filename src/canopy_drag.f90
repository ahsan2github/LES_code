    subroutine canopy_drag(Fdx,Fdy,Fdz,LAD_uvp,LAD_w,u,v,w,u_w,v_w, &
    mag_w,w_uvp,mag_uvp)

    use globals
    use canopyModule
    use wallBoundaryConditions
    implicit none

    integer*4 :: i,j,k
    real*8,dimension(:,:,:) :: Fdx,Fdy,Fdz,LAD_w,LAD_uvp,u,v,w,u_w, &
    v_w,mag_w,w_uvp,mag_uvp

!$$$cccc  Interpolate to find LAD on u,v  nodes
!$$$      if(t.eq.1)then
!$$$
!$$$         do k=2,Nzb+1
!$$$            do j=1,Nyb
!$$$               do i=1,Nx
!$$$                  if(LAD_w(i,j,k+1).eq.0.)then
!$$$                     LAD_uvp(i,j,k)=0.0
!$$$                  else
!$$$                     LAD_uvp(i,j,k) =
!$$$     +                    0.5d0*(LAD_w(i,j,k+1)+LAD_w(i,j,k))
!$$$                  endif
!$$$               enddo
!$$$            enddo
!$$$         enddo
!$$$
!$$$      endif

! cc  Interpolate to find LAD on w  nodes
    if(t == 1)then

        do k=2,Nzb+1
            do j=1,Nyb
                do i=1,Nx
                    if(LAD_uvp(i,j,k) == 0.)then
                        LAD_w(i,j,k)=0.0
                    else
                        LAD_w(i,j,k) = &
                        &                     0.5d0*(LAD_uvp(i,j,k)+LAD_uvp(i,j,k-1))
                    endif
                enddo
            enddo
        enddo

    endif

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

            ! cc           Interpolate to find u,v,w on uvp and w nodes
                if(k == nzb+1 .AND. vfact == vprocs)then
                    u_w(i,j,k) = 0.5d0*(u(i,j,k)+u(i,j,k-1))
                    v_w(i,j,k) = 0.5d0*(v(i,j,k)+v(i,j,k-1))
                    w_uvp(i,j,k) = w_uvp(i,j,k-1)
                elseif(k == 2 .AND. vfact == 0)then
                    u_w(i,j,k) = 0.d0
                    v_w(i,j,k) = 0.d0
                    w_uvp(i,j,k) = 0.5d0*(w(i,j,k+1)+w(i,j,k))
                else
                    u_w(i,j,k) = 0.5d0*(u(i,j,k)+u(i,j,k-1))
                    v_w(i,j,k) = 0.5d0*(v(i,j,k)+v(i,j,k-1))
                    w_uvp(i,j,k) = 0.5d0*(w(i,j,k+1)+w(i,j,k))
                endif

            ! cc           Velocity Magnitude on uvp and w nodes
                mag_uvp(i,j,k)=sqrt(u(i,j,k)**2+v(i,j,k)**2+ &
                w_uvp(i,j,k)**2)
                mag_w(i,j,k)=sqrt(u_w(i,j,k)**2+v_w(i,j,k)**2+ &
                w(i,j,k)**2)



            ! cc           Drag Force

                Fdx(i,j,k-1) = -Cd*LAD_uvp(i,j,k)*mag_uvp(i,j,k)* &
                u(i,j,k)
                Fdy(i,j,k-1) = -Cd*LAD_uvp(i,j,k)*mag_uvp(i,j,k)* &
                v(i,j,k)
                Fdz(i,j,k-1) = -Cd*LAD_w(i,j,k)*mag_w(i,j,k)*w(i,j,k)

            enddo
        enddo
    enddo

    return
    end subroutine canopy_drag


                   
