    subroutine rmsdiv(rms,du,dv,dw)

    use globals
    implicit none

    integer*4 :: i,j,k,kend
    real*8,dimension(:,:,:) :: du,dv,dw
    real*8 :: rms
          
    if (vfact == vprocs-1) then
        kend=nzb
    else
        kend=nzb+1
    end if
    rms=0.0

    do k=2,kend
        do j=1,nyb
            do i=1,ny
                rms = rms+abs(du(i,j,k)+dv(i,j,k)+dw(i,j,k))
            enddo
        enddo
    end do

    rms=rms/(Nx*Ny*Nz)

    return
    end subroutine rmsdiv
