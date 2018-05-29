    subroutine ddy_upwind(dfdy,f,v)

    use globals
    implicit none

    interface
    include './interfaces/updateHorz.f90'
    end interface

    integer*4 :: i,j,k
    real*8,dimension(:,:,:) :: dfdy,f,v

! cc differencing scheme

    do k=2,nzb+1
        do j=2,nyb+1
            do i=1,nx

                if(v(i,j-1,k) >= 0)then !use j-1 point
                    dfdy(i,j,k)=(f(i,j,k)-f(i,j,k))*idy
                else !use j+1 point
                    dfdy(i,j,k)=(f(i,j+1,k)-f(i,j,k))*idy
                endif

            enddo
        enddo
    enddo

    return
    end subroutine ddy_upwind
