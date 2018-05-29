    subroutine ddx_upwind(dfdx,f,u)

    use globals
    implicit none

    integer*4 :: i,j,k
    real*8,dimension(:,:,:) :: dfdx,f,u

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

                if(u(i,j,k) >= 0)then !use i-1 point
                    if(i == 1)then
                        dfdx(i,j,k)= &
                        (f(i,j,k)-f(Nx-1,j,k))*idx
                    else
                        dfdx(i,j,k)= &
                        (f(i,j,k)-f(i-1,j,k))*idx
                    endif
                else !use i+1 point
                    if(i == Nx)then
                        dfdx(i,j,k)= &
                        (f(2,j,k)-f(i,j,k))*idx
                    else
                        dfdx(i,j,k)= &
                        (f(i+1,j,k)-f(i,j,k))*idx
                    endif
                endif

            enddo
        enddo
    enddo

    return
    end subroutine ddx_upwind
