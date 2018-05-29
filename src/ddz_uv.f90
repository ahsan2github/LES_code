    subroutine ddz_uv (DFDZ,F)

!  First deriv in z direction for boundary layer (2nd order numerics)
!...F is on UVP nodes and dFdz is on w nodes
    use globals
    implicit none
         
    integer*4 :: i,j,k
    real*8,dimension(:,:,:)::F,dfdz
!...2*pi*aspect is the nondimensional vertical depth of the domain (Aspect is in dimen.h)
!...Skip wall level, which will be done later with M.O.
    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                dfdz(i,j,k)=(f(i,j,k)-f(i,j,k-1))*idz
            end do
        end do
    end do

    if(vfact == 0 .AND. verticalBC == 0) then
        do j=1,nyb
            do i=1,nx
                dfdz(i,j,2)=0.
            end do
        end do
    end if
    return
    end subroutine ddz_uv
