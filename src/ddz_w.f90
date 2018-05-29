    subroutine ddz_w (DFDZ,F)

!...F is on w nodes and dFdz is on uvp nodes
    use globals
    implicit none
         
    integer*4 :: i,j,k
    real*8,dimension(:,:,:)::F,dfdz

!...2*pi*aspect is the nondimensional vertical depth of the domain
!...(Aspect is in dimen.h)
        
    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                dfdz(i,j,k)=(f(i,j,k+1)-f(i,j,k))*idz
            end do
        end do
    end do

!...  Stress free lid
    	
    if (vfact==vprocs-1 .AND. verticalBC == 0) then
        dfdz(:,:,Nzb+1)=0.d0
    end if

    return
    end subroutine ddz_w
