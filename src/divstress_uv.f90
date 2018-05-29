    subroutine divstress_uv (divt, tx, ty, tz)

    use globals
    use scalars
    implicit none

    interface
    include './interfaces/ddx.f90'
    include './interfaces/ddy.f90'
    include './interfaces/ddz_w.f90'
    end interface

    integer*4 :: i,j,k
    real*8,dimension(:,:,:):: divt,tx,ty,tz
    real*8,dimension(size(tx,1),size(tx,2),size(tx,3))::tdx,tdy,tdz

!...  Compute stress gradients
    if(scalarCount >= 1)then
        call ddx(tdx,tx)
        call ddy(tdy,ty)
    else
        call ddx(tdx,tx)
        call ddy(tdy,ty)
    endif

    call ddz_w(tdz, tz)

    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                divt(i,j,k)=tdx(i,j,k)+ tdy(i,j,k)+ tdz(i,j,k)
            end do
        end do
    end do

    return
    end subroutine divstress_uv
