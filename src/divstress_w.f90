    subroutine divstress_w (divt, tx, ty, tdz)

    use globals
    implicit none

    interface
    include './interfaces/ddx.f90'
    include './interfaces/ddy.f90'
    include './interfaces/filt_da.f90'
    end interface

    integer*4 :: i,j,k
    real*8,dimension(:,:,:):: divt,tx,ty,tdz
    real*8,dimension(size(tx,1),size(tx,2),size(tx,3))::tdx,tdy
          
!...  Compute stress gradients
    call ddx(tdx,tx)
    call ddy(tdy,ty)

!...  Note ddz_uv does not return values for k=1 level.
    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                divt(i,j,k)=tdx(i,j,k)+ tdy(i,j,k)+ tdz(i,j,k)
            end do
        end do

    !...  At wall we have to assume that tdz(tzz)=0.0.  Any better ideas?
        if (vfact == 0 .AND. k == 2 .AND. verticalBC == 0)then
            do j=1,Nyb
                do i=1,Nx
                    divt(i,j,k)=tdx(i,j,k)+ tdy(i,j,k)
                end do
            end do
        end if
    end do
    return
    end subroutine divstress_w
















