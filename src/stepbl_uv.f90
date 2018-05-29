    subroutine STEPBL_uv (ui, RHSi, RHSi_f, force)
!...  For staggered grid where U and V are not forced to 0 at k=1

    use globals
    implicit none

    integer*4 :: i,j,k
    real*8,dimension(:,:,:)::ui,RHSi,RHSi_f
    real*8,dimension(:)::force

    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                ui(i,j,k)= ui(i,j,k)+DT*(1.5d0*RHSi(i,j,k)- &
                &               0.5d0*RHSi_f(i,j,k)+force(k))
            end do
        end do
    end do
          
!...No-stress top
          
    if(vfact == vprocs-1 .AND. verticalBC == 0) then
        do j=1,Nyb
            do i=1,Nx
                ui(i,j,Nzb+1)=ui(i,j,Nzb)
            end do
        end do
    elseif(verticalBC == 1)then
        if(vfact == 0)then
            call MPI_SEND(ui(1,1,2),size(ui,1)*size(ui,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me), &
            (nprocs-1)-(hprocs-1-me),nall,ierr )
        elseif(vfact == vprocs-1)then
            call MPI_RECV(ui(1,1,nzb+1),size(ui,1)*size(ui,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me,nall, &
            MPI_STATUS_IGNORE,ierr )
        endif
    end if
          
    return
    end subroutine STEPBL_uv
