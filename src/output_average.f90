! c NOTE ave_x is expected to be of size anx,nz2 not nz ccc
! c This routine does a reduce to me=0 for final stat   ccc
! c output.                                             ccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine output_average(ave_x,fnum)

    use globals
    implicit none

    integer*4 :: fnum,ii,nn
          
    real*8,dimension(:,:) :: ave_x
    real*8,dimension(size(ave_x,1),size(ave_x,2))::bar
    real*8,dimension(size(ave_x,1),(size(ave_x,2)-2)*vprocs) :: x_bar

    nn=size(ave_x,1)

! cc  me-hfact=0 collects and sums values from its level
!      call MPI_ALLREDUCE(ave_x,bar,nn*nz2,MPI_DOUBLE_PRECISION,MPI_SUM,
!     +     MPI_COMM_LEVEL,ierr)
    call MPI_REDUCE(ave_x,bar,nn*nz2,MPI_DOUBLE_PRECISION,MPI_SUM, &
    &      0,MPI_COMM_LEVEL,ierr)

! cc me-hfact=0 sends to me=0 for output

    if(me == 0)then

        x_bar(:,1:nzb)=bar(:,2:nzb+1)

        do ii=1,vprocs-1
            call MPI_RECV(x_bar(1,nzb*ii+1),nn*nzb,MPI_DOUBLE_PRECISION, &
            me+ii*hprocs,me+ii*hprocs,nall,MPI_STATUS_IGNORE,ierr)

        enddo

        write(fnum)x_bar
        call flush(fnum)

    elseif(me-hfact == 0)then

        call MPI_SEND(bar(1,2),nn*nzb,MPI_DOUBLE_PRECISION,0,me,nall, &
        ierr)

    endif


    return
    end subroutine output_average

