    subroutine spanwise_avg(XX,XX_avg)
    use globals
    implicit none

    integer*4 :: j
    real*8,dimension(:,:,:) :: XX
    real*8,dimension(:,:) :: XX_avg

    real*8,dimension(size(XX,1),size(XX,3)) :: XXX

    if(hprocs > 1)then
                 
        XXX(:,:)=0.d0
        do j=1,nyb
            XXX(:,:)=XXX(:,:)+XX(:,j,:)
        enddo

        call MPI_ALLREDUCE(XXX,XX_avg,size(XXX),MPI_DOUBLE_PRECISION, &
        MPI_SUM,MPI_COMM_LEVEL,ierr)

        XX_avg=XX_avg/Ny

    else

        XXX(:,:)=0.d0
        do j=1,nyb
            XXX(:,:)=XXX(:,:)+XX(:,j,:)
        enddo

        XX_avg=XX_avg/Ny

    endif

    call MPI_BARRIER(MPI_COMM_COLUMN,ierr)

    return
    end subroutine spanwise_avg
