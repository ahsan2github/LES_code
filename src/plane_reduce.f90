! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     PLANE_REDUCE indicially sums the XX arrays for processors on    c
!     the same horizontal level, and distributes the sum value back   c
!     out to each processor                                           c
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine plane_reduce(XX)

    use globals
    implicit none
          
    real*8 :: XX,XXX

    if(hprocs > 1)then

        call MPI_ALLREDUCE(XX,XXX,1,MPI_DOUBLE_PRECISION, &
        MPI_SUM,MPI_COMM_LEVEL,ierr)

        call MPI_BARRIER(MPI_COMM_LEVEL,ierr)

        XX=XXX

    endif

    return
    end subroutine plane_reduce
