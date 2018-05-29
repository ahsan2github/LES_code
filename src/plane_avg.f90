    subroutine plane_avg(XX,XX_avg)

    use globals
    implicit none
          
    integer*4 :: k
    real*8,dimension(:,:,:) :: XX
    real*8,dimension(:) :: XX_avg

    real*8,dimension(size(XX,3)) :: XXX
         
    IF(hprocs > 1)THEN
              
        do k=1,size(XX,3)
            XXX(k) = sum(XX(:,:,k))
        enddo
          
        call MPI_ALLREDUCE(XXX,XX_avg,size(XX,3),MPI_DOUBLE_PRECISION, &
        MPI_SUM,MPI_COMM_LEVEL,ierr)

        call MPI_BARRIER(MPI_COMM_LEVEL,ierr)

        XX_avg=XX_avg*inxny

    ELSE

        do k=1,size(XX,3)
            XX_avg(k) = sum(XX(:,:,k))
        enddo

        XX_avg=XX_avg*inxny
                 
    ENDIF

    return
    end subroutine plane_avg
