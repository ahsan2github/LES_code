!     ... MPI UPDATE GHOSTLAYERS for one variable

!     ... TRICKS for MPI: SEND   .... TO   PROCESSOR
!     ...                 RECEIVE ... FROM PROCESSOR
!     ... For example, the following indicates the information
!     ... from the processor me==0 is sent to processor me+1!

    subroutine update3(u1,u2,u3)

    use globals
    implicit none

    real*8, dimension(:,:,:):: u1,u2,u3

    IF (vprocs > 1 .AND. vfact == 0) then
        call MPI_SEND(u1(1,1,nzb+1),size(u1,1)*size(u1,2), &
        MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
        call MPI_SEND(u2(1,1,nzb+1),size(u1,1)*size(u1,2), &
        MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
        call MPI_SEND(u3(1,1,nzb+1),size(u1,1)*size(u1,2), &
        MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )

        call MPI_RECV(u1(1,1,nzb+2),size(u1,1)*size(u1,2), &
        MPI_DOUBLE_PRECISION,me+hprocs,me,nall,MPI_STATUS_IGNORE, &
        ierr )
        call MPI_RECV(u2(1,1,nzb+2),size(u1,1)*size(u1,2), &
        MPI_DOUBLE_PRECISION,me+hprocs,me,nall,MPI_STATUS_IGNORE, &
        ierr )
        call MPI_RECV(u3(1,1,nzb+2),size(u1,1)*size(u1,2), &
        MPI_DOUBLE_PRECISION,me+hprocs,me,nall,MPI_STATUS_IGNORE, &
        ierr )

    ELSEIF (vfact <= vprocs-2) then

    !.... FOR EVEN LEVELS ....
        if(mod(vfact,2) == 0)then

        !.....mpi send upper layer
            call MPI_SEND(u1(1,1,nzb+1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
            call MPI_SEND(u2(1,1,nzb+1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
            call MPI_SEND(u3(1,1,nzb+1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
        !     mpi recieve upper layers
            call MPI_RECV(u1(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
        !     mpi send lower layer
            call MPI_SEND(u1(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u2(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u3(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
        !     mpi RECV lower ghostlayer
            call MPI_RECV(u1(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )

        else
        !.... FOR ODD LEVELS ....

        !     mpi RECV lower ghostlayer
            call MPI_RECV(u1(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
        !     mpi send lower layer
            call MPI_SEND(u1(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u2(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u3(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
        !     mpi recieve upper layers
            call MPI_RECV(u1(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
        !.....mpi send upper layer
            call MPI_SEND(u1(1,1,nzb+1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
            call MPI_SEND(u2(1,1,nzb+1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )
            call MPI_SEND(u3(1,1,nzb+1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me+hprocs,me+hprocs,nall,ierr )

        end if

    ELSE
    !     mpi	upper node
             
    !.... IF vprocs is even ....
        if(mod(vprocs,2) == 0)then
            call MPI_RECV(u1(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_SEND(u1(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u2(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u3(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
        else
            call MPI_SEND(u1(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u2(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_SEND(u3(1,1,2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me-hprocs,nall,ierr )
            call MPI_RECV(u1(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hprocs,me,nall, &
            MPI_STATUS_IGNORE,ierr )
        endif
                    
    ENDIF

!...  FOR PEROIDIC VERTICAL BOUNDARY CONDITION
!     top and bottom trade at the end
    IF(verticalBC == 1)THEN

        if(vfact == 0)then

            call MPI_SEND(u1(1,1,3),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me), &
            (nprocs-1)-(hprocs-1-me),nall,ierr )
            call MPI_SEND(u2(1,1,3),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me), &
            (nprocs-1)-(hprocs-1-me),nall,ierr )
            call MPI_SEND(u3(1,1,3),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me), &
            (nprocs-1)-(hprocs-1-me),nall,ierr )

            call MPI_RECV(u1(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me),me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me),me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,1),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me),me,nall, &
            MPI_STATUS_IGNORE,ierr )
              
        elseif(vfact == vprocs-1)then
              
            call MPI_RECV(u1(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u2(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me,nall, &
            MPI_STATUS_IGNORE,ierr )
            call MPI_RECV(u3(1,1,nzb+2),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me,nall, &
            MPI_STATUS_IGNORE,ierr )

            call MPI_SEND(u1(1,1,nzb),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me-hfact,nall,ierr )
            call MPI_SEND(u2(1,1,nzb),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me-hfact,nall,ierr )
            call MPI_SEND(u3(1,1,nzb),size(u1,1)*size(u1,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me-hfact,nall,ierr )

        endif

    ENDIF

    return
          
    end subroutine update3
