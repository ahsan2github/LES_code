!     ... MPI UPDATE horizontal sides for one variable

    subroutine updateHorz(u)

    use globals
    implicit none

    real*8, dimension(:,:,:):: u
    real*8,allocatable,dimension(:,:) :: usend,urecv
    real*8,allocatable,dimension(:,:,:) :: uu

    if(size(u,3) > 1)then
        allocate(uu(nx,nyb+2,nzb),usend(nx,nzb),urecv(nx,nzb))
        uu=u(:,:,2:nzb+1)
    else
        allocate(uu(nx,nyb+2,1),usend(nx,1),urecv(nx,1))
        uu=u
    endif

    IF (hprocs == 1) THEN

        u(:,1,:) = u(:,ny+1,:)
        u(:,ny+2,:) = u(:,2,:)

    ELSEIF (me-hfact == 0) THEN

        usend=uu(:,nyb+1,:)

        call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
        me+1,me+1,nall,ierr )

        call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
        me+1,me,nall,MPI_STATUS_IGNORE,ierr )
        uu(:,nyb+2,:)=urecv

        usend=uu(:,2,:)
        call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
        me+hprocs-1,me+hprocs-1,nall,ierr )
                
        call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
        me+hprocs-1,me,nall,MPI_STATUS_IGNORE,ierr )
        uu(:,1,:)=urecv
                 
    ELSEIF (me-hfact <= hprocs-2) then

    !.... FOR EVENS ....
        if(mod(me-hfact,2) == 0)then
                        
        !.....mpi send upper layer
            usend=uu(:,nyb+1,:)

            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me+1,me+1,nall,ierr )
        !     mpi recieve upper layers
            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me+1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,nyb+2,:)=urecv
        !     mpi send lower layer
            usend=uu(:,2,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me-1,nall,ierr )
        !     mpi RECV lower ghostlayer
            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,1,:)=urecv
                        
        else
        !.... FOR ODDS ....
                 
        !     mpi RECV lower ghostlayer
            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,1,:)=urecv
        !     mpi send lower layer
            usend=uu(:,2,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me-1,nall,ierr )
        !     mpi recieve upper layers
            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me+1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,nyb+2,:)=urecv
        !.....mpi send upper layer
            usend=uu(:,nyb+1,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me+1,me+1,nall,ierr )

        endif
            
    ELSE
                
    !     mpi	upper node
                 
    !.... IF hprocs is even ....
        if(mod(hprocs,2) == 0)then

            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,1,:)=urecv
            usend=uu(:,2,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me-1,nall,ierr)

            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-hprocs+1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,nyb+2,:)=urecv
            usend=uu(:,nyb+1,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-hprocs+1,me-hprocs+1,nall,ierr)
                        
        else

            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-hprocs+1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,nyb+2,:)=urecv
            usend=uu(:,nyb+1,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-hprocs+1,me-hprocs+1,nall,ierr)

            usend=uu(:,2,:)
            call MPI_SEND(usend(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me-1,nall,ierr)
            call MPI_RECV(urecv(1,1),nx*size(uu,3),MPI_DOUBLE_PRECISION, &
            me-1,me,nall,MPI_STATUS_IGNORE,ierr )
            uu(:,1,:)=urecv
                        
        endif
                 
    ENDIF

    if(size(u,3) > 1)then
        u(:,:,2:nzb+1)=uu
        deallocate(uu,usend,urecv)
    else
        u=uu
        deallocate(uu,usend,urecv)
    endif
     
    return
          
    end subroutine updateHorz
          
