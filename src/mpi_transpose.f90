!     ... MPI transpose for one variable
! c strategy: processor me trades with modulo(j-me,hprocs) where
! c hprocs is the # of horizontal procs and j goes from 1->hprocs.
! c The higher me always recieves first in the paring.
! c Note: transposes for both real and complex arrays exist
! c usage: call mpi_transpose_xy(in,out) for real x-y transform
! c y-x simply replaces xy by yx in the subroutine name and for complex
! c add _cmplx after xy

! c forward from x-y transpose for complex numbers

    subroutine mpi_transpose_xy_cmplx(ux,uy)
    use globals
    implicit none
    integer*4 :: i,ii,j,k,nnxb,nnyb

    double complex,dimension(:,:):: ux
    double complex,dimension(:,:):: uy

    double complex,dimension((size(ux,1)-1)/hprocs,size(ux,2)) &
    ::usend,urecv
    double complex, dimension((size(ux,1)-1)/hprocs+1,size(ux,2)) &
    :: usendl,urecvl

    nnxb = (size(ux,1)-1)*2/hprocs
    nnyb = size(ux,2)

    do i=1,hprocs
                 
        ii=modulo(i-(me-hfact),hprocs)
         
        if(me-hfact == ii)then
                        
            uy(1:nnxb/2+ip,(ii)*nnyb+1:(ii+1)*nnyb) = &
            ux(ii*nnxb/2+1:(ii+1)*nnxb/2+ip,:)

        elseif(me-hfact > ii)then
                  
            if(me-hfact == hprocs-1)then

                call MPI_RECV(urecvl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,me,nall, &
                MPI_STATUS_IGNORE,ierr )
                              
                uy(1:nnxb/2+1,ii*nnyb+1:(ii+1)*nnyb) = urecvl(:,:)
                               
            else

                call MPI_RECV(urecv(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,me,nall, &
                MPI_STATUS_IGNORE,ierr )

                uy(1:nnxb/2,ii*nnyb+1:(ii+1)*nnyb) = urecv(:,:)

            endif

            if(ii == hprocs-1)then

                usendl(:,:) = ux(ii*nnxb/2+1:(ii+1)*nnxb/2+1,:)

                call MPI_SEND(usendl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,ii+hfact,nall,ierr )

            else

                usend(:,:) = ux(ii*nnxb/2+1:(ii+1)*nnxb/2,:)

                call MPI_SEND(usend(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,ii+hfact,nall,ierr )

            endif

        else

            if(ii == hprocs-1)then

                usendl(:,:) = ux(ii*nnxb/2+1:(ii+1)*nnxb/2+1,:)
                              
                call MPI_SEND(usendl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact, &
                ii+hfact,nall,ierr )

            else

                usend(:,:) = ux(ii*nnxb/2+1:(ii+1)*nnxb/2,:)

                call MPI_SEND(usend(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact, &
                ii+hfact,nall,ierr )
            endif

            if(me-hfact == hprocs-1)then

                call MPI_RECV(urecvl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact, &
                me,nall,MPI_STATUS_IGNORE,ierr )

                uy(1:nnxb/2+1,ii*nnyb+1:(ii+1)*nnyb) = urecvl(:,:)

            else
                            
                call MPI_RECV(urecv(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact, &
                me,nall,MPI_STATUS_IGNORE,ierr )
                            
                uy(1:nnxb/2,ii*nnyb+1:(ii+1)*nnyb) = urecv(:,:)

            endif

        endif
                 
    enddo

    return
          
    end subroutine mpi_transpose_xy_cmplx

! c going back from y-x for complex numbers ccc
    subroutine mpi_transpose_yx_cmplx(uy,ux)
    use globals
    implicit none
    integer*4 :: i,ii,j,k,nnxb,nnyb
        
    double complex,dimension(:,:):: uy
    double complex,dimension(:,:) :: ux

    double complex,dimension((size(ux,1)-1)/hprocs,size(ux,2)) &
    ::usend,urecv
    double complex, dimension((size(ux,1)-1)/hprocs+1,size(ux,2)) &
    :: usendl,urecvl

    nnxb = (size(ux,1)-1)*2/hprocs
    nnyb = size(ux,2)

    do i=1,hprocs
                 
        ii=modulo(i-(me-hfact),hprocs)

        if(me-hfact == ii)then

            ux(ii*nnxb/2+1:(ii+1)*nnxb/2+ip,:)= &
            uy(1:nnxb/2+ip,ii*nnyb+1:(ii+1)*nnyb)

        elseif(me-hfact > ii)then
                        
            if(ii == hprocs-1)then

                call MPI_RECV(urecvl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,me,nall, &
                MPI_STATUS_IGNORE,ierr )
                               
                ux(ii*nnxb/2+1:(ii+1)*nnxb/2+1,:) = urecvl(:,:)

            else

                call MPI_RECV(urecv(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,me,nall, &
                MPI_STATUS_IGNORE,ierr )

                ux(ii*nnxb/2+1:(ii+1)*nnxb/2,:) = urecv(:,:)

            endif

            if(me-hfact == hprocs-1)then

                usendl(:,:) = uy(:,ii*nnyb+1:(ii+1)*nnyb)

                call MPI_SEND(usendl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,ii+hfact,nall,ierr )

            else

                usend(:,:) = uy(:,ii*nnyb+1:(ii+1)*nnyb)

                call MPI_SEND(usend(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,ii+hfact,nall,ierr )

            endif

        else

            if(me-hfact == hprocs-1)then

                usendl(:,:) = uy(:,ii*nnyb+1:(ii+1)*nnyb)

                call MPI_SEND(usendl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,ii+hfact,nall,ierr )

            else

                usend(:,:) = uy(:,ii*nnyb+1:(ii+1)*nnyb)

                call MPI_SEND(usend(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,ii+hfact,nall,ierr )

            endif

            if(ii == hprocs-1)then
                               
                call MPI_RECV(urecvl(1,1),(nnxb/2+1)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,me,nall, &
                MPI_STATUS_IGNORE,ierr )
                               
                ux(ii*nnxb/2+1:(ii+1)*nnxb/2+1,:) = urecvl(:,:)

            else

                call MPI_RECV(urecv(1,1),(nnxb/2)*nnyb, &
                MPI_DOUBLE_COMPLEX,ii+hfact,me,nall, &
                MPI_STATUS_IGNORE,ierr )
                               
                ux(ii*nnxb/2+1:(ii+1)*nnxb/2,:) = urecv(:,:)

            endif

        endif
                 
    enddo
          
    return
          
    end subroutine mpi_transpose_yx_cmplx

! c forward from z-x transpose for complex numbers
    subroutine mpi_transpose_zx_cmplx(uxz,uxy)
    use globals
    implicit none

    integer*4 :: i,ii

    double complex,dimension(:,:,:) :: uxz
    double complex,dimension(:,:,:) :: uxy

    double complex, dimension(size(uxz,1),size(uxz,2)/vprocs, &
    size(uxz,3)-2) :: usend
    double complex, dimension(size(uxz,1),size(uxz,2)/vprocs, &
    size(uxz,3)-2) :: urecv

    do i=1,vprocs
                 
        ii=modulo(i-vfact,vprocs)

        if(me == ii*hprocs+(me-hfact))then
                        
            uxy(:,1:ny/vprocs,ii*nzb+1:(ii+1)*nzb) = &
            uxz(:,ii*ny/vprocs+1:(ii+1)*ny/vprocs,2:nzb+1)

        elseif(me > ii*hprocs+(me-hfact))then
                        
            usend(:,:,1:nzb) = uxz(:,ii*ny/vprocs+1:(ii+1)*ny/vprocs, &
            &               2:nzb+1)

            call MPI_RECV(urecv(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact),me, &
            nall,MPI_STATUS_IGNORE,ierr )

            uxy(:,:,ii*nzb+1:(ii+1)*nzb) = urecv(:,:,1:nzb)
                           
            call MPI_SEND(usend(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact), &
            ii*hprocs+(me-hfact),nall,ierr )

        else

            usend(:,:,1:nzb) = uxz(:,ii*ny/vprocs+1:(ii+1)*ny/vprocs, &
            &               2:nzb+1)

            call MPI_SEND(usend(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact), &
            ii*hprocs+(me-hfact),nall,ierr )

            call MPI_RECV(urecv(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact),me, &
            nall,MPI_STATUS_IGNORE,ierr )
                        
            uxy(:,:,ii*nzb+1:(ii+1)*nzb) = urecv(:,:,1:nzb)

        endif
                 
    enddo

    return
          
    end subroutine mpi_transpose_zx_cmplx

! c forward from x-z transpose for complex numbers
    subroutine mpi_transpose_xz_cmplx(uxy,uxz)
    use globals
    implicit none

    integer*4 :: i,ii

    double complex,dimension(:,:,:) :: uxz
    double complex,dimension(:,:,:) :: uxy

    double complex, dimension(size(uxz,1),size(uxz,2)/vprocs, &
    size(uxz,3)-2) :: usend
    double complex, dimension(size(uxz,1),size(uxz,2)/vprocs, &
    size(uxz,3)-2) :: urecv

    do i=1,vprocs
                 
        ii=modulo(i-vfact,vprocs)
                          
        if(me == ii*hprocs+(me-hfact))then
                        
            uxz(:,ii*ny/vprocs+1:(ii+1)*ny/vprocs,2:nzb+1) = &
            uxy(:,1:ny/vprocs,ii*nzb+1:(ii+1)*nzb)

        elseif(me > ii*hprocs+(me-hfact))then
                        
            usend(:,:,1:nzb) = uxy(:,1:ny/vprocs,ii*nzb+1:(ii+1)*nzb)

            call MPI_RECV(urecv(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact),me, &
            nall,MPI_STATUS_IGNORE,ierr )

            uxz(:,ii*ny/vprocs+1:(ii+1)*ny/vprocs,2:nzb+1) = &
            urecv(:,:,1:nzb)
                           
            call MPI_SEND(usend(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact), &
            ii*hprocs+(me-hfact),nall,ierr )

        else

            usend(:,:,1:nzb) = uxy(:,:,ii*nzb+1:(ii+1)*nzb)

            call MPI_SEND(usend(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact), &
            ii*hprocs+(me-hfact),nall,ierr )

            call MPI_RECV(urecv(1,1,1),(nxb/2+ip)*ny/vprocs*nzb, &
            MPI_DOUBLE_COMPLEX,ii*hprocs+(me-hfact), &
            me,nall,MPI_STATUS_IGNORE,ierr )
                        
            uxz(:,ii*ny/vprocs+1:(ii+1)*ny/vprocs,2:nzb+1) = &
            urecv(:,:,1:nzb)

        endif
                 
    enddo

    return
          
    end subroutine mpi_transpose_xz_cmplx


