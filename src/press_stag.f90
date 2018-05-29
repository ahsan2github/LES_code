    Subroutine Press_stag(P, RHSx, RHSy, RHSz, RHSx_f, RHSy_f, &
    RHSz_f,u,v,w,DFDX,DFDY,divtz)
          
    use globals
    implicit none
    include 'fftw3.f90'

    interface
    include './interfaces/mpi_transpose_xy_cmplx.f90'
    include './interfaces/mpi_transpose_yx_cmplx.f90'
    include './interfaces/mpi_transpose_xz_cmplx.f90'
    include './interfaces/mpi_transpose_zx_cmplx.f90'
    include './interfaces/tridag.f90'
    include './interfaces/cyclic.f90'
    end interface

    integer*4 :: i,j,k,ii,jj,sz,kstart
    integer*8 :: plan_Hx_x,plan_Hy_x,plan_Hz_x,plan_Hx_y,plan_Hy_y, &
    plan_Hz_y,plan_tz_x,plan_tz_y,plan_bP_y,plan_bdpdx_y, &
    plan_bdpdy_y,plan_bfx,plan_DFDXx,plan_DFDYx
    real*8 :: alpha,beta
          
    real*8,dimension(:,:,:):: RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f, &
    u,v,w,dfdx,dfdy,divtz,P
    real*8,volatile,dimension(size(RHSx,1))::Hx_1Dxf,Hy_1Dxf,Hz_1Dxf
    double complex,volatile,dimension(size(P,1)/2+1)::Hx_hat_1D, &
    Hy_hat_1D,Hz_hat_1D
    double complex,volatile,dimension(size(P,1)/2+1,size(P,2)):: &
    Hx_hat_2D,Hy_hat_2D,Hz_hat_2D
    double complex,dimension(size(P,1)/(2*hprocs)+ip,size(P,2)*hprocs) &
    ::bottomw,topw
    double complex,volatile,dimension(size(P,1)/(2*hprocs)+ip, &
    size(P,2)*hprocs) :: Hx_hat_2D_htrans,Hy_hat_2D_htrans, &
    Hz_hat_2D_htrans,F_hat_2D_trans,DFDX_hat_2D_trans, &
    DFDY_hat_2D_trans
    double complex,volatile,dimension(size(P,2)*hprocs)::Hx_1Dy, &
    Hy_1Dy,Hz_1Dy
    double complex,volatile,dimension(size(P,1)/(2*hprocs)+ip, &
    size(P,2)*hprocs,size(P,3))::Hx_hat_3D_htrans, &
    Hy_hat_3D_htrans,Hz_hat_3D_htrans
    double complex,volatile,dimension(size(P,1)/(2*hprocs)+ip, &
    size(P,2)*hprocs/vprocs,(size(P,3)-2)*vprocs) :: &
    Hx_hat_3D_vtrans,Hy_hat_3D_vtrans,Hz_hat_3D_vtrans, &
    P_hat_vtrans
    real*8,dimension((size(P,3)-2)*vprocs+1)::RHS_col,a1,b1,c1, &
    p_colr,p_coli
    double complex,dimension(size(P,1)/(2*hprocs)+ip,size(P,2)*hprocs, &
    size(P,3)) :: p_hat
    double complex,volatile,dimension(size(P,1)/2+1)::f_1Dxb,DFDX_1Dx, &
    DFDY_1Dx
    real*8,volatile,dimension(size(P,1))::F_1D,DFDX_1D,DFDY_1D

    if(t == 1)then

        call dfftw_plan_dft_r2c_1d(plan_Hx_x,Nx,Hx_1Dxf,Hx_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_r2c_1d(plan_Hy_x,Nx,Hy_1Dxf,Hy_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_r2c_1d(plan_Hz_x,Nx,Hz_1Dxf,Hz_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_Hx_y,Ny,Hx_1Dy,Hx_1Dy,FFTW_FORWARD, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_Hy_y,Ny,Hy_1Dy,Hy_1Dy,FFTW_FORWARD, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_Hz_y,Ny,Hz_1Dy,Hz_1Dy,FFTW_FORWARD, &
        FFT_FLAG)
        call dfftw_plan_dft_r2c_1d(plan_tz_x,Nx,Hx_1Dxf,Hx_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_tz_y,Ny,Hx_1Dy,Hx_1Dy,FFTW_FORWARD, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bp_y,Ny,Hz_1Dy,Hz_1Dy, &
        FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bdpdx_y,Ny,Hx_1Dy,Hx_1Dy, &
        FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bdpdy_y,Ny,Hy_1Dy,Hy_1Dy, &
        FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_bfx,Nx,F_1Dxb,F_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_DFDXx,Nx,DFDX_1Dx,DFDX_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_DFDYx,Nx,DFDY_1Dx,DFDY_1D, &
        FFT_FLAG)

    else

        call dfftw_plan_dft_r2c_1d(plan_Hx_x,Nx,Hx_1Dxf,Hx_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_r2c_1d(plan_Hy_x,Nx,Hy_1Dxf,Hy_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_r2c_1d(plan_Hz_x,Nx,Hz_1Dxf,Hz_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_Hx_y,Ny,Hx_1Dy,Hx_1Dy,FFTW_FORWARD, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_Hy_y,Ny,Hy_1Dy,Hy_1Dy,FFTW_FORWARD, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_Hz_y,Ny,Hz_1Dy,Hz_1Dy,FFTW_FORWARD, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_r2c_1d(plan_tz_x,Nx,Hx_1Dxf,Hx_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_tz_y,Ny,Hx_1Dy,Hx_1Dy,FFTW_FORWARD, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bp_y,Ny,Hz_1Dy,Hz_1Dy, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bdpdx_y,Ny,Hx_1Dy,Hx_1Dy, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bdpdy_y,Ny,Hy_1Dy,Hy_1Dy, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_bfx,Nx,F_1Dxb,F_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_DFDXx,Nx,DFDX_1Dx,DFDX_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_DFDYx,Nx,DFDY_1Dx,DFDY_1D, &
        FFTW_WISDOM_ONLY)

    endif

    if(verticalBC == 0)then
        kstart=2
    elseif(verticalBC == 1)then
        kstart=2
    endif
     
    do k=2,Nzb+1
                 
                
    !.....Forward transform the planes of RHS matrices in x-dir

        do j=1,Nyb

            Hx_1Dxf = RHSx(:,j,k)-(1.d0/3.d0)*RHSx_f(:,j,k)+ &
            u(:,j,k)*(2.d0/3.d0)/DT
            Hy_1Dxf = RHSy(:,j,k)-(1.d0/3.d0)*RHSy_f(:,j,k)+ &
            v(:,j,k)*(2.d0/3.d0)/DT
            Hz_1Dxf = RHSz(:,j,k)-(1.d0/3.d0)*RHSz_f(:,j,k)+ &
            w(:,j,k)*(2.d0/3.d0)/DT

            call dfftw_execute(plan_Hx_x)
            call dfftw_execute(plan_Hy_x)
            call dfftw_execute(plan_Hz_x)

            Hx_hat_2D(:,j) = Hx_hat_1D
            Hy_hat_2D(:,j) = Hy_hat_1D
            Hz_hat_2D(:,j) = Hz_hat_1D

        end do

    !......Forward xy transpose

        call mpi_transpose_xy_cmplx(Hx_hat_2D,Hx_hat_2D_htrans)
        call mpi_transpose_xy_cmplx(Hy_hat_2D,Hy_hat_2D_htrans)
        call mpi_transpose_xy_cmplx(Hz_hat_2D,Hz_hat_2D_htrans)

    !......Forward transform the planes of RHS matrices in y-dir
        do i=1,Nxb/2+ip

            Hx_1Dy = Hx_hat_2D_htrans(i,:)
            Hy_1Dy = Hy_hat_2D_htrans(i,:)
            Hz_1Dy = Hz_hat_2D_htrans(i,:)

            call dfftw_execute(plan_Hx_y)
            call dfftw_execute(plan_Hy_y)
            call dfftw_execute(plan_Hz_y)

            Hx_hat_2D_htrans(i,:) = Hx_1Dy
            Hy_hat_2D_htrans(i,:) = Hy_1Dy
            Hz_hat_2D_htrans(i,:) = Hz_1Dy

        end do

    !......Normalize Fourier coefficients
        Hx_hat_3D_htrans(:,:,k) = Hx_hat_2D_htrans*inxny
        Hy_hat_3D_htrans(:,:,k) = Hy_hat_2D_htrans*inxny
        Hz_hat_3D_htrans(:,:,k) = Hz_hat_2D_htrans*inxny

    end do
                           
!.......Transform bottom boundary condition and distribute
    if(verticalBC == 0)then
        IF (vfact == 0) then   !Bottom BC
            do j=1,Nyb
                               
                Hx_1Dxf=divtz(:,j,2)
                call dfftw_execute(plan_tz_x)
                Hx_hat_2D(:,j) = Hx_hat_1D
                               
            end do
                        
            call mpi_transpose_xy_cmplx(Hx_hat_2D,Hx_hat_2D_htrans)
                        
            do i=1,Nxb/2+ip
                               
                Hx_1Dy = Hx_hat_2D_htrans(i,:)
                call dfftw_execute(plan_tz_y)
                Hx_hat_2D_htrans(i,:) = Hx_1Dy
                               
            end do
                        
            bottomw=Hx_hat_2D_htrans*inxny
                        
        END IF

    ! c mpi send bottomw

        IF(nprocs > 1 .AND. vprocs > 1)THEN
            if(vfact > 0)then
                               
                call MPI_RECV(bottomw(1,1),(Nxb/2+ip)*Ny/vprocs, &
                MPI_DOUBLE_COMPLEX,me-vfact*hprocs,me, &
                nall,MPI_STATUS_IGNORE,ierr )
                               
            else
                do j=1,vprocs-1
                                      
                    call MPI_SEND(bottomw(1,j*Ny/vprocs+1), &
                    (Nxb/2+ip)*Ny/vprocs,MPI_DOUBLE_COMPLEX, &
                    j*hprocs+me,j*hprocs+me,nall,ierr)
                                      
                end do
            end if
        ENDIF

    !......Transform top boundary condition and distribute
        IF (vfact == vprocs-1) then !Top BC

            do j=1,Nyb
                               
                Hx_1Dxf=divtz(:,j,nzb+1)
                call dfftw_execute(plan_tz_x)
                Hx_hat_2D(:,j) = Hx_hat_1D
                               
            end do
                        
            call mpi_transpose_xy_cmplx(Hx_hat_2D,Hx_hat_2D_htrans)
                        
            do i=1,Nxb/2+ip

                Hx_1Dy = Hx_hat_2D_htrans(i,:)
                call dfftw_execute(plan_tz_y)
                Hx_hat_2D_htrans(i,:) = Hx_1Dy

            end do
                        
            topw=Hx_hat_2D_htrans*inxny
                        
        END IF

    ! c mpi send topw

        IF(nprocs > 1 .AND. vprocs > 1)THEN
            if(vfact < vprocs-1)then
                              
                call MPI_RECV(topw(1,1),(Nxb/2+ip)*Ny/vprocs, &
                MPI_DOUBLE_COMPLEX,me+(vprocs-vfact-1)*hprocs,me, &
                nall,MPI_STATUS_IGNORE,ierr )
                              
            else
                do j=1,vprocs-1
                                     
                    call MPI_SEND(topw(1,(j-1)*Ny/vprocs+1),(Nxb/2+ip)*Ny &
                    /vprocs,MPI_DOUBLE_COMPLEX, &
                    (j-1)*hprocs+(me-hfact),(j-1)*hprocs+(me-hfact), &
                    nall,ierr)

                end do
            !.............copy the portion of topw that I need back to the beginning of the               array
                topw(:,1:Ny/vprocs) = topw(:,Ny-Ny/vprocs+1:Ny)
            end if
        ENDIF
    endif

!......Forward zx transpose

    call mpi_transpose_zx_cmplx(Hx_hat_3D_htrans,Hx_hat_3D_vtrans)
    call mpi_transpose_zx_cmplx(Hy_hat_3D_htrans,Hy_hat_3D_vtrans)
    call mpi_transpose_zx_cmplx(Hz_hat_3D_htrans,Hz_hat_3D_vtrans)

!.....Set wave numbers for non-continuous x-y domain

    do j=1,Ny/vprocs

        jj=(j-1)+vfact*Ny/vprocs
        if(jj >= (Ny/2)) jj=jj-Ny
        jj=jj*l_r

        do i=1,Nxb/2+ip

            ii=(i-1)+(me-hfact)*Nxb/2
                            
        !...  Zero wave number problem is non-unique for matrix solution
        !.... So just integrate up from wall with arbitrary P=const starting point.

            if ((ii == 0) .AND. (jj == 0)) then

                if(verticalBC == 0)then
                    p_hat_vtrans(i,j,1)= &
                    dcmplx(0.d0-DZ*dreal(bottomw(i,j)),0.d0)
                elseif(verticalBC == 1)then
                    p_hat_vtrans(i,j,1)=dcmplx(0.d0,0.d0)
                endif
                 
                do k=2,Nz
                    p_hat_vtrans(i,j,k)=dcmplx(dreal(p_hat_vtrans(i,j, &
                    k-1))+dreal(Hz_hat_3D_vtrans(i,j,k))*DZ,0.d0)
                end do

            else

            !...  Near Wall Nodes
                if(verticalBC == 0)then
                    a1(1)=0.d0
                    b1(1)=-1.d0
                    c1(1)=1.d0
                    RHS_col(1)=-1.d0*DZ*dreal(bottomw(i,j))
                endif
            !...  Interior nodes

                do k=kstart,Nz
                    RHS_col(k)=-ii*dimag(Hx_hat_3D_vtrans(i,j,k-1)) - &
                    jj*dimag(Hy_hat_3D_vtrans(i,j,k-1)) + &
                    (dreal(Hz_hat_3D_vtrans(i,j,k))- &
                    dreal(Hz_hat_3D_vtrans(i,j,k-1)))/DZ
                    a1(k)=1.d0/(DZ**2)
                    b1(k)=(-ii*ii-jj*jj-2.d0/(DZ**2))
                    c1(k)=1.d0/(DZ**2)
                end do

            !..   top nodes
                if(verticalBC == 0)then
                    RHS_col(Nz+1)=-1.d0*dreal(topw(i,j))*DZ
                    b1(Nz+1)=  1.d0
                    a1(nz+1)= -1.d0
                    c1(nz+1)=  0.d0
                endif

                if(verticalBC == 0)then
                    call tridag (a1,b1,c1,RHS_col,P_colr)
                elseif(verticalBC == 1)then
                    alpha=1.d0/(DZ**2)
                    beta=1.d0/(DZ**2)
                    call cyclic(a1(2:Nz),b1(2:Nz),c1(2:Nz),alpha, &
                    beta,RHS_col(2:Nz),P_colr(2:Nz))
                    P_colr(Nz+1)=P_colr(2)
                endif

            !..   Imag PART .............................................
            !...  Near Wall Nodes
                if(verticalBC == 0)then
                    a1(1)=  0.d0
                    b1(1)= -1.d0
                    c1(1)=  1.d0
                    RHS_col(1)=-1.d0*DZ*dimag(bottomw(i,j))
                endif
            !...  Interior Nodes
                do k=kstart,Nz
                    RHS_col(k)=ii*dreal(Hx_hat_3D_vtrans(i,j,k-1)) + &
                    jj*dreal(Hy_hat_3D_vtrans(i,j,k-1)) + &
                    (dimag(Hz_hat_3D_vtrans(i,j,k))- &
                    dimag(Hz_hat_3D_vtrans(i,j,k-1)))/DZ
                    a1(k)=1.d0/(DZ**2)
                    b1(k)=(-ii*ii-jj*jj-2.d0/(DZ**2))
                    c1(k)=1.d0/(DZ**2)
                end do
            !...  Top nodes
                if(verticalBC == 0)then
                    RHS_col(Nz+1)=-1.d0*dimag(topw(i,j)*DZ)
                    c1(Nz+1)=  0.d0
                    b1(Nz+1)=  1.d0
                    a1(Nz+1)= -1.d0
                endif

                if(verticalBC == 0)then
                    call tridag (a1,b1,c1,RHS_col,P_coli)
                elseif(verticalBC == 1)then
                    alpha=1.d0/(DZ**2)
                    beta=1.d0/(DZ**2)
                    call cyclic(a1(2:Nz),b1(2:Nz),c1(2:Nz),alpha, &
                    beta,RHS_col(2:Nz),P_coli(2:Nz))
                    P_coli(Nz+1)=P_coli(2)
                endif

            !...  Put Pressure amplitudes in Matrix
                do k=1,Nz
                    p_hat_vtrans(i,j,k)=dcmplx(p_colr(k+1),p_coli(k+1))
                !......Cut the Nyquist
                    if ((ii == Nx/2) .OR. (abs(jj) == l_r*Ny/2))then
                        p_hat_vtrans(i,j,k)=p_hat_vtrans(i,j,k)*0.0
                    end if
                end do

            end if
                           
        end do
    end do
         
!...  Inverse xz transpose P_hat

    call mpi_transpose_xz_cmplx(p_hat_vtrans,p_hat)

    do k=2,Nzb+1
    !.....Differentiate along  x and y directions

        Hz_hat_2D_htrans = p_hat(:,:,k)

        do j=1,Ny

            jj=j-1
            if (jj > (Ny/2)) 	jj=jj-Ny
            jj=jj*l_r

            do i=1,Nxb/2+ip
                	       
                ii=(i-1)+(me-hfact)*Nxb/2
                	       
                Hx_hat_2D_htrans(i,j)=dcmplx(dimag(Hz_hat_2D_htrans(i,j))* &
                (-1.d0),dreal(Hz_hat_2D_htrans(i,j)))*ii
                Hy_hat_2D_htrans(i,j)=dcmplx(dimag(Hz_hat_2D_htrans(i,j))* &
                (-1.d0),dreal(Hz_hat_2D_htrans(i,j)))*jj
                		
            end do
        end do
                 
    ! c Inverse Transform along y
        do i=1,Nxb/2+ip
            Hz_1Dy = Hz_hat_2D_htrans(i,:)
            Hx_1Dy = Hx_hat_2D_htrans(i,:)
            Hy_1Dy = Hy_hat_2D_htrans(i,:)
            call dfftw_execute(plan_bP_y)
            call dfftw_execute(plan_bdpdx_y)
            call dfftw_execute(plan_bdpdy_y)
            Hz_hat_2D_htrans(i,:) = Hz_1Dy
            Hx_hat_2D_htrans(i,:) = Hx_1Dy
            Hy_hat_2D_htrans(i,:) = Hy_1Dy

        end do

    !.....perform inverse transposes
        call mpi_transpose_yx_cmplx(Hz_hat_2D_htrans,Hz_hat_2D)
        call mpi_transpose_yx_cmplx(Hx_hat_2D_htrans,Hx_hat_2D)
        call mpi_transpose_yx_cmplx(Hy_hat_2D_htrans,Hy_hat_2D)

    !.....inverse transform along x
        do j=1,Nyb
            f_1Dxb = Hz_hat_2D(:,j)
            DFDX_1Dx = Hx_hat_2D(:,j)
            DFDY_1Dx = Hy_hat_2D(:,j)
            call dfftw_execute(plan_bFx)
            call dfftw_execute(plan_DFDXx)
            call dfftw_execute(plan_DFDYx)
            P(:,j,k) = f_1D
            DFDX(:,j,k) = DFDX_1D
            DFDY(:,j,k) = DFDY_1D

        end do
    end do

    call dfftw_destroy_plan(plan_Hx_x)
    call dfftw_destroy_plan(plan_Hy_x)
    call dfftw_destroy_plan(plan_Hz_x)
    call dfftw_destroy_plan(plan_Hx_y)
    call dfftw_destroy_plan(plan_Hy_y)
    call dfftw_destroy_plan(plan_Hz_y)
    call dfftw_destroy_plan(plan_tz_x)
    call dfftw_destroy_plan(plan_tz_y)
    call dfftw_destroy_plan(plan_bP_y)
    call dfftw_destroy_plan(plan_bdpdx_y)
    call dfftw_destroy_plan(plan_bdpdy_y)
    call dfftw_destroy_plan(plan_bFx)
    call dfftw_destroy_plan(plan_DFDXx)
    call dfftw_destroy_plan(plan_DFDYx)

    return
    end Subroutine Press_stag

