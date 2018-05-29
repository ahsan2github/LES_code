    Subroutine FILT_DA (F,DFDX,DFDY)

    use globals
    implicit none
    include 'fftw3.f90'

    interface
    include './interfaces/mpi_transpose_xy_cmplx.f90'
    include './interfaces/mpi_transpose_yx_cmplx.f90'
    end interface
               
    real*8,dimension(:,:,:) :: F,DFDX,DFDY

    integer*4 :: i,j,k,ii,jj
    integer*8 :: plan_fx,plan_fy,plan_DFDXx,plan_DFDXy,plan_DFDYx, &
    plan_DFDYy,plan_bFx,plan_bFy
          
    real*8,volatile,dimension(size(F,1)) ::  F_1Dxf,F_1D,DFDX_1D, &
    DFDY_1D
    double complex,volatile,dimension(size(F,1)/(2*hprocs)+ip, &
    size(F,2)*hprocs) ::  F_hat_2D_trans, DFDX_hat_2D_trans, &
    DFDY_hat_2D_trans
    double complex,volatile,dimension(size(F,1)/2+1,size(F,2)) &
    :: F_hat_2D,DFDX_hat_2D, DFDY_hat_2D
    double complex,volatile,dimension(size(F,2)*hprocs) :: &
    F_1Dy,f_hat_1Dy,DFDX_hat_1Dy, DFDY_hat_1Dy
    double complex,volatile,dimension(size(F,1)/2+1) :: &
    F_1Dxb,f_hat_1D,DFDX_1Dx, DFDY_1Dx

    if(t == 1)then            !Create Wisdom

        call  dfftw_plan_dft_r2c_1d(plan_fx,Nx,F_1Dxf,F_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_fy,Ny,F_1Dy,F_hat_1Dy, &
        FFTW_FORWARD,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_DFDXx,Nx,DFDX_1Dx, &
        DFDX_1D,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_DFDYx,Nx,DFDY_1Dx, &
        DFDY_1D,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_bFx,Nx,F_1Dxb, &
        F_1D,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_DFDXy,Ny,DFDX_hat_1Dy, &
        DFDX_hat_1Dy,FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_DFDYy,Ny,DFDY_hat_1Dy, &
        DFDY_hat_1Dy,FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bFy,Ny,F_hat_1Dy,F_hat_1Dy, &
        FFTW_BACKWARD,FFT_FLAG)

    else

        call  dfftw_plan_dft_r2c_1d(plan_fx,Nx,F_1Dxf,F_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_fy,Ny,F_1Dy,F_hat_1Dy, &
        FFTW_FORWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_DFDXx,Nx,DFDX_1Dx, &
        DFDX_1D,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_DFDYx,Nx,DFDY_1Dx, &
        DFDY_1D,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_bFx,Nx,F_1Dxb, &
        F_1D,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_DFDXy,Ny,DFDX_hat_1Dy, &
        DFDX_hat_1Dy,FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_DFDYy,Ny,DFDY_hat_1Dy, &
        DFDY_hat_1Dy,FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bFy,Ny,F_hat_1Dy,F_hat_1Dy, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)

    endif

!.................loop through z-layers
    do k=2,Nzb+1

    !.....loop through y-rows and execute FFTW in x-dir
        do j=1,Nyb

            F_1Dxf = F(:,j,k)

            call dfftw_execute(plan_fx)
        !^^^^^^^^^^^^^inputs F_1Dxf(F) and outputs F_hat_1D

            F_hat_2D(:,j)=F_hat_1D

        end do

    !.....perform transpose
        call mpi_transpose_xy_cmplx(F_hat_2D,F_hat_2D_trans)

    !.....loop through x-rows and execute FFTW in y-dir
        do i=1,Nxb/2+ip

            F_1Dy=F_hat_2D_trans(i,:)

            call dfftw_execute(plan_fy)
        !^^^^^^^^^^^inputs F_1Dy(F_hat_2D_trans) and outputs F_hat_1Dy

            F_hat_2D_trans(i,:) = F_hat_1Dy
        end do

    !.....Normalize the Fourier coefficients
        F_hat_2D_trans=F_hat_2D_trans*inxny

    !.....filter
        do j=1,Ny
            jj=j-1
            if(jj > nint(Ny/2.)) jj=jj-Ny
            jj=jj*l_r
                     

            do i=1,Nxb/2+ip
                ii=(i-1)+(me-hfact)*Nxb/2
                if(ii >= nint(Nx/(2.*fgr)) .OR. abs(jj) >= nint(Ny/(2.*fgr))) &
                then
                    F_hat_2D_trans(i,j) = 0.d0
                           
                endif

            !......differentiate along y
                             
                DFDY_hat_2D_trans(i,j)=dcmplx(dimag(F_hat_2D_trans(i,j))* &
                (-1.d0), dreal(F_hat_2D_trans(i,j)))*jj

            !......differentiate along x
                     
                DFDX_hat_2D_trans(i,j)=dcmplx(dimag(F_hat_2D_trans(i,j))* &
                (-1.d0),dreal(F_hat_2D_trans(i,j)))*ii

            end do
        end do

    !.....inverse transforms along y
        do i=1,Nxb/2+ip
            f_hat_1Dy = f_hat_2D_trans(i,:)
            DFDX_hat_1Dy = DFDX_hat_2D_trans(i,:)
            DFDY_hat_1Dy = DFDY_hat_2D_trans(i,:)
            call dfftw_execute(plan_bFy)
            call dfftw_execute(plan_DFDXy)
            call dfftw_execute(plan_DFDYy)
            f_hat_2D_trans(i,:) = f_hat_1Dy
            DFDX_hat_2D_trans(i,:) = DFDX_hat_1Dy
            DFDY_hat_2D_trans(i,:) = DFDY_hat_1Dy
        end do

    !.....perform inverse transposes
        call mpi_transpose_yx_cmplx(F_hat_2D_trans,f_hat_2D)
        call mpi_transpose_yx_cmplx(DFDX_hat_2D_trans,DFDX_hat_2D)
        call mpi_transpose_yx_cmplx(DFDY_hat_2D_trans,DFDY_hat_2D)

    !.....inverse transform along x
        do j=1,Nyb
            f_1Dxb = f_hat_2D(:,j)
            DFDX_1Dx = DFDX_hat_2D(:,j)
            DFDY_1Dx = DFDY_hat_2D(:,j)
            call dfftw_execute(plan_bFx)
            call dfftw_execute(plan_DFDXx)
            call dfftw_execute(plan_DFDYx)
                     
            F(:,j,k) = f_1D
            DFDX(:,j,k) = DFDX_1D
            DFDY(:,j,k) = DFDY_1D

        end do
    end do

!.....Destroy plans and exit
    call dfftw_destroy_plan(plan_fx)
    call dfftw_destroy_plan(plan_fy)
    call dfftw_destroy_plan(plan_DFDXx)
    call dfftw_destroy_plan(plan_DFDXy)
    call dfftw_destroy_plan(plan_DFDYx)
    call dfftw_destroy_plan(plan_DFDYy)
    call dfftw_destroy_plan(plan_bFx)
    call dfftw_destroy_plan(plan_bFy)
     
    return
    end Subroutine FILT_DA
         
          

            
          
