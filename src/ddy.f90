    Subroutine DDY (DFDY,F)

    use globals
    implicit none
    include 'fftw3.f90'

    interface
    include './interfaces/mpi_transpose_xy_cmplx.f90'
    include './interfaces/mpi_transpose_yx_cmplx.f90'
    end interface

    real*8, dimension(:,:,:) :: F,DFDY

    integer*4 :: i,j,k,ii,jj
    integer*8 :: plan_fx,plan_fy,plan_DFDYx,plan_DFDYy

    real*8,volatile,dimension(size(F,1)) ::  F_1Dxf,DFDY_1D
    double complex,volatile,dimension(size(F,1)/(2*hprocs)+ip, &
    size(F,2)*hprocs) &
    ::  F_hat_2D_trans, DFDY_hat_2D_trans
    double complex,volatile,dimension(size(F,1)/2+1,size(F,2)) &
    :: F_hat_2D,DFDY_hat_2D
    double complex,volatile,dimension(size(F,2)*hprocs) :: &
    F_1Dy,f_hat_1Dy,DFDY_hat_1Dy
    double complex,volatile,dimension(size(F,1)/2+1) :: &
    f_hat_1D,DFDY_1Dx
          
    if(t == 1)then

        call  dfftw_plan_dft_r2c_1d(plan_fx,Nx,F_1Dxf,F_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_fy,Ny,F_1Dy,F_hat_1Dy, &
        FFTW_FORWARD,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_DFDYx,Nx,DFDY_1Dx, &
        DFDY_1D,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_DFDYy,Ny,DFDY_hat_1Dy, &
        DFDY_hat_1Dy,FFTW_BACKWARD,FFT_FLAG)

    else

        call  dfftw_plan_dft_r2c_1d(plan_fx,Nx,F_1Dxf,F_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_fy,Ny,F_1Dy,F_hat_1Dy, &
        FFTW_FORWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_DFDYx,Nx,DFDY_1Dx, &
        DFDY_1D,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_DFDYy,Ny,DFDY_hat_1Dy, &
        DFDY_hat_1Dy,FFTW_BACKWARD,FFTW_WISDOM_ONLY)

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

        do j=1,Ny
            jj=j-1
            if(jj > nint(Ny/2.)) jj=jj-Ny
            jj=jj*l_r
                     

            do i=1,Nxb/2+ip
                ii=(i-1)+(me-hfact)*Nxb/2

            !......filter and differentiate along y
                                    
                if ((ii == nint(Nx/2.)) .OR. (abs(jj) >= nint(Ny/2.)))then
                    DFDY_hat_2D_trans(i,j)= 0.d0
                else
                    DFDY_hat_2D_trans(i,j)=dcmplx(dimag(F_hat_2D_trans(i,j))* &
                    (-1.d0),dreal(F_hat_2D_trans(i,j)))*jj
                end if

            end do
        end do

    !.....inverse transforms along y
        do i=1,Nxb/2+ip
            DFDY_hat_1Dy = DFDY_hat_2D_trans(i,:)
            call dfftw_execute(plan_DFDYy)
            DFDY_hat_2D_trans(i,:) = DFDY_hat_1Dy
        end do

    !.....perform inverse transposes

        call mpi_transpose_yx_cmplx(DFDY_hat_2D_trans,DFDY_hat_2D)

    !.....inverse transform along x
        do j=1,Nyb

            DFDY_1Dx = DFDY_hat_2D(:,j)
            call dfftw_execute(plan_DFDYx)
            DFDY(:,j,k) = DFDY_1D

        end do
    end do

!.....Destroy plans and exit
    call dfftw_destroy_plan(plan_fx)
    call dfftw_destroy_plan(plan_fy)
    call dfftw_destroy_plan(plan_DFDYx)
    call dfftw_destroy_plan(plan_DFDYy)

    return
    end Subroutine DDY
