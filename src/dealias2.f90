    Subroutine Dealias2 (u2,u2_m)

    use globals
    implicit none
    include 'fftw3.f90'

    interface
    include './interfaces/mpi_transpose_xy_cmplx.f90'
    include './interfaces/mpi_transpose_yx_cmplx.f90'
    end interface

    real*8, dimension(:,:,:) :: u2,u2_m
    integer*4 :: i,j,k
    integer*8 :: plan_ux,plan_uy,plan_bUx,plan_bUy

    real*8,volatile,dimension(size(u2_m,1)) :: u_1Dxf_m
    double complex,volatile,dimension(size(u2_m,1)/2+1) :: u_hat_1D_m
    double complex,volatile,dimension(size(u2_m,1)/2+1,size(u2_m,2)) &
    :: u_hat_2D_m
    double complex,volatile,dimension(size(u2_m,1)/(2*hprocs)+ip, &
    size(u2_m,2)*hprocs) :: u_hat_2D_trans_m
    double complex,volatile,dimension(size(u2_m,1)/(2*hprocs)+ip, &
    size(u2,2)*hprocs) :: u_hat_2D_trans_mx
    double complex,volatile,dimension(size(u2_m,2)*hprocs) :: u_1Dy_m
    double complex,volatile,dimension(size(u2,1)/(2*hprocs)+ip, &
    size(u2,2)*hprocs):: u_hat_2D_trans
    double complex,volatile,dimension(size(u2,2)*hprocs) :: u_hat_1Dy
    double complex,volatile,dimension(size(u2,1)/2+1,size(u2,2)) &
    :: u_hat_2D
    double complex,volatile,dimension(size(u2_m,1)/2+1,size(u2,2)) &
    :: u_hat_2D_mx
    double complex,volatile,dimension(size(u2,1)/2+1) :: u_1Dxb
    real*8,volatile,dimension(size(u2,1)) :: u_1D
         
    if(t == 1)then

        call  dfftw_plan_dft_r2c_1d(plan_ux,Nx2,u_1Dxf_m, &
        u_hat_1D_m,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_uy,Ny2,u_1Dy_m,u_1Dy_m, &
        FFTW_FORWARD,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bUy,Ny,u_hat_1Dy,u_hat_1Dy, &
        FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_bUx,Nx,u_1Dxb,u_1D, &
        FFT_FLAG)

    else

        call  dfftw_plan_dft_r2c_1d(plan_ux,Nx2,u_1Dxf_m, &
        u_hat_1D_m,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_uy,Ny2,u_1Dy_m,u_1Dy_m, &
        FFTW_FORWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bUy,Ny,u_hat_1Dy,u_hat_1Dy, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_bUx,Nx,u_1Dxb,u_1D, &
        FFTW_WISDOM_ONLY)

    endif

    do k=2,nzb+1

    !.......perform forward FFT in x-dir

        do j=1,Nyb2

            u_1Dxf_m=u2_m(:,j,k)
            call dfftw_execute(plan_ux)
            u_hat_2D_m(:,j)=u_hat_1D_m

        enddo

    !.......perform transpose

        call mpi_transpose_xy_cmplx(u_hat_2D_m,u_hat_2D_trans_m)

    !.......perform FFT in y-dir
        do i=1,Nxb2/2+ip

            u_1Dy_m=u_hat_2D_trans_m(i,:)

            call dfftw_execute(plan_uy)

            u_hat_2D_trans_m(i,:)=u_1Dy_m
        end do

    !.......normalize Fourier coefficients

        u_hat_2D_trans_m=u_hat_2D_trans_m*inx2ny2

    !.......remove padding in y-dir

        u_hat_2D_trans_mx(:,:)=0.d0

        u_hat_2D_trans_mx(:,1:ny/2)=u_hat_2D_trans_m(:,1:ny/2)

        u_hat_2D_trans_mx(:,ny/2+1:ny)=u_hat_2D_trans_m(:,ny+1:ny2)
              
    !.......perform inverse transforms along y
        do i=1,Nxb2/2+ip

            u_hat_1Dy = u_hat_2D_trans_mx(i,:)
            call dfftw_execute(plan_bUy)
            u_hat_2D_trans_mx(i,:) = u_hat_1Dy

        end do

    !.......perform backward transpose

        call mpi_transpose_yx_cmplx(u_hat_2D_trans_mx,u_hat_2D_mx)

    !.......remove padding in x-dir

        u_hat_2D(:,:) = u_hat_2D_mx(1:Nx/2+1,:)

    !.......perform inverse transforms along x
        do j=1,Nyb

            u_1Dxb = u_hat_2D(:,j)
            call dfftw_execute(plan_bUx)
            u2(:,j,k)=u_1D

        end do

    end do

    call dfftw_destroy_plan(plan_ux)
    call dfftw_destroy_plan(plan_uy)
    call dfftw_destroy_plan(plan_bUx)
    call dfftw_destroy_plan(plan_bUy)

    return

    end Subroutine Dealias2
