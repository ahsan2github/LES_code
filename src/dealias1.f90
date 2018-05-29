    Subroutine Dealias1(u1,u1_m)

    use globals
    implicit none
    include 'fftw3.f90'

    interface
    include './interfaces/mpi_transpose_xy_cmplx.f90'
    include './interfaces/mpi_transpose_yx_cmplx.f90'
    end interface

    real*8, dimension(:,:,:) :: u1,u1_m
         
    integer*4 :: i,j,k
    integer*8 :: plan_ux,plan_uy,plan_bUx,plan_bUy

    real*8,volatile,dimension(size(u1,1)) :: u_1Dxf
    double complex,volatile,dimension(size(u1,2)*hprocs) :: u_1Dy
    double complex,volatile,dimension(size(u1,1)/2+1) :: u_hat_1D
    double complex,volatile,dimension(size(u1,1)/2+1,size(u1,2)) &
    :: u_hat_2D
    double complex,volatile,dimension(size(u1,1)/(2*hprocs)+ip, &
    size(u1,2)*hprocs):: u_hat_2D_trans
    double complex,volatile,dimension(size(u1,1)/(2*hprocs)+ip, &
    size(u1_m,2)*hprocs) :: u_hat_2D_trans_my
    double complex,volatile,dimension(size(u1_m,2)*hprocs) :: &
    u_hat_1Dy_m
    double complex,volatile,dimension(size(u1,1)/2+1,size(u1_m,2)) &
    :: u_hat_2D_my
    double complex,volatile,dimension(size(u1_m,1)/2+1,size(u1_m,2)) &
    :: u_hat_2D_m
    double complex,volatile,dimension(size(u1_m,1)/2+1) :: u_1Dxb_m
    real*8,volatile,dimension(size(u1_m,1)) :: u_1D_m

    if(t == 1)then

        call  dfftw_plan_dft_r2c_1d(plan_ux,Nx,u_1Dxf,u_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_uy,Ny,u_1Dy,u_1Dy,FFTW_FORWARD, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bUy,Ny2,u_hat_1Dy_m,u_hat_1Dy_m, &
        FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_bUx,Nx2,u_1Dxb_m,u_1D_m, &
        FFT_FLAG)

    else

        call  dfftw_plan_dft_r2c_1d(plan_ux,Nx,u_1Dxf,u_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_uy,Ny,u_1Dy,u_1Dy,FFTW_FORWARD, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bUy,Ny2,u_hat_1Dy_m,u_hat_1Dy_m, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_bUx,Nx2,u_1Dxb_m,u_1D_m, &
        FFTW_WISDOM_ONLY)

    endif
     
    do k=2,nzb+1
    !.......perform forward FFT in x-dir

        do j=1,Nyb

            u_1Dxf=u1(:,j,k)
            call dfftw_execute(plan_ux)
            u_hat_2D(:,j)=u_hat_1D

        enddo

    !.......perform transpose
        call mpi_transpose_xy_cmplx(u_hat_2D,u_hat_2D_trans)

    !.......perform FFT in y-dir
        do i=1,Nxb/2+ip

            u_1Dy=u_hat_2D_trans(i,:)
            call dfftw_execute(plan_uy)
            u_hat_2D_trans(i,:)=u_1Dy

        end do

    !.......normalize Fourier coefficients
        u_hat_2D_trans=u_hat_2D_trans*inxny

    ! ......add padding in y-dir

        u_hat_2D_trans_my(:,:)=0.d0

        u_hat_2D_trans_my(:,1:ny/2)=u_hat_2D_trans(:,1:ny/2)

        u_hat_2D_trans_my(:,ny+1:ny2)=u_hat_2D_trans(:,ny/2+1:ny)

    !.......perform inverse transforms along y
        do i=1,Nxb/2+ip
            u_hat_1Dy_m = u_hat_2D_trans_my(i,:)
            call dfftw_execute(plan_bUy)
            u_hat_2D_trans_my(i,:) = u_hat_1Dy_m
        end do

    !.......perform backward transpose

        call mpi_transpose_yx_cmplx(u_hat_2D_trans_my,u_hat_2D_my)

    !.......add padding in x-dir

        u_hat_2D_m(:,:) = 0.d0

        u_hat_2D_m(1:Nx/2+1,:) = u_hat_2D_my(:,:)

    !.......perform inverse transforms along x
        do j=1,Nyb2

            u_1Dxb_m = u_hat_2D_m(:,j)
            call dfftw_execute(plan_bUx)
            u1_m(:,j,k)=u_1D_m

        end do

    end do

    call dfftw_destroy_plan(plan_ux)
    call dfftw_destroy_plan(plan_uy)
    call dfftw_destroy_plan(plan_bUx)
    call dfftw_destroy_plan(plan_bUy)

    return

    end Subroutine Dealias1
