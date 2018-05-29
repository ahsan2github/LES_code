    Subroutine FILTER_2Laa (F_hat,F_hatd,F)

    use globals
    implicit none
    include 'fftw3.f90'

    interface
    include './interfaces/mpi_transpose_xy_cmplx.f90'
    include './interfaces/mpi_transpose_yx_cmplx.f90'
    end interface

    real*8, dimension(:,:,:) :: F,F_hat,F_hatd
       
    integer*8 :: i,j,k,ii,jj
    integer*8 :: plan_fx,plan_fy,plan_bfx,plan_bfx_d,plan_bfy,plan_bfy_d

    real*8,volatile,dimension(size(F,1)) :: F_1Dxf,F_1D,F_1D_d
    double complex,volatile,dimension(size(F,2)*hprocs) :: F_1Dy, &
    F_hat_1Dy,f_hat_1Dy_d
    double complex,volatile,dimension(size(F,1)/2+1) :: F_hat_1D
    double complex,volatile,dimension(size(F,1)/2+1,size(F,2)) &
    :: F_hat_2D,F_hat_2D_d
    double complex,volatile,dimension(size(F,1)/(2*hprocs)+ip, &
    size(F,2)*hprocs) :: F_hat_2D_trans,F_hat_2D_trans_d
    double complex,volatile,dimension(size(F,1)/2+1) :: f_1Dxb, &
    f_1Dxb_d

    if(t == 1)then

        call  dfftw_plan_dft_r2c_1d(plan_fx,Nx,F_1Dxf,F_hat_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_fy,Ny,F_1Dy,F_1Dy,FFTW_FORWARD, &
        FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bfy,Ny,F_hat_1Dy,F_hat_1Dy, &
        FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_1d(plan_bfy_d,Ny,F_hat_1Dy_d, &
        F_hat_1Dy_d,FFTW_BACKWARD,FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_bfx,Nx,F_1Dxb,F_1D, &
        FFT_FLAG)
        call dfftw_plan_dft_c2r_1d(plan_bfx_d,Nx,F_1Dxb_d, F_1D_d, &
        FFT_FLAG)

    else

        call  dfftw_plan_dft_r2c_1d(plan_fx,Nx,F_1Dxf,F_hat_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_fy,Ny,F_1Dy,F_1Dy,FFTW_FORWARD, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bfy,Ny,F_hat_1Dy,F_hat_1Dy, &
        FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_1d(plan_bfy_d,Ny,F_hat_1Dy_d, &
        F_hat_1Dy_d,FFTW_BACKWARD,FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_bfx,Nx,F_1Dxb,F_1D, &
        FFTW_WISDOM_ONLY)
        call dfftw_plan_dft_c2r_1d(plan_bfx_d,Nx,F_1Dxb_d, F_1D_d, &
        FFTW_WISDOM_ONLY)
    endif

    do k=2,nzb+1

    !.......perform forward FFT in x-dir

        do j=1,Nyb

            F_1Dxf=F(:,j,k)
            call dfftw_execute(plan_fx)
            F_hat_2D(:,j)=F_hat_1D

        enddo

    !.......perform transpose
        call mpi_transpose_xy_cmplx(F_hat_2D,F_hat_2D_trans)

    !.......perform FFT in y-dir
        do i=1,Nxb/2+ip

            F_1Dy=F_hat_2D_trans(i,:)

            call dfftw_execute(plan_fy)

            F_hat_2D_trans(i,:)=F_1Dy
        end do

    !.......normalize Fourier coefficients
        F_hat_2D_trans=F_hat_2D_trans*inxny
                   
        F_hat_2D_trans_d=F_hat_2D_trans

        do j=1,Ny

            jj=j-1
            if(jj > nint(Ny/2.)) jj=jj-Ny
            jj=jj*l_r

            do i=1,Nxb/2+ip
                ii=(i-1)+(me-hfact)*Nxb/2

            !%%%%%%%%%%%%%%%%%% Circular filter %%%%%%%%%%%%%%%%%%%%%%%%%%
            !$$$               if(ii.ge.nint(nx/(2.0*fgr*tfr*tfr)).and. &
            !$$$                   abs(jj).ge.nint(l_r*Ny/(2.0*fgr*tfr*tfr)))then
            !$$$                  ii=1000
            !$$$                  jj=1000
            !$$$               end if
            !$$$               d=+sqrt((1.*ii/Nx)**2.+(1.*jj/Ny)**2.)
            !$$$               if (d.ge.(1.0/(2.0*fgr*tfr*tfr))) then
            !$$$                  F_hat_2D_trans_d(i,j)= &
            !$$$                 F_hat_2D_trans_d(i,j)*0.0
            !$$$               end if
            !$$$
            !$$$               if(ii.ge.nint(nx/(2.0*fgr*tfr)).and. &
            !$$$                   abs(jj).ge.nint(l_r*Ny/(2.0*fgr*tfr)))then
            !$$$                  ii=1000
            !$$$                  jj=1000
            !$$$               end if
            !$$$               d=+sqrt((1.*ii/Nx)**2.+(1.*jj/Ny)**2.)
            !$$$               if (d.ge.(1.0/(2.0*fgr*tfr))) then
            !$$$                  F_hat_2D_trans(i,j)=F_hat_2D_trans(i,j)*0.0
            !$$$               end if
            !%%%%%%%%%%%%%%%%%% Square filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (ii >= nint(Nx/(2.0*fgr*tfr)))then
                    F_hat_2D_trans(i,j)=F_hat_2D_trans(i,j)*0.0
                elseif(abs(jj) >= nint(l_r*Ny/(2.0*fgr*tfr)))then
                    F_hat_2D_trans(i,j)=F_hat_2D_trans(i,j)*0.0
                end if

                if (ii >= nint(Nx/(2.0*fgr*tfr*tfr)))then
                    F_hat_2D_trans_d(i,j)= &
                    F_hat_2D_trans_d(i,j)*0.0
                elseif(abs(jj) >= nint(l_r*Ny/(2.0*fgr*tfr*tfr)))then
                    F_hat_2D_trans_d(i,j)= &
                    F_hat_2D_trans_d(i,j)*0.0
                end if
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end do
        end do

    !.......perform inverse transforms along y
        do i=1,Nxb/2+ip

            f_hat_1Dy = f_hat_2D_trans(i,:)
            f_hat_1Dy_d = f_hat_2D_trans_d(i,:)
            call dfftw_execute(plan_bFy)
            call dfftw_execute(plan_bFy_d)
            f_hat_2D_trans(i,:) = f_hat_1Dy
            f_hat_2D_trans_d(i,:) = f_hat_1Dy_d

        end do

    !.......perform backward transpose
        call mpi_transpose_yx_cmplx(F_hat_2D_trans,f_hat_2D)
        call mpi_transpose_yx_cmplx(F_hat_2D_trans_d,F_hat_2D_d)

    !.......perform inverse transforms along x
        do j=1,Nyb

            f_1Dxb = f_hat_2D(:,j)
            f_1Dxb_d = F_hat_2D_d(:,j)
            call dfftw_execute(plan_bFx)
            call dfftw_execute(plan_bFx_d)
            F_hat(:,j,k)=f_1D
            F_hatd(:,j,k)=f_1D_d

        end do
    end do
     
    call dfftw_destroy_plan(plan_fx)
    call dfftw_destroy_plan(plan_fy)
    call dfftw_destroy_plan(plan_bfx)
    call dfftw_destroy_plan(plan_bfx_d)
    call dfftw_destroy_plan(plan_bfy)
    call dfftw_destroy_plan(plan_bfy_d)

    return

    end Subroutine FILTER_2Laa
