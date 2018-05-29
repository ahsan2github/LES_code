    subroutine spectrum(uu,kk,spec)

    use globals
    implicit none

    include 'fftw3.f90'

    integer*4 :: i,j,k,ii,iii,im,kk,kko
    integer*8 :: plan_fspec

    real*8,dimension(:,:,:) :: uu
    real*8,dimension(size(uu,1)):: vel_c
    real*8,dimension(size(uu,1)/2+1):: C_2,sumP,Pgm
    real*8,dimension(:,:):: Spec
    double complex, dimension(size(uu,1)/2+1) :: vel_hat
          
    if(t == p_count)then
        call dfftw_plan_dft_r2c_1d(plan_fspec,Nx,vel_c,vel_hat, &
        FFT_FLAG)
    else
        call dfftw_plan_dft_r2c_1d(plan_fspec,Nx,vel_c,vel_hat, &
        FFTW_WISDOM_ONLY)
    endif

    sumP=0.d0
          
    do j=1,Nyb

        vel_c=uu(:,j,kk)
        call dfftw_execute(plan_fspec)
                 
        do im=1,Nx/2+1
            C_2(im)=dreal(vel_hat(im))**2+dimag(vel_hat(im))**2
        end do
                 
        do ii=1,Nx/2+1
            Pgm(ii)=1./(Nx**2.)*(2.d0*C_2(ii))
            Pgm(1)=1./(Nx**2.)*C_2(1)
            Pgm(Nx/2+1)=1./(Nx**2.)*C_2(Nx/2+1)
            sumP(ii)=sumP(ii)+Pgm(ii)
        end do
                 
    end do

    IF(nprocs == 1)THEN
        kko=kk-1
    else
        kko=kk
    endif

    do iii=1,Nx/2+1
        spec(iii,kko)=sumP(iii)/(1.d0*Ny)
    end do

    call dfftw_destroy_plan(plan_fspec)

    return
    end subroutine spectrum


