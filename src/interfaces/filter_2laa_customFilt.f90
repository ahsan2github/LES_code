Subroutine FILTER_2Laa_customFilt (F_hat,F, customFilt)

use globals
implicit none
include 'fftw3.f90'

interface
include './interfaces/mpi_transpose_xy_cmplx.f90'
include './interfaces/mpi_transpose_yx_cmplx.f90'
end interface

real*8, dimension(:,:,:) :: F,F_hat
    
integer*8 :: i,j,k,ii,jj
integer*8 :: plan_fx,plan_fy,plan_bfx,plan_bfy
real(kind = 8) :: customFilt 
real*8,volatile,dimension(size(F,1)) :: F_1Dxf,F_1D
double complex,volatile,dimension(size(F,2)*hprocs) :: F_1Dy, &
F_hat_1Dy
double complex,volatile,dimension(size(F,1)/2+1) :: F_hat_1D
double complex,volatile,dimension(size(F,1)/2+1,size(F,2)) &
:: F_hat_2D
double complex,volatile,dimension(size(F,1)/(2*hprocs)+ip, &
size(F,2)*hprocs) :: F_hat_2D_trans
double complex,volatile,dimension(size(F,1)/2+1) :: f_1Dxb
end Subroutine FILTER_2Laa_customFilt