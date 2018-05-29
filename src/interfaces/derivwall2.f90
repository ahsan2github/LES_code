    subroutine derivwall2 (dtdz,dudz,dvdz,u,v,fi,fi_h,t_flux, &
    ustar,M)
    integer*2 :: N,l
    real*8,dimension(:,:,:,:)::dtdz
    real*8,dimension(:,:,:)::dudz,dvdz,u,v,fi,fi_h,t_flux
    real*8,dimension(:,:):: ustar,M
    end subroutine derivwall2
