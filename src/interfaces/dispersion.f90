    subroutine dispersion(particle,u,v,w,ESGS,zo,ustar,LAD_w, &
    us_f,ESGS_f,ESGS_part_f,tau_f,idum,fx,fz,txx,txy,txz, &
    tyy,tyz,tzz,divtx,divty,divtz,TL,rogue_count, &
    trajectory_updates,rogue_flag)
    integer*4 :: idum,rogue_count,trajectory_updates
    integer*4,dimension(:) :: rogue_flag
    real*8,dimension(:) :: fx,fz,ESGS_part_f
    real*8,dimension(:,:) :: zo,ustar
    real*8,dimension(:,:,:) :: particle,u,v,w,ESGS,LAD_w, &
    txx,txy,txz,tyy,tyz,tzz,divtx,divty,divtz,us_f,TL, &
    ESGS_f
    real*8,dimension(:,:,:,:) :: tau_f
    end subroutine dispersion
