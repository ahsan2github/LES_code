    subroutine galileanGroundShift(gndScalars,soilHeatFlux, &
    ustar,scalarFlux,coolrate,zo,Psi,Psi0,fi,fi_H,ilow,jlow)
    integer*4, dimension(:) :: ilow,jlow
    real*8, dimension(:,:) :: zo,ustar,soilHeatFlux
    real*8, dimension(:,:,:) :: scalarFlux,coolrate,Psi,Psi0,fi,fi_H
    real*8, dimension(:,:,:,:) :: gndScalars
    end subroutine galileanGroundShift
