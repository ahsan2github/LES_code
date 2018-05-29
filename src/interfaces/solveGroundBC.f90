    subroutine solveGroundBC(Uref,scalarRef,scalar,u,v, &
    gndScalars,ustar,scalarFlux,soilHeatFlux,porosity, &
    satPotential,satHydrCond,soilExponent,heatCapSoil,albedo, &
    minAlbedo,measRad,netRad,ix,jy,Psi,Psi0,fi,fiH,zo, &
    obukhovLength,zGnd)
    integer*4 :: ix,jy
    real*8 :: Uref,ustar,soilHeatFlux,albedo,minAlbedo, &
    netRad,zo,obukhovLength
    real*8,dimension(:) :: scalarRef,scalarFlux,measRad,Psi,Psi0,fi, &
    fiH,zGnd,porosity,satPotential,satHydrCond,soilExponent, &
    heatCapSoil
    real*8,dimension(:,:) :: gndScalars
    real*8,dimension(:,:,:) :: u,v
    real*8,dimension(:,:,:,:) :: scalar
    end subroutine solveGroundBC
