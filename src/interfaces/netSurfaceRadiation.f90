    subroutine netSurfaceRadiation(surfaceTemperature,albedo, &
    minAlbedo,porosity,refTemp,aq,temperature,surfaceMoisture, &
    measRad,netRad,ix,jy,flag)
    integer*4 :: ix,jy,flag
    real*8 :: surfaceTemperature,refTemp,surfaceMoisture,netRad, &
    albedo,minAlbedo
    real*8,dimension(:) :: aq,temperature,measRad,porosity
    end subroutine netSurfaceRadiation
