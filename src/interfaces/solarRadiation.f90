    subroutine solarRadiation(shortNet,surfaceMoisture,ix,jy,albedo, &
    minAlbedo,porosity)
    integer*4 :: ix,iy
    real*8 :: shortNet,surfaceMoisture,albedo,minAlbedo
    real*8,dimension(:) :: porosity
    end subroutine solarRadiation
