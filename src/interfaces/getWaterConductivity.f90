    subroutine getWaterConductivity(moisture,diffCond,hydrCond,ix,jy, &
    porosity,satPotential,satHydrCond,soilExponent)
    integer*4 :: ix,jy
    real*8, dimension(:):: moisture,diffCond,hydrCond,porosity, &
    satPotential,satHydrCond,soilExponent
    end subroutine getWaterConductivity
