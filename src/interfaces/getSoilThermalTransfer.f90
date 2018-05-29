    subroutine getSoilThermalTransfer(moisture,thermalTransfer,flag, &
    ix,jy,porosity,satPotential,soilExponent,heatCapSoil)
    integer*4 :: flag,ix,jy
    real*8, dimension(:):: moisture,thermalTransfer,porosity, &
    satPotential,soilExponent,heatCapSoil
    end subroutine getSoilThermalTransfer
