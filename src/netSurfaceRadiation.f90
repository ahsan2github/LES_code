    subroutine netSurfaceRadiation(surfaceTemperature,albedo, &
    minAlbedo,porosity,refTemp,aq,temperature,surfaceMoisture, &
    measRad,netRad,ix,jy,flag)
! computes the net surface radiation
!  = incoming long - outgoing long + incoming short - outgoing short
! Incoming short wave radiation is all due to solar radiation and is computed by solarRadiation fcn
! Outgoing short is due to the reflection of the solar radiation
!     the amount of short wave radiation that is reflected depends on the
!     surface albedo.  This is included in the solarRadiation fcn since
!     the albedo depends on angle.
! Incoming long wave radiation is due to the radiation emitted by the atmosphere
!     and the environment and is currently neglected.
! Outgoing long wave radiation is due to the emittance of the surface.

! Outgoing long wave radiation emitted by the surface is computed
!     in this function using the Stefan-Boltzmann law and the surface emissivity
!     The result is added to the net short wave
!     radiation to obtain the net surface radiation.

! inputs:
!     UTC  - Coordinated Universal Time (UTC) of day (hours)
!     surfaceTemperature - temperature of surface (K)
! outputs:
!     netSurfaceRadiation - net radiation flux at the surface.
!            (negative is upward flux for all radiation terms)

    use globals
    use wallBoundaryConditions
    use scalars
    use SEBmodule
    implicit none

    interface
    include './interfaces/solarRadiation.f90'
    end interface

    integer*4 :: ix,jy,flag
    real*8 :: surfaceTemperature,refTemp,surfaceMoisture,netRad, &
    albedo,minAlbedo
    real*8,dimension(:) :: aq,temperature,measRad,porosity

    integer*4 :: i,j
    real*8 ::  longOut, longIn, shortIn, shortNet, longNet, solar, &
    precWater, totalPrecWater, z, &
    pressCorrection, deltaP, &
    press0,press1, vaporPress

!     if using measured radiation, interpolate and return
    if( radiationFlag > 0 )then

        i = mod(ttt-1,stepsPerRadVal)
        j = int( (ttt - i)/stepsPerRadVal ) + 1

        netRad = measRad(j) + i*( measRad(j+1) - measRad(j) ) &
        / stepsPerRadVal

        return
    endif
! precWater will be input from LES (precipitable water is the total amount of water in the atmos above when condensed)
!      totalPrecWater = 0
!      do i = 1,size(aq)
!         z = z_i*dz/2 + (i-1)*dz*z_i
!         if (i==1)then
!            press0 = pressureScale*100
!         else
!            press0 =100*pressureScale*(1 - lapseRate*(z-z_i*dz/2)
!     >           /(temperature(i-1)*scalarScales(1)))
!     >           **(1/(airGasConst/Cp_air))
!         endif
!         press1 =100*pressureScale*(1 - lapseRate*z
!     >        /(temperature(i)*scalarScales(1)))
!     >        **(1/(airGasConst/Cp_air))
!!         deltaP    = 100*pressureScale*(pressCorrection-1) ! in kPa
!         deltaP = press1 - press0
!         precWater = -(1.0/(g_hat*((u_star**2.0)/z_i)))*
!     >        aq(i)*deltaP
!         totalPrecWater = totalPrecWater + precWater ! 0 for a completely dry atmosphere
!      enddo
    totalPrecWater = 80.d0

! precipitable water computed using Prata (1996) (documented in Niemela 2001)
!      totalPrecWater = 0
!      do i = 1,size(aq)
!         vaporPress     = aq(i)*pressureScale/(0.622+aq(i))!same units as pressureScale(mb = hPa)
!         precWater = 465*vaporPress/(temperature(i)*scalarScales(1) ! vapor Pressure in hPa
!     >        /(1.0+0.61*aq(i)))
!         totalPrecWater = totalPrecWater + precWater
!      enddo

    longOut = emissivity * SB_Constant * surfaceTemperature**4
!      longOut = 0
! incoming longwave radiation model
!!!!!!!!!
!      Dilley and O'Brien (1998) model
    longIn =59.38d0+113.7d0*(refTemp/(273.16d0/scalarScales(1)))**6+ &
    &      96.96d0*sqrt(totalPrecWater/(25.d0)) ! *** COMPUTE precipitable water ****
    longIn  = longIn/(Cp_air*densityAir*scalarScales(1)*u_star) ! non-dimensionalize

!     Prata (1996) model
!      longIn  = (1-(1+totalPrecWater)*
!     >     exp(-(1.2+3.0*totalPrecWater)**(1.0/2.0)))
!     >     *SB_Constant*temperature(1)**4


    longNet = -0.04d0/(scalarScales(1)*u_star) ! NET longwave

!!!!!!!!!
!      rapid and accurate radiative transfer model (RRTM) (ref: Mlawer 1997)
!     ADD RRTM for better approximation of incoming long wave radiation from atmosphere

          
!      netSurfaceRadiation = (longIn - longOut) +
!     >     solarRadiation(UTC)

    call solarRadiation(shortNet,surfaceMoisture,ix,jy, &
    albedo,minAlbedo,porosity)
!      shortNet = (1-albedo)*shortIn
    netRad = shortNet + longNet

!      return
    end subroutine netSurfaceRadiation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine solarRadiation(shortNet,surfaceMoisture,ix,jy, &
    albedo,minAlbedo,porosity)
! solarRadiation approximates the incoming solar radiation
! based on the latitude, longitude, day of year and time of day
! neglects attenuation due to cloud cover (ie assumes clear sky)
! neglects attenuation due to moist air or air pollution
! assumes slope of ground surface is zero and
! earths rotation around the sun is circular

! inputs:
!     lat  - latitude of location to compute radiation (degrees)
!     long - longitude of location to compute radiation (degrees)
!     day  - julian day of year (1-365)
!     UTC  - Coordinated Universal Time (UTC) of day (hours)
! outputs:
!     radiation - solar short wave radiation flux at the
!                 earths surface given the above inputs (negative is upward, positive is downward).

! model reference: Stull, Roland B., An Introduction to Boundary Layer Meteorology.
!                         pg. 257

    use globals
    use wallBoundaryConditions
    use SEBmodule
    implicit none

    integer*4 :: ix,jy
    real*8 :: shortNet,surfaceMoisture,albedo,minAlbedo
    real*8,dimension(:) :: porosity

    real*8 :: Az, As, A, E, Z, declination, sinElevation, &
    transmissivity

    declination = 23.45d0*(pi/180.d0)*cos(2.d0*pi*(day-173)/365.25d0)
         
    sinElevation = sin(lat)*sin(declination) - &
    cos(lat)*cos(declination)* &
    cos( (2*pi*UTC/(24.d0*3600.d0/(z_i/u_star))) - long )


    if(sinElevation > 0)then
        transmissivity = (0.6d0 + 0.2d0*sinElevation)
    ! compute albedo based on material, moisture content and sinElevation
    !     ADD computation for albedo based on moisture and angle
        if( albedoFlag == 1 )then
            E  = asin(sinElevation)
            Z  = (pi/2.d0) - E
            Az = 0.01d0*(exp(0.003286d0*Z**1.5d0) - 1.d0)
            if( surfaceMoisture/porosity(1) <= 0.5)then
                As = albedo* &
                (1 - surfaceMoisture/porosity(1))
            else
                As = minAlbedo
            endif
            A  = As + Az
        ! or gables3 (from PILPS paper, online say albedo = 0.23 ??)
        !            A = albedo - minAlbedo*sinElevation
        else
            A = albedo
        endif

        shortNet = (1.0 - A)*solarIrradiance* &
        transmissivity*sinElevation
    else ! sinElevation <= 0
        shortNet = 0
    endif

    end subroutine solarRadiation
