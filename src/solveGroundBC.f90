    subroutine solveGroundBC(Uref,scalarRef,scalar,u,v, &
    gndScalars,ustar,scalarFlux,soilHeatFlux,porosity, &
    satPotential,satHydrCond,soilExponent,heatCapSoil,albedo, &
    minAlbedo,measRad,netRad,ix,jy,Psi,Psi0,fi,fiH,zo, &
    obukhovLength,zGnd)
     
    use globals
    use wallBoundaryConditions
    use scalars
    use SEBmodule
    implicit none

    interface
    include './interfaces/getStabilityCorrections.f90'
    include './interfaces/getSurfaceMixingRatio.f90'
    include './interfaces/integrateSoilDiffusion.f90'
    include './interfaces/getSoilThermalTransfer.f90'
    include './interfaces/getWaterConductivity.f90'
    include './interfaces/netSurfaceRadiation.f90'
    end interface

    integer*4 :: ix,jy
    real*8 :: Uref,ustar,soilHeatFlux,albedo,minAlbedo, &
    netRad,zo,obukhovLength
    real*8,dimension(:) :: scalarRef,scalarFlux,measRad,Psi,Psi0,fi, &
    fiH,zGnd,porosity,satPotential,satHydrCond,soilExponent, &
    heatCapSoil
    real*8,dimension(:,:) :: gndScalars
    real*8,dimension(:,:,:) :: u,v
    real*8,dimension(:,:,:,:) :: scalar

    integer*4 :: bisectFlag,iterateFlux,iterateTemp,posCnt, &
    negCnt,i,k,tempConvFlag ,sepFlag,computeLH, &
    skipSEBFlag,skipIntegrationFlag
    real*8, dimension(size(scalar,3)) :: aq,at
    real*8, dimension(size(Psi))::PsiH, PsiH0,moistPotential, &
    conductivity,diffusConduct,hydrConduct,lastSurfScalars, &
    lastScalarFlux
    real*8 :: denom,denomH,x, y,z,deta_dz,h,partialPressure,qs,q_gnd, &
    specHum,specHum_gnd,shortIn,shortNet,longNet,longIn,longOut, &
    soilMoistureFlux,mFluxConvergence,TsConvergence, &
    tFluxConvergence,measResidual,measPress,lastTemperature, &
    lastSoilMoistureFlux,SEB,SEBlast,dSEB_dT,lastSoilFlux, &
    lastTempFlux,lowT,highT,dTs_Old,dTs,dFluxM_dT, dFluxT_dT
    real*8, dimension(size(gndScalars,1)) :: tempT, tempM
           
    tempConvFlag = 0
    bisectFlag   = 0
    sepFlag      = 1
    computeLH    = 1

    lastSurfScalars(:) = gndScalars(1,:)
    lastTemperature    = gndScalars(1,temperatureIndex)

    if ( UTC == (startUTC + dt) .OR. t==1 )then ! first timeStep
        lowT    = gndScalars(1,temperatureIndex) - 5/scalarScales(1)  ! set bounds on temperature
        highT   = gndScalars(1,temperatureIndex) + 5/scalarScales(1)  ! for convergence
    else
    !     if solution converged 1st time step, it shouldn't change much there after
        lowT    = gndScalars(1,temperatureIndex) - 0.25/scalarScales(1) ! set bounds on temperature
        highT   = gndScalars(1,temperatureIndex) + 0.25/scalarScales(1) ! for convergence
    endif
    dTs_Old = highT - lowT
    dTs     = dTs_Old

    skipSEBFlag = 0
    if( ttt > endConstSEB .AND. &
    mod(ttt-endConstSEB,updateFreqSEB) /= 0 )then
        skipSEBFlag = 1
    endif

    if(skipSEBFlag == 1 .AND. ix == 1 .AND. jy == 1 .AND. me == 0)then
        write(*,*)'skipping SEB, t=',t
    endif
          
    if ( UTC == (startUTC + dt) .OR. t==1 .OR. skipSEBFlag == 1 )then

        Psi=0.d0
        Psi0=0.d0
        PsiH=0.d0
        PsiH0=0.d0

        do i=1,4
            denom = dlog( (dz/2.d0) / zo ) + Psi(1) - Psi0(1)
            ustar = Uref*vonk/denom
                        
        ! scalar flux
            denomH = dlog( (dz/2.d0) / (zt/z_i) ) + PsiH(1) - PsiH0(1)
                               
            scalarFlux(temperatureIndex) = &
            ( gndScalars( 1,temperatureIndex ) &
            - scalarRef(temperatureIndex) ) * ustar*vonk/denomH
                        
            call getSurfaceMixingRatio(gndScalars,q_gnd,measPress, &
            ix,jy,porosity,satPotential,soilExponent)
            if(computeLH == 1)then
                scalarFlux(moistureIndex) = ( q_gnd &
                - scalarRef(moistureIndex) )*ustar*vonk/denomH

            endif

        ! compute Psi and fi values for momentum and scalars from computed flux

            call getStabilityCorrections(ustar, &
            scalarRef(temperatureIndex), &
            scalarRef(moistureIndex),scalarFlux, &
            Psi,Psi0,fi,PsiH,PsiH0,fiH)
                          
        enddo
    endif

! solve ground budget iteratively for temperature and moisture
!     with solution of PsiH and PsiH0, compute flux of other scalars (if any).
    if(skipSEBFlag == 0)then
        do iterateFlux = 1,maxFluxIterations
            if( computeLH==1 ) then
            ! assume surface fluxes are constant and converge to surface temperature
            ! using Newton-Raphson, where SEB = f(Ts)
                do iterateTemp = 1, maxTempIterations

                ! SURFACE ENERGY BUDGET
                ! solve surface energy budget

                    call netSurfaceRadiation(gndScalars(1,temperatureIndex), &
                    albedo,minAlbedo,porosity,scalarRef(temperatureIndex), &
                    aq,at,gndScalars(1,moistureIndex),measRad, &
                    netRad, ix,jy, iterateFlux*iterateTemp)

                    call getSoilThermalTransfer(gndScalars(1:2,moistureIndex), &
                    conductivity,0,ix,jy,porosity,satPotential, &
                    soilExponent,heatCapSoil)

                    soilHeatFlux = ( gndScalars(1,temperatureIndex) - &
                    gndScalars(2,temperatureIndex) )* &
                    (sum(conductivity)/2.d0)/ (zGnd(2) - zGnd(1))


                    SEB = soilHeatFlux + scalarFlux(temperatureIndex) &
                    + scalarFlux(moistureIndex)*latentHeatWater &
                    - netRad !- shortNet - longNet
                !            netRad = shortNet + longNet

                ! compute dSEB_dT
                    dFluxT_dT = ustar*vonk/denomH

                    dSEB_dT = 4.d0*emissivity*SB_Constant* &
                    gndScalars(1,temperatureIndex)**3 &
                    + (sum(conductivity)/2)/zGnd(2) + dFluxT_dT

                    dTs = SEB / dSEB_dT

                    gndScalars(1,temperatureIndex) = &
                    gndScalars(1,temperatureIndex) - dTs

                    TsConvergence = abs(dTs)/gndScalars(1,temperatureIndex)

                    if ( TsConvergence < temperatureCriteria )then
                        exit
                    endif
                                
                    SEBlast = SEB
                    lastScalarFlux = scalarFlux
                    lastSoilFlux   = soilHeatFlux
                enddo ! iterateTemp
                      
            else
                         
                call getSoilThermalTransfer(gndScalars(1:2,moistureIndex), &
                conductivity,0,ix,jy,porosity,satPotential,soilExponent, &
                heatCapSoil)

                soilHeatFlux = ( gndScalars(1,temperatureIndex) - &
                gndScalars(2,temperatureIndex) )*(sum(conductivity)/2.d0) &
                / (zGnd(2) - zGnd(1))
            endif


            lastTempFlux = scalarFlux(temperatureIndex)
                                 
            scalarFlux(temperatureIndex) = ( gndScalars( 1,temperatureIndex ) &
            - scalarRef(temperatureIndex) ) * ustar*vonk/denomH
            tFluxConvergence = abs( scalarFlux(temperatureIndex)-lastTempFlux)

            if( sepFlag == 1 )then
                mFluxConvergence = 0
            endif

            if( ((mFluxConvergence < moistureCriteria) .AND. &
            (tFluxConvergence < tempFluxCriteria)) .OR. &
            (iterateFlux >= maxFluxIterations) )then
                exit
            endif

            gndScalars(1,temperatureIndex) = (lastTemperature+ &
            gndScalars(1,temperatureIndex))/2.d0

            lastTemperature = gndScalars(1,temperatureIndex)

            do i=1,4

                denom = dlog( (dz/2.d0) / zo ) + Psi(1) - Psi0(1)
                ustar = Uref*vonk/denom
                         
            ! scalar flux
                denomH = dlog( (dz/2.d0) / (zt/z_i) ) + PsiH(1) - PsiH0(1)
                         
                scalarFlux(temperatureIndex) = ( gndScalars(1,temperatureIndex) &
                - scalarRef(temperatureIndex) ) * ustar*vonk/denomH
                         
                call getSurfaceMixingRatio(gndScalars,q_gnd,measPress, &
                ix,jy,porosity,satPotential,soilExponent)
                if(computeLH == 1 )then
                    scalarFlux(moistureIndex) = ( q_gnd &
                    - scalarRef(moistureIndex) )*ustar*vonk/denomH
                endif
                         
            !     compute Psi and fi values for momentum and scalars from computed flux

                call getStabilityCorrections(ustar, &
                scalarRef(temperatureIndex),scalarRef(moistureIndex), &
                scalarFlux,Psi,Psi0,fi,PsiH,PsiH0,fiH)
                           
                        
            enddo

        enddo                     ! iterate=1,maxFluxIterations
         
    !    compute first soil moisture flux
        call getWaterConductivity(gndScalars(1:2,moistureIndex), &
        diffusConduct,hydrConduct,ix,jy,porosity,satPotential, &
        satHydrCond,soilExponent)

        soilMoistureFlux = densityWater*(sum(diffusConduct)/2.d0) &
        *(gndScalars(1,moistureIndex)-gndScalars(2,moistureIndex)) &
        /(zGnd(2)-zGnd(1)) + densityWater*sum(hydrConduct)/2.d0
             
        if( computeLH == 1 )then
            if( sepFlag == 1 )then
                do i = 1,200!maxTempIterations

                    call getSurfaceMixingRatio(gndScalars,q_gnd,measPress, &
                    ix,jy,porosity,satPotential,soilExponent)
                    if(computeLH == 1 )then
                        scalarFlux(moistureIndex) = ( q_gnd &
                        - scalarRef(moistureIndex) )*ustar*vonk/denomH
                    endif

                    call getWaterConductivity(gndScalars(1:2,moistureIndex), &
                    diffusConduct,hydrConduct,ix,jy,porosity,satPotential, &
                    satHydrCond,soilExponent)
                                
                    lastSoilMoistureFlux = soilMoistureFlux
                    soilMoistureFlux = convFactor*lastSoilMoistureFlux - &
                    (1.d0-convFactor)*scalarFlux(moistureIndex)
                                
                    moistPotential(2) = satPotential(2) &
                    *(porosity(2) &
                    /gndScalars(2,moistureIndex))**soilExponent(2)
                                
                    moistPotential(1) = moistPotential(2) + (zGnd(2)-zGnd(1)) &
                    *((soilMoistureFlux/(densityWater*sum(hydrConduct) &
                    /2.d0)) - 1.d0)

                    gndScalars(1,moistureIndex) = porosity(1)* &
                    (moistPotential(1)/satPotential(1)) &
                    **(-1.d0/soilExponent(1))

                    mFluxConvergence = abs( &
                    (scalarFlux(moistureIndex) &
                    - (-soilMoistureFlux)) &
                    /(scalarFlux(moistureIndex)))

                    if( mFluxConvergence < moistureCriteria )then
                        exit
                    endif
                enddo
            endif
        endif
                 
    endif ! if( skipSEBFlag == 0 )

    skipIntegrationFlag = 0
    if( ttt > endConstSEB .AND. &
    mod(ttt-endConstSEB,integrateSoilDiffFreq) /= 0 )then
        skipIntegrationFlag = 1
    endif

!  use surface fluxes and soil profiles to integrate
!    heat diffusion equation and soil moisture equation in time
    if( skipIntegrationFlag == 0 )then
        if (temperatureIndex /= 0)then
            call integrateSoilDiffusion(gndScalars, &
            lastSurfScalars,1,ix,jy,zGnd,porosity,satPotential, &
            soilExponent,heatCapSoil,satHydrCond)
        endif
                 
        if (moistureIndex /= 0)then
            call integrateSoilDiffusion(gndScalars, &
            lastSurfScalars,2,ix,jy,zGnd,porosity,satPotential, &
            soilExponent,heatCapSoil,satHydrCond)
        endif
    endif
    end subroutine solveGroundBC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine getStabilityCorrections(ustar,atmTemp,mixingRatio, &
    scalarFlux,Psi,Psi0,fi,PsiH,PsiH0,fiH)
    use globals
    use wallBoundaryConditions
    use SEBmodule
    implicit none

    real*8 :: ustar,atmTemp,mixingRatio
    real*8, dimension(:) :: scalarFlux,Psi,Psi0,fi,PsiH,PsiH0,fiH

    real*8 :: tempFlux,obukhovLength,x,y,z
    integer*4 :: negCnt, posCnt, k

    posCnt = 0
    negCnt = 0

    tempFlux = atmTemp*0.61d0*scalarFlux(2) &
    + scalarFlux(1)*(1.d0+0.61d0*mixingRatio)

    obukhovLength = -ustar**3*atmTemp / &
    ( vonk*g_hat*tempFlux )
          
    do k = 1,2
        z = k*0.5d0*dz
        if ( z/obukhovLength > 5.0 )then
            obukhovLength = z/5.0
        elseif( z/obukhovLength < -5.0)then
            obukhovLength = -z/5.0
        endif
        if( tempFlux > 0.0 )then
            if(k == 1)posCnt = posCnt+1
            x=(1.d0-(15.d0*z/obukhovLength))**0.25d0
            y=(1.d0-(15.d0*(zt/z_i)/obukhovLength))**0.25d0 ! zo changed to zt
            Psi(k) = -2.d0*dlog(0.5d0*(1+x)) - &
            dlog(0.5d0*(1.d0+x**2.d0))+2.d0*atan(x)-pi/2.d0
            Psi0(k) = -2.d0*dlog(0.5d0*(1+y)) - &
            dlog(0.5d0*(1.d0+y**2.d0)) + 2.d0*atan(y) - pi/2.d0
            fi(k) = 1.d0/x
        !     eqtn 11.9 and 11.14 from Arya
            PsiH(k)  = -2.d0*dlog(0.5d0*(1.d0+x**2.d0))
            PsiH0(k) = -2.d0*dlog(0.5d0*(1+y**2.d0))
            fiH(k)   = fi(k)**2.d0
        elseif( tempFlux < -0.0 )then
            if(k == 1)negCnt = negCnt + 1
        !     Psi(k)   = 4.7*z/obukhovLength
        ! for GABLS 5*z/L usually 4.7*z/L
            Psi(k)   = 5.d0*z/obukhovLength
            Psi0(k)  = 5.d0*(zt/z_i)/obukhovLength !zo changed to zt
            fi(k)    = 1.d0 + Psi(k)
            PsiH(k)  = 5.d0*z/obukhovLength
            PsiH0(k) = 5.d0*(zt/z_i)/obukhovLength ! zo changed to zt
            fiH(k)   = 0.74d0 + PsiH(k)
        else
            Psi(k)   = 0.d0
            Psi0(k)  = 0.d0
            fi(k)    = 1.d0
            PsiH(k)  = Psi(k)
            PsiH0(k) = Psi0(k)
            fiH(k)   = 0.74d0
        endif
    enddo
          
    end subroutine getStabilityCorrections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine getSurfaceMixingRatio(gndScalars,q_gnd,measPress, &
    ix,jy,porosity,satPotential,soilExponent)
    use globals
    use scalars
    use SEBmodule
    implicit none

    integer*4 :: ix,jy
    real*8 :: q_gnd,measPress
    real*8,dimension(:) :: porosity,satPotential,soilExponent
    real*8,dimension(:,:) :: gndScalars

    real*8 :: moistPotential(2), h, partialPressure, satHum, specHum_gnd

!     from McCumber (documented in Pielke, Mesoscale Meteorological Modeling (page 420)
    moistPotential(1) = satPotential(1)* &
    (porosity(1)/gndScalars(1,moistureIndex))**soilExponent(1)

    h = exp(g_hat*moistPotential(1) &
    / (waterGasConst*gndScalars(1,temperatureIndex)))
    partialPressure = 6.1078d0*exp(17.269d0* &
    (gndScalars(1,temperatureIndex) - 273.16d0/scalarScales(1)) &
    / (gndScalars(1,temperatureIndex) - 35.86d0/scalarScales(1)))

    satHum = 0.622d0*( partialPressure / &
    ( pressureScale - 0.378d0*partialPressure) )
    specHum_gnd = h*satHum/scalarScales(2)

!     convert from specific humidity to mixing ratio
    q_gnd = specHum_gnd/(1-specHum_gnd)
          
    end subroutine getSurfaceMixingRatio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine integrateSoilDiffusion(gndScalars,lastSurfScalars,ind, &
    ix,jy,zGnd,porosity,satPotential,soilExponent,heatCapSoil, &
    satHydrCond)
    use globals
    use SEBmodule
    implicit none

    interface
    include './interfaces/getSoilThermalTransfer.f90'
    include './interfaces/getWaterConductivity.f90'
    include './interfaces/solveTridiagonalSystem.f90'
    end interface

    integer*4 :: ix,jy,ind
    real*8,dimension(:) :: lastSurfScalars,zGnd,porosity,satPotential, &
    soilExponent,heatCapSoil,satHydrCond
    real*8,dimension(:,:) :: gndScalars

    integer :: i
    real*8 :: dt_
    real*8 :: D_mid(size(gndScalars,1)-1),D(size(gndScalars,1)), &
    K(size(gndScalars,1)),dKdz
    real*8 :: heatCapacity(size(gndScalars,1)), &
    k_mid(size(gndScalars,1)-1),z_mid(size(gndScalars,1)-1)
!     k_mid and z_mid is the value midway between nodes
          
! for moisture diffusion
!     D is diffuse conductivity and K is hydraulic conductivity
    real*8,dimension(size(gndScalars,1)-1) :: b,e,f,g
!     use e, f, and g to store diagonal 'columns' of data rather than a matrix M
!     this elimates the unnecessary storage of zeros in M
!     where M would be a tridiagonal matrix, e, f, and g are:
!     matrix is of form, M = [f(1), g(1), 0...       ...0;
!                             e(2), f(2), g(2), 0... ...0;
!                             ...                     ...;
!                             0...        ...0, e(end), f(end)]


! actual integration time,
! not == to dt when soil conditions are not updated every timestep

! set by integrateSoilDiffFreq
    if( (t+nrsub) > endConstSEB )then
        dt_ = dt*integrateSoilDiffFreq
    else
        dt_ = dt
    endif

    if( ind == 1)then
    !    compute soil conductivity based on moisture content
        call getSoilThermalTransfer(gndScalars(:,2),k,1,ix,jy,porosity, &
        satPotential,soilExponent,heatCapSoil)

    ! interpolate zGnd and k to get values at mid-levels
    ! for z do this once and save z_mid (never changes)
        do i = 1,soilLevels-1
            k_mid(i) = sum(k(i:i+1))/2.d0
            z_mid(i) = sum(zGnd(i:i+1))/2.d0
        enddo

    else

    ! compute diffusive and hydraulic water conductivity of soil
    !     given the soil properties from getGroundParams.f
        call getWaterConductivity(gndScalars(:,2),D,K,ix,jy,porosity, &
        satPotential,satHydrCond,soilExponent)
                  
        do i = 1,soilLevels-1
            z_mid(i) = sum(zGnd(i:i+1))/2.d0
            k_mid(i) = sum(D(i:i+1))/2.d0
        enddo

    endif


! set coefficients for node 2, the first row of M (top node, surface is solved )
! node 2 and soilLevels coefficients have a different form of implicit finite diff.
!   because of the boundary conditions.
    f(1)   = ( dt_/(2.d0*(z_mid(2)-z_mid(1))) ) * &
    ( k_mid(1)/(zGnd(2)-zGnd(1)) + k_mid(2)/(zGnd(3)-zGnd(2)))+ &
    &      1.d0
    g(1)   = - dt_*k_mid(2)/(2.d0*(z_mid(2)-z_mid(1))*(zGnd(3)- &
    zGnd(2)))
          

    b(1)   = (dt_/(2.d0*(z_mid(2)-z_mid(1)))) &
    * ( (k_mid(2)*(gndScalars(3,ind) - gndScalars(2,ind)) &
    / (zGnd(3)-zGnd(2))) - ( k_mid(1) &
    * (gndScalars(2,ind) - lastSurfScalars(ind)) &
    /(zGnd(2)-zGnd(1))) ) &
    + gndScalars(2,ind) + gndScalars(1,ind)*dt_*k_mid(1) &
    / (2.d0*(z_mid(2)-z_mid(1))*(zGnd(2)-zGnd(1)))

    if( ind == 2 )then
        dKdz = dt_*(K(1)*(zGnd(2)-zGnd(3)) &
        / ((zGnd(1)-zGnd(2))*(zGnd(1)-zGnd(3))) &
        + K(2)*(2.d0*zGnd(2)-zGnd(1)-zGnd(3)) &
        / ((zGnd(2)-zGnd(1))*(zGnd(2)-zGnd(3))) &
        + K(3)*(zGnd(2)-zGnd(1)) &
        / ((zGnd(3)-zGnd(1))*(zGnd(3)-zGnd(2))))

        b(1) = b(1) + dKdz
    endif

!  set coefficients for node soilLevels (last ground node)
    e(soilLevels-1) = -dt_*k_mid(soilLevels-1) &
    / (2.d0*(zGnd(soilLevels)-zGnd(soilLevels-1))**2)
    f(soilLevels-1) = 1.d0 + dt_*k_mid(soilLevels-1) &
    /(2.d0*(zGnd(soilLevels)-zGnd(soilLevels-1))**2)
    b(soilLevels-1) = gndScalars(soilLevels,ind) &
    - ( gndScalars(soilLevels,ind) &
    - gndScalars(soilLevels-1,ind) )*dt_*k_mid(soilLevels-1) &
    /(2.d0*(zGnd(soilLevels)-zGnd(soilLevels-1))**2)

    if( ind == 2 )then
        dKdz = dt_*(K(soilLevels)-K(soilLevels-1)) &
        / (zGnd(soilLevels)-zGnd(soilLevels-1))
                 
        b(soilLevels-1) = b(soilLevels-1) + dKdz
    endif

! set matrix coefficients for nodes 3 through soilLevels-1
! form of these coefficients is identical, there are no boundary conditions
    if( soilLevels > 3)then
        do i = 3,soilLevels-1
            e(i-1) = -dt_*k_mid(i-1)/(2.d0*(z_mid(i)-z_mid(i-1))* &
            (zGnd(i)-zGnd(i-1))) !M(i-1,i-2)
                        
            f(i-1) = 1.d0 + dt_*k_mid(i-1) &
            / (2.d0*(z_mid(i)-z_mid(i-1))*(zGnd(i)-zGnd(i-1))) &
            +dt_*k_mid(i) &
            / (2.d0*(z_mid(i)-z_mid(i-1))*(zGnd(i+1)-zGnd(i))) !M(i-1,i-1)
                        
            g(i-1)   = -dt_*k_mid(i) / &
            (2.d0*(z_mid(i)-z_mid(i-1))*(zGnd(i+1)-zGnd(i))) !M(i-1,i)
                        
            b(i-1)     = gndScalars(i,ind) &
            + (dt_/(2.d0*(z_mid(i)-z_mid(i-1))))* &
            (k_mid(i)*(gndScalars(i+1,ind) &
            -gndScalars(i,ind))/(zGnd(i+1)-zGnd(i)) &
            -k_mid(i-1)*(gndScalars(i,ind)-gndScalars(i-1,ind)) &
            / (zGnd(i)-zGnd(i-1)))

            if( ind == 2)then
            !     explicit finite difference for unevenly spaced data
                dKdz = dt_*(K(i-1)*(zGnd(i)-zGnd(i+1)) &
                / ((zGnd(i-1)-zGnd(i))*(zGnd(i-1)-zGnd(i+1))) &
                + K(i)*(2.d0*zGnd(i)-zGnd(i-1)-zGnd(i+1)) &
                / ((zGnd(i)-zGnd(i-1))*(zGnd(i)-zGnd(i+1))) &
                + K(i+1)*(zGnd(i)-zGnd(i-1)) &
                / ((zGnd(i+1)-zGnd(i-1))*(zGnd(i+1)-zGnd(i))))
                               
                b(i-1) = b(i-1) + dKdz

            endif
                        
        enddo
    endif

!      call solveTridiagonalSystem(e,f,g,b,
!     >     gndScalars(2:size(gndScalars,1),ind))
    call tridag(e,f,g,b,gndScalars(2:soilLevels,ind),soilLevels-1)

    end subroutine integrateSoilDiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine solveTridiagonalSystem(eDiag,fDiag,gDiag,bDiag,x)
!     given e, f, and g components of tridiagonal matrix (as described below)
!     and b (the RHS of a system of equations)
!     solveTridiagonalSystem() uses the Thomas algorithm to solve the tridiagonal system of eqtns

!     inputs:
!           e, f, and g to store diagonal 'columns' of data rather than matrix M
!           this elimates the unnecessary storage of zeros in M
!           where M would be a tridiagonal matrix, e, f, and g are:
!           matrix is of form, M = [f(1), g(1), 0...       ...0;
!                                  e(2), f(2), g(2), 0... ...0;
!                                  ...                     ...;
!                                  0...        ...0, e(end), f(end)]
!           b - RHS of system of equations where Mx=b, where x are unknowns

!     ouputs:
!           x - array of size soilLevels-1 = size(b), where Mx=b
!               x is the array of unknowns given a tridiagonal system of eqtns
    implicit none

    real*8, dimension(:):: bDiag, eDiag, fDiag, gDiag, x

    integer :: i, length

    length = size(bDiag)

!     decomposition
    do i=2,length
        eDiag(i) = eDiag(i)/fDiag(i-1)
        fDiag(i) = fDiag(i) - eDiag(i)*gDiag(i-1)
    enddo
          
!     forward substitution
    do i=2,length
        bDiag(i) = bDiag(i) - eDiag(i)*bDiag(i-1)
    enddo

!     back substitution
    x = bDiag/fDiag
    do i = length, 1, -1
        x(i) = ( bDiag(i) - gDiag(i)*x(i+1) )/fDiag(i)
    enddo
          
    end subroutine solveTridiagonalSystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine getSoilThermalTransfer(moisture,thermalTransfer,flag, &
    ix,jy,porosity,satPotential,soilExponent,heatCapSoil)
    use wallBoundaryConditions
    use SEBmodule
    implicit none

    integer*4 :: flag,ix,jy
    real*8, dimension(:):: moisture,thermalTransfer,porosity, &
    satPotential,soilExponent,heatCapSoil

    real*8, dimension(size(moisture)):: Pf, heatCap
    integer :: i
!  flag = 0 if thermal conductivity is requested (for solving Q = v dT/dz )
!  flag = 1 if thermal diffusivity is requested (for solving dT/dt = d( k dT/dz)dz

!     compute soil conductivity from emperical formula (McCumber 1980)
    Pf = log10( abs(100.d0*z_i*satPotential(:)* &
    (porosity(:)/moisture)**soilExponent(:) ) )
    do i=1,size(Pf)
        if ( Pf(i) <= (5.1) )then
            thermalTransfer(i) = 418.46d0*exp( -( Pf(i)+2.7d0 ) )
        else
            thermalTransfer(i) = 0.172d0
        endif
    enddo
!   emperical relationship for thermal conductivity is in [ J/(m*s*K) ]
!   non-dimensionalize the conductivity
    if(flag==1)then
        heatCap = (1-porosity(:))*heatCapSoil(:) &
        + moisture*heatCapWater
        thermalTransfer = thermalTransfer &
        /(heatCap*densityAir*Cp_air*z_i*u_star)
    else
        thermalTransfer = thermalTransfer &
        /(densityAir*Cp_air*z_i*u_star)
    endif

    end subroutine getSoilThermalTransfer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine getWaterConductivity(moisture,diffCond,hydrCond,ix,jy, &
    porosity,satPotential,satHydrCond,soilExponent)

    implicit none

    integer*4 :: ix,jy
    real*8, dimension(:):: moisture,diffCond,hydrCond,porosity, &
    satPotential,satHydrCond,soilExponent

!     Emperical relationships from Clapp and Hornberger (1978)
!     need to be non-dimensionalized
    diffCond = -(soilExponent(:)*satHydrCond(:)* &
    satPotential(:)/moisture(:)) &
    *(moisture(:)/porosity(:))**(soilExponent(:)+3.d0)

    hydrCond  = satHydrCond(:)*(moisture(:)/porosity(:))** &
    (2.d0*soilExponent(:) + 3.d0)
   
    end subroutine getWaterConductivity

