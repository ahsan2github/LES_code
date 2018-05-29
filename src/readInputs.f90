    subroutine readInputs()

    use globals
    use wallBoundaryConditions
    use scalars
    use sgs
    use mainModule
    use particleModule
    use canopyModule
    use press_force
    use SEBmodule
    use frameModule
    use visualizationModule
    implicit none
    include 'fftw3.f90'

    integer*4 :: i,ll, ind
    real*8 :: nx2Mult, ny2Mult, nscalar2Mult
    character(50) :: charIN
          
!      READ in parameters
                                                           
    open(unit=1,file='input/LESinputs.txt',status='old')
    do i = 1,6
        read(1,*)
    enddo
    read(1,*) vonk
    read(1,*) sc
    read(1,*) Co
    read(1,*) nnn
    read(1,*) pi
    do i=1,3
        read(1,*)
    enddo
    read(1,*) model
    read(1,*) averaging
    read(1,*) fgr
    read(1,*) tfr
    read(1,*) mom_nodes
    read(1,*) scl_nodes
    read(1,*) FFT_FLAG
    do i = 1,3
        read(1,*)
    enddo
    read(1,*) nsteps
    read(1,*) hprocs
    read(1,*) dt
    read(1,*) Nx
    read(1,*) Ny
    read(1,*) Nz
    read(1,*) aNx
    read(1,*) l_r
    read(1,*) z_i
    read(1,*) u_star
    read(1,*) l_z
    read(1,*) verticalBC
    read(1,*) Ugal
    read(1,*) Vgal
    read(1,*) Ug
    read(1,*) Vg
    read(1,*) press_cor
    read(1,*) f_p
    read(1,*) press_step
    read(1,*) u_avg
    read(1,*) theta_mean_wind
    read(1,*) f_c
    read(1,*) nu
    read(1,*) M_advec
    do i = 1,3
        read(1,*)
    enddo

! READ scalar parameters                                                            write(*,*)'reading scalar parameters'
    read(1,*) scalarCount

    if(scalarCount > 0)then
        ALLOCATE( scalarFlags(scalarCount), &
        surfaceFlags(scalarCount), &
        surfaceFluxes(scalarCount), &
        scalarScales(scalarCount), &
        dsdtHomogeneous(scalarCount), &
        S_advec(scalarCount), &
        inversion(scalarCount))
        do i=1,scalarCount
            read(1,*) scalarFlags(i)
        enddo
        do i=1,scalarCount
            read(1,*) surfaceFlags(i)
        enddo
        do i=1,scalarCount
            read(1,*) surfaceFluxes(i)
        enddo
        do i=1,scalarCount
            read(1,*) scalarScales(i)
        enddo
        do i = 1,scalarCount
            read(1,*) dsdtHomogeneous(i)
        enddo
        do i = 1,scalarCount
            read(1,*) S_advec(i)
        enddo
        do i = 1,scalarCount
            read(1,*) inversion(i)
        enddo
    else
        do i=1,7
            read(1,*)
        enddo
    endif
    read(1,*) theta_0
    read(1,*) sponge
    read(1,*) z_d
    read(1,*) rlx_time
    read(1,*) Ri_flag
    do i = 1,3
        read(1,*)
    enddo

! read computation and print frequency for stats and coefficients
    read(1,*) c_count
    read(1,*) p_count
    read(1,*) cs_count
    read(1,*)
!     comment
    read(1,*) nframe
    read(1,*) sframe
    read(1,*) framestep
    read(1,*) frameh
    read(1,*) framefiles
    do i=1,framefiles
        read(1,'(A)') charIN
        ind=index(charIN,'frame')+4
        framenames(i)=charIN(1:ind)
    end do
    do i = 1,3
        read(1,*)
    enddo

!     READ resubmission parameters
    read(1,*) initu
    read(1,*) inits
    read(1,*) ruler
    do i = 1,3
        read(1,*)
    enddo

!     READ particle dispersion parameters

    read(1,*) npart
    read(1,*) part_model
    read(1,*) start_release
    read(1,*) partstep
    read(1,*) skip_step
    read(1,*) freq
    read(1,*) nr
    read(1,*) Cop
    read(1,*) es_min
    read(1,*) det_min
    read(1,*) deposition
    read(1,*) nodep_copy
    read(1,*) wd
    read(1,*) Lv
    do i = 1,3
        read(1,*)
    enddo

!     READ canopy parameters

    read(1,*) c_flag
    read(1,*) Cd
    read(1,*) h_canopy
    do i = 1,3
        read(1,*)
    enddo

!     READ soil type parameters

    read(1,*) soilLevels
    read(1,*) zt
    read(1,*) pressureScale
    read(1,*) densityAir
    read(1,*) Cp_air
    read(1,*) densityWater
    read(1,*) latentHeatWater
    read(1,*) heatCapWater
    read(1,*) waterGasConst
    read(1,*) moistureCriteria
    read(1,*) temperatureCriteria
    read(1,*) tempFluxCriteria
    read(1,*) maxFluxIterations
    read(1,*) maxTempIterations
    read(1,*) convFactor
    read(1,*) endConstSEB
    read(1,*) updateFreqSEB
    read(1,*) integrateSoilDiffFreq
    read(1,*) albedoFlag
    do i = 1,3
        read(1,*)
    enddo

!     READ radiation parameters

    read(1,*) radiationFlag
    read(1,*) stepsPerRadVal
    read(1,*) SB_constant
    read(1,*) solarIrradiance
    read(1,*) lat
    read(1,*) long
    read(1,*) day
    read(1,*) emissivity
    close(1)

!     calculated grid parameters
          
    vprocs = nprocs/hprocs
    nx2 = nx*3/2
    nxb = Nx/hprocs
    nxb2 = nxb*3/2
    ny2 = ny*3/2
    nyb = Ny/hprocs
    nyb2 = nyb*3/2
    nzb = nz/vprocs
    nz2 = nzb+2
    dtl = dt*cs_count
    dx=2.0d0*Pi/Nx
    dy=2.0d0*Pi/Ny/l_r
    dz=(L_z/z_i)/(Nz-1)
    delta=fgr*(dx*dy*dz)**(1./3.)
    idz=1.d0/dz
    idx=1.d0/dx
    idy=1.d0/dy
    inxny=1.d0/(Nx*Ny)
    inx2ny2=1.d0/(Nx2*Ny2)
    hfact = floor(dble(me)/dble(hprocs))*hprocs
    vfact = hfact/hprocs
    if(Ugal /= 0.d0)then
        mhfx = ceiling(dx/(dt*Ugal)/2)
        swtx = nint(10*dmod((mhfx*dt*Ugal),dx)/dx)
        stepsPerDx = dx/(abs(Ugal)*dt)
    endif
    if(Vgal /= 0.d0)then
        mhfy = ceiling(dy/(dt*Vgal)/2)
        swty = nint(10*dmod((mhfy*dt*Vgal),dy)/dy)
        stepsPerDy = dy/(abs(Vgal)*dt)
    endif
    startUTC = 0.d0*3600/(z_i/u_star)

!     scalar parameters

    temperatureIndex=0
    moistureIndex=0
    if(scalarCount > 0)then
        do ll=1,scalarCount
            if (scalarFlags(ll) == 1)then
                temperatureIndex = ll
            elseif (scalarFlags(ll) == 2)then
                moistureIndex = ll
            endif
        enddo
    endif

!     nondimensionalizations

    dt=dt*u_star/z_i
    dtl=dtl*u_star/z_i
    f_p=f_p/u_star**2*z_i
    fc=f_c*10.**(-4)*z_i/u_star
    g_hat=9.81d0*(z_i/(u_star**2))
    Ugal = Ugal/u_star
    Vgal = Vgal/u_star
    Ug   = Ug/u_star
    Vg   = Vg/u_star
    u_avg = u_avg/u_star
    if(scalarCount >= 1)then
        do ll=1,scalarCount
            inversion(ll)=inversion(ll)*z_i/scalarScales(ll)
            surfaceFluxes(ll)= &
            surfaceFluxes(ll)/u_star/scalarScales(ll)*10000
        enddo
        if(temperatureIndex /= 0)then
            theta_0=theta_0/scalarScales(temperatureIndex)
        endif
    endif
    if(soilLevels > 0)then
        densityWater = densityWater/densityAir
        latentHeatWater = &
        latentHeatWater/Cp_air*scalarScales(temperatureIndex)
        heatCapWater = heatCapWater/(densityAir*Cp_air)
        waterGasConst = &
        waterGasConst*scalarScales(temperatureIndex)/u_star**2
        SB_constant=SB_constant*scalarScales(temperatureIndex)**3/ &
        (Cp_air*densityAir*u_star)
        solarIrradiance = &
        solarIrradiance/(scalarScales(temperatureIndex)*u_star)
        lat = lat*pi/180.d0
        long = long*pi/180.d0
    endif
    es_min=es_min/u_star**2
    if(deposition == 1)then
        wd=wd/u_star
        Lv=Lv/z_i
    endif

!     test for grid compatability

    IF(me == 0)THEN

        if(mod(nprocs,hprocs) /= 0)then
            write(*,*)'ERROR:nprocs must be evenly divisible by hprocs'
            stop
        elseif(mod(nx,2) /= 0)then
            write(*,*)'ERROR:Nx must be evenly divisible by 2'
            stop
        elseif(mod(nx/2,hprocs) /= 0)then
            write(*,*)'ERROR:Nx/2 must be evenly divisible by hprocs'
            stop
        elseif(mod(ny,hprocs) /= 0)then
            write(*,*)'ERROR:Ny must be evenly divisible by hprocs'
            stop
        elseif(mod(nz,vprocs) /= 0)then
            write(*,*)'ERROR:Nz must be evenly divisible by vprocs'
            stop
        elseif(mod(nx/hprocs,2) /= 0)then
            write(*,*)'ERROR:Nx/hprocs must be even'
            stop
        elseif(mod(ny/hprocs,2) /= 0)then
            write(*,*)'ERROR:Ny/hprocs must be even'
            stop
        elseif(mod(ny,vprocs) /= 0)then
            write(*,*)'ERROR:Ny must be evenly divisible by vprocs'
            stop
        endif

    ENDIF

    if(me-hfact == hprocs-1)then
        ip=1
    else
        ip=0
    end if

    if(FFT_FLAG == 1)then
        FFT_FLAG=FFTW_ESTIMATE
    elseif(FFT_FLAG == 2)then
        FFT_FLAG=FFTW_PATIENT
    elseif(FFT_FLAG == 3)then
        FFT_FLAG=FFTW_MEASURE
    elseif(FFT_FLAG == 4)then
        FFT_FLAG=FFTW_EXHAUSTIVE
    endif

    return
    end subroutine readInputs
