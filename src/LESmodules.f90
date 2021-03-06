!     LESmodules.f contains all the required modules for les_conv_scalar.f
!     Each module declares the variables required for a specific set of
!         subroutines and les_conv_scalar.f shall use all of the modules.
!         All subroutines use globals
!     The modules in this file include:
!              globals
!              momentum
!              statistics
!              wallBoundaryConditions
!              scalars
!              sgs
!              fileUnits
!              mainModule

    MODULE globals
!    use mpi
!     ALL program units will 'use globals'  (except [s_]zeroslice)
    implicit none
    include 'mpif.h'
    integer*4 :: Nx,fflag,Ny,Nz,Nx2,Nxb2,Ny2,Nyb2,aNx,l_r
          
    real*8 :: dt, inxny
    real*8 :: Co,fgr,tfr
    integer*4 :: t,ttt,me,nall,ierr
    integer*4 :: initu,inits
    integer*4 :: p_count,c_count
    integer*4 :: nsteps,nrsub
    integer*4 :: verticalBC,mom_nodes,scl_nodes,averaging
    integer*4 :: nxb,nyb,nzb,nz2,nprocs
    integer*4 :: vprocs,hprocs,hfact,vfact,ip
    integer*4 :: FFT_FLAG,MPI_COMM_LEVEL,MPI_COMM_COLUMN
    real*8 :: dx,dy,dz,delta,l_z
    real*8 :: Pi,Sc,Ugal,Vgal,g_hat
    real*8 :: UTC,startUTC,UTC_hrs
          
! lagrng_sd_s  lagrng_sd   ddz_w  ddz_uv_p   ddz_uv
    real*8 :: idz
! dealias2
    real*8 :: inx2ny2
! lagrng_sd_s  lagrng_sd
    real*8 :: dtl, idx, idy

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE wallBoundaryConditions

! surf_flux wallstress2
    real*8 :: z_i

! surf_flux sgs_stag wallstress2 derivwall2
    real*8 :: vonk

! surf_flux
    real*8 :: u_star,mhfx,mhfy,swtx,swty,stepsPerDx,stepsPerDy


    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE scalars

    integer*4 :: sponge,scalarCount,Ri_flag,temperatureIndex, &
    moistureIndex
    real*8 ::    theta_0,z_d,Ug,Vg,rlx_time
    integer*4,allocatable,dimension(:)::scalarFlags, &
    surfaceFlags,S_advec
    real*8,allocatable,dimension(:)::surfaceFluxes,scalarScales, &
    dsdtHomogeneous,inversion

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE sgs

    integer*4 :: model, cs_count
    real*8 ::    nnn, nu

    END MODULE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE mainModule

    integer*4 :: fileview_2D,fileview_3D
    integer*4 :: ruler, M_advec

    real*8 :: fc, f_c

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE particleModule

    integer*4 :: part_model,start_release,partstep,freq,nr, &
    skip_step,deposition,nodep_copy,fileview_part
    integer*8 :: npart,ipart
    real*8 :: Cop,es_min,det_min,wd,Lv

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE canopyModule

    integer*4 :: c_flag
    real*8 :: Cd,h_canopy

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
    MODULE press_force

    integer*4 :: press_cor,press_step
    real*8 :: u_avg,f_p,theta_mean_wind

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE SEBmodule

    integer*4 :: soilLevels,maxFluxIterations,maxTempIterations, &
    endConstSEB,updateFreqSEB,integrateSoilDiffFreq, &
    radiationFlag,stepsPerRadVal,day,albedoFlag,fileview_soil
    real*8 :: zt,pressureScale,densityAir,Cp_air,densityWater, &
    latentHeatWater,heatCapWater,waterGasConst,moistureCriteria, &
    temperatureCriteria,tempFluxCriteria,convFactor,SB_constant, &
    solarIrradiance,lat,long,emissivity

    END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    MODULE frameModule
          
    integer*4 :: nframe, sframe, framestep, frameh, max_frames, &
    framefiles, frame_cnt, frame_possibilities, &
    MPI_COMM_FRAMES

    parameter(frame_possibilities=10) !NOTE: this corresponds to how many variables are listed in the header as possible frame choices

    character(50) :: framenames(frame_possibilities)
    logical,dimension(frame_possibilities) :: output_frames

    END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    MODULE visualizationModule
    integer(kind = 4) :: filtLevelsLambda
    real(kind = 8), allocatable, dimension(:) :: filtsLambda    
    END MODULE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
