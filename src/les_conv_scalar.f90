!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c    This program is developed for parallel            c
!     c    computation of scalar transport in ABL            c
!     c                                                      c
!     c                 July, 2002                           c
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Program LES_CONV_SCALAR

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
interface
include 'interface_main.f90'
end interface

!!!!!!Variable Declarations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: idum
integer*4 :: b,c,i,j,k,l,mm,cst,pt,vt,kk,t01,t02,ii,iii,il,ssi, &
temperatureFlag,mfst,kstart,resubnum,ios,oldsteps,loopCount, &
fh,scalarIter,test_u,test_s,tstat,NNz, &
fh1,fh2,fh3,fh4,fh5,fh6,fh7,fh8,fh9,fh10,fh11,fh12,fh13,fh14, &
rogue_count,trajectory_updates, looper

integer (kind=MPI_OFFSET_KIND) :: view_disp
logical :: dir_ext

real*8 :: sum_s,z,tt,rmsdivvel,rmsdivve,kee,ke,wgx,oldtime, &
ztemp,uave,tave,ubar,wvar,Ugeo_g,Uadv_d,Vadv_d,Sadv_d,Vgeo_g, &
Ptime,ssr,Uref,netRad

real :: ran1,randn

integer*4,allocatable,dimension(:)::ilow,jlow,nowTime, &
dsdxUnits,dsdyUnits,dsdzUnits, &
scalarMeanUnits,scalar2Units,scalar3Units,xScalarFluxUnits, &
yScalarFluxUnits,zScalarFluxUnits,awsUnits,ausUnits,avsUnits, &
scalarSpectraUnits,scalarPrandtlUnits,cs2prUnits,beta2Units, &
ETUnits,surfaceFluxUnits, &
obukovUnits,scalarStarUnits, rogue_flag

real*8,allocatable,dimension(:)::force_x,force_y,q_bar,ke_bar, &
Sadv,Uadv,Vadv,Ugeo,Vgeo,awtT,awtT_global,zGnd,measRad,fx,fz, &
ESGS_part_f

real*8,allocatable,dimension(:,:)::zo,ustar,M,atxz_s,au,av,aw,u2, &
v2,w2,w3,p2,auw,avw,auv,ap,adudz,adudx,atxz,atyy,atxx,atyz, &
atzz,atxy,e,aCs2,aCs,adwdz,adwdx,abeta1,advdz,aESGS,aDSGS,specu, &
specv,specw,specp,release_pos,temp,albedo,minAlbedo, &
soilHeatFlux,obukhovL,u_hat,v_hat,tempFrameA2, &
particle_tmp,xz_bar, aTL

real*8,allocatable,dimension(:,:,:)::u,v,w,P,cx,cy,cz,RHSx,RHSy, &
RHSz,RHSx_f,RHSy_f,RHSz_f,dudx,dudy,dudz,dvdx,dvdy,dvdz, &
dwdx,dwdy,dwdz,txx,txy,txz,tyy,tyz,tzz,divtx,divty,divtz, &
dpdx,dpdy,dpdz,Beta,ddtzz,u_m,v_m,w_m,scalarMean,t2,t3, &
asgs_t1,asgs_t2,asgs_t3,aut,avt,awt,adsdx,adsdy, &
adsdz,aPr,aCs2Pr,abeta2,aqz_s,aET,surfaceScalar,scalarFlux, &
coolrate,thetaFactor,spect,a1_old,b1_old,c1_old,d1_old, &
e1_old,a2_old,b2_old,c2_old,d2_old,e2_old,tmp_old,Cs2,Cs2_m, &
beta1,fi,fi_h,Psi,Psi0,rdmp,ESGS3D,DSGS3D,u_global,v_global, &
w_global,LAD_uvp,LAD_w,Fdx,Fdy,Fdz,u_w,v_w,mag_w,w_uvp, &
mag_uvp,porosity,satPotential,satHydrCond,soilExponent, &
heatCapSoil,scalar_hat,particle,TL,us_f,ESGS_f, lambda2, hold4now

real*8,allocatable,dimension(:,:,:,:):: scalar,qx,qy,qz,scalarRHS, &
scalarRHSf,S,dsdx,dsdy,dsdz,Pr2,Pr2_m,beta2, &
a4_old,b4_old,c4_old,d4_old,e4_old, &
a8_old,b8_old,c8_old,d8_old,e8_old,ESCL3D,gndScalars,tau_f

character(50) :: fileStr1,fileStr2,fileStr3,fileStr4,fileStr5,fileStr6
character(100) :: frameStr, fileStr21, fileStr22, fileStr23, fileStr24
character(50) :: fileStr25, fileStr31, fileStr32, fileStr33, fileStr34
character(50) :: fileStr35

!!!!!!Initialize MPI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, me, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
call MPI_COMM_DUP( MPI_COMM_WORLD, nall, ierr)

!!!!!!Read Input Values!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call readInputs()

call MPI_BARRIER(nall,ierr)  !synchronize after reading inputs

!!!!!!Define MPI Communicator for Each LEVEL!!!!!!!!!!!!!

call MPI_COMM_SPLIT(nall,vfact,vfact,MPI_COMM_LEVEL,ierr)

!!!!!!Define MPI Communicator for Each Column!!!!!!!!!!!

call MPI_COMM_SPLIT(nall,me-hfact,me-hfact,MPI_COMM_COLUMN,ierr)

!!!!!!Define MPI Communicator for Frame Output!!!!!!!!!!

if((vfact+1)*nzb <= frameh)then
    i=1
else
    i=0
endif

call MPI_COMM_SPLIT(nall,i,i,MPI_COMM_FRAMES,ierr)

!!!!!!Allocate Variable Sizes!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(ilow(nx),jlow(nyb),nowtime(3))

allocate(force_x(nz2),force_y(nz2),q_bar(nz2),ke_bar(nz2), &
Sadv(nz2),Uadv(nz2),Vadv(nz2),Ugeo(nz2),Vgeo(nz2),awtT(nzb), &
awtT_global(nz))

allocate(atxz_s(nx,nyb),au(anx,nz2),av(anx,nz2),aw(anx,nz2), &
u2(anx,nz2),v2(anx,nz2),w2(anx,nz2),w3(anx,nz2),p2(anx,nz2), &
auw(anx,nz2),avw(anx,nz2),auv(anx,nz2),ap(anx,nz2), &
adudz(anx,nz2),adudx(anx,nz2),atxz(anx,nz2),atyy(anx,nz2), &
atxx(anx,nz2),atyz(anx,nz2),atzz(anx,nz2),atxy(anx,nz2), &
e(anx,nz2),aCs2(anx,nz2),aCs(anx,nz2),adwdz(anx,nz2), &
adwdx(anx,nz2),abeta1(anx,nz2),advdz(anx,nz2),aESGS(anx,nz2), &
specu(Nx/2+1,nz2),specv(Nx/2+1,nz2),specw(Nx/2+1,nz2), &
specp(Nx/2+1,nz2),xz_bar(nx,nz2), aDSGS(anx, nz2), aTL(anx, nz2))

allocate(u(nx,nyb,nz2),v(nx,nyb,nz2),w(nx,nyb,nz2),P(nx,nyb,nz2), &
cx(nx,nyb,nz2),cy(nx,nyb,nz2),cz(nx,nyb,nz2),RHSx(nx,nyb,nz2), &
RHSy(nx,nyb,nz2),RHSz(nx,nyb,nz2),RHSx_f(nx,nyb,nz2), &
RHSy_f(nx,nyb,nz2),RHSz_f(nx,nyb,nz2),dudx(nx,nyb,nz2), &
dudy(nx,nyb,nz2),dudz(nx,nyb,nz2),dvdx(nx,nyb,nz2), &
dvdy(nx,nyb,nz2),dvdz(nx,nyb,nz2),dwdx(nx,nyb,nz2), &
dwdy(nx,nyb,nz2),dwdz(nx,nyb,nz2),txx(nx,nyb,nz2), &
txy(nx,nyb,nz2),txz(nx,nyb,nz2),tyy(nx,nyb,nz2), &
tyz(nx,nyb,nz2),tzz(nx,nyb,nz2),divtx(nx,nyb,nz2), &
divty(nx,nyb,nz2),divtz(nx,nyb,nz2),dpdx(nx,nyb,nz2), &
dpdy(nx,nyb,nz2),dpdz(nx,nyb,nz2),Beta(nx,nyb,nz2), &
ddtzz(nx,nyb,nz2),u_m(nx2,nyb2,nz2),v_m(nx2,nyb2,nz2), &
w_m(nx2,nyb2,nz2),tempFrameA2(nx,nyb), &
Cs2(nx,nyb,nz2),Cs2_m(nx2,nyb2,nz2),beta1(nx,nyb,nz2), &
rdmp(nx,nyb,nz2),ESGS3D(nx,nyb,nz2),DSGS3D(nx,nyb,nz2), &
u_global(nx,ny,nz),v_global(nx,ny,nz),w_global(nx,ny,nz), &
TL(nx,nyb,nz2),lambda2(nx, nyb, nz2), hold4now(nx,nyb,nz2))

if(verticalBC == 0)then
    allocate(zo(nx,nyb),ustar(nx,nyb),M(nx,nyb))
endif
! nly bottom procs need allocate surface variables
if(vfact == 0 .AND. verticalBC == 0)then
    allocate(Psi(nx,nyb,2),Psi0(nx,nyb,2),fi(nx,nyb,2), &
    fi_h(nx,nyb,2))
endif

if(model == 2)then
    allocate(a1_old(nx,nyb+2,nz2),b1_old(nx,nyb+2,nz2), &
    tmp_old(nx,nyb,nzb))
elseif(model == 3 .AND. averaging == 1)then
    allocate(a1_old(nx,nyb+2,nz2),b1_old(nx,nyb+2,nz2), &
    c1_old(nx,nyb+2,nz2),d1_old(nx,nyb+2,nz2), &
    e1_old(nx,nyb+2,nz2),a2_old(nx,nyb+2,nz2), &
    b2_old(nx,nyb+2,nz2),c2_old(nx,nyb+2,nz2), &
    d2_old(nx,nyb+2,nz2),e2_old(nx,nyb+2,nz2), &
    tmp_old(nx,nyb,nzb))
endif

! cc  Scalar Variables
IF(scalarCount >= 1)THEN
    allocate(scalar(nx,nyb,nz2,scalarCount), &
    scalarRHS(nx,nyb,nz2,scalarCount), &
    scalarRHSf(nx,nyb,nz2,scalarCount), &
    S(nx,nyb,nz2,scalarCount),qx(nx,nyb,nz2,scalarCount), &
    qy(nx,nyb,nz2,scalarCount),qz(nx,nyb,nz2,scalarCount), &
    dsdx(nx,nyb,nz2,scalarCount),dsdy(nx,nyb,nz2,scalarCount), &
    dsdz(nx,nyb,nz2,scalarCount),Pr2(nx,nyb,nz2,scalarCount), &
    Pr2_m(nx2,nyb2,nz2,scalarCount), &
    beta2(nx,nyb,nz2,scalarCount), &
    surfaceScalar(nx,nyb,scalarCount), &
    scalarFlux(nx,nyb,scalarCount), &
    coolrate(nx,nyb,scalarCount),thetaFactor(nx,nyb,nz2), &
    spect(nx/2+1,nz,scalarCount), &
    ESCL3D(nx,nyb,nz2,scalarCount))
    allocate(dsdxUnits(scalarCount),dsdyUnits(scalarCount), &
    dsdzUnits(scalarCount), &
    scalarMeanUnits(scalarCount),scalar2Units(scalarCount), &
    scalar3Units(scalarCount),xScalarFluxUnits(scalarCount), &
    yScalarFluxUnits(scalarCount), &
    zScalarFluxUnits(scalarCount), &
    awsUnits(scalarCount),ausUnits(scalarCount), &
    avsUnits(scalarCount),scalarSpectraUnits(scalarCount), &
    scalarPrandtlUnits(scalarCount),cs2prUnits(scalarCount), &
    beta2Units(scalarCount),ETUnits(scalarCount), &
    surfaceFluxUnits(scalarCount),obukovUnits(scalarCount), &
    scalarStarUnits(scalarCount))
    allocate(scalarMean(anx,nz2,scalarCount), &
    t2(anx,nz2,scalarCount),t3(anx,nz2,scalarCount), &
    asgs_t1(anx,nz2,scalarCount),asgs_t2(anx,nz2,scalarCount), &
    asgs_t3(anx,nz2,scalarCount),aut(anx,nz2,scalarCount), &
    avt(anx,nz2,scalarCount),awt(anx,nz2,scalarCount), &
    adsdx(anx,nz2,scalarCount),adsdy(anx,nz2,scalarCount), &
    adsdz(anx,nz2,scalarCount),aPr(anx,nz2,scalarCount), &
    aCs2Pr(anx,nz2,scalarCount),abeta2(anx,nz2,scalarCount), &
    aET(anx,nz2,scalarCount),aqz_s(nx,ny,scalarCount))

    if(model == 2)then
        allocate(a4_old(nx,nyb+2,nz2,scalarCount), &
        b4_old(nx,nyb+2,nz2,scalarCount))
    elseif(model == 3 .AND. averaging == 1)then
        allocate(a4_old(nx,nyb+2,nz2,scalarCount), &
        b4_old(nx,nyb+2,nz2,scalarCount), &
        c4_old(nx,nyb+2,nz2,scalarCount), &
        d4_old(nx,nyb+2,nz2,scalarCount), &
        e4_old(nx,nyb+2,nz2,scalarCount), &
        a8_old(nx,nyb+2,nz2,scalarCount), &
        b8_old(nx,nyb+2,nz2,scalarCount), &
        c8_old(nx,nyb+2,nz2,scalarCount), &
        d8_old(nx,nyb+2,nz2,scalarCount), &
        e8_old(nx,nyb+2,nz2,scalarCount))
    endif
ENDIF

! cc  Particle Model Variables
IF(npart > 0)THEN
    allocate(release_pos(nr,3),ESGS_f(nx,nyb,nz2),us_f(npart,3,2), &
    particle_tmp(npart,6),ESGS_part_f(npart), &
    rogue_flag(npart))
    if(part_model == 1)then
        allocate(particle(npart,6,1))
    elseif(part_model == 2)then
        allocate(particle(npart,6,2))
    elseif(part_model == 3)then
        allocate(particle(npart,6,3),tau_f(nx,nyb,nz2,6))
    endif
    if(deposition == 1)then
        allocate(fx(nz),fz(nz))
    endif
ENDIF

! cc  Canopy Model Variables
IF(c_flag == 1)THEN
    allocate(LAD_uvp(nx,nyb,nz2),LAD_w(nx,nyb,nz2),Fdx(nx,nyb,nzb), &
    Fdy(nx,nyb,nzb),Fdz(nx,nyb,nzb),u_w(nx,nyb,nz2), &
    v_w(nx,nyb,nz2),mag_w(nx,nyb,nz2),w_uvp(nx,nyb,nz2), &
    mag_uvp(nx,nyb,nz2))
ENDIF

! cc  Soil Model Variables
IF(scalarCount > 0)THEN
    IF(surfaceFlags(1) == 2)THEN
        allocate(porosity(nx,nyb,soilLevels), &
        satPotential(nx,nyb,soilLevels), &
        satHydrCond(nx,nyb,soilLevels), &
        soilExponent(nx,nyb,soilLevels), &
        heatCapSoil(nx,nyb,soilLevels), &
        gndScalars(nx,nyb,soilLevels,2),albedo(nx,nyb), &
        minAlbedo(nx,nyb),soilHeatFlux(nx,nyb),obukhovL(nx,nyb), &
        zGnd(soilLevels),u_hat(nx,nyb),v_hat(nx,nyb), &
        scalar_hat(nx,nyb,scalarCount), &
        measRad(nsteps/stepsPerRadVal+1))
    ENDIF
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(me == 0)then
    call itime(nowtime)
    write(*,*)'time = ',nowtime(1),':',nowtime(2),':',nowtime(3)
endif

scalarIter = 0
! cccc  File unit variables                     ccccc

!       ***** initialize unit array variables for scalar files
!     starting at 600 and going to 860 making it possible for there
!     to be a maximum of 10 scalars (10 scalars will not be computed anytime in
!     the near future)
!  Unit=888 is used for another file ( don't go above this without changing)
do l=0,scalarCount-1

! output files
    dsdxUnits(l+1)                      = 680 + l
    dsdyUnits(l+1)                      = 690 + l
    dsdzUnits(l+1)                      = 700 + l
    scalarMeanUnits(l+1)                = 710 + l
    scalar2Units(l+1)                   = 720 + l
    scalar3Units(l+1)                   = 730 + l
    xScalarFluxUnits(l+1)               = 740 + l
    yScalarFluxUnits(l+1)               = 750 + l
    zScalarFluxUnits(l+1)               = 760 + l
    awsUnits(l+1)                       = 770 + l
    ausUnits(l+1)                       = 780 + l
    avsUnits(l+1)                       = 790 + l
    scalarSpectraUnits(l+1)             = 800 + l
    scalarPrandtlUnits(l+1)             = 810 + l
    cs2prUnits(l+1)                     = 820 + l
    beta2Units(l+1)                     = 830 + l
    ETUnits(l+1)                        = 840 + l
    obukovUnits(l+1)                    = 850 + l
    scalarStarUnits(l+1)                = 860 + l
    surfaceFluxUnits(l+1)               = 870 + l
enddo
!.....

temperatureFlag = 0
moistureIndex   = 0
do l=1,scalarCount
    if (scalarFlags(l) == 1)then
        temperatureFlag = 1
    elseif (scalarFlags(l) == 2)then
        moistureIndex = l
    endif
enddo
if (scalarCount >= 1)then
    loopCount = scalarCount
else
    loopCount = 1
endif


frame_cnt=0
cst=0

! cc Create Subarrays for MPI I/O

! momentum
call MPI_TYPE_CREATE_SUBARRAY(2,(/nx,ny/),(/nx,nyb/), &
(/0,(me-hfact)*nyb/),MPI_ORDER_FORTRAN, &
MPI_DOUBLE_PRECISION,fileview_2D,ierr)

call MPI_TYPE_CREATE_SUBARRAY(3,(/nx,ny,nz/), &
(/nx,nyb,nzb/),(/0,(me-hfact)*nyb,vfact*nzb/), &
MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,fileview_3D, &
ierr)
call MPI_TYPE_COMMIT(fileview_2D,ierr)
call MPI_TYPE_COMMIT(fileview_3D,ierr)

! soil
if(scalarCount > 0)then
    if(surfaceFlags(1) == 2 )then
        call MPI_TYPE_CREATE_SUBARRAY(3,(/nx,ny,soilLevels/), &
        (/nx,nyb,soilLevels/), &
        (/0,(me-hfact)*nyb,0/), &
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, &
        fileview_soil,ierr)
        call MPI_TYPE_COMMIT(fileview_soil,ierr)
    endif
endif

call MPI_BARRIER(nall,ierr)  !syncronize before I/O

IF(initu == 0)THEN

    nrsub = 0
    frame_cnt = 0
    tstat = 0

    if(me == nprocs-1)then
        open(unit=334,file='checkpoints/break',form='formatted')
        write(334,*)0
        write(334,*)0
        write(334,*)0
        close(334)
    endif
    call MPI_BARRIER(nall,ierr)

! cc  Read INITIAL u,v,w Velocity Fields

! cc  Read INITIAL u,v,w Velocity Fields

    view_disp=0
    if(me == 0)write(*,*)'start reading u,v,w'
    call read3D(u,1,1,2,'input/u.ini',fileview_3D)
    call read3D(v,1,1,2,'input/v.ini',fileview_3D)
    call read3D(w,1,1,2,'input/w.ini',fileview_3D)

! on-dimensionalize
    u=u/u_star
    v=v/u_star
    w=w/u_star

! allilean transformation
    u=u-Ugal
    v=v-Vgal

! cc  Read INITIAL Momentum Surface Files
    IF(vfact == 0 .AND. verticalBC == 0)THEN

        call MPI_FILE_OPEN(MPI_COMM_LEVEL,'input/zo.ini', &
        MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
        fileview_2D,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_READ_ALL(fh,zo(1,1), &
        Nx*Nyb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)
    ! on-dimensionalize
        zo=zo/z_i

    ENDIF

! cc  Read INITIAL Particle Release Positions
    if(npart > 0)then
        call MPI_FILE_OPEN(nall,'./input/release_pos.ini', &
        MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_READ_ALL(fh, &
        release_pos(1,1),nr*3,MPI_DOUBLE_PRECISION, &
        MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)
        release_pos=release_pos/z_i
    endif

! cc  Read Canopy Field
    if(c_flag == 1)then
        call read3D(LAD_uvp,1,1,2,'input/PlantDensity.ini',fileview_3D)

    ! on-Dimensionalize leaf density
        LAD_uvp = LAD_uvp*z_i

        call update1(LAD_uvp)

    ! if the canopy is not in my domain, don't calculate drag
        if(sum(LAD_uvp) == 0.d0)c_flag=0

        if(npart > 0 .AND. deposition == 1 .AND. me == 0)then

        ! lant area fractions
            call MPI_FILE_OPEN(MPI_COMM_SELF,'input/fx.ini', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_READ(fh,fx(1),Nz,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
            call MPI_FILE_OPEN(MPI_COMM_SELF,'input/fz.ini', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_READ(fh,fz(1),Nz,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

        endif

    endif
    if(me == 0) then
        write(*, *) "nx: ", Nx
        write(*, *) "nyb: ", Nyb
        write(*, *) "nzb: ", Nzb
        write(*, *) "nz2: ", Nz2
    end if
ELSEIF(initu == 1)THEN

    if(me == 0)then
        open(unit=111,file='checkpoints/break',form='formatted')
        read(111,*)nrsub
        read(111,*)frame_cnt
        read(111,*)tstat
        close(unit=111)
        write(*,*)'Reading Break Files'
    endif
    call MPI_BCAST(nrsub,1,MPI_INTEGER,0,nall,ierr)
    call MPI_BCAST(frame_cnt,1,MPI_INTEGER,0,nall,ierr)
    call MPI_BCAST(tstat,1,MPI_INTEGER,0,nall,ierr)
    
    call MPI_BARRIER(nall,ierr)


! cc  Read BREAK u,v,w Velocity Fields


    view_disp=0
    if(me == 0) write(*,*)'reading break files:'
    call read3D(u,1,1,2,'checkpoints/u.break',fileview_3D)
    call read3D(v,1,1,2,'checkpoints/v.break',fileview_3D)
    call read3D(w,1,1,2,'checkpoints/w.break',fileview_3D)



! non-dimensionalize
    u=u/u_star
    v=v/u_star
    w=w/u_star

! gallilean transformation
    u=u-Ugal
    v=v-Vgal

! cc  Read INITIAL Momentum Surface Files
    IF(vfact == 0 .AND. verticalBC == 0)THEN

        call MPI_FILE_OPEN(MPI_COMM_LEVEL,'input/zo.ini', &
        MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
        fileview_2D,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_READ_ALL(fh,zo(1,1), &
        Nx*Nyb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)

    ! non-dimensionalize
        zo=zo/z_i

    ENDIF

! cc  Read INITIAL Particle Release Positions
    if(npart > 0)then
        call MPI_FILE_OPEN(nall,'./input/release_pos.ini', &
        MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_READ_ALL(fh, &
        release_pos(1,1),nr*3,MPI_DOUBLE_PRECISION, &
        MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)
        release_pos=release_pos/z_i
    endif

! cc  Read BREAK Canopy Field
    if(c_flag == 1)then
    ! ote that PlantDensity.ini is on w-nodes
        call read3D(LAD_uvp,1,1,2,'input/PlantDensity.ini',fileview_3D)

    ! on-Dimensionalize leaf density
        LAD_uvp = LAD_uvp*z_i

        call update1(LAD_uvp)

    ! f the canopy is not in my domain, don't calculate drag
        if(sum(LAD_uvp) == 0.d0)c_flag=0

        if(npart > 0 .AND. deposition == 1 .AND. me == 0)then
        ! lant area fractions
            call MPI_FILE_OPEN(MPI_COMM_SELF,'input/fx.ini', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_READ(fh,fx(1),Nz,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
            call MPI_FILE_OPEN(MPI_COMM_SELF,'input/fz.ini', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_READ(fh,fz(1),Nz,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
        endif

    endif

! cc  Read BREAK RHS Momentum
    view_disp = 0
    call read3D(RHSx,1,1,2,'checkpoints/RHS_momX.break', &
    fileview_3D)
    call read3D(RHSy,1,1,2,'checkpoints/RHS_momY.break', &
    fileview_3D)
    call read3D(RHSz,1,1,2,'checkpoints/RHS_momZ.break', &
    fileview_3D)


    if(model == 2)then
    ! cc      Read BREAK a1_old-b1_old

        call read3D(tmp_old,1,1,1,'checkpoints/a1.break', &
        fileview_3D)
        a1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/b1.break', &
        fileview_3D)
        b1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old


    elseif(model == 3 .AND. averaging == 1)then
    !!      Read BREAK a1_old-e1_old and a2_old-e2_old

        call read3D(tmp_old,1,1,1,'checkpoints/a1.break', &
        fileview_3D)
        a1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old       
        call read3D(tmp_old,1,1,1,'checkpoints/b1.break', &
        fileview_3D)
        b1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/c1.break', &
        fileview_3D)
        c1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/d1.break', &
        fileview_3D)
        d1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/e1.break', &
        fileview_3D)
        e1_old(:, 2:nyb+1, 2:nzb+1) = tmp_old

        call read3D(tmp_old,1,1,1,'checkpoints/a2.break', &
        fileview_3D)
        a2_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/b2.break', &
        fileview_3D)
        b2_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/c2.break', &
        fileview_3D)
        c2_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/d2.break', &
        fileview_3D)
        d2_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
        call read3D(tmp_old,1,1,1,'checkpoints/e2.break', &
        fileview_3D)
        e2_old(:, 2:nyb+1, 2:nzb+1) = tmp_old
    endif

    if(npart > 0)then
    ! cc  Read BREAK Particle Positions

        open(unit=390,file= &
        'checkpoints/ParticlePositions_nomodel.break', &
        status='old')
        read(390,*)ipart
        do i=1,ipart
            read(390,*)particle(i,1,1),particle(i,2,1), &
            particle(i,3,1),particle(i,4,1), &
            particle(i,5,1),particle(i,6,1)
        enddo
        close(unit=390)
        if(part_model == 2)then
            open(unit=390,file= &
            'checkpoints/ParticlePositions_weiliso.break', &
            status='old')
            read(390,*)ipart
            do i=1,ipart
                read(390,*)particle(i,1,2),particle(i,2,2), &
                particle(i,3,2),particle(i,4,2), &
                particle(i,5,2),particle(i,6,2)
            enddo
            close(unit=390)
        endif

        if(part_model >= 2)then
        ! cc  Read BREAK SGS Particle Velocity
            call MPI_FILE_OPEN(nall, &
            'checkpoints/SGS_velocity.break', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_READ_ALL(fh,us_f(1,1,1),npart*3*2, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

        ! cc  Read BREAK Particle ESGS
            call MPI_FILE_OPEN(nall,'checkpoints/ESGS_f.break', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_READ_ALL(fh,ESGS_f(1,1,2),nx*nyb*nzb, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
        endif

    endif


ENDIF                    !end IF(initu == 1)

! cc  Scalar Inputs Loop

IF (scalarCount >= 1) then
    scalarIter = 0
    do l=1,scalarCount
        if (scalarFlags(l) == 1)then ! scalar is temperature
            write(fileStr1, '(A)') &
            'checkpoints/temperature.break'
            write(fileStr21, '(A)') &
            'checkpoints/a4_temperature.break'
            write(fileStr22, '(A)') &
            'checkpoints/b4_temperature.break'
            write(fileStr23, '(A)') &
            'checkpoints/c4_temperature.break'
            write(fileStr24, '(A)') &
            'checkpoints/d4_temperature.break'
            write(fileStr25, '(A)') &
            'checkpoints/e4_temperature.break'
            write(fileStr31, '(A)') &
            'checkpoints/a8_temperature.break'
            write(fileStr32, '(A)') &
            'checkpoints/b8_temperature.break'
            write(fileStr33, '(A)') &
            'checkpoints/c8_temperature.break'
            write(fileStr34, '(A)') &
            'checkpoints/d8_temperature.break'
            write(fileStr35, '(A)') &
            'checkpoints/e8_temperature.break'
            write(fileStr4, '(A)') &
            'checkpoints/RHS_temperature.break'
            write(fileStr5, '(A)')'checkpoints/coolrate.break'
            write(fileStr6, '(A)') &
            'checkpoints/surfaceTemperature.break'
        elseif ( scalarFlags(l) == 2)then ! scalar is moisture
            write(fileStr1, '(A)')'checkpoints/moisture.break'
            write(fileStr21, '(A)') &
            'checkpoints/a4_moisture.break'
            write(fileStr22, '(A)') &
            'checkpoints/b4_moisture.break'
            write(fileStr23, '(A)') &
            'checkpoints/c4_moisture.break'
            write(fileStr24, '(A)') &
            'checkpoints/d4_moisture.break'
            write(fileStr25, '(A)') &
            'checkpoints/e4_moisture.break'
            write(fileStr31, '(A)') &
            'checkpoints/a8_moisture.break'
            write(fileStr32, '(A)') &
            'checkpoints/b8_moisture.break'
            write(fileStr33, '(A)') &
            'checkpoints/c8_moisture.break'
            write(fileStr34, '(A)') &
            'checkpoints/d8_moisture.break'
            write(fileStr35, '(A)') &
            'checkpoints/e8_moisture.break'
            write(fileStr4, '(A)') &
            'checkpoints/RHS_moisture.break'
            write(fileStr5, '(A)') &
            'checkpoints/dsdt_Moisture.break'
            write(fileStr6, '(A)') &
            'checkpoints/surfaceMoisture.break'
        else             ! = 0, scalar is passive
            scalarIter = scalarIter + 1
            write(fileStr1, '(A,I1,A)') &
            'checkpoints/scalar.break',scalarIter,'.break'
            write(fileStr21, '(A,I1,A)') &
            'checkpoints/a4_scalar.break',scalarIter,'.break'
            write(fileStr22, '(A,I1,A)') &
            'checkpoints/b4_scalar.break',scalarIter,'.break'
            write(fileStr23, '(A,I1,A)') &
            'checkpoints/c4_scalar.break',scalarIter,'.break'
            write(fileStr24, '(A,I1,A)') &
            'checkpoints/d4_scalar.break',scalarIter,'.break'
            write(fileStr25, '(A,I1,A)') &
            'checkpoints/e4_scalar.break',scalarIter,'.break'
            write(fileStr31, '(A,I1,A)') &
            'checkpoints/a8_scalar.break',scalarIter,'.break'
            write(fileStr32, '(A,I1,A)') &
            'checkpoints/b8_scalar.break',scalarIter,'.break'
            write(fileStr33, '(A,I1,A)') &
            'checkpoints/c8_scalar.break',scalarIter,'.break'
            write(fileStr34, '(A,I1,A)') &
            'checkpoints/d8_scalar.break',scalarIter,'.break'
            write(fileStr35, '(A,I1,A)') &
            'checkpoints/e8_scalar.break',scalarIter,'.break'
            write(fileStr4, '(A,I1,A)') &
            'checkpoints/RHS_scalar.break',scalarIter,'.break'
            write(fileStr5, '(A,I1,A)') &
            'checkpoints/dsdt_scalar.break',scalarIter,'.break'
            write(fileStr6, '(A,I1,A)') &
            'checkpoints/surfaceScalar',scalarIter,'.break'
        endif

    ! cc  Read Scalar(s)

        call MPI_FILE_OPEN(nall,fileStr1,MPI_MODE_RDONLY, &
        MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_SET_VIEW(fh,view_disp, &
        MPI_DOUBLE_PRECISION,fileview_3D,"native", &
        MPI_INFO_NULL,ierr)
        call MPI_FILE_READ_ALL(fh,scalar(1,1,2,l), &
        Nx*Nyb*Nzb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
        ierr)
        call MPI_FILE_CLOSE(fh,ierr)

    ! cc  Read Scalar Surface Files

        if(vfact == 0 .AND. surfaceFlags(l) == 1)then
            call MPI_FILE_OPEN(MPI_COMM_LEVEL,fileStr5, &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp, &
            MPI_DOUBLE_PRECISION,fileview_2D, &
            "native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,coolrate(1,1,l),Nx*Nyb, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

            call MPI_FILE_OPEN(MPI_COMM_LEVEL,fileStr6, &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp, &
            MPI_DOUBLE_PRECISION,fileview_2D, &
            "native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,surfaceScalar(1,1,l),Nx*Nyb, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
        endif

        if(inits == 1)then

            if(model == 2)then
            ! cc            Read a4-b4

                call read3D(tmp_old,1,1,1,fileStr21,fileview_3D)
                a4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr22,fileview_3D)
                b4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old


            elseif(model == 3 .AND. averaging == 1)then
            ! cc            Read a4-e4 and a8-e8

                call read3D(tmp_old,1,1,1,fileStr21,fileview_3D)
                a4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr22,fileview_3D)
                b4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr23,fileview_3D)
                c4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr24,fileview_3D)
                d4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr25,fileview_3D)
                e4_old(:,2:nyb+1,2:nzb+1,l)=tmp_old

                call read3D(tmp_old,1,1,1,fileStr31,fileview_3D)
                a8_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr32,fileview_3D)
                b8_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr33,fileview_3D)
                c8_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr34,fileview_3D)
                d8_old(:,2:nyb+1,2:nzb+1,l)=tmp_old
                call read3D(tmp_old,1,1,1,fileStr35,fileview_3D)
                e8_old(:,2:nyb+1,2:nzb+1,l)=tmp_old

            endif

        ! cc  Read Scalar RHS

            call MPI_FILE_OPEN(nall,filestr4,MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp, &
            MPI_DOUBLE_PRECISION,fileview_3D,"native", &
            MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,scalarRHS(1,1,2,l),Nx*Nyb*Nzb, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

        endif
    enddo
endif

! cc Read INITIAL SEB files
call MPI_BARRIER(nall,ierr)

if(scalarCount > 0)then
    if(surfaceFlags(1) == 2 )then
        IF(vfact == 0)THEN

        ! ote: using porosity as temporary variable
            call MPI_FILE_OPEN(MPI_COMM_LEVEL, &
            'input/soilTemperature.ini',MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,porosity(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
            gndScalars(:,:,:,1) = porosity

            call MPI_FILE_OPEN(MPI_COMM_LEVEL, &
            'input/soilMoisture.ini',MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,porosity(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
            gndScalars(:,:,:,2) = porosity

            call MPI_FILE_OPEN(MPI_COMM_LEVEL, &
            'input/soilTypeParams.ini',MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,porosity(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            view_disp=8*Nx*Ny*soilLevels
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,satPotential(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            view_disp=2*8*Nx*Ny*soilLevels
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,satHydrCond(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            view_disp=3*8*Nx*Ny*soilLevels
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,soilExponent(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            view_disp=4*8*Nx*Ny*soilLevels
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_soil,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,heatCapSoil(1,1,1), &
            Nx*Nyb*soilLevels,MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
            view_disp=0

            call MPI_FILE_OPEN(MPI_COMM_LEVEL, &
            'input/albedo.ini',MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fh,ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
            fileview_2D,"native",MPI_INFO_NULL,ierr)
            call MPI_FILE_READ_ALL(fh,albedo(1,1), &
            Nx*Nyb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

            if(albedoFlag == 1)then
                call MPI_FILE_OPEN(MPI_COMM_LEVEL, &
                'input/minAlbedo.ini',MPI_MODE_RDONLY, &
                MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_SET_VIEW(fh,view_disp, &
                MPI_DOUBLE_PRECISION,fileview_2D,"native", &
                MPI_INFO_NULL,ierr)
                call MPI_FILE_READ_ALL(fh,minAlbedo(1,1), &
                Nx*Nyb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                ierr)
                call MPI_FILE_CLOSE(fh,ierr)
            endif

        ENDIF

        call MPI_FILE_OPEN(nall,'input/zGnd.ini',MPI_MODE_RDONLY, &
        MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_READ_ALL(fh,zGnd(1),soilLevels, &
        MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)
        zGnd=zGnd/z_i

    endif
endif

! cc me=0 Opens Statistic Output Files
if(me == 0)then

    call openfiles(dsdxUnits,dsdyUnits,dsdzUnits, &
    scalarMeanUnits,scalar2Units,scalar3Units, &
    xScalarFluxUnits,yScalarFluxUnits,zScalarFluxUnits, &
    awsUnits,ausUnits,avsUnits,scalarSpectraUnits, &
    scalarPrandtlUnits,cs2prUnits,beta2Units,ETUnits, &
    obukovUnits,scalarStarUnits,surfaceFluxUnits)

endif

! cc  All procs have to run through this initialization for statistics


call zeroslice(au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy, &
atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz,adwdx, &
e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL)

if(scalarCount >= 1)then
    call s_zeroslice(scalarMean,t2,t3,asgs_t1,asgs_t2,asgs_t3,aut, &
    avt,awt,adsdx,adsdy,adsdz,aPr,aCs2Pr,abeta2,aqz_s,aET)
endif

if(me == nprocs-1) write(*,*) 'Previous timesteps = ',nrsub
if(me == nprocs-1) write(*,*) 'Previous frames = ',frame_cnt
if(me == nprocs-1) write(*,*) 'Last stat output = ',tstat

! cc Initialize Counter Variables
mfst=1
vt=0

call MPI_BARRIER(nall,ierr)
if(me == 1) write(*, *) 'dx:', dx, 'dy:',dy,'dz:',dz
test_u=1
test_s=0
!!!!INPUTS CHECK!!!!!!!!!!!!!!!!!!!!!!!
if(me == 0 .AND. test_u == 1)then
    write(*,*)me,sum(u(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'u_o'
    write(*,*)me,sum(abs(v(:,:,2:nzb+1)))/(Nx*Nyb*Nzb),'v_o'
    write(*,*)me,sum(abs(w(:,:,2:nzb+1)))/(Nx*Nyb*Nzb),'w_o'
    if(verticalBC == 0)then
        write(*,*)me,sum(zo(:,:)/(Nx*Nyb))*z_i,'zo_o'
    endif
endif
if(me == 0 .AND. scalarcount > 0)then
    write(*,*)me,sum(scalar(:,:,2:nzb+1,1))/(Nx*Nyb*Nzb),'scalar_o'
    write(*,*)me,sum(coolrate(:,:,1))/(Nx*Nyb),'coolrate_o'
    write(*,*)me,sum(surfaceScalar(:,:,1))/(Nx*Nyb), &
    'surfaceTemp_o'
endif

!!!!INPUTS CHECK!!!!!!!!!!!!!!!!!!!!!!!

!     ================================================S
!     ================================================S
!     ================================================S
!     ================================================S
!     ================================================S
!     ================================================S
!     ... Start Time Loop      =======================S
!     ... Start Time Loop      =======================S
!     ... Start Time Loop      =======================S
!     ================================================S
!     ================================================S
!     ================================================S
!     ================================================S
!     ================================================S
!     ================================================S

call MPI_BARRIER(nall,ierr)  !syncronize before time loop

timeStep:   DO  ttt=nrsub+1,Nsteps

    t=ttt-nrsub
    IF(me == 0 .AND. t == 2) t01=MPI_WTIME()

    if(me == 0)write(*,*)'t=',ttt

!...  Truncate wavenumbers and calculate d_dx & d_dy

    call FILT_DA(u,dudx,dudy)

    call FILT_DA(v,dvdx,dvdy)

    if(scalarCount >= 1)then
        call FILT_DA(w,dwdx,dwdy)
        if(t == 1)then
            do l=1,scalarCount
                call FILT_DA(scalar(:,:,:,l),dsdx(:,:,:,l), &
                dsdy(:,:,:,l))
            enddo
        endif
    else
        call FILT_DA(w,dwdx,dwdy)
    endif

!     mpi update ghostlayers

    call update3(u,v,w)

    if(t == 1)then
        if(scalarCount >= 1)then
            do l=1,scalarCount
                call update1(scalar(:,:,:,l))
            enddo
        endif
    endif

!     ... Save previous time's right-hand-sides for Adams-Bashforth
!     ... Integration (In subroutine "STEP" use first order time
!     ... advancement on first time step)

    RHSx_f=RHSx
    RHSy_f=RHSy
    RHSz_f=RHSz

!     ... Compute spatial derivatives from u,v,w (phys. space) and
!     ... place in phys space.

    call DDZ_UV_p(dudz,u)
    call DDZ_UV_p(dvdz,v)
    call DDZ_W(dwdz,w)

    IF (scalarCount >= 1) then
        do l=1,scalarCount
            call ddz_uv_p(dsdz(:,:,:,l),scalar(:,:,:,l))
        enddo
    ENDIF

!     ... MPI UPDATE THE GHOSTLAYERS for the derivatives
!     ... (this is important for the SGS model)

    call update9(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)

    if(scalarCount >= 1)then
        do l=1,scalarCount
            call update3(dsdx(:,:,:,l),dsdy(:,:,:,l),dsdz(:,:,:,l))
        enddo
    endif

!!!!GRADIENTS CHECK!!!!!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum(dudx(:,:,3:nzb+1))/(Nx*Nyb*Nzb), &
        sum(dudy(:,:,3:nzb+1))/(Nx*Nyb*Nzb),'horz gradients'
        write(*,*)sum(dudz(:,:,3:nzb+1))/(Nx*Nyb*Nzb), &
        sum(dwdz(:,:,3:nzb+1))/(Nx*Nyb*Nzb), &
        'vert gradients'
        if(test_s == 1)then
            do l=1,scalarCount
                write(*,*)sum(dsdx(:,:,3:nzb+1,l))/(Nx*Nyb*Nzb), &
                sum(dsdz(:,:,3:nzb+1,l))/(Nx*Nyb*Nzb), &
                'scalar gradients'
            enddo
        endif
    endif
!!!!GRADIENTS CHECK!!!!!!!!!!!!!!!!!!!!!!

!     Compute Surface Boundary Conditions*************************
    IF (vfact == 0 .AND. verticalBC == 0) THEN
        if(scalarCount > 0)then
            if( surfaceFlags(1) == 2 )then
                UTC = startUTC + (t + nrsub)*dt

                call Filter_2dsl(u_hat,u(:,:,2))
                call Filter_2dsl(v_hat,v(:,:,2))
                do l = 1,scalarCount
                    call Filter_2dsl(scalar_hat(:,:,l), &
                    scalar(:,:,2,l))
                enddo

                call galileanGroundShift(gndScalars,soilHeatFlux, &
                ustar,scalarFlux,coolrate,zo, &
                Psi,Psi0,fi,fi_H,ilow,jlow)

                do j = 1,Nyb
                    do i = 1,Nx

                        Uref    = ((u_hat(i,j)+Ugal)**2 &
                        + (v_hat(i,j)+Vgal)**2)**0.5d0
                        M(i,j) = Uref

                        call solveGroundBC(Uref,scalar_hat(i,j,:), &
                        scalar,u,v,gndScalars(i,j,:,:),ustar(i,j), &
                        scalarFlux(i,j,:),soilHeatFlux(i,j), &
                        porosity(i,j,:),satPotential(i,j,:), &
                        satHydrCond(i,j,:),soilExponent(i,j,:), &
                        heatCapSoil(i,j,:),albedo(i,j), &
                        minAlbedo(i,j),measRad,netRad,i,j, &
                        Psi(i,j,:),Psi0(i,j,:),fi(i,j,:), &
                        fi_H(i,j,:),zo(i,j),obukhovL(i,j),zGnd)

                    enddo
                enddo

            else

                do l=1,scalarCount
                    call surf_flux(scalar(:,:,:,l),thetaFactor,u,v, &
                    scalarFlux(:,:,l),Psi,Psi0,fi,fi_H,zo, &
                    surfaceScalar(:,:,l),coolrate(:,:,l),l)
                enddo

            endif

        else

            if (mod(t,c_count) == 0) then
                ssr = floor((dble(t+nrsub)*dt*Ugal)/dx)
                ssi = int(ssr)
                iii = mod(ssi,nx)

                do i=1,nx
                    if(i-iii < 1)then
                        ilow(i)=i-iii+nx
                    else
                        ilow(i)=i-iii
                    endif
                enddo
            endif

            Psi(:,:,1:2)=0.
            Psi0(:,:,1:2)=0.
            fi(:,:,1:2)=1.
            fi_H(:,:,1:2)=1.

        endif

        call wallstress2(txz,tyz,u,v,Psi,Psi0,zo,ustar,M)

        call derivwall2 (dsdz,dudz,dvdz,u,v,fi,fi_H,scalarFlux, &
        ustar,M)

    !!!!  SURFACE CONDITIONS CHECK!!!!!!!!!!!!!!!!!!!!!
        if(me == 0 .AND. test_u == 1)then
            write(*,*)sum(txz(:,:,2))/(Nx*Nyb), &
            sum(tyz(:,:,2))/(Nx*Nyb),'wall stress'
            write(*,*)sum(dudz(:,:,2))/(Nx*Nyb), &
            sum(dvdz(:,:,2))/(Nx*Nyb), &
            'wall gradients'
            write(*,*)sum(ustar(:,:))/(Nx*Nyb),'ustar'
        endif
    !              if(me.eq.0.and.test_s.eq.1)then
    !                 do l=1,scalarCount
    !                    write(*,*)sum(dsdz(:,:,2,l))/(Nx*Nyb),
    !     +                   'scalar wall gradients'
    !                 enddo
    !              endif
    !!!!SURFACE CONDITIONS CHECK!!!!!!!!!!!!!!!!!!!!!

    ENDIF

!         Compute Statistics every c_count ****************************
    IF (mod(ttt,c_count) == 0) then

        rmsdivve=0.0
        kee=0.0
        ke=0.0

        if(vfact == 0) then
            kstart=3
        else
            kstart=2
        end if

        do k=kstart,Nzb+1
            do j=1,Nyb
                do i=1,Nx
                    ke=ke+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+ &
                    w(i,j,k)*w(i,j,k)
                end do
            end do
        end do

        ke=ke/(1.d0*(nx*ny*(nz-1)))

        call rmsdiv (rmsdivvel,dudx,dvdy,dwdz)

        call MPI_BARRIER(nall,ierr)
        call MPI_REDUCE(rmsdivvel,rmsdivve,1,MPI_DOUBLE_PRECISION, &
        MPI_SUM,0,nall,ierr)
        call MPI_REDUCE(ke,kee,1,MPI_DOUBLE_PRECISION, &
        MPI_SUM,0,nall,ierr)

        if(me == 0)then
            write (11,*) ttt*dt,dt,kee,rmsdivve
        endif

        call MPI_Bcast(ilow,nx,MPI_INTEGER,0,nall,ierr)


        call avgslice(u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz, &
        dudx,dvdz,dwdz,dwdx,Cs2,beta1,ESGS3D,DSGS3D,TL,au,av, &
        aw,ap,u2,v2,w2,p2,w3,atxx,atxz,atyy,atyz,atzz,atxy,auw, &
        avw,auv,adudz,adudx,advdz,adwdz,adwdx,e,aCs2,aCs, &
        abeta1,atxz_s,aESGS,aDSGS,aTL,ilow,wgx)

        If (scalarCount >= 1) then
            do l=1,scalarCount
                call update1(Pr2(:,:,:,l))
                call scalar_slice(u,v,w,scalar(:,:,:,l), &
                qx(:,:,:,l),qy(:,:,:,l),qz(:,:,:,l), &
                dsdx(:,:,:,l),dsdy(:,:,:,l),dsdz(:,:,:,l), &
                Pr2(:,:,:,l), Cs2, beta2(:,:,:,l), &
                scalarMean(:,:,l),t2(:,:,l),t3(:,:,l), &
                asgs_t1(:,:,l),asgs_t2(:,:,l),asgs_t3(:,:,l), &
                aut(:,:,l),avt(:,:,l),awt(:,:,l), &
                adsdx(:,:,l),adsdy(:,:,l),adsdz(:,:,l), &
                aPr(:,:,l),aCs2Pr(:,:,l),abeta2(:,:,l), &
                aqz_s(:,:,l),ESCL3D(:,:,:,l),aET(:,:,l), &
                ilow,wgx)
            enddo

        ENDIF

    ENDIF

!     ... Compute Convective term in physical space.

    call dealias1(u,u_m)
    call dealias1(v,v_m)
    call dealias1(w,w_m)

    call update3(u_m,v_m,w_m)

    Call CONVEC(Cx,Cy,Cz,u_m,v_m,w_m,dudy,dudz,dvdx,dvdz,dwdx,dwdy)

!!!!CONVEC CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum(Cx(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(Cy(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(Cz(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'convection'
    endif
!!!!CONVEC CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!

    Call SGS_STAG(qx,qy,qz,u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz, &
    dwdx,dwdy,dwdz,dsdx,dsdy,dsdz,pt,Cs2,Cs2_m,Pr2,Pr2_m, &
    cst,txx,txy,txz,tyy,tyz,tzz,scalar,a1_old,b1_old, &
    c1_old,d1_old,e1_old,a2_old,b2_old,c2_old,d2_old,e2_old, &
    a4_old,b4_old,c4_old,d4_old,e4_old,a8_old, b8_old, &
    c8_old,d8_old,e8_old,beta1,beta2,ESGS3D,DSGS3D,ESCL3D, &
    TL)

!!!!SGS CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum(txx(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(txy(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(txz(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'tx'
        write(*,*)sum(tyy(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(tyz(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        'ty'
        write(*,*)sum(tzz(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'tz'
    endif
    if(me == 0 .AND. test_s == 1)then
        write(*,*)sum(Pr2(:,:,2:nzb+1,:)),'Pr2'
    endif
!!!!SGS CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ... Compute Scalar RHS

    IF (scalarCount >= 1) then
        do l=1,scalarCount
            call update1(qz(:,:,:,l))
        enddo

        scalarRHSf=scalarRHS

        do l=1,scalarCount
            call Scalar_RHS(scalar(:,:,:,l),u_m,v_m,w_m, &
            qx(:,:,:,l),qy(:,:,:,l),qz(:,:,:,l),txz,tyz, &
            dsdx(:,:,:,l),dsdy(:,:,:,l),dsdz(:,:,:,l), &
            scalarRHS(:,:,:,l),scalarFlux(:,:,l))
        enddo

        IF (t == 1 .AND. inits == 0) scalarRHSf=scalarRHS

    ! cc  Scalar Sponge Layer

        if(sponge /= 0)then
            call sponge_layer(scalar,u,v,w,scalarRHS,RHSx,RHSy,RHSz, &
            rdmp,2)
        endif

    ! cc  Scalar Advection

        call scalar_advection(scalarRHS)

    ! cc  Time Advance Scalar

        do l=1,scalarCount
            call Step_S (scalar(:,:,:,l), scalarRHS(:,:,:,l), &
            scalarRHSf(:,:,:,l),l)

            call FILT_DA(scalar(:,:,:,l),dsdx(:,:,:,l), &
            dsdy(:,:,:,l))

            call update1(scalar(:,:,:,l))

        enddo

    ! cc  Calculate Buoyancy Term

        call calcbeta (scalar,beta)

    ENDIF

!     ... Compute Divergence of SGS shear stresses

!     mpi  possible to calculate ddxytxx, txy in sgs_stag

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    call ddz_uv(ddtzz,tzz)
    call update1(ddtzz)
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    call divstress_UV (divtx, txx, txy, txz)
    call divstress_UV (divty, txy, tyy, tyz)
    call divstress_W (divtz, txz, tyz, ddtzz)

!!!!DIVSTRESS CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum(divtx(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'divtx'
        write(*,*)sum(divtz(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'divtz'
    endif
!!!!DIVSTRESS CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ... Calculate the proper Ugeo for forcing
    if(press_cor == 1)then
        Ugeo(:)=Ug
        Vgeo(:)=Vg
    elseif(press_cor == 2)then
        Ptime = startUTC+dble(ttt)*dt*z_i/u_star/3600.d0
        if(Ptime < 3.d0)then
            Ugeo_g=-6.5d0/u_star+(Ptime+1.d0)/4.d0* &
            (-5.d0+6.5d0)/u_star
            Vgeo_g=4.5d0/u_star+(Ptime+1.d0)/4.d0* &
            (4.5d0-4.5d0)/u_star
        elseif(Ptime < 6.d0)then
            Ugeo_g=-5.d0/u_star+dmod(Ptime,3.d0)/3.d0* &
            (-5.d0+5.d0)/u_star
            Vgeo_g=4.5d0/u_star+dmod(Ptime,3.d0)/3.d0* &
            (4.5d0-4.5d0)/u_star
        else
            Ugeo_g=-5.d0/u_star+dmod(Ptime,6.d0)/6.d0* &
            (-6.5d0+5.d0)/u_star
            Vgeo_g=4.5d0/u_star+dmod(Ptime,6.d0)/6.d0* &
            (2.5d0-4.5d0)/u_star
        endif

        do k=2,nzb+1
            ztemp=(vfact*nzb+k-1)*dz*z_i-dz*z_i/2.d0
            Ugeo(k)=Ugeo_g*((Ug/Ugeo_g-1.d0)* &
            ztemp/2000.d0+1.d0)
            Vgeo(k)=Vgeo_g*((Vg/Vgeo_g-1.d0)* &
            ztemp/2000.d0+1.d0)
        enddo
    endif

!     ... Calculate dynamic tendencies for momentum
    if(M_advec == 0)then
        Uadv=0.d0
        Vadv=0.d0
    else
        Ptime = startUTC+dble(ttt)*dt*z_i/u_star/3600.d0
        if(Ptime < 3)then
            Uadv_d = 0.0005/(u_star*u_star/z_i)
            Vadv_d = 0.0000/(u_star*u_star/z_i)
        else
            Uadv_d = 0.0000/(u_star*u_star/z_i)
            Vadv_d = 0.0000/(u_star*u_star/z_i)
        endif

        do k=2,nzb+1
            ztemp=(vfact*nzb+k-1)*dz*z_i-dz*z_i/2.d0
            if(ztemp >= 200)then
                Uadv(k)=Uadv_d
                Vadv(k)=Vadv_d
            else
                Uadv(k)=Uadv_d*ztemp/200.d0
                Vadv(k)=Vadv_d*ztemp/200.d0
            endif
        enddo
    endif

!     ... Compute preliminary RHS matrices for pressure calculation

    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                if (press_cor == 0 .OR. press_cor == 3) then
                    RHSx(i,j,k)=-Cx(i,j,k)-divtx(i,j,k)+Uadv(k)
                    RHSy(i,j,k)=-Cy(i,j,k)-divty(i,j,k)+Vadv(k)
                elseif(press_cor == 1 .OR. press_cor == 2)then
                    RHSx(i,j,k)=-Cx(i,j,k)-divtx(i,j,k) &
                    -fc*(Vgeo(k)-Vgal-v(i,j,k))+Uadv(k)
                    RHSy(i,j,k)=-Cy(i,j,k)-divty(i,j,k) &
                    +fc*(Ugeo(k)-Ugal-u(i,j,k))+Vadv(k)
                end if
                if (temperatureFlag == 1) then! temperature is computed
                    RHSz(i,j,k)=-Cz(i,j,k)-divtz(i,j,k)+beta(i,j,k)
                else
                    RHSz(i,j,k)=-Cz(i,j,k)-divtz(i,j,k)
                end if
            end do
        end do
    end do

! cc  Momentum Sponge Layer

    if(sponge /= 0 .AND. scalarCount >= 1)then
        call sponge_layer(scalar,u,v,w,scalarRHS,RHSx,RHSy,RHSz, &
        rdmp,2)
    endif


!!!!RHS CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum(RHSx(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(RHSy(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(RHSz(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'RHS'
    endif
!!!!RHS CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ... Add drag term for canopy

!$$$         if(vfact.eq.1)then
!$$$         if(me-hfact.eq.0)then
!$$$            do k=1,nzb
!$$$               write(*,*)me,k+vfact*nzb,LAD_uvp(5,5,k+1),Fdx(5,5,k),
!$$$     +              RHSx(5,5,k+1),u(5,5,k+1),v(5,5,k+1)
!$$$            enddo
!$$$         endif
!$$$         endif
!$$$         !write(*,*)me,'stopping'
!$$$         call MPI_BARRIER(nall,ierr)
!$$$         stop

    if(c_flag == 1)then

        call canopy_drag(Fdx,Fdy,Fdz,LAD_uvp,LAD_w,u,v,w,u_w,v_w, &
        mag_w,w_uvp,mag_uvp)

        RHSx(:,:,2:nzb+1) = RHSx(:,:,2:nzb+1) + Fdx
        RHSy(:,:,2:nzb+1) = RHSy(:,:,2:nzb+1) + Fdy
        RHSz(:,:,2:nzb+1) = RHSz(:,:,2:nzb+1) + Fdz

    endif

!     ... Solve Poisson Equation for pressure.

    call Press_STAG(P,RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f, &
    u,v,w,dpdx,dpdy,divtz)
    call update1(P)

!     ... Calculate pressure gradients here

    Call DDZ_UV (dpdz,P)

!!!!PRESSURE CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum(P(:,:,2:nzb+1))/(Nx*Nyb*Nzb),'Pressure'
        write(*,*)sum(dpdx(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(dpdy(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        sum(dpdz(:,:,2:nzb+1))/(Nx*Nyb*Nzb), &
        'Pressure Gradients'
    endif
!!!!PRESSURE CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!

!...  Finalize right hand side for explicit time advancement

    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                RHSx(i,j,k)=  RHSx(i,j,k) - dpdx(i,j,k)
                RHSy(i,j,k)=  RHSy(i,j,k) - dpdy(i,j,k)
                RHSz(i,j,k)=  RHSz(i,j,k) - dpdz(i,j,k)
            end do
        end do
    end do

!     ... Start first step using new C as former C

    IF (t == 1 .AND. initu == 0) then
        RHSx_f=RHSx
        RHSy_f=RHSy
        RHSz_f=RHSz
    ENDIF

    if(press_cor == 0)then

        force_x(2:nzb+1)= f_p*cos(theta_mean_wind*pi/180.d0)
        force_y(2:nzb+1)= f_p*sin(theta_mean_wind*pi/180.d0)


    elseif(press_cor == 1 .OR. press_cor == 2)then

        force_x(2:nzb+1)=0.d0
        force_y(2:nzb+1)=0.d0

    elseif(press_cor == 3)then

    !     ... Compute required mean forcing to keep total x momentum constant.

    !    if(mod(t,press_step) == 0)then

        call pressure_forcing(force_x,force_y,u,v,RHSx,RHSy, &
        RHSx_f,RHSy_f)

        !   endif

    endif

!     ... Step forward in time

    call stepbl_UV(u,RHSx,RHSx_f,force_x)

    Call stepbl_UV (v,RHSy,RHSy_f,force_y)

    Call stepbl_W (w, RHSz,RHSz_f)

!!!!TIME ADVANCEMENT CHECK!!!!!!!!!!!!!!!!!!
    if(me == 0 .AND. test_u == 1)then
        write(*,*)sum((u(:,:,2:nzb+1)+Ugal)*u_star)/(Nx*Nyb*Nzb), &
        'u adv'
        write(*,*)sum((abs(v(:,:,2:nzb+1))+Vgal)*u_star) &
        /(Nx*Nyb*Nzb),'v adv'
        write(*,*)sum(abs(w(:,:,2:nzb+1))*u_star)/(Nx*Nyb*Nzb), &
        'w adv'
        if(scalarCount.ge.1)then
            write(*,*)sum(scalar(:,:,2:nzb+1,1))/(Nx*Nyb*Nzb), &
                'scalar adv'
        endif
                write(*,*)ttt,'time steps completed'
    endif
!!!!TIME ADVANCEMENT CHECK!!!!!!!!!!!!!!!!!!

! cc  Check for NaN
    if(mod(t,p_count) == 0)then
        if(u(5,5,5) /= u(5,5,5))then
            write(*,*)'NaN error'
            write(*,*)'t=',ttt
            write(*,*)me,sum(u(:,:,2:nzb+1)),'u'
            write(*,*)me,sum(v(:,:,2:nzb+1)),'v'
            write(*,*)me,sum(w(:,:,2:nzb+1)),'w'
            stop
        endif
    endif

    IF(npart > 0)THEN

        if(ttt == 1)then
            us_f=0.d0
            ESGS_f=0.d0
            if(part_model == 3)then
                tau_f=0.d0
            endif
            rogue_count=0
            trajectory_updates=0
            rogue_flag(:)=0
        endif
        if(t == 1)then
            idum=-me*10
        endif

    ! release particles
        if(ttt >= start_release .AND. me == 0)then
            if(mod(ttt-start_release,partstep) == 0)then
                if(ipart < npart)then
                    if(t == start_release)ipart=0
                    call point_release(particle,release_pos,zo,idum)
                endif
            endif
            write(*,*)'ipart=',ipart
        endif


    ! update particle positions
        if(ttt >= start_release)then
            if(mod(ttt-start_release,skip_step) == 0)then
                call dispersion(particle,u,v,w,ESGS3D,zo,ustar, &
                LAD_w,us_f,ESGS_f,ESGS_part_f,tau_f,idum,fx,fz, &
                txx,txy,txz,tyy,tyz,tzz,divtx,divty,divtz,TL, &
                rogue_count,trajectory_updates,rogue_flag)
            endif
            if(mod(ttt,10) == 0 .AND. me == 0)then
                write(*,*)'rogue trajectory 1 out of every', &
                dble(trajectory_updates)/dble(rogue_count), &
                'updates'
            endif
        endif


    ENDIF

!     ... Print stats to file


    IF (mod(ttt,p_count) == 0) then

        if(initu == 0 .OR. t > tstat-nrsub)then

            if(me == nprocs-1)then
                open(unit=334,file='checkpoints/break', &
                form='formatted')
                read(334,*)i
                read(334,*)j
                close(334)
                open(unit=334,file='checkpoints/break', &
                form='formatted',status='replace')
                write(334,*)i
                write(334,*)j
                write(334,*)ttt
                close(334)
            endif

            do k=2,nzb+1

                call spectrum (u,k,specu)
                call spectrum (v,k,specv)
                call spectrum (p,k,specp)
                IF (scalarCount >= 1) then
                    do l=1,scalarCount
                        call spectrum (scalar(:,:,:,l),k,spect(:,:,l))
                    enddo
                ENDIF
                call spectrum (w,k,specw)
            end do

            call output_average(specu,299)
            call output_average(specv,298)
            call output_average(specw,297)
            call output_average(specp,295)

            IF (scalarCount >= 1) then
                do l=1,scalarCount
                    call output_average(spect(:,:,l), &
                    scalarSpectraUnits(l))
                enddo
            ENDIF

        ! D averages
            do k=2,nzb+1
                do i=1,nx
                    xz_bar(i,k)=sum(TL(i,:,k))/nyb
                enddo
            enddo
            call output_average(xz_bar/hprocs,320)
            do k=2,nzb+1
                do i=1,nx
                    xz_bar(i,k)=sum(ESGS3D(i,:,k))/nyb
                enddo
            enddo
            call output_average(xz_bar/hprocs,321)
            do k=2,nzb+1
                do i=1,nx
                    xz_bar(i,k)=sum((u(i,:,k)-sum(u(i,:,k))/nyb)**2)/nyb
                enddo
            enddo
            call output_average(xz_bar/hprocs,322)
            do k=2,nzb+1
                do i=1,nx
                    xz_bar(i,k)=sum((v(i,:,k)-sum(v(i,:,k))/nyb)**2)/nyb
                enddo
            enddo
            call output_average(xz_bar/hprocs,323)
            do k=2,nzb+1
                do i=1,nx
                    xz_bar(i,k)=sum((w(i,:,k)-sum(w(i,:,k))/nyb)**2)/nyb
                enddo
            enddo
            call output_average(xz_bar/hprocs,324)

            call output_average(au,71)
            call output_average(av,72)
            call output_average(aw,73)
            call output_average(ap,85)
            call output_average(u2,74)
            call output_average(v2,75)
            call output_average(w2,76)
            call output_average(w3,40)
            call output_average(atxx,57)
            call output_average(atxz,78)
            call output_average(atyy,79)
            call output_average(atyz,80)
            call output_average(atzz,87)
            call output_average(atxy,90)
            call output_average(p2,81)
            call output_average(auv,77)
            call output_average(auw,82)
            call output_average(avw,83)
            call output_average(adudz,88)
            call output_average(adudx,91)
            call output_average(advdz,95)
            call output_average(adwdz,93)
            call output_average(adwdx,92)
            call output_average(e,39)
            call output_average(aCs2,398)
            call output_average(aCs,399)
            call output_average(abeta1,978)
            call output_average(aESGS,888)
            call output_average(aDSGS,321)
            call output_average(aTL,320)

            if(scalarCount >= 1)then
                do l=1,scalarCount
                    call output_average(scalarMean(:,:,l), &
                    scalarMeanUnits(l))
                    call output_average(t2(:,:,l), &
                    scalar2Units(l))
                    call output_average(t3(:,:,l), &
                    scalar3Units(l))
                    call output_average(asgs_t1(:,:,l), &
                    xScalarFluxUnits(l))
                    call output_average(asgs_t2(:,:,l), &
                    yScalarFluxUnits(l))
                    call output_average(asgs_t3(:,:,l), &
                    zScalarFluxUnits(l))
                    call output_average(aut(:,:,l), &
                    ausUnits(l))
                    call output_average(avt(:,:,l), &
                    avsUnits(l))
                    call output_average(awt(:,:,l), &
                    awsUnits(l))
                    call output_average(adsdx(:,:,l), &
                    dsdxUnits(l))
                    call output_average(adsdy(:,:,l), &
                    dsdyUnits(l))
                    call output_average(adsdz(:,:,l), &
                    dsdzUnits(l))
                    call output_average(aPr(:,:,l), &
                    scalarPrandtlUnits(l))
                    call output_average(aCs2Pr(:,:,l), &
                    cs2prUnits(l))
                    call output_average(abeta2(:,:,l), &
                    beta2Units(l))
                    call output_average(aqz_s(:,:,l), &
                    surfaceFluxUnits(l))
                    call output_average(aET(:,:,l), &
                    ETUnits(l))

                enddo
            endif

        endif

        call zeroslice(au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy, &
        atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz, &
        adwdx,e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL)

        if(scalarCount >= 1)then
            call s_zeroslice(scalarMean,t2,t3,asgs_t1,asgs_t2, &
            asgs_t3,aut,avt,awt,adsdx,adsdy,adsdz,aPr,aCs2Pr, &
            abeta2,aqz_s,aET)
        endif

    ENDIF


!     =========== C
!     =========== C

!   Write frames

        if(nframe>0)then

            if(t>1 .and. output_frames(8).and. &
                ttt.ge.sframe.and.ttt.le.sframe+(nframe-1)*framestep)then
                if(ttt.eq.sframe.or.mod(ttt-sframe,framestep).eq.0)then
                    call calc_lambda2(lambda2,dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                        dwdx,dwdy,dwdz)
                end if
            end if
        
            call write_frames(u,v,w,p,ESGS3D,DSGS3D,TL,lambda2,particle, &
                scalar(:,:,:,1))            
                    
        end if !! nframe
! cc  Write Break Files
    IF (mod(ttt,ruler) == 0) then

        if(me == nprocs-1)then
            open(unit=334,file='checkpoints/break',form='formatted')
            read(334,*)
            read(334,*)
            read(334,*)i
            close(334)
            open(unit=334,file='checkpoints/break',form='formatted', &
            status='replace')
            write(334,*)ttt
            write(334,*)frame_cnt
            write(334,*)i
            close(334)
        endif
        call MPI_BARRIER(nall,ierr)

    ! undo Gallilean transformation
        u=u+Ugal
        v=v+Vgal
    ! dimensionalize
        u=u*u_star
        v=v*u_star
        w=w*u_star

        view_disp=0
        call MPI_FILE_OPEN(nall,'checkpoints/u.break', &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
        fh,ierr)
        call write3D(fh,u,1,1,2,fileview_3D,view_disp)
        call MPI_FILE_CLOSE(fh,ierr)
        call MPI_FILE_OPEN(nall,'checkpoints/v.break', &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
        fh,ierr)
        call write3D(fh,v,1,1,2,fileview_3D,view_disp)
        call MPI_FILE_CLOSE(fh,ierr)
        call MPI_FILE_OPEN(nall,'checkpoints/w.break', &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
        fh,ierr)
        call write3D(fh,w,1,1,2,fileview_3D,view_disp)
        call MPI_FILE_CLOSE(fh,ierr)


    ! non-dimensionalize
        u=u/u_star
        v=v/u_star
        w=w/u_star

    ! gallilean transformation
        u=u-Ugal
        v=v-Vgal

        IF(npart > 0 .AND. me == 0)THEN
        ! cc  Write Particle Positions

            open(unit=390,file= &
            'checkpoints/ParticlePositions_nomodel.break')
            write(390,*)ipart
            do i=1,ipart
                write(390,*)particle(i,1,1),particle(i,2,1), &
                particle(i,3,1),particle(i,4,1), &
                particle(i,5,1),particle(i,6,1)
            enddo
            close(unit=390)
            if(part_model == 2)then
                open(unit=390,file= &
                'checkpoints/ParticlePositions_weiliso.break', &
                form='formatted')
                write(390,*)ipart
                do i=1,ipart
                    write(390,*)particle(i,1,2),particle(i,2,2), &
                    particle(i,3,2),particle(i,4,2), &
                    particle(i,5,2),particle(i,6,1)
                enddo
                close(unit=390)
            endif

            if(part_model >= 2)then
            ! cc  Write SGS Particle Velocity
                call MPI_FILE_OPEN(MPI_COMM_SELF, &
                'checkpoints/SGS_velocity.break',MPI_MODE_WRONLY &
                +MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_WRITE(fh,us_f(1,1,1),npart*3*2, &
                MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                call MPI_FILE_CLOSE(fh,ierr)

            ! cc  Write Particle ESGS
                call MPI_FILE_OPEN(MPI_COMM_SELF, &
                'checkpoints/ESGS_f.break',MPI_MODE_WRONLY &
                +MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_WRITE(fh,ESGS_f(1,1,2),nx*nyb*nzb, &
                MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                call MPI_FILE_CLOSE(fh,ierr)
            endif

        ENDIF

    ! cc Write RHS Momentum

        call MPI_FILE_OPEN(nall,'checkpoints/RHS_momX.break', &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
        fh,ierr)
        call write3D(fh,RHSx,1,1,2,fileview_3D,view_disp)
        call MPI_FILE_CLOSE(fh,ierr)
        call MPI_FILE_OPEN(nall,'checkpoints/RHS_momY.break', &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
        fh,ierr)
        call write3D(fh,RHSy,1,1,2,fileview_3D,view_disp)
        call MPI_FILE_CLOSE(fh,ierr)
        call MPI_FILE_OPEN(nall,'checkpoints/RHS_momZ.break', &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
        fh,ierr)
        call write3D(fh,RHSz,1,1,2,fileview_3D,view_disp)
        call MPI_FILE_CLOSE(fh,ierr)


        if(model == 2)then
        ! cc     Write a1-b1

            tmp_old=a1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/a1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=b1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/b1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)


        elseif(model == 3 .AND. averaging == 1)then
        ! cc     Write a1_old-e1_old and a2_old-e2_old

            tmp_old=a1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/a1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=b1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/b1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=c1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/c1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=d1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/d1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=e1_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/e1.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=a2_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/a2.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=b2_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/b2.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=c2_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/c2.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=d2_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/d2.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=e2_old(:,2:nyb+1,2:nzb+1)
            call MPI_FILE_OPEN(nall,'checkpoints/e2.break', &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)



        endif

    ! cc  Write Scalar(s)

        IF (scalarCount >= 1) then
            scalarIter = 0
            do l=1,scalarCount
                if (scalarFlags(l) == 1)then ! scalar is temperature
                    write(fileStr1, '(A)')'checkpoints/temperature.break'
                    write(fileStr21, '(A)') &
                    'checkpoints/a4_temperature.break'
                    write(fileStr22, '(A)') &
                    'checkpoints/b4_temperature.break'
                    write(fileStr23, '(A)') &
                    'checkpoints/c4_temperature.break'
                    write(fileStr24, '(A)') &
                    'checkpoints/d4_temperature.break'
                    write(fileStr25, '(A)') &
                    'checkpoints/e4_temperature.break'
                    write(fileStr31, '(A)') &
                    'checkpoints/a8_temperature.break'
                    write(fileStr32, '(A)') &
                    'checkpoints/b8_temperature.break'
                    write(fileStr33, '(A)') &
                    'checkpoints/c8_temperature.break'
                    write(fileStr34, '(A)') &
                    'checkpoints/d8_temperature.break'
                    write(fileStr35, '(A)') &
                    'checkpoints/e8_temperature.break'
                    write(fileStr4, '(A)') &
                    'checkpoints/RHS_temperature.break'
                    write(fileStr5, '(A)')'checkpoints/coolrate.break'
                    write(fileStr6, '(A)') &
                    'checkpoints/surfaceTemperature.break'
                elseif ( scalarFlags(l) == 2)then ! scalar is moisture
                    write(fileStr1, '(A)')'checkpoints/moisture.break'
                    write(fileStr21, '(A)') &
                    'checkpoints/a4_moisture.break'
                    write(fileStr22, '(A)') &
                    'checkpoints/b4_moisture.break'
                    write(fileStr23, '(A)') &
                    'checkpoints/c4_moisture.break'
                    write(fileStr24, '(A)') &
                    'checkpoints/d4_moisture.break'
                    write(fileStr25, '(A)') &
                    'checkpoints/e4_moisture.break'
                    write(fileStr3, '(A)') &
                    'checkpoints/a8_e8_moisture.break'
                    write(fileStr31, '(A)') &
                    'checkpoints/a8_moisture.break'
                    write(fileStr32, '(A)') &
                    'checkpoints/b8_moisture.break'
                    write(fileStr33, '(A)') &
                    'checkpoints/c8_moisture.break'
                    write(fileStr34, '(A)') &
                    'checkpoints/d8_moisture.break'
                    write(fileStr35, '(A)') &
                    'checkpoints/e8_moisture.break'
                    write(fileStr4, '(A)')'checkpoints/RHS_moisture.break'
                    write(fileStr5, '(A)') &
                    'checkpoints/dsdt_Moisture.break'
                    write(fileStr6, '(A)') &
                    'checkpoints/surfaceMoisture.break'
                else             ! = 0, scalar is passive
                    scalarIter = scalarIter + 1
                    write(fileStr1, '(A,I1,A)') &
                    'checkpoints/scalar.break',scalarIter,'.break'
                    write(fileStr21, '(A,I1,A)') &
                    'checkpoints/a4_scalar.break',scalarIter,'.break'
                    write(fileStr22, '(A,I1,A)') &
                    'checkpoints/b4_scalar.break',scalarIter,'.break'
                    write(fileStr23, '(A,I1,A)') &
                    'checkpoints/c4_scalar.break',scalarIter,'.break'
                    write(fileStr24, '(A,I1,A)') &
                    'checkpoints/d4_scalar.break',scalarIter,'.break'
                    write(fileStr25, '(A,I1,A)') &
                    'checkpoints/e4_scalar.break',scalarIter,'.break'
                    write(fileStr31, '(A,I1,A)') &
                    'checkpoints/a8_scalar.break',scalarIter,'.break'
                    write(fileStr32, '(A,I1,A)') &
                    'checkpoints/b8_scalar.break',scalarIter,'.break'
                    write(fileStr33, '(A,I1,A)') &
                    'checkpoints/c8_scalar.break',scalarIter,'.break'
                    write(fileStr34, '(A,I1,A)') &
                    'checkpoints/d8_scalar.break',scalarIter,'.break'
                    write(fileStr35, '(A,I1,A)') &
                    'checkpoints/e8_scalar.break',scalarIter,'.break'
                    write(fileStr4, '(A,I1,A)') &
                    'checkpoints/RHS_scalar.break',scalarIter,'.break'
                    write(fileStr5, '(A,I1,A)') &
                    'checkpoints/dsdt_scalar.break',scalarIter,'.break'
                    write(fileStr6, '(A,I1,A)') &
                    'checkpoints/surfaceScalar',scalarIter,'.break'
                endif

                call MPI_FILE_OPEN(nall,fileStr1,MPI_MODE_WRONLY+ &
                MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                fileview_3D,"native",MPI_INFO_NULL,ierr)
                call MPI_FILE_WRITE_ALL(fh,scalar(1,1,2,l),Nx*Nyb*Nzb, &
                MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                call MPI_FILE_CLOSE(fh,ierr)

                if(model == 2)then
                ! cc            Write a4-b4

                    call MPI_FILE_OPEN(nall,filestr21,MPI_MODE_WRONLY+ &
                    MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
                    tmp_old=a4_old(:,2:nyb+1,2:nzb+1,l)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)

                    call MPI_FILE_OPEN(nall,filestr22,MPI_MODE_WRONLY+ &
                    MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
                    tmp_old=b4_old(:,2:nyb+1,2:nzb+1,l)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)


                elseif(model == 3 .AND. averaging == 1)then
                ! cc           Write a4-e4 and a8-e8


                    tmp_old=a4_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,filestr21, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=b4_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,filestr22, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=c4_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,filestr23, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=d4_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,filestr24, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=e4_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,filestr25, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    view_disp=0

                    tmp_old=a8_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,fileStr31, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=b8_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,fileStr32, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=c8_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,fileStr33, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=d8_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,fileStr34, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    tmp_old=e8_old(:,2:nyb+1,2:nzb+1,l)
                    call MPI_FILE_OPEN(nall,fileStr35, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
                    call MPI_FILE_CLOSE(fh,ierr)
                    view_disp=0

                endif

            ! cc  Write Scalar RHS Momentum

                call MPI_FILE_OPEN(nall,filestr4,MPI_MODE_WRONLY+ &
                MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
                call MPI_FILE_SET_VIEW(fh,view_disp, &
                MPI_DOUBLE_PRECISION,fileview_3D,"native", &
                MPI_INFO_NULL,ierr)
                call MPI_FILE_WRITE_ALL(fh,scalarRHS(1,1,2,l),Nx*Nyb*Nzb, &
                MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                call MPI_FILE_CLOSE(fh,ierr)

            ! cc  Write Scalar Surface Files

                if(vfact == 0 .AND. surfaceFlags(l) == 1)then
                    call MPI_FILE_OPEN(MPI_COMM_LEVEL,fileStr5, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh, &
                    ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_2D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,coolrate(1,1,l),Nx*Nyb, &
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                    call MPI_FILE_CLOSE(fh,ierr)

                    call MPI_FILE_OPEN(MPI_COMM_LEVEL,fileStr6, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh, &
                    ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_2D,"native",MPI_INFO_NULL,ierr)
                    tempFrameA2(:,:)=surfaceScalar(:,:,l)*scalarScales(l)
                    call MPI_FILE_WRITE_ALL(fh,tempFrameA2(1,1),Nx*Nyb, &
                    MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                    call MPI_FILE_CLOSE(fh,ierr)

                endif               !end if(vfact == 0)

            enddo  !! l=1, scalarCount
        ENDIF      !!end IF(scalarCount >= 1)
    ENDIF          !!end IF(mod(t,ruler) == 0)

END DO timeStep

!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ... End Time Loop  ========================E
!     ... End Time Loop  ========================E
!     ... End Time Loop  ========================E
!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ===========================================E
!     ===========================================E

if(me == 0)then
    t02=MPI_WTIME()
    call closefiles(dsdxUnits,dsdyUnits,dsdzUnits, &
    scalarMeanUnits,scalar2Units,scalar3Units,xScalarFluxUnits, &
    yScalarFluxUnits,zScalarFluxUnits,awsUnits,ausUnits,avsUnits, &
    scalarSpectraUnits,scalarPrandtlUnits,cs2prUnits,beta2Units, &
    ETUnits,obukovUnits,scalarStarUnits,surfaceFluxUnits)
endif

! cc  Write Velocity Fields

! undo Gallilean transformation
u=u+Ugal
v=v+Vgal

! dimensionalize
u=u*u_star
v=v*u_star
w=w*u_star

view_disp=0
call MPI_FILE_OPEN(nall,'checkpoints/u.break', &
MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
fh,ierr)
call write3D(fh,u,1,1,2,fileview_3D,view_disp)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_FILE_OPEN(nall,'checkpoints/v.break', &
MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
fh,ierr)
call write3D(fh,v,1,1,2,fileview_3D,view_disp)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_FILE_OPEN(nall,'checkpoints/w.break', &
MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
fh,ierr)
call write3D(fh,w,1,1,2,fileview_3D,view_disp)
call MPI_FILE_CLOSE(fh,ierr)


IF(npart > 0 .AND. me == 0)THEN
! cc  Write Particle Positions

    open(unit=390, &
    file='checkpoints/ParticlePositions_nomodel.break', &
    form='formatted')
    write(390,*)ipart
    do i=1,ipart
        write(390,*)particle(i,:,1)
    enddo
    close(unit=390)
    if(part_model == 2)then
        open(unit=390, &
        file='checkpoints/ParticlePositions_weiliso.break', &
        form='formatted')
        write(390,*)ipart
        do i=1,ipart
            write(390,*)particle(i,:,2)
        enddo
        close(unit=390)
    endif

! cc  Write SGS Particle Velocity

    if(part_model >= 2)then
    ! cc  Write SGS Particle Velocity
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
        'checkpoints/SGS_velocity.break',MPI_MODE_WRONLY &
        +MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_WRITE(fh,us_f(1,1,1),npart*3*2, &
        MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)

    ! cc  Write Particle ESGS
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
        'checkpoints/ESGS_f.break',MPI_MODE_WRONLY &
        +MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_WRITE(fh,ESGS_f(1,1,2),nx*nyb*nzb, &
        MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)
    endif

ENDIF

!! Write RHS Momentum

view_disp=0

call MPI_FILE_OPEN(nall,'checkpoints/RHS_momX.break', &
MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
fh,ierr)
call write3D(fh,RHSx,1,1,2,fileview_3D,view_disp)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_FILE_OPEN(nall,'checkpoints/RHS_momY.break', &
MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
fh,ierr)
call write3D(fh,RHSy,1,1,2,fileview_3D,view_disp)
call MPI_FILE_CLOSE(fh,ierr)
call MPI_FILE_OPEN(nall,'checkpoints/RHS_momZ.break', &
MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
fh,ierr)
call write3D(fh,RHSz,1,1,2,fileview_3D,view_disp)
call MPI_FILE_CLOSE(fh,ierr)

view_disp=0

if(model == 2)then
! cc     Write a1-b1

    tmp_old=a1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/a1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=b1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/b1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)


elseif(model == 3 .AND. averaging == 1)then
! cc     Write a1_old-e1_old and a2_old-e2_old

    tmp_old=a1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/a1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=b1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/b1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=c1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/c1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=d1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/d1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=e1_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/e1.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)

    tmp_old=a2_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/a2.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=b2_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/b2.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=c2_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/c2.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=d2_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/d2.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)
    tmp_old=e2_old(:,2:nyb+1,2:nzb+1)
    call MPI_FILE_OPEN(nall,'checkpoints/e2.break', &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
    fh,ierr)
    call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
    call MPI_FILE_CLOSE(fh,ierr)


endif

! cc  Write Scalar(s)

IF (scalarCount >= 1) then
    scalarIter = 0
    do l=1,scalarCount
        if (scalarFlags(l) == 1)then ! scalar is temperature
            write(fileStr1, '(A)')'checkpoints/temperature.break'
            write(fileStr2, '(A)') &
            'checkpoints/a4_e4_temperature.break'
            write(fileStr3, '(A)') &
            'checkpoints/a8_e8_temperature.break'
            write(fileStr4, '(A)') &
            'checkpoints/RHS_temperature.break'
        elseif ( scalarFlags(l) == 2)then ! scalar is moisture
            write(fileStr1, '(A)')'checkpoints/moisture.break'
            write(fileStr2, '(A)') &
            'checkpoints/a4_e4_moisture.break'
            write(fileStr3, '(A)') &
            'checkpoints/a8_e8_moisture.break'
            write(fileStr4, '(A)')'checkpoints/RHS_moisture.break'
        else             ! = 0, scalar is passive
            scalarIter = scalarIter + 1
            write(fileStr1, '(A,I1,A)')'checkpoints/scalar', &
            scalarIter,'.break'
            write(fileStr2, '(A,I1,A)')'checkpoints/a4_e4_scalar', &
            scalarIter,'.break'
            write(fileStr3, '(A,I1,A)')'checkpoints/a8_e8_scalar', &
            scalarIter,'.break'
            write(fileStr4, '(A,I1,A)')'checkpoints/RHS_scalar', &
            scalarIter,'.break'
        endif

        call MPI_FILE_OPEN(nall,fileStr1,MPI_MODE_WRONLY+ &
        MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
        fileview_3D,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_WRITE_ALL(fh,scalar(1,1,2,l),Nx*Nyb*Nzb, &
        MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)

    ! cc            Write a4-b4
        if(model == 2)then

            tmp_old=a4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr21, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=b4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr22, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)

        ! cc        Write a4-e4 and a8-e8
        elseif(model == 3 .AND. averaging == 1)then

            tmp_old=a4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr21, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=b4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr22, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=c4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr23, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=d4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr24, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=e4_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr25, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
    
            tmp_old=a8_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr31, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=b8_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr32, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=c8_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr33, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=d8_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr34, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)
            tmp_old=e8_old(:,2:nyb+1,2:nzb+1,l)
            call MPI_FILE_OPEN(nall,fileStr35, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
            fh,ierr)
            call write3D(fh,tmp_old,1,1,1,fileview_3D,view_disp)
            call MPI_FILE_CLOSE(fh,ierr)


        endif

    ! cc  Write Scalar RHS Momentum

        call MPI_FILE_OPEN(nall,filestr4,MPI_MODE_WRONLY+ &
        MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
        call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
        fileview_3D,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_WRITE_ALL(fh,scalarRHS(1,1,2,l),Nx*Nyb*Nzb, &
        MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)

    ! cc  Write Scalar Surface Files
        if (scalarFlags(l) == 1)then !temperature scalar
            write(fileStr1, '(A)') 'checkpoints/coolrate.break'
            write(fileStr2, '(A)') &
            'checkpoints/surfaceTemperature.break'
        elseif (scalarFlags(l) == 2)then !moisture scalar
            write(fileStr1, '(A)') &
            'checkpoints/dsdt_Moisture.break'
            write(fileStr2, '(A)') &
            'checkpoints/surfaceMoisture.break'
        else
            scalarIter = scalarIter + 1
            write(fileStr1, '(A,I1,A)') &
            'checkpoints/dsdt_scalar',scalarIter,'.break'
            write(fileStr2, '(A,I1,A)') &
            'checkpoints/surfaceScalar',scalarIter,'.break'
        endif

        if(vfact == 0 .AND. surfaceFlags(l) == 1)then
            call MPI_FILE_OPEN(MPI_COMM_LEVEL,fileStr1, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh, &
            ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp, &
            MPI_DOUBLE_PRECISION,fileview_2D,"native", &
            MPI_INFO_NULL,ierr)
            call MPI_FILE_WRITE_ALL(fh,coolrate(1,1,l),Nx*Nyb, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

            call MPI_FILE_OPEN(MPI_COMM_LEVEL,fileStr2, &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh, &
            ierr)
            call MPI_FILE_SET_VIEW(fh,view_disp, &
            MPI_DOUBLE_PRECISION,fileview_2D,"native", &
            MPI_INFO_NULL,ierr)
            tempFrameA2(:,:)=surfaceScalar(:,:,l)*scalarScales(l)
            call MPI_FILE_WRITE_ALL(fh,tempFrameA2(1,1),Nx*Nyb, &
            MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

        endif               !end if(vfact == 0)
    enddo                  !end do l=1,scalarCount
endif                     !end if(scalarCount >= 1)

if (me==0) then
    write(*,*) 'Total time = ', t02-t01
end if

call MPI_BARRIER(nall,ierr) !synchronize before exit

999 call MPI_FINALIZE(ierr)

END PROGRAM

!=================================================================
