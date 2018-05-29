    subroutine dispersion(particle,u,v,w,ESGS,zo,ustar,LAD_w, &
    us_f,ESGS_f,ESGS_part_f,tau_f,idum,fx,fz,txx,txy,txz, &
    tyy,tyz,tzz,divtx,divty,divtz,TL,rogue_count, &
    trajectory_updates,rogue_flag)
    use globals
    use wallBoundaryConditions
    use mainModule
    use particleModule
    implicit none

    interface
    include 'interfaces/reshape_global_array.f90'
    include 'interfaces/ddx.f90'
    include 'interfaces/ddy.f90'
    include 'interfaces/update1.f90'
    include 'interfaces/spanwise_avg.f90'
    include 'interfaces/interp3D_particle.f90'
    include 'interfaces/ddx_upwind.f90'
    include 'interfaces/ddy_upwind.f90'
    end interface

    integer :: idum
    integer*4 :: i,j,k,ii,jj,kk,p,M,hp,vp,positive_T1s,rogue_count, &
    trajectory_updates,second_sign
    real :: randn,ran1
    real*8 :: dtp,iw,jw,kw,ur_,vr_,wr_, &
    us,vs,ws,us_1,vs_1,ws_1,us_2,vs_2,ws_2,us_3,vs_3,ws_3, &
    desdx_,desdy_,desdz_,es_, &
    LAD_,Ez,Gv,Gg,ustar_,zo_,mag,fx_,fz_, &
    txx_,txy_,txz_,tyy_,tyz_,tzz_,divtx_,divty_,divtz_, &
    lxx_,lxy_,lxz_,lyy_,lyz_,lzz_,t_1,t_2, &
    ddt_txx_,ddt_txy_,ddt_txz_,ddt_tyy_,ddt_tyz_,ddt_tzz_, &
    ddx_txx_,ddx_txy_,ddx_txz_,ddy_txy_,ddy_tyy_,ddy_tyz_, &
    ddz_txz_,ddz_tyz_,ddz_tzz_,TL_,ddt_es_, &
    DDt_max,var_gamma,T1,total_deriv,dz2
    integer*4,dimension(:) :: rogue_flag
    real*8,dimension(:) :: fx,fz,ESGS_part_f
    real*8,dimension(:,:) :: zo,ustar
    real*8,dimension(:,:,:) :: particle,u,v,w,ESGS,LAD_w, &
    txx,txy,txz,tyy,tyz,tzz,divtx,divty,divtz,us_f,ESGS_f, &
    TL
    real*8,dimension(:,:,:,:) :: tau_f
     
    real*8,dimension(size(u,1),size(u,2),size(u,3)) :: desdx,desdy, &
    desdz,u2,v2,w2,es,ddx_txx,ddx_txy,ddx_tyy,ddx_txz,ddy_txy, &
    ddy_tyy,ddy_tyz,ddt_esgs,ddt_txx,ddt_txy,ddt_txz,ddt_tyy, &
    ddt_tyz,ddt_tzz

    real*8,dimension(size(u,1)*size(u,2)*hprocs*(size(u,3)-2)*vprocs) &
    :: rbuf
    real*8,allocatable,dimension(:,:) :: SGS_timescale_vars
    real*8,allocatable,dimension(:,:,:) :: u_t,v_t,w_t,LAD_t,esgs_t, &
    desdx_t,desdy_t,desdz_t,zo_t,ustar_t, &
    txx_t,txy_t,txz_t,tyy_t,tyz_t,tzz_t,divtx_t,divty_t,divtz_t, &
    ddx_txz_t,ddy_tyz_t,ddx_txx_t,ddx_txy_t, &
    ddy_txy_t,ddy_tyy_t,TL_t,ddt_esgs_t, &
    ddt_txx_t,ddt_txy_t,ddt_txz_t,ddt_tyy_t,ddt_tyz_t,ddt_tzz_t

    character(LEN=60) :: filestr

    if(me == 0)then
        allocate(u_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
        v_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
        w_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
        LAD_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
        zo_t(size(u,1),size(u,2)*hprocs,2), &
        ustar_t(size(u,1),size(u,2)*hprocs,2), &
        SGS_timescale_vars(npart,6))
        if(part_model >= 2)then
            allocate( &
            esgs_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            ddt_esgs_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            desdx_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            desdy_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            desdz_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            TL_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs))
        endif
        if(part_model == 3)then
            allocate( &
            txx_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            txy_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            txz_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            tyy_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            tyz_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            tzz_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            divtx_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            divty_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            divtz_t(size(u,1),size(u,2)*hprocs,(size(u,3)-2)*vprocs), &
            ddx_txx_t(size(u,1),size(u,2)*hprocs,2), &
            ddx_txy_t(size(u,1),size(u,2)*hprocs,2), &
            ddx_txz_t(size(u,1),size(u,2)*hprocs,2), &
            ddy_txy_t(size(u,1),size(u,2)*hprocs,2), &
            ddy_tyy_t(size(u,1),size(u,2)*hprocs,2), &
            ddy_tyz_t(size(u,1),size(u,2)*hprocs,2), &
            ddt_txx_t(size(u,1),size(u,2)*hprocs,2), &
            ddt_txy_t(size(u,1),size(u,2)*hprocs,2), &
            ddt_txz_t(size(u,1),size(u,2)*hprocs,2), &
            ddt_tyy_t(size(u,1),size(u,2)*hprocs,2), &
            ddt_tyz_t(size(u,1),size(u,2)*hprocs,2), &
            ddt_tzz_t(size(u,1),size(u,2)*hprocs,2))
        endif
    endif

    dz2=0.5d0*dz

    dtp=dt*skip_step

    positive_T1s=0
    second_sign=0

! cc  me=0 Gather u,v,w,LAD
    call MPI_BARRIER(nall,ierr)
    call MPI_GATHER(u(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
    rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
    call reshape_global_array(rbuf,u_t,3)
    call MPI_GATHER(v(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
    rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
    call reshape_global_array(rbuf,v_t,3)
    call MPI_GATHER(w(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
    rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
    call reshape_global_array(rbuf,w_t,3)
    if(me == 0)w_t(:,:,1)=0.d0
    if(deposition == 1)then
        call MPI_GATHER(LAD_w(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,LAD_t,3)
    endif

! cc  me=0 Gather zo,ustar from vfact=0
    if(vfact == 0)then
        call MPI_GATHER(zo(1,1),nx*nyb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL,ierr)
        call reshape_global_array(rbuf(1:nx*ny),zo_t,2)
        if(me == 0)then
            zo_t(:,:,2)=zo_t(:,:,1)
        endif
        call MPI_GATHER(ustar(1,1),nx*nyb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL,ierr)
        call reshape_global_array(rbuf(1:nx*ny),ustar_t,2)
        if(me == 0)then
            ustar_t(:,:,2)=ustar_t(:,:,1)
        endif
    endif

    IF(part_model >= 2)THEN

        ESGS=ESGS*2.d0/3.d0    !change esgs to sigma^2_s

        if(vfact == 0)then
            ESGS(:,:,2)=0.d0
        endif

    ! cc  Calculate Spatial Derivs of SGS TKE on Eulerian Grid

        call ddx_upwind(desdx,ESGS,u)
        call ddy_upwind(desdy,ESGS,v)

    ! all ddx(desdx,ESGS)
    ! all ddy(desdy,ESGS)
                 
        call update1(ESGS)
        do k=2,nzb+1
            desdz(:,:,k)=(ESGS(:,:,k+1)-ESGS(:,:,k))/dz
        enddo
        if(vfact == vprocs-1)then
            desdz(:,:,Nzb+1)=0.d0
        elseif(vfact == 0)then
            desdz(:,:,2)=0.d0
        endif

        call MPI_GATHER(esgs(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,esgs_t,3)
        call MPI_GATHER(desdx(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,desdx_t,3)
        call MPI_GATHER(desdy(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,desdy_t,3)
        call MPI_GATHER(desdz(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,desdz_t,3)

        call MPI_GATHER(TL(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,TL_t,3)


    ! cc Calculate time derivative of SGS TKE on Eulerian Grid

        if(sum(abs(ESGS_f)) == 0)then
            print*,me,'assuming ddt_esgs=0'
            ddt_esgs=0.d0
        else
            ddt_esgs=(ESGS-ESGS_f)/dtp
        endif
        ESGS_f=ESGS

        call MPI_GATHER(ddt_esgs(1,1,2),nx*nyb*nzb, &
        MPI_DOUBLE_PRECISION,rbuf(1),nx*nyb*nzb, &
        MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,ddt_esgs_t,3)

    ENDIF
    IF(part_model == 3)THEN

    ! cc  U-node horizontal gradients at first node
    ! TODO this is extremely wasteful, work on making sure
    ! his is only calculated on first level
        call ddx(ddx_txx,txx)
        call ddx(ddx_txy,txy)
        call ddx(ddx_txz,txz)
        call ddy(ddy_txy,txy)
        call ddy(ddy_tyy,tyy)
        call ddy(ddy_tyz,tyz)

    ! cc  me=0 Gather tij, divtx, divty, divtz

        call MPI_GATHER(txx(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,txx_t,3)
        call MPI_GATHER(txy(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,txy_t,3)
        call MPI_GATHER(txz(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,txz_t,3)
        call MPI_GATHER(tyy(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,tyy_t,3)
        call MPI_GATHER(tyz(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,tyz_t,3)
        call MPI_GATHER(tzz(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,tzz_t,3)
        call MPI_GATHER(divtx(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,divtx_t,3)
        call MPI_GATHER(divty(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,divty_t,3)
        call MPI_GATHER(divtz(1,1,2),nx*nyb*nzb,MPI_DOUBLE_PRECISION, &
        rbuf(1),nx*nyb*nzb,MPI_DOUBLE_PRECISION,0,nall,ierr)
        call reshape_global_array(rbuf,divtz_t,3)

    ! cc  me=0 Gather ddx_t

        if(vfact == 0)then
            call MPI_GATHER(ddx_txx(1,1,2),nx*nyb,MPI_DOUBLE_PRECISION, &
            rbuf(1),nx*nyb,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL, &
            ierr)
            call reshape_global_array(rbuf(1:nx*ny),ddx_txx_t,2)
            if(me == 0)then
                ddx_txx_t(:,:,2)=0.d0
            endif
            call MPI_GATHER(ddx_txy(1,1,2),nx*nyb,MPI_DOUBLE_PRECISION, &
            rbuf(1),nx*nyb,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL, &
            ierr)
            call reshape_global_array(rbuf(1:nx*ny),ddx_txy_t,2)
            if(me == 0)then
                ddx_txy_t(:,:,2)=0.d0
            endif
            call MPI_GATHER(ddx_txz(1,1,2),nx*nyb*2,MPI_DOUBLE_PRECISION, &
            rbuf(1),nx*nyb*2,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL, &
            ierr)
            call reshape_global_array(rbuf(1:nx*ny*2),ddx_txz_t,3)
            call MPI_GATHER(ddy_txy(1,1,2),nx*nyb,MPI_DOUBLE_PRECISION, &
            rbuf(1),nx*nyb,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL, &
            ierr)
            call reshape_global_array(rbuf(1:nx*ny),ddy_txy_t,2)
            if(me == 0)then
                ddx_txy_t(:,:,2)=0.d0
            endif
            call MPI_GATHER(ddy_tyy(1,1,2),nx*nyb,MPI_DOUBLE_PRECISION, &
            rbuf(1),nx*nyb,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL, &
            ierr)
            call reshape_global_array(rbuf(1:nx*ny),ddy_tyy_t,2)
            if(me == 0)then
                ddy_tyy_t(:,:,2)=0.d0
            endif
            call MPI_GATHER(ddy_tyz(1,1,2),nx*nyb*2,MPI_DOUBLE_PRECISION, &
            rbuf(1),nx*nyb*2,MPI_DOUBLE_PRECISION,0,MPI_COMM_LEVEL, &
            ierr)
            call reshape_global_array(rbuf(1:nx*ny*2),ddy_tyz_t,3)
        endif

    ENDIF

    call MPI_BARRIER(nall,ierr)

! cc  me=0 loops through particles and does calculations

    if(me == 0)then

        do M=1,size(particle,3)
            do p=1,ipart
                if(particle(p,4,M) == 0)then

                ! lobal indices
                    ii = floor(particle(p,1,M)/dx)+1
                    jj = floor(particle(p,2,M)/dy)+1
                    kk = floor(particle(p,3,M)/dz)+1
                                      
                    iw=(particle(p,1,M)/dx-floor(particle(p,1,M)/dx))
                    jw=(particle(p,2,M)/dy-floor(particle(p,2,M)/dy))
                    kw=(particle(p,3,M)/dz-floor(particle(p,3,M)/dz))

                ! cc  Boundary Conditions
                    if(ii >= Nx)then !particle is beyond x-boundaries
                        ii=ii-(Nx-1)*(floor(real(ii-1)/(Nx-1)))
                    elseif(ii < 1)then
                        ii=ii+(Nx-1)*abs(floor(real(ii-1)/(Nx-1)))
                    endif
                    if(ii < 1 .OR. ii >= Nx)then
                        ii=1
                        particle(p,4,M)=1
                    endif
                                      
                    if(jj >= Ny)then !particle is beyond y-boundaries
                        jj=jj-(Ny-1)*(floor(real(jj-1)/(Ny-1)))
                    elseif(jj < 1)then
                        jj=jj+(Ny-1)*abs(floor(real(jj-1)/(Ny-1)))
                    endif
                    if(jj < 1 .OR. jj >= Ny)then
                        jj=1
                        particle(p,4,M)=1
                    endif
                                      
                    if(kk < 1 .OR. kk >= Nz)then
                        write(*,*)'warning: particle out of bounds'
                        kk=1
                    endif

                ! cc  Do Interpolations

                ! ariables on w-nodes
                    call interp3D_particle(zo_,iw,jw,0.d0, &
                    zo_t(ii:ii+1,jj:jj+1,:))
                    zo_=zo_/z_i !nondimensionalize
                    if(deposition == 1)then
                        call interp3D_particle(LAD_,iw,jw,kw, &
                        LAD_t(ii:ii+1,jj:jj+1,kk:kk+1))
                    endif
                    if(particle(p,3,M) >= dz)then !particle above first grid node
                        call interp3D_particle(wr_,iw,jw,kw, &
                        w_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        if(M == 3)then
                            call interp3D_particle(txz_,iw,jw,kw, &
                            txz_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(tyz_,iw,jw,kw, &
                            tyz_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(tzz_,iw,jw,kw, &
                            tzz_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(divtz_,iw,jw,kw, &
                            divtz_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        endif
                        if(M > 1)then
                            call interp3D_particle(desdz_,iw,jw,kw, &
                            desdz_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        endif
                    else !particle below first grid node
                        call interp3D_particle(wr_,iw,jw,0.d0, &
                        w_t(ii:ii+1,jj:jj+1,2:3))
                    ! r_=(particle(p,3,M)-zo_)/(dz-zo_)*wr_
                    ! r_=dlog(particle(p,3,M)/zo_)/dlog(dz/zo_)*wr_
                        if(M == 3)then
                            call interp3D_particle(txz_,iw,jw,0.d0, &
                            txz_t(ii:ii+1,jj:jj+1,2:3))
                        ! xz_=(particle(p,3,M)-zo_)/(dz-zo_)*txz_
                            call interp3D_particle(tyz_,iw,jw,0.d0, &
                            tyz_t(ii:ii+1,jj:jj+1,2:3))
                        ! yz_=(particle(p,3,M)-zo_)/(dz-zo_)*tyz_
                            call interp3D_particle(tzz_,iw,jw,0.d0, &
                            tzz_t(ii:ii+1,jj:jj+1,2:3))
                        ! zz_=(particle(p,3,M)-zo_)/(dz-zo_)*tzz_
                            call interp3D_particle(divtz_,iw,jw,0.d0, &
                            divtz_t(ii:ii+1,jj:jj+1,2:3))
                        ! ivtz_=(particle(p,3,M)-zo_)/(dz-zo_)*divtz_
                        endif
                        if(M > 1)then
                            call interp3D_particle(desdz_,iw,jw,0.d0, &
                            desdz_t(ii:ii+1,jj:jj+1,2:3))
                        ! esdz_=(particle(p,3,M)-zo_)/(dz-zo_)*desdz_
                        endif
                    endif
                                   
                ! hift indices to u,v nodes
                    if(kw < 0.5)then
                        kw=kw+0.5
                    elseif(kw >= 0.5)then
                        kw=kw-0.5
                    endif

                    if(particle(p,3,M) >= dz2)then
                        call interp3D_particle(ur_,iw,jw,kw, &
                        u_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        call interp3D_particle(vr_,iw,jw,kw, &
                        v_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        if(M > 1)then
                            call interp3D_particle(es_,iw,jw,kw, &
                            esgs_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(ddt_es_,iw,jw,kw, &
                            ddt_esgs_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(desdx_,iw,jw,kw, &
                            desdx_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(desdy_,iw,jw,kw, &
                            desdy_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(TL_,iw,jw,kw, &
                            TL_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        endif
                        if(M == 3)then
                            call interp3D_particle(txx_,iw,jw,kw, &
                            txx_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(txy_,iw,jw,kw, &
                            txy_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(tyy_,iw,jw,kw, &
                            tyy_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(divtx_,iw,jw,kw, &
                            divtx_t(ii:ii+1,jj:jj+1,kk:kk+1))
                            call interp3D_particle(divty_,iw,jw,kw, &
                            divty_t(ii:ii+1,jj:jj+1,kk:kk+1))
                        endif
                    else !particles below lowest node

                    !u,v use log-law and partition
                        call interp3D_particle(ustar_,iw,jw,0.d0, &
                        ustar_t(ii:ii+1,jj:jj+1,:))
                        call interp3D_particle(ur_,iw,jw,0.d0, &
                        u_t(ii:ii+1,jj:jj+1,1:2))
                        call interp3D_particle(vr_,iw,jw,0.d0, &
                        v_t(ii:ii+1,jj:jj+1,1:2))

                    ! ag=sqrt(ur_**2+vr_**2)
                        ur_=dlog(particle(p,3,M)/zo_)/dlog(dz2/zo_)*ur_
                        vr_=dlog(particle(p,3,M)/zo_)/dlog(dz2/zo_)*vr_
                    ! r_=dlog(particle(p,3,M)/zo_)*ur_/mag
                    ! r_=dlog(particle(p,3,M)/zo_)*vr_/mag

                        if(particle(p,3,M) <= zo_)then
                            write(*,*)'WARNING: particle #',p,'is below zo,'
                            write(*,*)'resolved velocity will be infinite.'
                        endif

                        if(M > 1)then
                        ! s linear interpolate to 0 at z_o
                            call interp3D_particle(es_,iw,jw,0.d0, &
                            esgs_t(ii:ii+1,jj:jj+1,1:2))
                            es_=(particle(p,3,M)-zo_)/(dz2-zo_)*es_
                        ! dt_es linear interpolate to 0 at z_o
                            call interp3D_particle(ddt_es_,iw,jw,0.d0, &
                            ddt_esgs_t(ii:ii+1,jj:jj+1,1:2))
                            ddt_es_=(particle(p,3,M)-zo_)/(dz2-zo_)*ddt_es_
                        ! s gradients interpolate  to 0 at z_o
                            call interp3D_particle(desdx_,iw,jw,0.d0, &
                            desdx_t(ii:ii+1,jj:jj+1,1:2))
                            desdx_=(particle(p,3,M)-zo_)/(dz2-zo_)*desdx_
                            call interp3D_particle(desdy_,iw,jw,0.d0, &
                            desdy_t(ii:ii+1,jj:jj+1,1:2))
                            desdy_=(particle(p,3,M)-zo_)/(dz2-zo_)*desdy_
                        ! L linear interpolate to 0 at z_o
                            call interp3D_particle(TL_,iw,jw,0.d0, &
                            TL_t(ii:ii+1,jj:jj+1,1:2))
                            TL_=(particle(p,3,M)-zo_)/(dz2-zo_)*TL_
                        endif
                        if(M == 3)then
                        ! xx, txy, tyy linear to 0 at z_o
                            call interp3D_particle(txx_,iw,jw,0.d0, &
                            txx_t(ii:ii+1,jj:jj+1,1:2))
                            txx_=(particle(p,3,M)-zo_)/(dz2-zo_)*txx_
                            call interp3D_particle(txy_,iw,jw,0.d0, &
                            txy_t(ii:ii+1,jj:jj+1,1:2))
                            txy_=(particle(p,3,M)-zo_)/(dz2-zo_)*txy_
                            call interp3D_particle(tyy_,iw,jw,0.d0, &
                            tyy_t(ii:ii+1,jj:jj+1,1:2))
                            tyy_=(particle(p,3,M)-zo_)/(dz2-zo_)*tyy_
                        ! au gradients (divergence terms)
                        ! dx(txz) and ddy(tyz) interpolate btwn k=1,2
                            call interp3D_particle(ddx_txz_,iw,jw,kw-0.5, &
                            ddx_txz_t(ii:ii+1,jj:jj+1,1:2))
                            call interp3D_particle(ddy_tyz_,iw,jw,kw-0.5, &
                            ddy_tyz_t(ii:ii+1,jj:jj+1,1:2))
                        ! dz(txz) and ddz(tyz) assume grad is const
                            call interp3D_particle(t_1,iw,jw,0.d0, &
                            txz_t(ii:ii+1,jj:jj+1,1:2))
                            call interp3D_particle(t_2,iw,jw,1.d0, &
                            txz_t(ii:ii+1,jj:jj+1,1:2))
                            ddz_txz_=(t_2-t_1)/dz
                            call interp3D_particle(t_1,iw,jw,0.d0, &
                            tyz_t(ii:ii+1,jj:jj+1,1:2))
                            call interp3D_particle(t_2,iw,jw,1.d0, &
                            tyz_t(ii:ii+1,jj:jj+1,1:2))
                            ddz_tyz_=(t_2-t_1)/dz
                            ddz_tzz_=0.d0
                        ! dx(txx), ddx(txy), ddx(txz), ddy(txy), ddy(tyy), ddy(tyz) linear to zero at z_o
                            call interp3D_particle(ddx_txx_,iw,jw,0.d0, &
                            ddx_txx_t(ii:ii+1,jj:jj+1,1:2))
                            ddx_txx_=(particle(p,3,M)-zo_)/(dz2-zo_)*ddx_txx_
                            call interp3D_particle(ddx_txy_,iw,jw,0.d0, &
                            ddx_txy_t(ii:ii+1,jj:jj+1,1:2))
                            ddx_txy_=(particle(p,3,M)-zo_)/(dz2-zo_)*ddx_txy_
                            call interp3D_particle(ddy_txy_,iw,jw,0.d0, &
                            ddy_txy_t(ii:ii+1,jj:jj+1,1:2))
                            ddy_txy_=(particle(p,3,M)-zo_)/(dz2-zo_)*ddy_txy_
                            call interp3D_particle(ddy_tyy_,iw,jw,0.d0, &
                            ddy_tyy_t(ii:ii+1,jj:jj+1,1:2))
                            ddy_tyy_=(particle(p,3,M)-zo_)/(dz2-zo_)*ddy_tyy_

                            divtx_=ddx_txx_+ddy_txy_+ddz_txz_
                            divty_=ddx_txy_+ddy_tyy_+ddz_tyz_
                            divtz_=ddx_txz_+ddy_tyz_+ddz_tzz_
                                                    
                        endif

                        if(particle(p,6,M) == 1)then
                            fx_=kw*(fx(kk+1)-fx(kk))+fx(kk)
                            fz_=kw*(fz(kk+1)-fz(kk))+fz(kk)
                        endif

                    endif

                ! cc              Calculate SGS velocities
                    if(M == 1)then
                        us=0.d0
                        vs=0.d0
                        ws=0.d0
                    elseif(M == 2)then

                    !$$$                     !d/dt ESGS
                    !$$$                     if(ESGS_part_f(p).eq.0.d0)then
                    !$$$                        ddt_es_=0.d0
                    !$$$                     else
                    !$$$                        ddt_es_=(es_-ESGS_part_f(p))/dtp
                    !$$$                     endif
                                             
                    ! alculate SGS velocty (Weil isotropic)
                        if(TL_ > 1E-14 .AND. es_ > 1E-14)then

                        !--------------- forward Euler scheme ------------------!

                        !$$$                        us=us_f(p,1,1)+(-us_f(p,1,1)/TL_+
                        !$$$     +                       0.5d0*((ddt_es_+(ur_+us_f(p,1,1))*desdx_+
                        !$$$     +                       (vr_+us_f(p,2,1))*desdy_+
                        !$$$     +                       (wr_+us_f(p,3,1))*desdz_)
                        !$$$     +                       *us_f(p,1,1)/es_+desdx_))*dtp+
                        !$$$     +                       sqrt(2.d0*es_/TL_)*randn(idum)*sqrt(dtp)
                        !$$$                        vs=us_f(p,2,1)+(-us_f(p,2,1)/TL_+
                        !$$$     +                       0.5d0*((ddt_es_+(ur_+us_f(p,1,1))*desdx_+
                        !$$$     +                       (vr_+us_f(p,2,1))*desdy_+
                        !$$$     +                       (wr_+us_f(p,3,1))*desdz_)
                        !$$$     +                       *us_f(p,2,1)/es_+desdy_))*dtp+
                        !$$$     +                       sqrt(2.d0*es_/TL_)*randn(idum)*sqrt(dtp)
                        !$$$                        ws=us_f(p,3,1)+(-us_f(p,3,1)/TL_+
                        !$$$     +                       0.5d0*((ddt_es_+(ur_+us_f(p,1,1))*desdx_+
                        !$$$     +                       (vr_+us_f(p,2,1))*desdy_+
                        !$$$     +                       (wr_+us_f(p,3,1))*desdz_)
                        !$$$     +                       *us_f(p,3,1)/es_+desdz_))*dtp+
                        !$$$     +                       sqrt(2.d0*es_/TL_)*randn(idum)*sqrt(dtp)

                        !--------------- implicit linear/explicit nonlinear (fractional step) ------------------!

                        ! irst sub-step (simple implicit)
                            total_deriv=ddt_es_+ur_*desdx_+vr_*desdy_+ &
                            wr_*desdz_
                            us_1=(us_f(p,1,1)+0.5d0*desdx_*dtp+ &
                            sqrt(2.d0*es_/TL_)*randn(idum)*sqrt(dtp))/ &
                            (1.d0+(1.d0/TL_- &
                            &                        0.5d0*total_deriv/es_)*dtp)
                            vs_1=(us_f(p,2,1)+0.5d0*desdy_*dtp+ &
                            sqrt(2.d0*es_/TL_)*randn(idum)*sqrt(dtp))/ &
                            (1.d0+(1.d0/TL_- &
                            &                        0.5d0*total_deriv/es_)*dtp)
                            ws_1=(us_f(p,3,1)+0.5d0*desdz_*dtp+ &
                            sqrt(2.d0*es_/TL_)*randn(idum)*sqrt(dtp))/ &
                            (1.d0+(1.d0/TL_- &
                            &                        0.5d0*total_deriv/es_)*dtp)
                                                    
                            if(abs(ur_+us_1) > dx/dtp &
                             .OR. abs(vr_+vs_1) > dy/dtp &
                             .OR. abs(wr_+ws_1) > dz/dtp)then
                                rogue_count=rogue_count+1
                                rogue_flag(p)=rogue_flag(p)+1
                            endif
                                                    
                            if(abs(ur_+us_1) > dx/dtp .OR. us /= us)then
                                us=us_f(p,1,1)
                            elseif(abs(vr_+vs_1) > dy/dtp .OR. vs /= vs)then
                                vs=us_f(p,2,1)
                            elseif(abs(wr_+ws_1) > dz/dtp .OR. ws /= ws)then
                                ws=us_f(p,3,1)
                            endif

                        ! econd sub-step (forward Euler)
                            total_deriv=desdx_*us_1+desdy_*vs_1+desdz_*ws_1
                            us=us_1+0.5d0*total_deriv*us_1/es_*dtp
                            vs=vs_1+0.5d0*total_deriv*vs_1/es_*dtp
                            ws=ws_1+0.5d0*total_deriv*ws_1/es_*dtp

                            if(abs(ur_+us) > dx/dtp &
                             .OR. abs(vr_+vs) > dy/dtp &
                             .OR. abs(wr_+ws) > dz/dtp)then
                                rogue_flag(p)=rogue_flag(p)+1
                            endif

                            if(abs(ur_+us) > dx/dtp .OR. us /= us)then
                                us=us_1
                            ! s=us_f(p,1,1)
                            endif
                            if(abs(vr_+vs) > dy/dtp .OR. vs /= vs)then
                                vs=vs_1
                            ! s=us_f(p,2,1)
                            endif
                            if(abs(wr_+ws) > dz/dtp .OR. ws /= ws)then
                                ws=ws_1
                            ! s=us_f(p,3,1)
                            endif
                            if(particle(p,3,2)+(wr_+ws)*dtp > l_z/z_i)then
                                ws=0.d0
                            elseif(particle(p,3,2)+(wr_+ws)*dtp < zo_)then
                                ws=0.d0
                            endif

                        else
                            us=0.d0
                            ws=0.d0
                            ws=0.d0
                        endif

                        trajectory_updates=trajectory_updates+1
                                          
                    ! tore 'old' SGS velocities
                        us_f(p,1,1)=us
                        us_f(p,2,1)=vs
                        us_f(p,3,1)=ws

                    elseif(M == 3)then

                    ! nvert SGS stress tensor
                        call invert_stress(lxx_,lxy_,lxz_,lyy_,lyz_,lzz_, &
                        txx_,txy_,txz_,tyy_,tyz_,tzz_)

                    ! alculate D/Dt(tau)
                    !$$$                     DDt_max=10.d0
                    !$$$                     if(t.eq.start_release)then
                    !$$$                        ddt_txx=0.d0
                    !$$$                        ddt_txy=0.d0
                    !$$$                        ddt_txz=0.d0
                    !$$$                        ddt_tyy=0.d0
                    !$$$                        ddt_tyz=0.d0
                    !$$$                        ddt_tzz=0.d0
                    !$$$                        tau_f(p,:,1)=(/txx_,txy_,txz_,tyy_,tyz_,tzz_/)
                    !$$$                     elseif(t.eq.start_release+skip_step)then
                    !$$$                        ddt_txx=(txx_-tau_f(p,1,1))/dtp
                    !$$$                        ddt_txy=(txy_-tau_f(p,2,1))/dtp
                    !$$$                        ddt_txz=(txz_-tau_f(p,3,1))/dtp
                    !$$$                        ddt_tyy=(tyy_-tau_f(p,4,1))/dtp
                    !$$$                        ddt_tyz=(tyz_-tau_f(p,5,1))/dtp
                    !$$$                        ddt_tzz=(tzz_-tau_f(p,6,1))/dtp
                    !$$$
                    !$$$                        tau_f(p,:,2)=tau_f(p,:,1)
                    !$$$                        tau_f(p,:,1)=(/txx_,txy_,txz_,tyy_,tyz_,tzz_/)
                    !$$$                     else
                    !$$$                        ddt_txx=(txx_-tau_f(p,1,1))/dtp
                    !$$$                        ddt_txy=(txy_-tau_f(p,2,1))/dtp
                    !$$$                        ddt_txz=(txz_-tau_f(p,3,1))/dtp
                    !$$$                        ddt_tyy=(tyy_-tau_f(p,4,1))/dtp
                    !$$$                        ddt_tyz=(tyz_-tau_f(p,5,1))/dtp
                    !$$$                        ddt_tzz=(tzz_-tau_f(p,6,1))/dtp
                    !$$$
                    !$$$                        if(abs(ddt_txx).gt.DDt_max)then
                    !$$$                           ddt_txx=sign(DDt_max,ddt_txx)
                    !$$$                        endif
                    !$$$                        if(abs(ddt_txy).gt.DDt_max)then
                    !$$$                           ddt_txy=sign(DDt_max,ddt_txy)
                    !$$$                        endif
                    !$$$                        if(abs(ddt_txz).gt.DDt_max)then
                    !$$$                           ddt_txz=sign(DDt_max,ddt_txz)
                    !$$$                        endif
                    !$$$                        if(abs(ddt_tyy).gt.DDt_max)then
                    !$$$                           ddt_tyy=sign(DDt_max,ddt_tyy)
                    !$$$                        endif
                    !$$$                        if(abs(ddt_tyz).gt.DDt_max)then
                    !$$$                           ddt_tyz=sign(DDt_max,ddt_tyz)
                    !$$$                        endif
                    !$$$                        if(abs(ddt_tzz).gt.DDt_max)then
                    !$$$                           ddt_tzz=sign(DDt_max,ddt_tzz)
                    !$$$                        endif
                    !$$$
                    !$$$                        tau_f(p,:,2)=tau_f(p,:,1)
                    !$$$                        tau_f(p,:,1)=(/txx_,txy_,txz_,tyy_,tyz_,tzz_/)
                    !$$$                     endif

                    ! alculate SGS velocty (Weil anisotropic)

                    !$$$         us=us_f(p,1,2)-0.5d0*fs_*Cop*eps_*
                    !$$$     +   (lxx_*us_f(p,1,2)+lxy_*us_f(p,2,2)+lxz_*us_f(p,3,2))*dtp+0.5d0*
                    !$$$     +   (lxx_*ddt_txx*us_f(p,1,2)+lxy_*ddt_txx*us_f(p,2,2)+
                    !$$$     +       lxz_*ddt_txx*us_f(p,3,2)+
                    !$$$     +    lxy_*ddt_txy*us_f(p,1,2)+lyy_*ddt_txy*us_f(p,2,2)+
                    !$$$     +       lyz_*ddt_txy*us_f(p,3,2)+
                    !$$$     +    lxz_*ddt_txz*us_f(p,1,2)+lyz_*ddt_txz*us_f(p,2,2)+
                    !$$$     +        lzz_*ddt_txz*us_f(p,3,2)+
                    !$$$     +    divtx_)*dtp+sqrt(fs_*Cop*eps_)*randn(idum)*dtp
                    !$$$
                    !$$$         vs=us_f(p,2,2)-0.5d0*fs_*Cop*eps_*
                    !$$$     +   (lxy_*us_f(p,1,2)+lyy_*us_f(p,2,2)+lyz_*us_f(p,3,2))*dtp+0.5d0*
                    !$$$     +   (lxx_*ddt_txy*us_f(p,1,2)+lxy_*ddt_txy*us_f(p,2,2)+
                    !$$$     +        lxz_*ddt_txy*us_f(p,3,2)+
                    !$$$     +    lxy_*ddt_tyy*us_f(p,1,2)+lyy_*ddt_tyy*us_f(p,2,2)+
                    !$$$     +        lyz_*ddt_tyy*us_f(p,3,2)+
                    !$$$     +    lxz_*ddt_tyz*us_f(p,1,2)+lyz_*ddt_tyz*us_f(p,2,2)+
                    !$$$     +        lzz_*ddt_tyz*us_f(p,3,2)+
                    !$$$     +    divty_)*dtp+sqrt(fs_*Cop*eps_)*randn(idum)*dtp
                    !$$$
                    !$$$         ws=us_f(p,3,2)-0.5d0*fs_*Cop*eps_*
                    !$$$     +   (lxz_*us_f(p,1,2)+lyz_*us_f(p,2,2)+lzz_*us_f(p,3,2))*dtp+0.5d0*
                    !$$$     +   (lxx_*ddt_txz*us_f(p,1,2)+lxy_*ddt_txz*us_f(p,2,2)+
                    !$$$     +       lxz_*ddt_txz*us_f(p,3,2)+
                    !$$$     +    lxy_*ddt_tyz*us_f(p,1,2)+lyy_*ddt_tyz*us_f(p,2,2)+
                    !$$$     +       lyz_*ddt_tyz*us_f(p,3,2)+
                    !$$$     +    lxz_*ddt_tzz*us_f(p,1,2)+lyz_*ddt_tzz*us_f(p,2,2)+
                    !$$$     +        lzz_*ddt_tzz*us_f(p,3,2)+
                    !$$$     +    divtz_)*dtp+sqrt(fs_*Cop*eps_)*randn(idum)*dtp


                    ! tore 'old' SGS velocities
                        us_f(p,1,2)=us
                        us_f(p,2,2)=vs
                        us_f(p,3,2)=ws

                    endif

                ! pdate particle positions

                    particle(p,1,M)=particle(p,1,M)+(ur_+us)*dtp
                    particle(p,2,M)=particle(p,2,M)+(vr_+vs)*dtp
                    particle(p,3,M)=particle(p,3,M)+(wr_+ws-wd)*dtp

                    if(particle(p,1,M) /= particle(p,1,M) .OR. &
                    particle(p,2,M) /= particle(p,2,M) .OR. &
                    particle(p,3,M) /= particle(p,3,M))then
                        write(*,*)'particle is NaN'
                        write(*,*)M,p,ur_,us,vr_,vs,wr_,ws
                        stop
                    endif

                    if(M == 2)then
                        if(abs(ur_) > 100.d0)ur_=0.d0
                        if(abs(us) > 100.d0)us=0.d0
                        if(abs(vr_) > 100.d0)vr_=0.d0
                        if(abs(vs) > 100.d0)vs=0.d0
                        if(abs(wr_) > 100.d0)wr_=0.d0
                        if(abs(ws) > 100.d0)ws=0.d0
                        SGS_timescale_vars(p,:)= &
                        (/ur_,us,vr_,vs,wr_,ws/)
                    end if
                                         
                    if(particle(p,6,M) == 1)then
                    ! oliage deposition model
                        Ez=0.86d0/ &
                        (1.d0+0.442d0*(abs(ur_+us)*wd/(Lv*g_hat)) &
                        **(-1.967d0))
                        Gv = (wd*fx_*LAD_+abs(ur_+us)*fz_*Ez*LAD_)*dtp
                        if(ran1(idum) < Gv)then
                            particle(p,4,M)=1
                        endif
                    endif
                                         
                ! round deposition model
                    if(particle(p,3,M) <= zo_)then
                        if(particle(p,6,M) == 0)then
                            Gg=0.d0
                        elseif(wr_+ws < -wd)then
                            Gg=2.d0*wd/(wd-(wr_+ws))
                        else
                            Gg=1.d0
                        endif
                        if(ran1(idum) < Gg)then
                            particle(p,4,M)=1
                        else
                        ! rite(*,*)'particle',p,'was reflected'
                        !z-position
                            particle(p,3,M)=particle(p,3,M)-(wr_+ws-wd)*dtp
                            particle(p,3,M)=abs(wr_+ws-wd)*dtp &
                            -(particle(p,3,M)-zo_)+zo_
                        endif
                    endif

                    if(particle(p,3,M) > L_z/z_i)then !particle is above z-boundary
                        particle(p,4,M)=1
                    elseif(particle(p,3,M) <= zo_)then !particle below zo
                        write(*,*)'warning: depositing particle',p, &
                        'on ground'
                        particle(p,4,M)=1
                    endif
                                      
                endif
            enddo
        enddo

    ! RITING LAGRANGIAN AUTOCORRELATION STUFF
    !$$$         if(ttt.ge.sframe .and. mod(ttt,framestep).eq.0)then
    !$$$
    !$$$            write(filestr,'(A,I7.7,A)')
    !$$$     +           'output/part_frame/SGS_timescale_vars',
    !$$$     +           ttt,'.out'
    !$$$            open(unit=1010,file=filestr,action='write')
    !$$$            do p=1,npart
    !$$$               write(1010,
    !$$$     +              '(F15.7,1X,F15.7,1X,F15.7,1X,F15.7,1X,
    !$$$     +              F15.7,1X,F15.7)',
    !$$$     +              advance='no')SGS_timescale_vars(p,:)
    !$$$               write(1010,*)
    !$$$            end do
    !$$$            close(1010)
    !$$$
    !$$$            open(unit=1010,file='output/part_frame/rogue_flag.out',
    !$$$     +           status='replace',action='write')
    !$$$            do p=1,npart
    !$$$               write(1010,*)rogue_flag(p)
    !$$$            end do
    !$$$            close(1010)
    !$$$
    !$$$         end if

    endif


    if(me == 0)then
    ! rite(*,*)dble(positive_T1s)/dble(ipart),'percent positive T1s'
    ! rite(*,*)dble(second_sign)/dble(ipart)/3.,'% second skips'
        deallocate(u_t,v_t,w_t,LAD_t,zo_t,ustar_t,SGS_timescale_vars)
        if(part_model >= 2)then
            deallocate(esgs_t,ddt_esgs_t,desdx_t,desdy_t,desdz_t,TL_t)
        endif
        if(part_model == 3)then
            deallocate(txx_t,txy_t,txz_t,tyy_t,tyz_t,tzz_t, &
            divtx_t,divty_t,divtz_t,ddx_txx_t,ddx_txy_t,ddx_txz_t, &
            ddy_txy_t,ddy_tyy_t,ddy_tyz_t,ddt_txx_t,ddt_txy_t, &
            ddt_txz_t,ddt_tyy_t,ddt_tyz_t,ddt_tzz_t)
        endif
    endif

    call MPI_BARRIER(nall,ierr)

    return
    end subroutine dispersion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reshape_global_array(rbuf,u,flag)
    use globals
    implicit none

    integer*4 :: flag,hp,vp,i,k
    real*8,dimension(:) :: rbuf
    real*8,dimension(:,:,:) :: u

    if(me == 0)then

        if(flag == 3)then
            if(size(u,3) == Nz)then
                do i=0,nprocs-1
                    hp=mod(i,hprocs)
                    vp=floor(dble(i)/dble(hprocs))
                    u(:,hp*nyb+1:(hp+1)*nyb,vp*nzb+1:(vp+1)*nzb)= &
                    reshape(rbuf(i*nx*nyb*nzb+1:(i+1)*nx*nyb*nzb), &
                    (/nx,nyb,nzb/))
                enddo
            else
                do k=0,1
                    do i=0,hprocs-1
                        hp=mod(i,hprocs)
                        u(:,hp*nyb+1:(hp+1)*nyb,k+1)= &
                        reshape(rbuf(i*nx*nyb+k*nx*ny+1: &
                        (i+1)*nx*nyb+k*nx*ny),(/nx,nyb/))
                    enddo
                enddo
            endif
        elseif(flag == 2)then
            do i=0,hprocs-1
                hp=mod(i,hprocs)
                u(:,hp*nyb+1:(hp+1)*nyb,1)= &
                reshape(rbuf(i*nx*nyb+1:(i+1)*nx*nyb),(/nx,nyb/))
            enddo
        endif
    endif

    return
    end subroutine reshape_global_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine invert_stress(lxx,lxy,lxz,lyy,lyz,lzz,txx,txy,txz,tyy, &
    tyz,tzz)
    use globals
    use particleModule
    implicit none

    real*8 :: D,lxx,lxy,lxz,lyy,lyz,lzz,txx,txy,txz,tyy, &
    tyz,tzz

! eterminant
    D = txx*(tyy*tzz-tyz**2) - txy*(txy*tzz-tyz*txz) + &
    txz*(txy*tyz-tyy*txz)

! heck for ~0 determinants
    if(D<det_min)then
        D=10.**10.
    endif
          
! alculate inverse
    lxx = (tyy*tzz-tyz**2)/D
    lxy = (tyz*txz-txy*tzz)/D
    lxz = (txy*tyz-tyy*txz)/D
    lyy = (txx*tzz-txz**2)/D
    lyz = (txy*txz-txx*tyz)/D
    lzz = (txx*tyy-txy**2)/D

    return
    end subroutine invert_stress

