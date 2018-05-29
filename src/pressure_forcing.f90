    subroutine pressure_forcing(force_x,force_y,u,v,RHSx,RHSy,RHSx_f, &
    RHSy_f)

    use globals
    use press_force
    implicit none

    real*8,dimension(:) :: force_x, force_y
    real*8,dimension(:,:,:) :: u,v,RHSx,RHSy,RHSx_f,RHSy_f

    integer*4 :: k, k_min
    real*8 :: u_bar,v_bar,u_bar_mag,RHSx_bar,RHSy_bar,RHSxf_bar,RHSyf_bar
    real*8,dimension(size(u,1),size(u,2),size(u,3)) :: u_pln,v_pln, &
    RHSx_pln,RHSy_pln,RHSxf_pln,RHSyf_pln
    real*8,dimension(size(u,3)-2) :: u_avg_pln,v_avg_pln,u_avg_mag, &
    RHSx_avg_pln,RHSy_avg_pln,RHSxf_avg_pln,RHSyf_avg_pln
    real*8,dimension((size(u,3)-2)*vprocs) :: u_avg_vert,v_avg_vert, &
    RHSx_avg_vert,RHSy_avg_vert,RHSxf_avg_vert,RHSyf_avg_vert
    real*8,dimension(size(force_x)) :: force

! cc  Find average u, RHSx, RHSx_f in vertical dir

! ake plane averages
    call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
    call MPI_ALLREDUCE(u,u_pln,nx*nyb*nz2, &
    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEVEL,ierr)
    call MPI_ALLREDUCE(v,v_pln,nx*nyb*nz2, &
    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEVEL,ierr)
    call MPI_ALLREDUCE(RHSx,RHSx_pln,nx*nyb*nz2, &
    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEVEL,ierr)
    call MPI_ALLREDUCE(RHSy,RHSy_pln,nx*nyb*nz2, &
    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEVEL,ierr)
    call MPI_ALLREDUCE(RHSx_f,RHSxf_pln,nx*nyb*nz2, &
    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEVEL,ierr)
    call MPI_ALLREDUCE(RHSy_f,RHSyf_pln,nx*nyb*nz2, &
    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEVEL,ierr)

    do k=2,nzb+1
        u_avg_pln(k-1)=sum(u_pln(:,:,k))*inxny
        v_avg_pln(k-1)=sum(v_pln(:,:,k))*inxny
        RHSx_avg_pln(k-1)=sum(RHSx_pln(:,:,k))*inxny
        RHSy_avg_pln(k-1)=sum(RHSy_pln(:,:,k))*inxny
        RHSxf_avg_pln(k-1)=sum(RHSxf_pln(:,:,k))*inxny
        RHSyf_avg_pln(k-1)=sum(RHSyf_pln(:,:,k))*inxny
    enddo

! pi_allgather averages
    call MPI_BARRIER(MPI_COMM_COLUMN,ierr)
    call MPI_ALLGATHER(u_avg_pln(1),nzb,MPI_DOUBLE_PRECISION, &
    u_avg_vert,nzb,MPI_DOUBLE_PRECISION,MPI_COMM_COLUMN,ierr)
    call MPI_ALLGATHER(v_avg_pln(1),nzb,MPI_DOUBLE_PRECISION, &
    v_avg_vert,nzb,MPI_DOUBLE_PRECISION,MPI_COMM_COLUMN,ierr)
    call MPI_ALLGATHER(RHSx_avg_pln(1),nzb,MPI_DOUBLE_PRECISION, &
    RHSx_avg_vert,nzb,MPI_DOUBLE_PRECISION,MPI_COMM_COLUMN, &
    ierr)
    call MPI_ALLGATHER(RHSy_avg_pln(1),nzb,MPI_DOUBLE_PRECISION, &
    RHSy_avg_vert,nzb,MPI_DOUBLE_PRECISION,MPI_COMM_COLUMN, &
    ierr)
    call MPI_ALLGATHER(RHSxf_avg_pln(1),nzb,MPI_DOUBLE_PRECISION, &
    RHSxf_avg_vert,nzb,MPI_DOUBLE_PRECISION,MPI_COMM_COLUMN, &
    ierr)
    call MPI_ALLGATHER(RHSyf_avg_pln(1),nzb,MPI_DOUBLE_PRECISION, &
    RHSyf_avg_vert,nzb,MPI_DOUBLE_PRECISION,MPI_COMM_COLUMN, &
    ierr)

! ertical average
    u_bar = sum(u_avg_vert)/nz
    v_bar = sum(v_avg_vert)/nz
    RHSx_bar = sum(RHSx_avg_vert)/nz
    RHSy_bar = sum(RHSy_avg_vert)/nz
    RHSxf_bar = sum(RHSxf_avg_vert)/nz
    RHSyf_bar = sum(RHSyf_avg_vert)/nz

! cc  Calculate required dpdx based on Adams-Bashforth scheme

    force_x(2:nzb+1) = (u_avg*cos(theta_mean_wind*pi/180.d0) &
    -u_bar)/dt - 1.5d0*RHSx_bar+0.5d0*RHSxf_bar
    force_y(2:nzb+1) = (u_avg*sin(theta_mean_wind*pi/180.d0) &
    -v_bar)/dt - 1.5d0*RHSy_bar+0.5d0*RHSyf_bar

    return
    end subroutine pressure_forcing
