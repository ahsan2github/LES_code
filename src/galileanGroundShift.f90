    subroutine galileanGroundShift(gndScalars,soilHeatFlux, &
    ustar,scalarFlux,coolrate,zo,Psi,Psi0,fi,fi_H,ilow,jlow)

    use globals
    use SEBmodule
    use scalars
    implicit none

    integer*4, dimension(:) :: ilow,jlow
    real*8, dimension(:,:) :: zo,ustar,soilHeatFlux
    real*8, dimension(:,:,:) :: scalarFlux,coolrate,Psi,Psi0,fi,fi_H
    real*8, dimension(:,:,:,:) :: gndScalars

    integer*4 :: i,j,ind
    integer*4 :: il,jl,swtx, swty, ddd,xAbsNodes, xRelNodes, yAbsNodes, &
    yRelNodes
    real*8, dimension(size(zo,1),size(zo,2)) :: zo_,soilHeatFlux_, &
    ustar_
    real*8, dimension(size(gndScalars,1),size(gndScalars,2), &
    size(gndScalars,3),2):: gndScalars_
    real*8, dimension(size(scalarFlux,1),size(scalarFlux,2), &
    size(scalarFlux,3))::scalarFlux_,coolrate_
    real*8, dimension(size(Psi,1),size(Psi,2),2)::Psi_,Psi0_,fi_,fi_H_
    real*8 :: mhfx, mhfy, xNodes, yNodes
    real*8 :: stepsPerDx, stepsPerDy
    integer*4 :: dxIntervals, lastDxIntervals, &
    dyIntervals, lastDyIntervals,s_row,r_row,send_to,recv_from
    real*8,dimension(size(gndScalars,1)) :: send1,recv1
    real*8,dimension(size(gndScalars,1),size(gndScalars,3)) :: send2, &
    recv2
    real*8,dimension(size(gndScalars,1),size(gndScalars,3),2) :: &
    send3,recv3
    save swtx, swty, stepsPerDx, stepsPerDy, &
    lastDxIntervals, lastDyIntervals

    xNodes = floor((dble(t+nrsub)*dt*Ugal)/dx)
    xAbsNodes = int(xNodes)
    xRelNodes = mod(xAbsNodes,nx) ! ?? what is this
    ddd       = nint(10*dmod((dble(t+nrsub)*dt*Ugal),dx)/dx)
    do i=1,nx
        il = i-xRelNodes
        if(il < 1)then
            il = il+nx
        endif
        ilow(i) = il
    enddo

    yNodes    = floor((dble(t+nrsub)*dt*Vgal)/dy)
    yAbsNodes = int(yNodes)
    yRelNodes = mod(yAbsNodes,ny)
    ddd       = nint(10*dmod((dble(t+nrsub)*dt*Vgal),dy)/dy)
    do j=1,ny
        jl = j - yRelNodes
        if(jl < 1)then
            jl = jl + ny
        endif
        jlow(j)=jl
    enddo

    lastDxIntervals = floor(((t+nrsub-1)*dt*abs(Ugal) - dx/2)/dx)
    lastDyIntervals = floor(((t+nrsub-1)*dt*abs(Vgal) - dy/2)/dy)

    dxIntervals = floor(((t+nrsub)*dt*abs(Ugal) - dx/2)/dx)
    dyIntervals = floor(((t+nrsub)*dt*abs(Vgal) - dy/2)/dy)

    if( dxIntervals > lastDxIntervals )then
        if(Ugal /= 0.)then
                        
        !     shift ground to account for galilean shift
        !     so that ground info is below correct atmosphere profile
        !     do for x then y

            do i=1,nx

                if( Ugal < 0 )then ! Ugal < 0 => shift ground positive direction
                    if( i == 1 )then
                        ind = nx
                    else
                        ind = i - 1
                    endif
                elseif( Ugal > 0 )then
                    if( i == nx )then
                        ind = 1
                    else
                        ind = i + 1
                    endif
                endif

                gndScalars_(i,:,:,:) = gndScalars(ind,:,:,:)
                soilHeatFlux_(i,:)   = soilHeatFlux(ind,:)
                ustar_(i,:)          = ustar(ind,:)
                scalarFlux_(i,:,:)   = scalarFlux(ind,:,:)
                coolrate_(i,:,:)     = coolrate(ind,:,:)
                zo_(i,:)             = zo(ind,:)
                Psi_(i,:,:)          = Psi(ind,:,:)
                Psi0_(i,:,:)         = Psi0(ind,:,:)
                fi_(i,:,:)           = fi(ind,:,:)
                fi_H_(i,:,:)         = fi_H(ind,:,:)
                             
            enddo

        endif

        gndScalars   = gndScalars_
        soilHeatFlux = soilHeatFlux_
        ustar        = ustar_
        scalarFlux   = scalarFlux_
        coolrate     = coolrate_
        zo           = zo_
        Psi          = Psi_
        Psi0         = Psi0_
        fi           = fi_
        fi_H         = fi_H_
                 
    endif

    if( dyIntervals > lastDyIntervals )then

        if(me == 0)write(*,*)'shifting y'

    ! apply y shift from Vgal

        if(hprocs > 1)then
            if(Vgal < 0)then
                s_row=nyb
                r_row=1
                if(me == 0)then
                    send_to = me+1
                    recv_from = hprocs-1
                elseif(me == hprocs-1)then
                    send_to = 0
                    recv_from = me-1
                else
                    send_to=me+1
                    recv_from = me-1
                endif
            elseif(Vgal > 0)then
                s_row=1
                r_row=nyb
                if(me == 0)then
                    send_to = hprocs-1
                    recv_from = me+1
                elseif(me == hprocs-1)then
                    send_to = me-1
                    recv_from = 0
                else
                    send_to=me-1
                    recv_from = me+1
                endif
            endif
                           
            send3=gndScalars(:,s_row,:,:)
            call MPI_SEND(send3(1,1,1),size(send3), &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv3(1,1,1),size(recv3), &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            gndScalars_(:,r_row,:,:)=recv3
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send1=soilHeatFlux(:,s_row)
            call MPI_SEND(send1(1),nx, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv1(1),nx, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            soilHeatFlux_(:,r_row)=recv1
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send1=ustar(:,s_row)
            call MPI_SEND(send1(1),nx, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv1(1),nx, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            ustar_(:,r_row)=recv1
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send2=scalarFlux(:,s_row,:)
            call MPI_SEND(send2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            scalarFlux_(:,r_row,:)=recv2
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send2=coolrate(:,s_row,:)
            call MPI_SEND(send2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            coolrate_(:,r_row,:)=recv2
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send1=zo(:,s_row)
            call MPI_SEND(send1(1),nx, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv1(1),nx, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            zo_(:,r_row)=recv1
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send2=psi(:,s_row,:)
            call MPI_SEND(send2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            psi_(:,r_row,:)=recv2
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send2=psi0(:,s_row,:)
            call MPI_SEND(send2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            psi0_(:,r_row,:)=recv2
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send2=fi(:,s_row,:)
            call MPI_SEND(send2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            fi_(:,r_row,:)=recv2
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
            send2=fi_H(:,s_row,:)
            call MPI_SEND(send2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,send_to,me,MPI_COMM_LEVEL, &
            MPI_STATUS_IGNORE,ierr)
            call MPI_RECV(recv2(1,1),nx*scalarCount, &
            MPI_DOUBLE_PRECISION,recv_from,recv_from, &
            MPI_COMM_LEVEL,MPI_STATUS_IGNORE,ierr)
            fi_H_(:,r_row,:)=recv2
            call MPI_BARRIER(MPI_COMM_LEVEL,ierr)
                           
        endif
                    
        do j=1,nyb
                           
            if( Vgal < 0)then ! Vgal < 0 => shift ground positive direction
                if( j == 1)then
                    ind = ny
                else
                    ind = j - 1
                endif

                if(hprocs == 1 .OR. j /= 1)then
                    gndScalars_(:,j,:,:) = gndScalars(:,ind,:,:)
                    soilHeatFlux_(:,j)   = soilHeatFlux(:,ind)
                    ustar_(:,j)          = ustar(:,ind)
                    scalarFlux_(:,j,:)   = scalarFlux(:,ind,:)
                    coolrate_(:,j,:)     = coolrate(:,ind,:)
                    zo_(:,j)             = zo(:,ind)
                    Psi_(:,j,:)          = Psi(:,ind,:)
                    Psi0_(:,j,:)         = Psi0(:,ind,:)
                    fi_(:,j,:)           = fi(:,ind,:)
                    fi_H_(:,j,:)         = fi_H(:,ind,:)
                endif
                                  
            elseif( Vgal > 0 )then
                if( j == ny)then
                    ind = 1
                else
                    ind = j + 1
                endif
                                  
                if(hprocs == 1 .OR. j /= nyb)then
                    gndScalars_(:,j,:,:) = gndScalars(:,ind,:,:)
                    soilHeatFlux_(:,j)   = soilHeatFlux(:,ind)
                    ustar_(:,j)          = ustar(:,ind)
                    scalarFlux_(:,j,:)   = scalarFlux(:,ind,:)
                    coolrate_(:,j,:)     = coolrate(:,ind,:)
                    zo_(:,j)             = zo(:,ind)
                    Psi_(:,j,:)          = Psi(:,ind,:)
                    Psi0_(:,j,:)         = Psi0(:,ind,:)
                    fi_(:,j,:)           = fi(:,ind,:)
                    fi_H_(:,j,:)         = fi_H(:,ind,:)
                endif
                                  
            endif
                           
        enddo

        gndScalars   = gndScalars_
        soilHeatFlux = soilHeatFlux_
        ustar        = ustar_
        scalarFlux   = scalarFlux_
        coolrate     = coolrate_
        zo           = zo_
        Psi          = Psi_
        Psi0         = Psi0_
        fi           = fi_
        fi_H         = fi_H_
                    
    endif

    lastDxIntervals = dxIntervals
    lastDyIntervals = dyIntervals
    end subroutine galileanGroundShift
