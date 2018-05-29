    subroutine sponge_layer(scalar,u,v,w,scalarRHS,RHSx,RHSy,RHSz, &
    rdmp,flag)
! lag=1 for scalar sponge layer, flag=2 for momenturm sponge layer
          
    use globals
    use scalars
    use wallBoundaryConditions
    implicit none

    interface
    include './interfaces/update1.f90'
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: flag
    real*8,dimension(:,:,:) :: rdmp,u,v,w,RHSx,RHSy,RHSz
    real*8,dimension(:,:,:,:) :: scalar,scalarRHS

    integer*4 :: i,j,k,l
    real*8 :: cfrdmp,ztemp
    real*8,dimension(size(scalar,3)) :: t_bar,u_bar,v_bar
          
    if(flag == 1)then

        cfrdmp = 1./(rlx_time*u_star/z_i)
        do k=2,Nzb+1
            ztemp=(vfact*Nzb+(k-2.))*dz*z_i
            if (ztemp >= z_d .AND. ztemp <= l_z) then
                rdmp(:,:,k)=cfrdmp*0.5* &
                (1.0-COS(pi*(ztemp-z_d)/(l_z-z_d)))
            else
                rdmp(:,:,k)=0.0
            end if
        end do
                 
        if (vfact == (vprocs-1)) then
            do j=1,Nyb
                do i=1,Nx
                    rdmp(i,j,Nzb+2)=rdmp(i,j,Nzb+1)
                end do
            end do
        end if
        call update1(rdmp)
                
        do l=1,scalarCount     !do for each scalar
                           
            call plane_avg(scalar(:,:,:,l),t_bar)
                           
            do j=1,Nyb
                do i=1,Nx
                    scalarRHS(i,j,k,l)=scalarRHS(i,j,k,l)- &
                    &                  0.5*(rdmp(i,j,k)+rdmp(i,j,k+1))* &
                    (scalar(i,j,k,l)-t_bar(k))
                end do
            end do
        end do
                 
    elseif(flag == 2)then

        if (sponge == 1) then
            do k=2,Nzb+1
                do j=1,Nyb
                    do i=1,Nx
                        RHSx(i,j,k)=RHSx(i,j,k)-0.5* &
                        (rdmp(i,j,k)+rdmp(i,j,k+1))* &
                        (u(i,j,k)-(Ug-Ugal))
                        RHSy(i,j,k)=RHSy(i,j,k)-0.5* &
                        (rdmp(i,j,k)+rdmp(i,j,k+1))* &
                        (v(i,j,k)-(Vg-Vgal))
                    end do
                end do
            end do
            do k=2,Nzb+1
                do j=1,Nyb
                    do i=1,Nx
                        RHSz(i,j,k)=RHSz(i,j,k)-rdmp(i,j,k)*w(i,j,k)
                    end do
                end do
            end do
        elseif (sponge == 2) then
            call plane_avg(u,u_bar)
            call plane_avg(v,v_bar)
            do k=2,Nzb+1
                do j=1,Nyb
                    do i=1,Nx
                        RHSx(i,j,k)=RHSx(i,j,k)-0.5* &
                        (rdmp(i,j,k)+rdmp(i,j,k+1))* &
                        (u(i,j,k)-u_bar(k))
                        RHSy(i,j,k)=RHSy(i,j,k)-0.5* &
                        (rdmp(i,j,k)+rdmp(i,j,k+1))* &
                        (v(i,j,k)-v_bar(k))
                    end do
                end do
            end do
            do k=2,Nzb+1
                do j=1,Nyb
                    do i=1,Nx
                        RHSz(i,j,k)=RHSz(i,j,k)-rdmp(i,j,k)*w(i,j,k)
                    end do
                end do
            end do
        end if


    endif
         
          

    return
    end subroutine sponge_layer
