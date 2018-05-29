    Subroutine Scalar_RHS(s,u_m,v_m,w_m,qx,qy,qz,txz,tyz,dsdx,dsdy, &
    dsdz,RHS,Surf_flux)

    use globals
    use scalars
    implicit none

    interface
    include './interfaces/dealias1.f90'
    include './interfaces/dealias2.f90'
    include './interfaces/update1.f90'
    include './interfaces/ddx.f90'
    include './interfaces/ddy.f90'
    include './interfaces/ddz_w.f90'
    end interface

    integer*4 :: i,j,k

    real*8, dimension(:,:,:):: dsdx,dsdy,dsdz,RHS,s,txz,tyz, &
    qx,qy,qz,u_m,v_m,w_m

    real*8, dimension(:,:)::Surf_flux

    real*8, dimension(size(s,1),size(s,2),size(s,3)):: dtemp,temp

    real*8, dimension(size(u_m,1),size(u_m,2),size(u_m,3)):: &
    dsdx_m,dsdy_m,dsdz_m,cc,temp_m,s_m

    real*8 :: S_Surf(size(s,1),size(s,2)),l(size(s,3))

    call dealias1(dsdx,dsdx_m)
    call dealias1(dsdy,dsdy_m)
    call dealias1(dsdz,dsdz_m)

    call update1(dsdz_m)
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    Do k=2,Nzb+1
        Do j=1,Nyb2
            Do i=1,Nx2
                               
                cc(i,j,k) = u_m(i,j,k)*dsdx_m(i,j,k)+ &
                v_m(i,j,k)*dsdy_m(i,j,k)+ &
                (w_m(i,j,k)*dsdz_m(i,j,k)+ &
                w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2.
                               
            end do
        end do
    end do

!     ... TOP LIP

    if (vfact == vprocs-1 .AND. verticalBC == 0) then
        do j=1,Nyb2
            do i=1,Nx2

                cc(i,j,Nzb+1)= u_m(i,j,Nzb+1)*dsdx_m(i,j,Nzb+1)+ &
                v_m(i,j,Nzb+1)*dsdy_m(i,j,Nzb+1)+ &
                w_m(i,j,Nzb+1)*dsdz_m(i,j,Nzb+1)

            end do
        end do
    endif

    call dealias2(RHS,cc)

!     ... Now building the SGS part of the RHS.
!     ... Note: Since we bring the Conective term to RHS its sign changes.
!     ... Below "Temp" is used for SGS flux; its divergence is added to RHS

!     ... XXXXXXXXXXXXXXXXXXXX

    Do k=2,Nzb+1
        Do j=1,Nyb
            Do i=1,Nx
                temp(i,j,k)=qx(i,j,k)
            end do
        end do
    end do
         
    Call DDX(dtemp,temp)

    Do k=2,Nzb+1
        Do j=1,Nyb
            Do i=1,Nx
                RHS(i,j,k) = -RHS(i,j,k)-dtemp(i,j,k)
            end do
        end do
    end do

!     ... YYYYYYYYYYYYYYYYYYYY

    do k=2,Nzb+1
        Do j=1,Nyb
            Do i=1,Nx
                temp(i,j,k)=qy(i,j,k)
            end do
        end do
    end do
            
    Call DDY (dtemp,temp)

!...  And Smagorinsky in interior
!.....Note dsdz on W nodes...

    Do k=2,Nzb+1
        Do j=1,Nyb
            Do i=1,Nx
                RHS(i,j,k)=RHS(i,j,k)-dtemp(i,j,k)
            end do
        end do
    end do

!     ... ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    Do k=2,Nzb+2
        Do j=1,Nyb
            Do i=1,Nx
                temp(i,j,k)=qz(i,j,k)
            end do
        end do
    end do

!...  Use MO flux at wall

    if (vfact == 0 .AND. verticalBC == 0) then
        Do j=1,Nyb
            Do i=1,Nx
                temp(i,j,2)= surf_flux(i,j)
                qz(i,j,2)  = surf_flux(i,j)
            end do
        end do
    endif

!     ... The SGS_z flux is on the W nodes,
!     ... but DDZ_W will put it back on UVP nodes

    call DDZ_w (dtemp, temp)

    Do k=2,Nzb+1
        Do j=1,Nyb
            Do i=1,Nx
                                   
                RHS(i,j,k) = RHS(i,j,k)-dtemp(i,j,k)

            end do
        end do
    end do

    return

    end Subroutine Scalar_RHS
