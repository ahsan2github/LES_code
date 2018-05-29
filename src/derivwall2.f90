    subroutine derivwall2 (dtdz,dudz,dvdz,u,v,fi,fi_h,t_flux, &
    ustar,M)

    use globals
    use wallBoundaryConditions
    use scalars
    implicit none

    interface
    include './interfaces/filter_2dsl.f90'
    end interface

    integer*4 :: i,j
    integer*2 :: N,l
    real*8,dimension(:,:,:,:):: dtdz
    real*8,dimension(:,:,:):: dudz,dvdz,u,v,fi,fi_h,t_flux
    real*8,dimension(:,:):: ustar,M
    real*8,dimension(size(M,1),size(M,2))::u_hat,v_hat
    	  
!     Option 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$$$
!$$$      do j=1,nyb
!$$$         do i=1,nx
!$$$
!$$$            dvdz(i,j,2)=+0.
!$$$     +           +0.5*(v(i,j,3)+v(i,j,2)-0.)/dz
!$$$            dudz(i,j,2)=2.*dudz(i,j,3)*fi(i,j,1)/fi(i,j,2)
!$$$            dtdz(i,j,2)=2.*dtdz(i,j,3)*fi_h(i,j,1)/fi_h(i,j,2)
!$$$
!$$$         end do
!$$$      end do
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         	  
    call Filter_2dsl(u_hat,u(:,:,2))
    call Filter_2dsl(v_hat,v(:,:,2))
     
!     Option 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    do j=1,nyb
        do i=1,nx
            dudz(i,j,2)=fi(i,j,1)*ustar(i,j)*(u_hat(i,j)+Ugal)/ &
            (M(i,j)*vonk*0.5*dz)
            dvdz(i,j,2)=fi(i,j,1)*ustar(i,j)*(v_hat(i,j)+Vgal)/ &
            (M(i,j)*vonk*0.5*dz)
                         
            if(scalarCount >= 1)then
                do l=1,scalarCount
                    dtdz(i,j,2,l)=fi_h(i,j,1)*(-t_flux(i,j,l)/ustar(i,j))/ &
                    (vonk*0.5*dz)
                enddo
            endif
        enddo
    end do

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return
    end subroutine derivwall2
