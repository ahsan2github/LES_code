    subroutine wallstress2(txz,tyz,u,v,psi,psi0,zo,ustar,M)

    use globals
    use wallBoundaryConditions
    implicit none

    interface
    include './interfaces/filter_2dsl.f90'
    include './interfaces/plane_reduce.f90'
    end interface

    integer*4 ::  i,j,nx_decorr
    	
    real*8,dimension(:,:,:) :: txz,tyz,u,v,psi,psi0
    real*8,dimension(:,:) :: zo,ustar,M

    real*8,dimension(size(txz,1),size(txz,2)) :: U_res,u_hat,v_hat

    real*8 :: denom,ang_deg,factor,u_used,v_used,u_res_sum

!		angle to account for inclined structure in velocity/shear
!        correlation
!        ang_deg=13.
!		 factor=(dz/2.)/(tan(ang_deg*Pi/180.)*(2.*Pi/Nx2))
    factor=0
    nx_decorr=0
         
    call Filter_2dsl(u_hat,u(:,:,2))
    call Filter_2dsl(v_hat,v(:,:,2))

!	  u_hat=u(:,:,2)
!	  v_hat=v(:,:,2)

!       Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	  U_res_sum = 0.
!	  do j=1,Nyb
!	     do i=1,Nx
!	      U_res_sum = U_res_sum + ((u_hat(i,j)+Ugal)**2.+
!    +                    (v_hat(i,j)+Vgal)**2.)**0.5
!	   end do
!	  end do
!         call plane_reduce(U_res_sum)
!	  do j=1,Nyb
!	   do i=1,Nx
!	      U_res(i,j)=U_res_sum*iNxNy
!	   end do
!	  end do
!	  M(:,:) = U_res(:,:)

!       Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    do j=1,Nyb
        do i=1,Nx
            M(i,j) = ((u_hat(i,j)+Ugal)**2.+(v_hat(i,j)+Vgal)**2.)**0.5
        enddo
    enddo
!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    do j=1,Nyb
        do i=1,Nx

            if(i < Nx)then
                u_used=factor*(u_hat(i+1,j)+Ugal)+(1.-factor)* &
                (u_hat(i,j)+Ugal)
                v_used=factor*(v_hat(i+1,j)+Vgal)+(1.-factor)* &
                (v_hat(i,j)+Vgal)
            else
                u_used=factor*(u_hat(1,j)+Ugal)+(1.-factor)* &
                (u_hat(Nx,j)+Ugal)
                v_used=factor*(v_hat(1,j)+Vgal)+(1.-factor)* &
                (v_hat(Nx,j)+Vgal)
            endif

            denom=dlog((dz/2.)/zo(i,j))+Psi(i,j,1)-Psi0(i,j,1)
            ustar(i,j)= M(i,j)*vonk/denom
                
            txz(i,j,2)=-ustar(i,j)**2.*u_used/M(i,j)
            tyz(i,j,2)=-ustar(i,j)**2.*v_used/M(i,j)
            	      
        enddo
    enddo

    return
    end subroutine wallstress2
    
