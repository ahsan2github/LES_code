    subroutine ddz_uv_p (DFDZ,F)

!  First deriv in z direction for boundary layer (2nd order numerics)
!...F is on UVP nodes and dFdz is on w nodes
    use globals
    implicit none

    integer*4 :: i,j,k
    real*8,dimension(:,:,:)::F,dfdz

!     real*8 zz(size(F,3)), d_avg
              
!      do k=1,Nzb+2
!         zz(k)=(k-1.+vfact*nzb-0.5)*DZ
!      end do

    do k=2,Nzb+1
    !         d_avg=0.
        do j=1,Nyb
            do i=1,Nx
                dfdz(i,j,k)=(f(i,j,k)-f(i,j,k-1))*idz
            !               d_avg=d_avg+dfdz(i,j,k)*inxny
            end do
        end do
    enddo

!$$$         call plane_avg(dfdz,d_avg)
!$$$
!$$$         if(vfact.eq.0)then
!$$$            do j=1,nyb
!$$$               do i=1,nx
!$$$                  dfdz(i,j,3)=d_avg*dz/((zz(2)+0.5*dz)*
!$$$     +                 dlog(zz(3)/zz(2)))+(dfdz(i,j,3)-d_avg(3))
!$$$               end do
!$$$            end do
!$$$         end if

    return
    end subroutine ddz_uv_p
