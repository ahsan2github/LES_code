    subroutine calc_lambda2(lambda2,dudx,dudy,dudz,dvdx,dvdy,dvdz, &
    dwdx,dwdy,dwdz)

    use globals
    implicit none

    integer*4 :: i,j,k
    real*8,dimension(:,:,:) :: lambda2,dudx,dudy,dudz,dvdx,dvdy,dvdz, &
    dwdx,dwdy,dwdz

    real*8 ::m,p1,det,q,p,r,phi,e1,e3
    real*8,dimension(size(lambda2,1),size(lambda2,2),size(lambda2,3)) &
    :: ux,uy,uz,vx,vy,vz,wx,wy,wz, &
    S11,S22,S33,S12,S13,S23, &
    O12,O21,O13,O31,O23,O32, &
    M11,M22,M33,M12,M13,M23

! interpolations to u-v-p nodes
    do k=2,nzb+1
        ux(:,:,k)=     dudx(:,:,k)
        uy(:,:,k)=     dudy(:,:,k)
        uz(:,:,k)=0.5d0*(dudz(:,:,k)+dudz(:,:,k+1))
        vx(:,:,k)=     dvdx(:,:,k)
        vy(:,:,k)=     dvdy(:,:,k)
        vz(:,:,k)=0.5d0*(dvdz(:,:,k)+dvdz(:,:,k+1))
        wx(:,:,k)=0.5d0*(dwdx(:,:,k)+dwdx(:,:,k+1))
        wy(:,:,k)=0.5d0*(dwdy(:,:,k)+dwdy(:,:,k+1))
        wz(:,:,k)=     dwdz(:,:,k)
    end do

! bottom boundary condition
    if (vfact == 0 .AND. verticalBC == 0) then
        ux(:,:,2)=     dudx(:,:,2)
        uy(:,:,2)=     dudy(:,:,2)
        uz(:,:,2)=     dudz(:,:,2)
        vx(:,:,2)=     dvdx(:,:,2)
        vy(:,:,2)=     dvdy(:,:,2)
        vz(:,:,2)=     dvdz(:,:,2)
        wx(:,:,2)=0.5d0*(dwdx(:,:,3)+dwdx(:,:,2))
        wy(:,:,2)=0.5d0*(dwdy(:,:,3)+dwdy(:,:,2))
        wz(:,:,2)=     dwdz(:,:,2)
    endif

! cc  Calculate S_ij
    S11(:,:,2:nzb+1)=0.5*(ux(:,:,2:nzb+1)+ux(:,:,2:nzb+1))
    S33(:,:,2:nzb+1)=0.5*(wz(:,:,2:nzb+1)+wz(:,:,2:nzb+1))
    S22(:,:,2:nzb+1)=0.5*(vy(:,:,2:nzb+1)+vy(:,:,2:nzb+1))
    S12(:,:,2:nzb+1)=0.5*(uy(:,:,2:nzb+1)+vx(:,:,2:nzb+1))
    S13(:,:,2:nzb+1)=0.5*(uz(:,:,2:nzb+1)+wx(:,:,2:nzb+1))
    S23(:,:,2:nzb+1)=0.5*(vz(:,:,2:nzb+1)+wy(:,:,2:nzb+1))

! cc  Calculate Omega_ij  (note antisymmetric with zero diagonal)
    O12(:,:,2:nzb+1)=0.5*(uy(:,:,2:nzb+1)-vx(:,:,2:nzb+1))
    O21=-O12
    O13(:,:,2:nzb+1)=0.5*(uz(:,:,2:nzb+1)-wx(:,:,2:nzb+1))
    O31=-O13
    O23(:,:,2:nzb+1)=0.5*(vz(:,:,2:nzb+1)-wy(:,:,2:nzb+1))
    O32=-O23

! cc  S^2+O^2
    M11(:,:,2:nzb+1)=S11(:,:,2:nzb+1)**2+S12(:,:,2:nzb+1)**2+ &
    S13(:,:,2:nzb+1)**2+ &
    O12(:,:,2:nzb+1)*O21(:,:,2:nzb+1)+ &
    O13(:,:,2:nzb+1)*O31(:,:,2:nzb+1)
    M22(:,:,2:nzb+1)=S12(:,:,2:nzb+1)**2+S22(:,:,2:nzb+1)**2+ &
    S23(:,:,2:nzb+1)**2+ &
    O12(:,:,2:nzb+1)*O21(:,:,2:nzb+1)+ &
    O23(:,:,2:nzb+1)*O32(:,:,2:nzb+1)
    M33(:,:,2:nzb+1)=S13(:,:,2:nzb+1)**2+S23(:,:,2:nzb+1)**2+ &
    S33(:,:,2:nzb+1)**2+ &
    O13(:,:,2:nzb+1)*O31(:,:,2:nzb+1)+ &
    O23(:,:,2:nzb+1)*O32(:,:,2:nzb+1)
    M12(:,:,2:nzb+1)=S11(:,:,2:nzb+1)*S12(:,:,2:nzb+1)+ &
    S12(:,:,2:nzb+1)*S22(:,:,2:nzb+1)+ &
    S13(:,:,2:nzb+1)*S23(:,:,2:nzb+1)+ &
    O13(:,:,2:nzb+1)*O32(:,:,2:nzb+1)
    M13(:,:,2:nzb+1)=S11(:,:,2:nzb+1)*S13(:,:,2:nzb+1)+ &
    S12(:,:,2:nzb+1)*S23(:,:,2:nzb+1)+ &
    S13(:,:,2:nzb+1)*S33(:,:,2:nzb+1)+ &
    O12(:,:,2:nzb+1)*O23(:,:,2:nzb+1)
    M23(:,:,2:nzb+1)=S12(:,:,2:nzb+1)*S13(:,:,2:nzb+1)+ &
    S22(:,:,2:nzb+1)*S23(:,:,2:nzb+1)+ &
    S23(:,:,2:nzb+1)*S33(:,:,2:nzb+1)+ &
    O21(:,:,2:nzb+1)*O13(:,:,2:nzb+1)

! cc eigenvector of M_ij
    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx
                m = (M11(i,j,k)+M22(i,j,k)+M33(i,j,k))/3;
                p1 = M12(i,j,k)**2 + M13(i,j,k)**2 + M23(i,j,k)**2
                if (p1 == 0)then
                    lambda2(i,j,k) = M22(i,j,k)
                else
                    det=(M11(i,j,k)-m)*(M22(i,j,k)-m)*(M33(i,j,k)-m)+ &
                    &                  2.d0*M12(i,j,k)*M23(i,j,k)*M13(i,j,k)- &
                    M13(i,j,k)**2*(M22(i,j,k)-m)- &
                    M12(i,j,k)**2*(M33(i,j,k)-m)- &
                    M23(i,j,k)**2*(M11(i,j,k)-m)
                    q = det / 2.d0;
                    p = sqrt((2.d0*p1+(M11(i,j,k)-m)**2+(M22(i,j,k)-m)**2+ &
                    (M33(i,j,k)-m)**2) / 6)
                    r = q / p**3;
                     
                ! theoretically for a real symmetric matrix, -1 <= r <= 1. But
                ! computation error can result in an r slightly outside that range.
                    if (r <= -1)then
                        phi = pi / 3.d0
                    elseif (r >= 1)then
                        phi = 0.d0
                    else
                        phi = acos(r) / 3.d0
                    endif
                     
                ! the eigenvalues satisfy eig3 <= eig2 <= eig1
                    e1 = m + 2.d0 * p * cos(phi)
                    e3 = m + 2.d0 * p * cos(phi + pi * (2.d0/3.d0))
                    lambda2(i,j,k) = 3.d0 * m - e1 - e3
                endif
            enddo
        enddo
    enddo


    return
    end subroutine calc_lambda2
