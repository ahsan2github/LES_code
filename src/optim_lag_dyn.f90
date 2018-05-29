    Subroutine OPTIM_LAG_DYN (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat, &
    u_,v_,w_,L,LM_old,MM_old)

    use globals
    use scalars
    implicit none

    interface
    include './interfaces/filter_la.f90'
    include './interfaces/updateHorz.f90'
    include './interfaces/lagrng_dyn.f90'
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: i,j,k
    real*8,dimension(:,:,:):: &
    S11,S22,S33,S12,S13,S23,S,S_hat,w_,u_,v_,LM_old,MM_old,Cs2

    real*8,dimension(size(Cs2,1),size(Cs2,2),size(Cs2,3)):: &
    S11_hat,S22_hat,S33_hat,S12_hat,S13_hat,S23_hat,u_hat,v_hat, &
    w_hat,uu_hat,uv_hat,uw_hat,vv_hat,vw_hat,ww_hat, &
    SS11_hat,SS12_hat,SS13_hat,SS22_hat,SS23_hat,SS33_hat, &
    uu,uv,uw,vv,vw,ww,SS11,SS22,SS33,SS12,SS13,SS23

    real*8,dimension(size(Cs2,1),size(Cs2,2),size(Cs2,3)):: LM,MM

    real*8 :: L(:),M11,M12,M13,M22,M23,M33, &
    L13,L23,L12,L22,L33,L11, &
    ll,trace,Cs_ave,uin,vin,win

    real*8,dimension(size(Cs2,3)):: &
    LM_avg,MM_avg,Z

    real*8 :: a1,b1,c1,d1,e1

    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                 
                S(i,j,k)=sqrt(2.*(S11(i,j,k)**2.+S22(i,j,k)**2.+ &
                S33(i,j,k)**2.+2.*S12(i,j,k)**2.+ &
                &               2.*S13(i,j,k)**2.+ &
                &               2.*S23(i,j,k)**2.))

                SS11(i,j,k)=S(i,j,k)*S11(i,j,k)
                SS33(i,j,k)=S(i,j,k)*S33(i,j,k)
                SS22(i,j,k)=S(i,j,k)*S22(i,j,k)
                SS12(i,j,k)=S(i,j,k)*S12(i,j,k)
                SS13(i,j,k)=S(i,j,k)*S13(i,j,k)
                SS23(i,j,k)=S(i,j,k)*S23(i,j,k)
                               
                uu(i,j,k)=u_(i,j,k)**2.
                vv(i,j,k)=v_(i,j,k)**2.
                ww(i,j,k)=w_(i,j,k)**2.
                uv(i,j,k)=u_(i,j,k)*v_(i,j,k)
                vw(i,j,k)=v_(i,j,k)*w_(i,j,k)
                uw(i,j,k)=u_(i,j,k)*w_(i,j,k)

            end do
        end do
    end do
          
!C...Filtering to get the _hat variables (coarser resolution)

    Call Filter_La(u_hat,u_)
    Call Filter_La(v_hat,v_)
    Call Filter_La(w_hat,w_)
    Call Filter_La(uu_hat,uu)
    Call Filter_La(vv_hat,vv)
    Call Filter_La(ww_hat,ww)
    Call Filter_La(uv_hat,uv)
    Call Filter_La(uw_hat,uw)
    Call Filter_La(vw_hat,vw)
    Call Filter_La(S11_hat,S11)
    Call Filter_La(S22_hat,S22)
    Call Filter_La(S33_hat,S33)
    Call Filter_La(S12_hat,S12)
    Call Filter_La(S13_hat,S13)
    Call Filter_La(S23_hat,S23)
    Call Filter_La(SS11_hat,SS11)
    Call Filter_La(SS22_hat,SS22)
    Call Filter_La(SS33_hat,SS33)
    Call Filter_La(SS12_hat,SS12)
    Call Filter_La(SS13_hat,SS13)
    Call Filter_La(SS23_hat,SS23)

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx
                S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+ &
                S22_hat(i,j,k)**2.+S33_hat(i,j,k)**2.+ &
                &               2.*S12_hat(i,j,k)**2.+2.*S13_hat(i,j,k)**2.+ &
                &               2.*S23_hat(i,j,k)**2.))
            end do
        end do
    end do
          
    do k=2,Nzb+1
        do j=1,nyb
            do i=1,nx
                L11=(uu_hat(i,j,k))-(u_hat(i,j,k))**2.
                L22=(vv_hat(i,j,k))-(v_hat(i,j,k))**2.
                L12=(uv_hat(i,j,k))-(u_hat(i,j,k)*v_hat(i,j,k))
                L13=(uw_hat(i,j,k))-(u_hat(i,j,k)*w_hat(i,j,k))
                L23=(vw_hat(i,j,k))-(v_hat(i,j,k)*w_hat(i,j,k))
                L33=(ww_hat(i,j,k))-(w_hat(i,j,k))**2.

                a1=-2.*l(k)**2.*tfr**2.*S_hat(i,j,k)* &
                (L11*S11_hat(i,j,k)+L22*S22_hat(i,j,k)+ &
                L33*S33_hat(i,j,k)+2.*(L12*S12_hat(i,j,k)+L13* &
                S13_hat(i,j,k)+L23*S23_hat(i,j,k)))
                               
                b1=-2.*l(k)**2.*(L11*SS11_hat(i,j,k)+L22* &
                SS22_hat(i,j,k)+L33*SS33_hat(i,j,k)+ &
                &               2.*(L12*SS12_hat(i,j,k)+L13*SS13_hat(i,j,k) &
                +L23*SS23_hat(i,j,k)))
                               
                c1=(2.*l(k)**2.)**2.* &
                (SS11_hat(i,j,k)**2.+SS22_hat(i,j,k)**2. &
                +SS33_hat(i,j,k)**2.+2.*(SS12_hat(i,j,k)**2.+ &
                SS13_hat(i,j,k)**2.+SS23_hat(i,j,k)**2.))
                               
                d1=(4.*l(k)**4.)*tfr**4.*S_hat(i,j,k)**2.* &
                (S11_hat(i,j,k)**2.+S22_hat(i,j,k)**2. &
                +S33_hat(i,j,k)**2.+2.*(S12_hat(i,j,k)**2.+ &
                S13_hat(i,j,k)**2.+S23_hat(i,j,k)**2.))
                               
                e1=(8.*l(k)**4.)*tfr**2.* &
                S_hat(i,j,k)*(S11_hat(i,j,k)* &
                SS11_hat(i,j,k)+S22_hat(i,j,k)* &
                SS22_hat(i,j,k)+S33_hat(i,j,k)* &
                SS33_hat(i,j,k)+2.*(S12_hat(i,j,k)* &
                SS12_hat(i,j,k)+S13_hat(i,j,k)* &
                SS13_hat(i,j,k)+S23_hat(i,j,k)* &
                SS23_hat(i,j,k)))

                LM(i,j,k) = a1 - b1
                MM(i,j,k) = c1 + d1 - e1

            enddo
        enddo
    enddo

!     **************initilization for lagrangian averaging***********
    if(t <= 50 .AND. INITU /= 1)then

        call plane_avg(LM,LM_avg)
        call plane_avg(MM,MM_avg)
                
        do k=2,nzb+1
            do j=1,Nyb
                do i=1,Nx
                                      
                    LM(i,j,k)=LM_avg(k)
                    MM(i,j,k)=MM_avg(k)

                    LM_old(i,j+1,k)=LM(i,j,k)
                    MM_old(i,j+1,k)=MM(i,j,k)

                enddo
            enddo
        enddo

    endif
!     **************end of the lagrangian initilization**************

    call updateHorz(LM_old)
    call updateHorz(MM_old)

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

                if(vfact == 0 .AND. k <= 2 .AND. verticalBC == 0)then
                    uin = u_(i,j,k)
                    vin = v_(i,j,k)
                    win = w_hat(i,j,k)
                else
                    uin = u_(i,j,k)
                    vin = v_(i,j,k)
                    win = w_hat(i,j,k)
                endif

                call lagrng_dyn(LM(i,j,k),MM(i,j,k),LM_old,MM_old, &
                uin,vin,win,i,j,k)

                if(LM(i,j,k) < 0 .OR. MM(i,j,k) < (1e-10))then
                                      
                    cs2(i,j,k)=0.0
                                      
                    LM(i,j,k)=0.
                !                 MM(i,j,k)=0.
                                         
                else
                                      
                    cs2(i,j,k)=LM(i,j,k)/MM(i,j,k)
                                     
                endif
                                  
            enddo
        enddo
    enddo

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

                LM_old(i,j+1,k)=LM(i,j,k)
                MM_old(i,j+1,k)=MM(i,j,k)
                               
            enddo
        enddo
    enddo
          
    return
          
    end Subroutine OPTIM_LAG_DYN
          














