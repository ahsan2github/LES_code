    Subroutine OPTIM_WONG (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd, &
    u_,v_,w_,L,betaa)

    use globals
    use scalars
    use sgs
    implicit none

    interface
    include './interfaces/filter_2laa.f90'
    include './interfaces/updateHorz.f90'
    include './interfaces/plane_reduce.f90'
    end interface

    integer*4 :: i,j,k,Mx,My,ii,jj
    real*8,dimension(:,:,:):: &
    S,S11,S22,S33,S12,S13,S23,S_hat,S_hatd,w_,u_,v_,Cs2,betaa

    real*8,dimension(size(Cs2,1),size(Cs2,2),size(Cs2,3)):: &
    S11_hatd,S22_hatd,S33_hatd,S12_hatd,S13_hatd,S23_hatd,u_hatd, &
    v_hatd,w_hatd,uu_hatd,uv_hatd,uw_hatd,vv_hatd,vw_hatd, &
    ww_hatd,SS11_hat,SS12_hat,SS13_hat,SS22_hat,SS23_hat, &
    SS33_hat,SS11_hatd,SS12_hatd,SS13_hatd,SS22_hatd,SS23_hatd, &
    SS33_hatd,uu,uv,uw,vv,vw,ww,SS11,SS22,SS33,SS12,SS13,SS23

    real*8,dimension(size(Cs2,1),size(Cs2,2)+2,size(Cs2,3)):: &
    S11_hat,S22_hat,S33_hat,S12_hat,S13_hat,S23_hat,u_hat,v_hat, &
    w_hat,uu_hat,uv_hat,uw_hat,vv_hat,vw_hat,ww_hat

    real*8 :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10, &
    aa,bb,cc,dd,ee,ff,betaa1D(size(Cs2,3)),rtnewt

    real*8 :: L(:),M11,M12,M13,M22,M23,M33, &
    LM,MM,L13,L23,L12,L22,L33,L11,Q13,Q23,Q12,Q22,Q33,Q11

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
          
!     C...Filtering to get the _hat variables (coarser resolution)
          
    Call Filter_2Laa(u_hat,u_hatd,u_)
    Call Filter_2Laa(v_hat,v_hatd,v_)
    Call Filter_2Laa(w_hat,w_hatd,w_)
    Call Filter_2Laa(uu_hat,uu_hatd,uu)
    Call Filter_2Laa(vv_hat,vv_hatd,vv)
    Call Filter_2Laa(ww_hat,ww_hatd,ww)
    Call Filter_2Laa(uv_hat,uv_hatd,uv)
    Call Filter_2Laa(uw_hat,uw_hatd,uw)
    Call Filter_2Laa(vw_hat,vw_hatd,vw)
    Call Filter_2Laa(S11_hat,S11_hatd,S11)
    Call Filter_2Laa(S22_hat,S22_hatd,S22)
    Call Filter_2Laa(S33_hat,S33_hatd,S33)
    Call Filter_2Laa(S12_hat,S12_hatd,S12)
    Call Filter_2Laa(S13_hat,S13_hatd,S13)
    Call Filter_2Laa(S23_hat,S23_hatd,S23)
    Call Filter_2Laa(SS11_hat,SS11_hatd,SS11)
    Call Filter_2Laa(SS22_hat,SS22_hatd,SS22)
    Call Filter_2Laa(SS33_hat,SS33_hatd,SS33)
    Call Filter_2Laa(SS12_hat,SS12_hatd,SS12)
    Call Filter_2Laa(SS13_hat,SS13_hatd,SS13)
    if(scalarCount >= 1)then
        Call Filter_2Laa(SS23_hat,SS23_hatd,SS23)
    else
        Call Filter_2Laa(SS23_hat,SS23_hatd,SS23)
    endif
          
    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx
                S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+ &
                S22_hat(i,j,k)**2.+S33_hat(i,j,k)**2.+ &
                &               2.*S12_hat(i,j,k)**2.+2.*S13_hat(i,j,k)**2.+ &
                &               2.*S23_hat(i,j,k)**2.))
                S_hatd(i,j,k)=sqrt(2.*(S11_hatd(i,j,k)**2.+ &
                S22_hatd(i,j,k)**2.+S33_hatd(i,j,k)**2.+ &
                &               2.*S12_hatd(i,j,k)**2.+2.*S13_hatd(i,j,k)**2.+ &
                &               2.*S23_hatd(i,j,k)**2.))
            end do
        end do
    end do
          
    do k=2,Nzb+1

        a1=0.
        a2=0.
        a3=0.
        a4=0.
        a5=0.
        a6=0.
        a7=0.
        a8=0.
        a9=0.
        a10=0.
                 
        do j=1,nyb
            do i=1,nx
                L11=(uu_hat(i,j,k))-(u_hat(i,j,k))**2.
                L22=(vv_hat(i,j,k))-(v_hat(i,j,k))**2.
                L33=(ww_hat(i,j,k))-(w_hat(i,j,k))**2.
            !               L33=ww_hat(i,j,k)-w_hat2(i,j,k)
                L12=(uv_hat(i,j,k))-(u_hat(i,j,k)*v_hat(i,j,k))
                L13=(uw_hat(i,j,k))-(u_hat(i,j,k)*w_hat(i,j,k))
                L23=(vw_hat(i,j,k))-(v_hat(i,j,k)*w_hat(i,j,k))
                               
                Q11=(uu_hatd(i,j,k))-(u_hatd(i,j,k))**2.
                Q22=(vv_hatd(i,j,k))-(v_hatd(i,j,k))**2.
                Q33=(ww_hatd(i,j,k))-(w_hatd(i,j,k))**2.
            !               Q33=ww_hatd(i,j,k)-w_hatd2(i,j,k)
                Q12=(uv_hatd(i,j,k))-(u_hatd(i,j,k)*v_hatd(i,j,k))
                Q13=(uw_hatd(i,j,k))-(u_hatd(i,j,k)*w_hatd(i,j,k))
                Q23=(vw_hatd(i,j,k))-(v_hatd(i,j,k)*w_hatd(i,j,k))
                       
                a1=a1+(Q11*S11_hatd(i,j,k)+ &
                Q22*S22_hatd(i,j,k)+ &
                Q33*S33_hatd(i,j,k)+ &
                &               2.*(Q12*S12_hatd(i,j,k)+ &
                Q23*S23_hatd(i,j,k)+ &
                Q13*S13_hatd(i,j,k)))

                a2=a1*(-tfr**(8./3.))

                a3=a3+(S11_hat(i,j,k)**2.+ &
                S22_hat(i,j,k)**2.+ &
                S33_hat(i,j,k)**2.+ &
                &               2.*(S12_hat(i,j,k)**2.+ &
                S23_hat(i,j,k)**2.+ &
                S13_hat(i,j,k)**2.))

                a4=(-2.*tfr**(4./3.))*a3
                     
                a5=(tfr**(8./3.))*a3

                a6=a6+(L11*S11_hat(i,j,k)+ &
                L22*S22_hat(i,j,k)+ &
                L33*S33_hat(i,j,k)+ &
                &               2.*(L12*S12_hat(i,j,k)+ &
                L23*S23_hat(i,j,k)+ &
                L13*S13_hat(i,j,k)))

                a7=a6*(-tfr**(4./3.))
                               
                a8=a8+(S11_hatd(i,j,k)**2.+ &
                S22_hatd(i,j,k)**2.+ &
                S33_hatd(i,j,k)**2.+ &
                &               2.*(S12_hatd(i,j,k)**2.+ &
                S23_hatd(i,j,k)**2.+ &
                S13_hatd(i,j,k)**2.))

                a9=(-2.*tfr**(8./3.))*a8

                a10=(tfr**(16./3.))*a8

            end do
        end do


        call plane_reduce(a1)
        call plane_reduce(a2)
        call plane_reduce(a3)
        call plane_reduce(a4)
        call plane_reduce(a5)
        call plane_reduce(a6)
        call plane_reduce(a7)
        call plane_reduce(a8)
        call plane_reduce(a9)
        call plane_reduce(a10)

        aa = a1*a3-a6*a8
        bb = a1*a4-a7*a8
        cc = a2*a3+a1*a5-a6*a9
        dd = a2*a4-a7*a9
        ee = a2*a5-a6*a10
        ff = -a7*a10

        if (t <= 1000) then
            call root8(rtnewt,aa,bb,cc,dd,ee,ff)
        else
            call newroots(rtnewt,aa,bb,cc,dd,ee,ff)
        end if

        betaa1D(k)=rtnewt
         
        if(betaa1D(k) <= 0 .OR. betaa1D(k) > 5.0) then
            betaa1D(k)=1.0
        end if

        if(model == 2)then
            betaa1D(k)=1.
        end if
    end do

!     Create and Update Ghost Layers*************

    do j=1,nyb
        jj=abs(j-nyb-1)
        u_hat(:,jj+1,:) = u_hat(:,jj,:)
        v_hat(:,jj+1,:) = v_hat(:,jj,:)
        w_hat(:,jj+1,:) = w_hat(:,jj,:)
        uu_hat(:,jj+1,:) = uu_hat(:,jj,:)
        vv_hat(:,jj+1,:) = vv_hat(:,jj,:)
        ww_hat(:,jj+1,:) = ww_hat(:,jj,:)
        uv_hat(:,jj+1,:) = uv_hat(:,jj,:)
        uw_hat(:,jj+1,:) = uw_hat(:,jj,:)
        vw_hat(:,jj+1,:) = vw_hat(:,jj,:)
        S11_hat(:,jj+1,:) = S11_hat(:,jj,:)
        S22_hat(:,jj+1,:) = S22_hat(:,jj,:)
        S33_hat(:,jj+1,:) = S33_hat(:,jj,:)
        S12_hat(:,jj+1,:) = S12_hat(:,jj,:)
        S13_hat(:,jj+1,:) = S13_hat(:,jj,:)
        S23_hat(:,jj+1,:) = S23_hat(:,jj,:)
    end do

    call updateHorz(u_hat)
    call updateHorz(v_hat)
    call updateHorz(w_hat)
    call updateHorz(uu_hat)
    call updateHorz(vv_hat)
    call updateHorz(ww_hat)
    call updateHorz(uv_hat)
    call updateHorz(uw_hat)
    call updateHorz(vw_hat)
    call updateHorz(S11_hat)
    call updateHorz(S22_hat)
    call updateHorz(S33_hat)
    call updateHorz(S12_hat)
    call updateHorz(S13_hat)
    call updateHorz(S23_hat)

          
              
!     Local Averaging**************************
    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx
                LM=0.
                MM=0.
                do mX=-1,1
                    do mY=-1,1
                        if(i+mX > Nx)then
                            ii=i+mX-Nx
                        elseif(i+mX < 1)then
                            ii=i+mX+Nx
                        else
                            ii=i+mX
                        endif

                        jj=j+mY+1

                        L11=(uu_hat(ii,jj,k))-(u_hat(ii,jj,k))**2.
                        L22=(vv_hat(ii,jj,k))-(v_hat(ii,jj,k))**2.
                        L33=(ww_hat(ii,jj,k))-(w_hat(ii,jj,k))**2.
                    !       L33=ww_hat(ii,jj,k)-w_hat2(ii,jj,k)
                        L12=(uv_hat(ii,jj,k))-(u_hat(ii,jj,k)*v_hat(ii,jj,k))
                        L13=(uw_hat(ii,jj,k))-(u_hat(ii,jj,k)*w_hat(ii,jj,k))
                        L23=(vw_hat(ii,jj,k))-(v_hat(ii,jj,k)*w_hat(ii,jj,k))

                        M11=2.*L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        S11_hat(ii,jj,k)
                        M22=2.*L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        S22_hat(ii,jj,k)
                        M33=2.*L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        S33_hat(ii,jj,k)
                        M12=2.*L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        S12_hat(ii,jj,k)
                        M13=2.*L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        S13_hat(ii,jj,k)
                        M23=2.*L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        S23_hat(ii,jj,k)

                        LM=LM + (L11*M11+L22*M22+L33*M33)+ &
                        &               2.*(L12*M12+L13*M13+L23*M23)
                                       
                        MM=MM + (M11*M11+M22*M22+M33*M33)+ &
                        &               2.*(M12*M12+M13*M13+M23*M23)
                                       
                    end do
                end do

                if(abs(MM) < (1e-10) .OR. (LM/MM) < 0 .OR. (LM/MM) > 10) then
                    Cs2(i,j,k) = 0.0
                else
                    Cs2(i,j,k) = LM/MM
                end if
                       
            end do
        end do
              
        betaa(:,:,k) = betaa1D(k)
         
    end do

    return
    end Subroutine OPTIM_WONG
          














