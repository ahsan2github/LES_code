    Subroutine OPTIM_PLN (Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd, &
    u_,v_,w_,L,betaa)

    use globals
    use scalars
    use sgs
    implicit none

    interface
    include './interfaces/filter_2laa.f90'
    include './interfaces/plane_reduce.f90'
    end interface

    integer*4 :: i,j,k,Mx,My,ii,jj
    real*8,dimension(:,:,:):: &
    S,S11,S22,S33,S12,S13,S23,S_hat,S_hatd,w_,u_,v_,Cs2,betaa

    real*8,dimension(size(Cs2,1),size(Cs2,2),size(Cs2,3)):: &
    S11_hat,S22_hat,S33_hat,S12_hat,S13_hat,S23_hat,u_hat,v_hat, &
    w_hat,uu_hat,uv_hat,uw_hat,vv_hat,vw_hat,ww_hat, &
    S11_hatd,S22_hatd,S33_hatd,S12_hatd,S13_hatd,S23_hatd,u_hatd, &
    v_hatd,w_hatd,uu_hatd,uv_hatd,uw_hatd,vv_hatd,vw_hatd, &
    ww_hatd,SS11_hat,SS12_hat,SS13_hat,SS22_hat,SS23_hat, &
    SS33_hat,SS11_hatd,SS12_hatd,SS13_hatd,SS22_hatd,SS23_hatd, &
    SS33_hatd,uu,uv,uw,vv,vw,ww,SS11,SS22,SS33,SS12,SS13,SS23
    real*8,dimension(size(Cs2,3)) :: betaa1D,Cs21D

    real*8 :: a1,b1,c1,d1,e1,a2,b2,c2,d2,e2,LM,MM,aa,bb,cc,dd,ee,ff, &
    b,rtnewt

    real*8 :: L(:),M11,M12,M13,M22,M23,M33, &
    L13,L23,L12,L22,L33,L11,Q13,Q23,Q12,Q22,Q33,Q11

    do k=2,Nzb+1
                       
        S(:,:,k)=sqrt(2.*(S11(:,:,k)**2.+S22(:,:,k)**2.+ &
        S33(:,:,k)**2.+2.*S12(:,:,k)**2.+ &
        &               2.*S13(:,:,k)**2.+ &
        &               2.*S23(:,:,k)**2.))
                       
        SS11(:,:,k)=S(:,:,k)*S11(:,:,k)
        SS33(:,:,k)=S(:,:,k)*S33(:,:,k)
        SS22(:,:,k)=S(:,:,k)*S22(:,:,k)
        SS12(:,:,k)=S(:,:,k)*S12(:,:,k)
        SS13(:,:,k)=S(:,:,k)*S13(:,:,k)
        SS23(:,:,k)=S(:,:,k)*S23(:,:,k)
                       
        uu(:,:,k)=u_(:,:,k)**2.
        vv(:,:,k)=v_(:,:,k)**2.
        ww(:,:,k)=w_(:,:,k)**2.
        uv(:,:,k)=u_(:,:,k)*v_(:,:,k)
        vw(:,:,k)=v_(:,:,k)*w_(:,:,k)
        uw(:,:,k)=u_(:,:,k)*w_(:,:,k)

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

        S_hat(:,:,k)=sqrt(2.*(S11_hat(:,:,k)**2.+ &
        S22_hat(:,:,k)**2.+S33_hat(:,:,k)**2.+ &
        &               2.*S12_hat(:,:,k)**2.+2.*S13_hat(:,:,k)**2.+ &
        &               2.*S23_hat(:,:,k)**2.))
        S_hatd(:,:,k)=sqrt(2.*(S11_hatd(:,:,k)**2.+ &
        S22_hatd(:,:,k)**2.+S33_hatd(:,:,k)**2.+ &
        &               2.*S12_hatd(:,:,k)**2.+2.*S13_hatd(:,:,k)**2.+ &
        &               2.*S23_hatd(:,:,k)**2.))

    end do
          
    do k=2,Nzb+1
                 
        a1=0.
        a2=0.
        b1=0.
        b2=0.
        c1=0.
        c2=0.
        d1=0.
        d2=0.
        e1=0.
        e2=0.
                 
        do j=1,nyb
            do i=1,nx
                L11=(uu_hat(i,j,k))-(u_hat(i,j,k))**2.
                L22=(vv_hat(i,j,k))-(v_hat(i,j,k))**2.
                L33=(ww_hat(i,j,k))-(w_hat(i,j,k))**2.
            !     L33=ww_hat(i,j,k)-w_hat2(i,j,k)
                L12=(uv_hat(i,j,k))-(u_hat(i,j,k)*v_hat(i,j,k))
                L13=(uw_hat(i,j,k))-(u_hat(i,j,k)*w_hat(i,j,k))
                L23=(vw_hat(i,j,k))-(v_hat(i,j,k)*w_hat(i,j,k))
                Q11=(uu_hatd(i,j,k))-(u_hatd(i,j,k))**2.
                Q22=(vv_hatd(i,j,k))-(v_hatd(i,j,k))**2.
                Q33=(ww_hatd(i,j,k))-(w_hatd(i,j,k))**2.
            !     Q33=ww_hatd(i,j,k)-w_hatd2(i,j,k)
                Q12=(uv_hatd(i,j,k))-(u_hatd(i,j,k)*v_hatd(i,j,k))
                Q13=(uw_hatd(i,j,k))-(u_hatd(i,j,k)*w_hatd(i,j,k))
                Q23=(vw_hatd(i,j,k))-(v_hatd(i,j,k)*w_hatd(i,j,k))
                               
                               
                a1=a1+2.*l(k)**2.* &
                (L11*SS11_hat(i,j,k)+ &
                L22*SS22_hat(i,j,k)+ &
                L33*SS33_hat(i,j,k)+ &
                &               2.*(L12*SS12_hat(i,j,k)+ &
                L13*SS13_hat(i,j,k)+ &
                L23*SS23_hat(i,j,k)))
                               
                a2=a2+2.*l(k)**2.* &
                (Q11*SS11_hatd(i,j,k)+ &
                Q22*SS22_hatd(i,j,k)+ &
                Q33*SS33_hatd(i,j,k)+ &
                &               2.*(Q12*SS12_hatd(i,j,k)+ &
                Q13*SS13_hatd(i,j,k)+ &
                Q23*SS23_hatd(i,j,k)))
                               
                b1=b1+2.*l(k)**2.*tfr**2.*S_hat(i,j,k)* &
                (L11*S11_hat(i,j,k)+ &
                L22*S22_hat(i,j,k)+ &
                L33*S33_hat(i,j,k)+ &
                &               2.*(L12*S12_hat(i,j,k)+ &
                L13*S13_hat(i,j,k)+ &
                L23*S23_hat(i,j,k)))
                               
                b2=b2+2.*l(k)**2.*tfr**4.*S_hatd(i,j,k)* &
                (Q11*S11_hatd(i,j,k)+ &
                Q22*S22_hatd(i,j,k)+ &
                Q33*S33_hatd(i,j,k)+ &
                &               2.*(Q12*S12_hatd(i,j,k)+ &
                Q13*S13_hatd(i,j,k)+ &
                Q23*S23_hatd(i,j,k)))
                               
                c1=c1+(2.*l(k)**2.)**2.* &
                (SS11_hat(i,j,k)**2.+ &
                SS22_hat(i,j,k)**2.+ &
                SS33_hat(i,j,k)**2.+ &
                &               2.*(SS12_hat(i,j,k)**2.+ &
                SS13_hat(i,j,k)**2.+ &
                SS23_hat(i,j,k)**2.))
                               
                c2=c2+(2.*l(k)**2.)**2.* &
                (SS11_hatd(i,j,k)**2.+ &
                SS22_hatd(i,j,k)**2.+ &
                SS33_hatd(i,j,k)**2.+ &
                &               2.*(SS12_hatd(i,j,k)**2.+ &
                SS13_hatd(i,j,k)**2.+ &
                SS23_hatd(i,j,k)**2.))
                               
                d1=d1+(4.*l(k)**4.)*tfr**4.*S_hat(i,j,k)**2.* &
                (S11_hat(i,j,k)**2.+ &
                S22_hat(i,j,k)**2.+ &
                S33_hat(i,j,k)**2.+ &
                &               2.*(S12_hat(i,j,k)**2.+ &
                S13_hat(i,j,k)**2.+ &
                S23_hat(i,j,k)**2.))
                               
                d2=d2+(4.*l(k)**4.)*tfr**8.*S_hatd(i,j,k)**2.* &
                (S11_hatd(i,j,k)**2.+ &
                S22_hatd(i,j,k)**2.+ &
                S33_hatd(i,j,k)**2.+ &
                &               2.*(S12_hatd(i,j,k)**2.+ &
                S13_hatd(i,j,k)**2.+ &
                S23_hatd(i,j,k)**2.))
                               
                e1=e1+(8.*l(k)**4.)*tfr**2.*S_hat(i,j,k)* &
                (S11_hat(i,j,k)*SS11_hat(i,j,k)+ &
                S22_hat(i,j,k)*SS22_hat(i,j,k)+ &
                S33_hat(i,j,k)*SS33_hat(i,j,k)+ &
                &               2.*(S12_hat(i,j,k)*SS12_hat(i,j,k)+ &
                S13_hat(i,j,k)*SS13_hat(i,j,k)+ &
                S23_hat(i,j,k)*SS23_hat(i,j,k)))
                               
                e2=e2+(8.*l(k)**4.)*tfr**4.*S_hatd(i,j,k)* &
                (S11_hatd(i,j,k)*SS11_hatd(i,j,k)+ &
                S22_hatd(i,j,k)*SS22_hatd(i,j,k)+ &
                S33_hatd(i,j,k)*SS33_hatd(i,j,k)+ &
                &               2.*(S12_hatd(i,j,k)*SS12_hatd(i,j,k)+ &
                S13_hatd(i,j,k)*SS13_hatd(i,j,k)+ &
                S23_hatd(i,j,k)*SS23_hatd(i,j,k)))
            end do
        end do

        call plane_reduce(a1)
        call plane_reduce(a2)
        call plane_reduce(b1)
        call plane_reduce(b2)
        call plane_reduce(c1)
        call plane_reduce(c2)
        call plane_reduce(d1)
        call plane_reduce(d2)
        call plane_reduce(e1)
        call plane_reduce(e2)

        aa = a1*c2-a2*c1
        bb = -b1*c2+a2*e1
        cc = -a1*e2+b2*c1-a2*d1
        dd = b1*e2-b2*e1
        ee = a1*d2+b2*d1
        ff = -b1*d2

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
          
!     Plane Averaging**************************
    do k=2,nzb+1
        LM = 0.
        MM = 0.

        do j=1,nyb
            do i=1,nx
                                     
                L11=(uu_hat(i,j,k))-(u_hat(i,j,k))**2.
                L22=(vv_hat(i,j,k))-(v_hat(i,j,k))**2.
                L33=(ww_hat(i,j,k))-(w_hat(i,j,k))**2.
            !     L33=ww_hat(i,j,k)-w_hat2(i,j,k)
                L12=(uv_hat(i,j,k))-(u_hat(i,j,k)*v_hat(i,j,k))
                L13=(uw_hat(i,j,k))-(u_hat(i,j,k)*w_hat(i,j,k))
                L23=(vw_hat(i,j,k))-(v_hat(i,j,k)*w_hat(i,j,k))
                       
                M11=2.*L(k)**2.*SS11_hat(i,j,k)- &
                &       2.*(tfr*L(k))**2.*betaa1D(k)*S_hat(i,j,k)*S11_hat(i,j,k)
                M22=2.*L(k)**2.*SS22_hat(i,j,k)- &
                &      2.*(tfr*L(k))**2.*betaa1D(k)*S_hat(i,j,k)*S22_hat(i,j,k)
                M33=2.*L(k)**2.*SS33_hat(i,j,k)- &
                &      2.*(tfr*L(k))**2.*betaa1D(k)*S_hat(i,j,k)*S33_hat(i,j,k)
                M12=2.*L(k)**2.*SS12_hat(i,j,k)- &
                &      2.*(tfr*L(k))**2.*betaa1D(k)*S_hat(i,j,k)*S12_hat(i,j,k)
                M13=2.*L(k)**2.*SS13_hat(i,j,k)- &
                &      2.*(tfr*L(k))**2.*betaa1D(k)*S_hat(i,j,k)*S13_hat(i,j,k)
                M23=2.*L(k)**2.*SS23_hat(i,j,k)- &
                &      2.*(tfr*L(k))**2.*betaa1D(k)*S_hat(i,j,k)*S23_hat(i,j,k)
                                              
                LM=LM + (L11*M11+L22*M22+L33*M33)+ &
                &               2.*(L12*M12+L13*M13+L23*M23)
                               
                MM=MM + (M11*M11+M22*M22+M33*M33)+ &
                &               2.*(M12*M12+M13*M13+M23*M23)


            end do
        end do

        call plane_reduce(LM)
        call plane_reduce(MM)

        if(abs(MM) < (1e-10) .OR. (LM/MM) < 0 .OR. (LM/MM) > 1.0) then
            Cs21D(k) = 0.0
        else
            Cs21D(k) = LM/MM
        end if

        betaa(:,:,k) = betaa1D(k)
        Cs2(:,:,k) = Cs21D(k)

    end do
          
    return
    end Subroutine OPTIM_PLN
          














