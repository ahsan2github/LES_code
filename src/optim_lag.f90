    Subroutine OPTIM_LAG(Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd, &
    u_,v_,w_,L,betaa,a1_old,b1_old,c1_old,d1_old,e1_old &
    ,a2_old,b2_old,c2_old,d2_old,e2_old,TL2)

    use globals
    use scalars
    implicit none

    interface
    include './interfaces/filter_2laa.f90'
    include './interfaces/updateHorz.f90'
    include './interfaces/lagrng_sd.f90'
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: i,j,k
    real*8,dimension(:,:,:):: &
    Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd,u_,v_,w_,betaa, &
    a1_old,b1_old,c1_old,d1_old,e1_old,a2_old,b2_old,c2_old, &
    d2_old,e2_old, TL2

    real*8,dimension(:) :: L

    real*8,dimension(size(Cs2,1),size(Cs2,2),size(Cs2,3)):: &
    S11_hat,S22_hat,S33_hat, &
    S12_hat,S13_hat,S23_hat,u_hat,v_hat,w_hat,uu_hat,uv_hat, &
    uw_hat,vv_hat,vw_hat,ww_hat,S11_hatd,S22_hatd, &
    S33_hatd,S12_hatd,S13_hatd,S23_hatd,u_hatd,v_hatd,w_hatd, &
    uu_hatd,uv_hatd,uw_hatd,vv_hatd,vw_hatd,ww_hatd, &
    SS11_hat,SS12_hat,SS13_hat,SS22_hat,SS23_hat,SS33_hat, &
    SS11_hatd,SS12_hatd,SS13_hatd,SS22_hatd,SS23_hatd,SS33_hatd, &
    uu,uv,uw,vv,vw,ww,SS11,SS22,SS33,SS12,SS13,SS23,TL4

    real*8,dimension(size(Cs2,1),size(Cs2,2),size(Cs2,3)):: &
    a1,b1,c1,d1,e1,a2,b2,c2,d2,e2,b

    real*8 :: M11,M12,M13,M22,M23,M33, &
    LM,MM,L13,L23,L12,L22,L33,L11,Q13,Q23,Q12,Q22,Q33,Q11, &
    ll,QN,NN,trace,Cs_ave,b_ave,uin,vin,win

    real*8,dimension(size(Cs2,3)):: &
    a1_avg,b1_avg,c1_avg,d1_avg,e1_avg,a2_avg,b2_avg,c2_avg, &
    d2_avg,e2_avg,Z

    real*8 :: rtnewt
        
    do k=2,Nzb+1
         
        S(:,:,k)=sqrt(2.d0*(S11(:,:,k)**2+S22(:,:,k)**2+ &
        S33(:,:,k)**2+2.d0*S12(:,:,k)**2+ &
        &               2.d0*S13(:,:,k)**2+ &
        &               2.d0*S23(:,:,k)**2))

        SS11(:,:,k)=S(:,:,k)*S11(:,:,k)
        SS33(:,:,k)=S(:,:,k)*S33(:,:,k)
        SS22(:,:,k)=S(:,:,k)*S22(:,:,k)
        SS12(:,:,k)=S(:,:,k)*S12(:,:,k)
        SS13(:,:,k)=S(:,:,k)*S13(:,:,k)
        SS23(:,:,k)=S(:,:,k)*S23(:,:,k)
                       
        uu(:,:,k)=u_(:,:,k)**2
        vv(:,:,k)=v_(:,:,k)**2
        ww(:,:,k)=w_(:,:,k)**2
        uv(:,:,k)=u_(:,:,k)*v_(:,:,k)
        vw(:,:,k)=v_(:,:,k)*w_(:,:,k)
        uw(:,:,k)=u_(:,:,k)*w_(:,:,k)

    end do
         
!C...Filtering to get the _hat variables (coarser resolution)

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
    Call Filter_2Laa(SS23_hat,SS23_hatd,SS23)

    do k=2,nzb+1

        S_hat(:,:,k)=sqrt(2.d0*(S11_hat(:,:,k)**2+ &
        S22_hat(:,:,k)**2.d0+S33_hat(:,:,k)**2+ &
        &               2.d0*S12_hat(:,:,k)**2+2.d0*S13_hat(:,:,k)**2+ &
        &               2.d0*S23_hat(:,:,k)**2))
        S_hatd(:,:,k)=sqrt(2.d0*(S11_hatd(:,:,k)**2+ &
        S22_hatd(:,:,k)**2+S33_hatd(:,:,k)**2+ &
        &               2.d0*S12_hatd(:,:,k)**2+2.d0*S13_hatd(:,:,k)**2+ &
        &               2.d0*S23_hatd(:,:,k)**2))

    end do

    do k=2,Nzb+1
        do j=1,nyb
            do i=1,nx
                L11=(uu_hat(i,j,k))-(u_hat(i,j,k))**2
                L22=(vv_hat(i,j,k))-(v_hat(i,j,k))**2
                L12=(uv_hat(i,j,k))-(u_hat(i,j,k)*v_hat(i,j,k))
                L13=(uw_hat(i,j,k))-(u_hat(i,j,k)*w_hat(i,j,k))
                L23=(vw_hat(i,j,k))-(v_hat(i,j,k)*w_hat(i,j,k))
                L33=(ww_hat(i,j,k))-(w_hat(i,j,k))**2
                Q11=(uu_hatd(i,j,k))-(u_hatd(i,j,k))**2
                Q22=(vv_hatd(i,j,k))-(v_hatd(i,j,k))**2
                Q12=(uv_hatd(i,j,k))-(u_hatd(i,j,k)*v_hatd(i,j,k))
                Q13=(uw_hatd(i,j,k))-(u_hatd(i,j,k)*w_hatd(i,j,k))
                Q23=(vw_hatd(i,j,k))-(v_hatd(i,j,k)*w_hatd(i,j,k))
                Q33=(ww_hatd(i,j,k))-(w_hatd(i,j,k))**2

                a1(i,j,k)=-2.d0*l(k)**2*tfr**2.d0*S_hat(i,j,k)* &
                (L11*S11_hat(i,j,k)+L22*S22_hat(i,j,k)+ &
                L33*S33_hat(i,j,k)+2.d0*(L12*S12_hat(i,j,k)+L13* &
                S13_hat(i,j,k)+L23*S23_hat(i,j,k)))
                               
                a2(i,j,k)=-2.d0*l(k)**2*tfr**4*S_hatd(i,j,k)* &
                (Q11*S11_hatd(i,j,k)+Q22*S22_hatd(i,j,k)+Q33* &
                S33_hatd(i,j,k)+2.d0*(Q12*S12_hatd(i,j,k) &
                +Q13*S13_hatd(i,j,k)+Q23*S23_hatd(i,j,k)))
                               
                b1(i,j,k)=-2.d0*l(k)**2*(L11*SS11_hat(i,j,k)+L22* &
                SS22_hat(i,j,k)+L33*SS33_hat(i,j,k)+ &
                &               2.d0*(L12*SS12_hat(i,j,k)+L13*SS13_hat(i,j,k) &
                +L23*SS23_hat(i,j,k)))
                               
                b2(i,j,k)=-2.d0*l(k)**2*(Q11*SS11_hatd(i,j,k)+ &
                Q22*SS22_hatd(i,j,k)+Q33*SS33_hatd(i,j,k)+2.d0* &
                (Q12*SS12_hatd(i,j,k)+Q13*SS13_hatd(i,j,k)+ &
                Q23*SS23_hatd(i,j,k)))
                               
                c1(i,j,k)=(2.d0*l(k)**2)**2* &
                (SS11_hat(i,j,k)**2+SS22_hat(i,j,k)**2 &
                +SS33_hat(i,j,k)**2+2.d0*(SS12_hat(i,j,k)**2+ &
                SS13_hat(i,j,k)**2+SS23_hat(i,j,k)**2))
                               
                c2(i,j,k)=(2.d0*l(k)**2)**2*(SS11_hatd(i,j,k)**2+ &
                SS22_hatd(i,j,k)**2+SS33_hatd(i,j,k)**2 &
                +2.d0*(SS12_hatd(i,j,k)**2+SS13_hatd(i,j,k)**2+ &
                SS23_hatd(i,j,k)**2.))
                               
                d1(i,j,k)=(4.d0*l(k)**4)*tfr**4*S_hat(i,j,k)**2* &
                (S11_hat(i,j,k)**2+S22_hat(i,j,k)**2 &
                +S33_hat(i,j,k)**2+2.d0*(S12_hat(i,j,k)**2+ &
                S13_hat(i,j,k)**2+S23_hat(i,j,k)**2))
                               
                d2(i,j,k)=(4.d0*l(k)**4)*tfr**8*S_hatd(i,j,k)**2 &
                *(S11_hatd(i,j,k)**2+S22_hatd(i,j,k)**2+ &
                S33_hatd(i,j,k)**2+2.d0*(S12_hatd(i,j,k)**2+ &
                S13_hatd(i,j,k)**2+S23_hatd(i,j,k)**2))
                               
                e1(i,j,k)=(8.d0*l(k)**4)*tfr**2* &
                S_hat(i,j,k)*(S11_hat(i,j,k)* &
                SS11_hat(i,j,k)+S22_hat(i,j,k)* &
                SS22_hat(i,j,k)+S33_hat(i,j,k)* &
                SS33_hat(i,j,k)+2.d0*(S12_hat(i,j,k)* &
                SS12_hat(i,j,k)+S13_hat(i,j,k)* &
                SS13_hat(i,j,k)+S23_hat(i,j,k)* &
                SS23_hat(i,j,k)))
                               
                e2(i,j,k)=(8.d0*l(k)**4)*tfr**4*S_hatd(i,j,k)* &
                (S11_hatd(i,j,k)*SS11_hatd(i,j,k)+ &
                S22_hatd(i,j,k)*SS22_hatd(i,j,k)+ &
                S33_hatd(i,j,k)*SS33_hatd(i,j,k)+2.d0 &
                *(S12_hatd(i,j,k)*SS12_hatd(i,j,k)+ &
                S13_hatd(i,j,k)*SS13_hatd(i,j,k)+ &
                S23_hatd(i,j,k)*SS23_hatd(i,j,k)))
                               
            enddo
        enddo
    enddo

!     **************initilization for lagrangian averaging***********
    if(t <= 50 .AND. INITU /= 1)then
         
        call plane_avg(a1,a1_avg)
        call plane_avg(a2,a2_avg)
        call plane_avg(b1,b1_avg)
        call plane_avg(b2,b2_avg)
        call plane_avg(c1,c1_avg)
        call plane_avg(c2,c2_avg)
        call plane_avg(d1,d1_avg)
        call plane_avg(d2,d2_avg)
        call plane_avg(e1,e1_avg)
        call plane_avg(e2,e2_avg)
                 
        do k=2,nzb+1
            do j=1,nyb
                do i=1,Nx
                                      
                    a1(:,j,k)=a1_avg(k)
                    b1(:,j,k)=b1_avg(k)
                    c1(:,j,k)=c1_avg(k)
                    d1(:,j,k)=d1_avg(k)
                    e1(:,j,k)=e1_avg(k)
                    a2(:,j,k)=a2_avg(k)
                    b2(:,j,k)=b2_avg(k)
                    c2(:,j,k)=c2_avg(k)
                    d2(:,j,k)=d2_avg(k)
                    e2(:,j,k)=e2_avg(k)
                                      
                    a1_old(:,j+1,k)=a1(:,j,k)
                    a2_old(:,j+1,k)=a2(:,j,k)
                    b1_old(:,j+1,k)=b1(:,j,k)
                    b2_old(:,j+1,k)=b2(:,j,k)
                    c1_old(:,j+1,k)=c1(:,j,k)
                    d1_old(:,j+1,k)=d1(:,j,k)
                    e1_old(:,j+1,k)=e1(:,j,k)
                    c2_old(:,j+1,k)=c2(:,j,k)
                    d2_old(:,j+1,k)=d2(:,j,k)
                    e2_old(:,j+1,k)=e2(:,j,k)
                                      
                    if(t <= 5)then
                        betaa(:,j,k)=1.
                    endif
                enddo
            enddo
        enddo
                 
    endif
!     **************end of the lagrangian initilization**************

    call updateHorz(a1_old)
    call updateHorz(b1_old)
    call updateHorz(c1_old)
    call updateHorz(d1_old)
    call updateHorz(e1_old)
    call updateHorz(a2_old)
    call updateHorz(b2_old)
    call updateHorz(c2_old)
    call updateHorz(d2_old)
    call updateHorz(e2_old)

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
                end if
                         
                call lagrng_sd(a1(i,j,k),b1(i,j,k),c1(i,j,k), &
                d1(i,j,k),e1(i,j,k),a1_old,b1_old,c1_old, &
                d1_old,e1_old,betaa(i,j,k),TL2(i,j,k),uin,vin,win, &
                i,j,k,0)
                               
                call lagrng_sd(a2(i,j,k),b2(i,j,k),c2(i,j,k), &
                d2(i,j,k),e2(i,j,k),a2_old,b2_old,c2_old, &
                d2_old,e2_old,betaa(i,j,k),TL4(i,j,k),uin,vin,win, &
                i,j,k,1)

                call root(rtnewt,a1(i,j,k),a2(i,j,k),b1(i,j,k), &
                b2(i,j,k),c1(i,j,k),c2(i,j,k),d1(i,j,k), &
                d2(i,j,k),e1(i,j,k),e2(i,j,k))
                              
                if(rtnewt <= 0 .OR. rtnewt > (100.0))then
                                                        
                    betaa(i,j,k)=0.0
                    cs2(i,j,k)=0.0
                                      
                    a1(i,j,k)=0.
                    b1(i,j,k)=0.
                    c1(i,j,k)=0.
                    d1(i,j,k)=0.
                    e1(i,j,k)=0.
                                      
                    a2(i,j,k)=0.
                    b2(i,j,k)=0.
                    c2(i,j,k)=0.
                    d2(i,j,k)=0.
                    e2(i,j,k)=0.
                                      
                else
                                      
                    if(t <= 1 .AND. INITU /= 1) then
                        betaa(i,j,k)=1.
                    else
                        betaa(i,j,k)=rtnewt
                    end if
                                      
                    LM=(a1(i,j,k)*betaa(i,j,k)-b1(i,j,k))
                    MM=(c1(i,j,k)+d1(i,j,k)*betaa(i,j,k)**2.- &
                    e1(i,j,k)*betaa(i,j,k))
                    QN=(a2(i,j,k)*betaa(i,j,k)**2.-b2(i,j,k))
                    NN=(c2(i,j,k)+d2(i,j,k)*betaa(i,j,k)**4.- &
                    e2(i,j,k)*betaa(i,j,k)**2.)

                    if(LM < 0 .OR. MM < (1e-10))then
                                          
                        cs2(i,j,k)=0.0
                                             
                        a1(i,j,k)=0.
                        b1(i,j,k)=0.
                    !                     c1(i,j,k)=0.
                    !                     d1(i,j,k)=0.
                    !                     e1(i,j,k)=0.
                                             
                        a2(i,j,k)=0.
                        b2(i,j,k)=0.
                    !                     c2(i,j,k)=0.
                    !                     d2(i,j,k)=0.
                    !                     e2(i,j,k)=0.
                                             
                    elseif(QN < (0.0) .OR. NN < (1e-10))then
                                             
                        cs2(i,j,k)=0.0
                                             
                        a1(i,j,k)=0.
                        b1(i,j,k)=0.
                    !                    c1(i,j,k)=0.
                    !                    d1(i,j,k)=0.
                    !                    e1(i,j,k)=0.
                                             
                        a2(i,j,k)=0.
                        b2(i,j,k)=0.
                    !                    c2(i,j,k)=0.
                    !                    d2(i,j,k)=0.
                    !                    e2(i,j,k)=0.
                                             
                    else
                                             
                        cs2(i,j,k)=LM/MM

                    endif
                endif
                         
            enddo
        enddo
    enddo

    do k=2,nzb+1
        do j=1,nyb
            a1_old(:,j+1,k)=a1(:,j,k)
            b1_old(:,j+1,k)=b1(:,j,k)
            c1_old(:,j+1,k)=c1(:,j,k)
            d1_old(:,j+1,k)=d1(:,j,k)
            e1_old(:,j+1,k)=e1(:,j,k)
            a2_old(:,j+1,k)=a2(:,j,k)
            b2_old(:,j+1,k)=b2(:,j,k)
            c2_old(:,j+1,k)=c2(:,j,k)
            d2_old(:,j+1,k)=d2(:,j,k)
            e2_old(:,j+1,k)=e2(:,j,k)
        end do
    enddo
          
    return
          
    end Subroutine OPTIM_LAG
          














