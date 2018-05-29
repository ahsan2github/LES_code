    Subroutine optim_scl_lag(Pr2,S11,S33,S22,S12,S13,S23,S_hat,S_hatd, &
    u_,v_,w_,L,t_,tx,ty,tz,betaa1,a2_old,b2_old,c2_old, &
    d2_old,e2_old,a4_old,b4_old,c4_old,d4_old,e4_old)

    use globals
    use scalars
    use sgs
    implicit none

    interface
    include './interfaces/filter_2laa.f90'
    include './interfaces/updateHorz.f90'
    include './interfaces/lagrng_sd_s.f90'
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: i,j,k,m
    real*8,dimension(:,:,:):: Pr2,tz,tx,ty,S11,S22,S33,S12,S13,S23, &
    S_hat,t_,S_hatd,w_,u_,v_,betaa1, &
    a2_old,b2_old,c2_old,d2_old,e2_old,a4_old,b4_old,c4_old, &
    d4_old,e4_old

    real*8,dimension(:) :: L

    real*8,dimension(size(Pr2,1),size(Pr2,2),size(Pr2,3)):: tx_hat, &
    ty_hat,tz_hat,S,S11_hat,S22_hat,S33_hat,S12_hat,S13_hat, &
    S23_hat,u_hat,v_hat,w_hat,t_hat,S11_hatd,S22_hatd,S33_hatd, &
    S12_hatd,S13_hatd,S23_hatd,u_hatd,v_hatd,w_hatd,s_tx,s_ty, &
    s_tz,ut,vt,wt,ut_hat,vt_hat,wt_hat,ut_hatd,vt_hatd,wt_hatd, &
    t_hatd,tx_hatd,ty_hatd,tz_hatd,s_tx_hat,s_ty_hat, &
    s_tz_hat,s_tx_hatd,s_ty_hatd,s_tz_hatd

    real*8,dimension(size(Pr2,1),size(Pr2,2),size(Pr2,3)):: a2,b2,c2, &
    d2,e2,a4,b4,c4,d4,e4

    real*8,dimension(size(Pr2,3))::a2_avg,b2_avg,c2_avg,d2_avg,e2_avg, &
    a4_avg,b4_avg,c4_avg,d4_avg,e4_avg,sig_t,avg_t

    real*8 :: up1,down1,aa,bb,cc,dd,ee,ff,root1,root2, &
    rtnewt,Tn,b_ave, &
    pr_ave,KX1,XX1,KX2,XX2,uin,vin,win,plane_ave
         
    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

                S(i,j,k)=sqrt(2.*(S11(i,j,k)**2.+S22(i,j,k)**2.+ &
                S33(i,j,k)**2.+2.*S12(i,j,k)**2.+ &
                &               2.*S13(i,j,k)**2.+2.*S23(i,j,k)**2.))

                S_tx(i,j,k)=S(i,j,k)*tx(i,j,k)
                S_ty(i,j,k)=S(i,j,k)*ty(i,j,k)
                S_tz(i,j,k)=S(i,j,k)*tz(i,j,k)
                ut(i,j,k)=u_(i,j,k)*t_(i,j,k)
                vt(i,j,k)=v_(i,j,k)*t_(i,j,k)
                wt(i,j,k)=w_(i,j,k)*t_(i,j,k)

            enddo
        enddo
    enddo

!C...Filtering to get the _hat variables (coarser resolution)

    Call Filter_2Laa(u_hat,u_hatd,u_)
    Call Filter_2Laa(v_hat,v_hatd,v_)
    Call Filter_2Laa(w_hat,w_hatd,w_)
    Call Filter_2Laa(t_hat,t_hatd,t_)
    Call Filter_2Laa(ut_hat,ut_hatd,ut)
    Call Filter_2Laa(vt_hat,vt_hatd,vt)
    Call Filter_2Laa(wt_hat,wt_hatd,wt)
    Call Filter_2Laa(tx_hat,tx_hatd,tx)
    Call Filter_2Laa(ty_hat,ty_hatd,ty)
    Call Filter_2Laa(tz_hat,tz_hatd,tz)
          
    if(scl_nodes /= mom_nodes)then

        Call Filter_2Laa(S11_hat,S11_hatd,S11)
        Call Filter_2Laa(S22_hat,S22_hatd,S22)
        Call Filter_2Laa(S33_hat,S33_hatd,S33)
        Call Filter_2Laa(S12_hat,S12_hatd,S12)
        Call Filter_2Laa(S13_hat,S13_hatd,S13)
        Call Filter_2Laa(S23_hat,S23_hatd,S23)

        do k=2,nzb+1
            do j=1,nyb
                do i=1,nx
                                      
                    S_hat(i,j,k)=sqrt(2.*(S11_hat(i,j,k)**2.+ &
                    S22_hat(i,j,k)**2.+ &
                    S33_hat(i,j,k)**2.+ &
                    &                  2.*S12_hat(i,j,k)**2.+ &
                    &                  2.*S13_hat(i,j,k)**2.+ &
                    &                  2.*S23_hat(i,j,k)**2.))

                    S_hatd(i,j,k)=sqrt(2.*(S11_hatd(i,j,k)**2. &
                    +S22_hatd(i,j,k)**2. &
                    +S33_hatd(i,j,k)**2. &
                    +2.*S12_hatd(i,j,k)**2. &
                    +2.*S13_hatd(i,j,k)**2. &
                    +2.*S23_hatd(i,j,k)**2.))

                end do
            end do
        end do

    endif

    Call Filter_2Laa(s_tx_hat,s_tx_hatd,s_tx)
    Call Filter_2Laa(s_ty_hat,s_ty_hatd,s_ty)
    Call Filter_2Laa(s_tz_hat,s_tz_hatd,s_tz)
    Call Filter_2Laa(s_tz_hat,s_tz_hatd,s_tz)

    do k=2,nzb+1
                 
        do j=1,nyb
            do i=1,nx

                a2(i,j,k)=+L(k)**2.* &
                ((ut_hat(i,j,k)-u_hat(i,j,k)*t_hat(i,j,k))* &
                s_tx_hat(i,j,k)+(vt_hat(i,j,k)-v_hat(i,j,k)* &
                t_hat(i,j,k))*s_ty_hat(i,j,k)+(wt_hat(i,j,k)- &
                w_hat(i,j,k)*t_hat(i,j,k))*s_tz_hat(i,j,k))

                b2(i,j,k)=-(tfr*L(k))**2.*((ut_hat(i,j,k)-u_hat(i,j,k)* &
                t_hat(i,j,k))*s_hat(i,j,k)*tx_hat(i,j,k) &
                +(vt_hat(i,j,k)-v_hat(i,j,k)*t_hat(i,j,k))* &
                s_hat(i,j,k)*ty_hat(i,j,k)+(wt_hat(i,j,k)- &
                w_hat(i,j,k)*t_hat(i,j,k))*s_hat(i,j,k)* &
                tz_hat(i,j,k))

                c2(i,j,k)=+L(k)**4.*((s_tx_hat(i,j,k))**2. &
                +(s_ty_hat(i,j,k))**2.+(s_tz_hat(i,j,k))**2.)

                d2(i,j,k)=-2.*L(k)**4.*(tfr**2.)* &
                (s_tx_hat(i,j,k)*s_hat(i,j,k)*tx_hat(i,j,k) &
                +s_ty_hat(i,j,k)*s_hat(i,j,k)*ty_hat(i,j,k) &
                +s_tz_hat(i,j,k)*s_hat(i,j,k)*tz_hat(i,j,k))

                e2(i,j,k)=+(tfr*L(k))**4.* &
                ((s_hat(i,j,k)*tx_hat(i,j,k))**2. &
                +(s_hat(i,j,k)*ty_hat(i,j,k))**2. &
                +(s_hat(i,j,k)*tz_hat(i,j,k))**2.)

                a4(i,j,k)=+L(k)**2.* &
                ((ut_hatd(i,j,k)-u_hatd(i,j,k)*t_hatd(i,j,k))* &
                s_tx_hatd(i,j,k)+(vt_hatd(i,j,k)-v_hatd(i,j,k)* &
                t_hatd(i,j,k))*s_ty_hatd(i,j,k)+(wt_hatd(i,j,k)- &
                w_hatd(i,j,k)*t_hatd(i,j,k))*s_tz_hatd(i,j,k))

                b4(i,j,k)=-(L(k))**2.*(tfr**4.)*((ut_hatd(i,j,k)- &
                u_hatd(i,j,k)*t_hatd(i,j,k))*s_hatd(i,j,k)* &
                tx_hatd(i,j,k)+(vt_hatd(i,j,k)-v_hatd(i,j,k)* &
                t_hatd(i,j,k))*s_hatd(i,j,k)*ty_hatd(i,j,k)+ &
                (wt_hatd(i,j,k)-w_hatd(i,j,k)*t_hatd(i,j,k))* &
                s_hatd(i,j,k)*tz_hatd(i,j,k))

                c4(i,j,k)=+L(k)**4.*((s_tx_hatd(i,j,k))**2. &
                +(s_ty_hatd(i,j,k))**2.+(s_tz_hatd(i,j,k))**2.)

                d4(i,j,k)=-2.*L(k)**4.*tfr**4.* &
                (s_tx_hatd(i,j,k)*s_hatd(i,j,k)*tx_hatd(i,j,k) &
                +s_ty_hatd(i,j,k)*s_hatd(i,j,k)*ty_hatd(i,j,k) &
                +s_tz_hatd(i,j,k)*s_hatd(i,j,k)*tz_hatd(i,j,k))

                e4(i,j,k)=+(L(k))**4.*(tfr**8.)* &
                ((s_hatd(i,j,k)*tx_hatd(i,j,k))**2. &
                +(s_hatd(i,j,k)*ty_hatd(i,j,k))**2. &
                +(s_hatd(i,j,k)*tz_hatd(i,j,k))**2.)

            end do
        end do
    end do

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     **************initilization for lagrangian averaging***********
    if(t <= 50 .AND. INITS == 0)then
                   
        call plane_avg(a2,a2_avg)
        call plane_avg(b2,b2_avg)
        call plane_avg(c2,c2_avg)
        call plane_avg(d2,d2_avg)
        call plane_avg(e2,e2_avg)
        call plane_avg(a4,a4_avg)
        call plane_avg(b4,b4_avg)
        call plane_avg(c4,c4_avg)
        call plane_avg(d4,d4_avg)
        call plane_avg(e4,e4_avg)
                 
        do k=2,nzb+1
            do j=1,Nyb
                do i=1,Nx
                                      
                    a2(i,j,k)=a2_avg(k)
                    b2(i,j,k)=b2_avg(k)
                    c2(i,j,k)=c2_avg(k)
                    d2(i,j,k)=d2_avg(k)
                    e2(i,j,k)=e2_avg(k)
                                      
                    a4(i,j,k)=a4_avg(k)
                    b4(i,j,k)=b4_avg(k)
                    c4(i,j,k)=c4_avg(k)
                    d4(i,j,k)=d4_avg(k)
                    e4(i,j,k)=e4_avg(k)
                                      
                    a2_old(i,j+1,k)=a2_avg(k)
                    b2_old(i,j+1,k)=b2_avg(k)
                    c2_old(i,j+1,k)=c2_avg(k)
                    d2_old(i,j+1,k)=d2_avg(k)
                    e2_old(i,j+1,k)=e2_avg(k)
                                      
                    a4_old(i,j+1,k)=a4_avg(k)
                    b4_old(i,j+1,k)=b4_avg(k)
                    c4_old(i,j+1,k)=c4_avg(k)
                    d4_old(i,j+1,k)=d4_avg(k)
                    e4_old(i,j+1,k)=e4_avg(k)
                ! ccccccccccccccccccccccccccccccccccc
                    if(t <= 5)then
                        betaa1(i,j,k)=1.
                    end if
                ! ccccccccccccccccccccccccccccccccccc
                                      
                end do
            end do
        end do
    endif
!     **************end of the lagrangian initilization**************

    call updateHorz(a2_old)
    call updateHorz(b2_old)
    call updateHorz(c2_old)
    call updateHorz(d2_old)
    call updateHorz(e2_old)
    call updateHorz(a4_old)
    call updateHorz(b4_old)
    call updateHorz(c4_old)
    call updateHorz(d4_old)
    call updateHorz(e4_old)
      
    do k=2,nzb+1
        avg_t(k)=0.0
        do j=1,nyb
            do i=1,nx
                avg_t(k)=avg_t(k)+t_(i,j,k)
            enddo
        enddo
        avg_t(k)=avg_t(k)/float(nx*nyb)
    enddo

    do k=2,nzb+1
        sig_t(k)=0.0
        do j=1,nyb
            do i=1,nx
                sig_t(k)=sig_t(k)+(t_(i,j,k)-avg_t(k))**2.
            enddo
        enddo
        sig_t(k)=(sig_t(k)/float(nx*nyb-1))**0.5
    enddo

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

                call lagrng_sd_s(a2(i,j,k),b2(i,j,k),c2(i,j,k), &
                d2(i,j,k),e2(i,j,k),a2_old,b2_old,c2_old, &
                d2_old,e2_old,betaa1(i,j,k),uin,vin,win, &
                i,j,k,0,sig_t(k))
                               
                call lagrng_sd_s(a4(i,j,k),b4(i,j,k),c4(i,j,k), &
                d4(i,j,k),e4(i,j,k),a4_old,b4_old,c4_old, &
                d4_old,e4_old,betaa1(i,j,k),uin,vin,win, &
                i,j,k,1,sig_t(k))

                aa= a2(i,j,k)*c4(i,j,k)-a4(i,j,k)*c2(i,j,k)
                bb=-a4(i,j,k)*d2(i,j,k)+b2(i,j,k)*c4(i,j,k)
                cc=-c2(i,j,k)*b4(i,j,k)+a2(i,j,k)*d4(i,j,k)- &
                a4(i,j,k)*e2(i,j,k)
                dd= b2(i,j,k)*d4(i,j,k)-b4(i,j,k)*d2(i,j,k)
                ee= a2(i,j,k)*e4(i,j,k)-b4(i,j,k)*e2(i,j,k)
                ff= b2(i,j,k)*e4(i,j,k)
                               
                call root8(rtnewt,aa,bb,cc,dd,ee,ff)
                 
                if(rtnewt <= 0 .OR. rtnewt > 100.)then

                    betaa1(i,j,k)=0.0
                    Pr2(i,j,k)=0.0

                    a2(i,j,k)=0.0
                    b2(i,j,k)=0.0
                    c2(i,j,k)=0.0
                    d2(i,j,k)=0.0
                    e2(i,j,k)=0.0

                    a4(i,j,k)=0.0
                    b4(i,j,k)=0.0
                    c4(i,j,k)=0.0
                    d4(i,j,k)=0.0
                    e4(i,j,k)=0.0

                else

                    if(t <= 1 .AND. INITS /= 1) then
                        betaa1(i,j,k)=1.
                    else
                        betaa1(i,j,k)=rtnewt
                    end if

                    KX1=a2(i,j,k)+betaa1(i,j,k)*b2(i,j,k)
                    XX1=c2(i,j,k)+betaa1(i,j,k)*d2(i,j,k)+ &
                    e2(i,j,k)*betaa1(i,j,k)**2

                    KX2=a4(i,j,k)+betaa1(i,j,k)**2.*b4(i,j,k)
                    XX2=c4(i,j,k)+betaa1(i,j,k)**2.*d4(i,j,k)+ &
                    e4(i,j,k)*betaa1(i,j,k)**4.

                    if(KX1 < 0 .OR. XX1 < (1e-10))then

                        Pr2(i,j,k)=0.0
                                             
                        a2(i,j,k)=0.0
                        b2(i,j,k)=0.0
                    !                     c2(i,j,k)=0.0
                    !                     d2(i,j,k)=0.0
                    !                     e2(i,j,k)=0.0
                                             
                        a4(i,j,k)=0.0
                        b4(i,j,k)=0.0
                    !                     c4(i,j,k)=0.0
                    !                     d4(i,j,k)=0.0
                    !                     e4(i,j,k)=0.0

                    elseif(KX2 < 0 .OR. XX2 < (1e-10))then
                                             
                        Pr2(i,j,k)=0.0
                                             
                        a2(i,j,k)=0.0
                        b2(i,j,k)=0.0
                    !                     c2(i,j,k)=0.0
                    !                     d2(i,j,k)=0.0
                    !                     e2(i,j,k)=0.0
                                             
                        a4(i,j,k)=0.0
                        b4(i,j,k)=0.0
                    !                     c4(i,j,k)=0.0
                    !                     d4(i,j,k)=0.0
                    !                     e4(i,j,k)=0.0
                                             
                    else

                    !     Note: Pr is actually Cs2/Prt

                        Pr2(i,j,k)=KX1/XX1

                    endif

                end if
                               
            enddo
        enddo
    enddo

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

                a4_old(i,j+1,k)=a4(i,j,k)
                b4_old(i,j+1,k)=b4(i,j,k)
                c4_old(i,j+1,k)=c4(i,j,k)
                d4_old(i,j+1,k)=d4(i,j,k)
                e4_old(i,j+1,k)=e4(i,j,k)
                               
                a2_old(i,j+1,k)=a2(i,j,k)
                b2_old(i,j+1,k)=b2(i,j,k)
                c2_old(i,j+1,k)=c2(i,j,k)
                d2_old(i,j+1,k)=d2(i,j,k)
                e2_old(i,j+1,k)=e2(i,j,k)

            enddo
        enddo
    enddo
           
!     **************end of betaa calculation*************************

    return
    end Subroutine optim_scl_lag











