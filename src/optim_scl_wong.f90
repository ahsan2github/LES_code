    Subroutine optim_scl_wong(Pr2,S11,S33,S22,S12,S13,S23,S_hat, &
    S_hatd,u_,v_,w_,L,t_,tx,ty,tz,betaa1)

    use globals
    use scalars
    use sgs
    implicit none

    interface
    include './interfaces/filter_2laa.f90'
    include './interfaces/updateHorz.f90'
    include './interfaces/plane_reduce.f90'
    end interface

    integer*4 :: i,j,k,ii,jj,mX,mY
    real*8,dimension(:,:,:):: S11,S22,S33,S12,S13,S23,S_hat,S_hatd, &
    u_,v_,w_,t_,tx,ty,tz,betaa1,Pr2

    real*8,dimension(size(Pr2,1),size(Pr2,2),size(Pr2,3)):: &
    S,S11_hat,S22_hat,S33_hat,S12_hat,S13_hat,S23_hat,S11_hatd, &
    S22_hatd,S33_hatd,S12_hatd,S13_hatd,S23_hatd,u_hatd,v_hatd, &
    w_hatd,s_tx,s_ty,s_tz,ut,vt,wt,ut_hatd,vt_hatd,wt_hatd, &
    t_hatd,tx_hatd,ty_hatd,tz_hatd,s_tx_hat,s_ty_hat,s_tz_hat, &
    s_tx_hatd,s_ty_hatd,s_tz_hatd

    real*8,dimension(size(Pr2,1),size(Pr2,2)+2,size(Pr2,3)):: &
    tx_hat,ty_hat,tz_hat,u_hat,v_hat,w_hat,t_hat,ut_hat,vt_hat, &
    wt_hat

    real*8 :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,betaa1D(size(Pr2,3))

    real*8 :: aa,bb,cc,dd,ee,ff,L1,L2,L3,M1,M2,M3,Q1,Q2,Q3, &
    LM,MM,rtnewt,L(:)
         
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

    do k=2,nzb+1

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

                L1=(ut_hat(i,j,k))-(u_hat(i,j,k))*(t_hat(i,j,k))
                L2=(vt_hat(i,j,k))-(v_hat(i,j,k))*(t_hat(i,j,k))
                L3=(wt_hat(i,j,k))-(w_hat(i,j,k))*(t_hat(i,j,k))
                               
                Q1=(ut_hatd(i,j,k))-(u_hatd(i,j,k))*(t_hatd(i,j,k))
                Q2=(vt_hatd(i,j,k))-(v_hatd(i,j,k))*(t_hatd(i,j,k))
                Q3=(wt_hatd(i,j,k))-(w_hatd(i,j,k))*(t_hatd(i,j,k))
                       
                a1=a1+(Q1*tx_hatd(i,j,k)+ &
                Q2*ty_hatd(i,j,k)+ &
                Q3*tz_hatd(i,j,k))
                 
                a2=a1*(-tfr**(8./3.))

                a3=a3+(tx_hat(i,j,k)**2.+ &
                ty_hat(i,j,k)**2.+ &
                tz_hat(i,j,k)**2.)

                a4=(-2.*tfr**(4./3.))*a3
                     
                a5=(tfr**(8./3.))*a3

                a6=a6+(L1*tx_hat(i,j,k)+ &
                L2*ty_hat(i,j,k)+ &
                L3*tz_hat(i,j,k))

                a7=a6*(-tfr**(4./3.))
                               
                a8=a8+(tx_hatd(i,j,k)**2.+ &
                ty_hatd(i,j,k)**2.+ &
                tz_hatd(i,j,k)**2.)

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

!     Create and Update Ghost Layers**************

    do j=1,nyb
        jj=abs(j-nyb-1)
        u_hat(:,jj+1,:) = u_hat(:,jj,:)
        v_hat(:,jj+1,:) = v_hat(:,jj,:)
        w_hat(:,jj+1,:) = w_hat(:,jj,:)
        t_hat(:,jj+1,:) = t_hat(:,jj,:)
        ut_hat(:,jj+1,:) = ut_hat(:,jj,:)
        vt_hat(:,jj+1,:) = vt_hat(:,jj,:)
        wt_hat(:,jj+1,:) = wt_hat(:,jj,:)
        tx_hat(:,jj+1,:) = tx_hat(:,jj,:)
        ty_hat(:,jj+1,:) = ty_hat(:,jj,:)
        tz_hat(:,jj+1,:) = tz_hat(:,jj,:)
    end do

    call updateHorz(u_hat)
    call updateHorz(v_hat)
    call updateHorz(w_hat)
    call updateHorz(t_hat)
    call updateHorz(ut_hat)
    call updateHorz(vt_hat)
    call updateHorz(wt_hat)
    call updateHorz(tx_hat)
    call updateHorz(ty_hat)
    call updateHorz(tz_hat)

!     Local Averaging****************************

    do k=2,nzb+1
        do j=1,nyb
            do i=1,nx

                LM = 0.
                MM = 0.
                               
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
                                          
                        L1=(ut_hat(ii,jj,k))-(u_hat(ii,jj,k))*(t_hat(ii,jj,k))
                        L2=(vt_hat(ii,jj,k))-(v_hat(ii,jj,k))*(t_hat(ii,jj,k))
                        L3=(wt_hat(ii,jj,k))-(w_hat(ii,jj,k))*(t_hat(ii,jj,k))
                                          
                        M1=L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        tx_hat(ii,jj,k)

                        M2=L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        ty_hat(ii,jj,k)

                        M3=L(k)**(4./3.)*(1.-tfr**(4./3.)*betaa1D(k))* &
                        tz_hat(ii,jj,k)
                                       
                        LM=LM + L1*M1+L2*M2+L3*M3
                        MM=MM + M1*M1+M2*M2+M3*M3
                                       
                    end do
                end do

                if(abs(MM) < (1e-10) .OR. (LM/MM) < 0 .OR. (LM/MM) > 10) then
                    Pr2(i,j,k) = 0.0
                else
                    Pr2(i,j,k) = LM/MM
                end if
                            
            end do
        end do

        betaa1(:,:,k) = betaa1D(k)
              
    end do
          
    return
    end Subroutine optim_scl_wong











