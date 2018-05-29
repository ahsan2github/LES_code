    subroutine scalar_slice(u,v,w,theta,sgs_t1,sgs_t2,sgs_t3,dtdx, &
    dtdy,dtdz,Pr2,Cs2,beta2,at,t2,t3,asgs_t1,asgs_t2,asgs_t3, &
    aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aqz_s, &
    ET3D,aET,ilow,wgx)

    use globals
    use scalars
    implicit none

    interface
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: i,j,k,l,kk,h
    real*8,dimension(:,:,:):: u,v,w,theta,sgs_t1,sgs_t2,sgs_t3, &
    dtdx,dtdy,dtdz,Pr2,Cs2,beta2,ET3D
    real*8,dimension(:,:) :: at,t2,t3,asgs_t1,asgs_t2,asgs_t3, &
    aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aET
    real*8,dimension(:,:) :: aqz_s
    integer*4,dimension(:) :: ilow
          
    real*8,dimension(size(at,1),size(at,2)) :: tu1,tt1,tt2,tt3,tsgst1, &
    tsgst2,tsgst3,tut,tvt,twt,tdtdx,tdtdy,tCs2,tdtdz,tCs2Pr, &
    tbeta2,tET,tw1
    real*8,dimension(size(at,2)) :: u_bar,v_bar,w_bar,t_bar
    real*8 :: arg1,arg2,fr,norm,ftn,wgx

    fr=(1./p_count)*c_count
    norm=1.d0/(Ny*(1+Nx-aNx))
    ftn=fr*norm

! c compute the plane averages of terms involved in products ccc

    call plane_avg(theta,t_bar)
    call plane_avg(u,u_bar)
    call plane_avg(v,v_bar)
    call plane_avg(w,w_bar)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    tu1=0.
    tt1=0.
    tt2=0.
    tt3=0.
    tsgst1=0.
    tsgst2=0.
    tsgst3=0.
    tut=0.
    tvt=0.
    twt=0.
    tdtdx=0.
    tdtdy=0.
    tdtdz=0.
    tCs2Pr=0.
    tCs2=0.
    tbeta2=0.
    tET=0.
          
    do k=2,Nzb+1
        do i=1,aNx

            do j=1,Nyb
                do l=1,1+Nx-aNx

                    tt1(i,k)=tt1(i,k)+theta(i-1+l,j,k)
                    tu1(i,k)=tu1(i,k)+u(i-1+l,j,k)+Ugal
                    tw1(i,k)=tw1(i,k)+w(i-1+l,j,k)
                    tsgst1(i,k)=tsgst1(i,k)+sgs_t1(i-1+l,j,k)
                    tsgst2(i,k)=tsgst2(i,k)+sgs_t2(i-1+l,j,k)
                    tsgst3(i,k)=tsgst3(i,k)+sgs_t3(i-1+l,j,k)
                    tdtdx(i,k)=tdtdx(i,k)+dtdx(i-1+l,j,k)
                    tdtdy(i,k)=tdtdy(i,k)+dtdy(i-1+l,j,k)
                    tdtdz(i,k)=tdtdz(i,k)+dtdz(i-1+l,j,k)
                    tCs2Pr(i,k)=tCs2Pr(i,k)+Pr2(i-1+l,j,k)
                    tCs2(i,k)=tCs2(i,k)+Cs2(i-1+l,j,k)
                    tbeta2(i,k)=tbeta2(i,k)+beta2(i-1+l,j,k)
                    tET(i,k)=tET(i,k)+ET3D(i-1+l,j,k)
                                      
                    tt2(i,k)=tt2(i,k)+(theta(i-1+l,j,k) - t_bar(k))**2
                    tt3(i,k)=tt3(i,k)+(theta(i-1+l,j,k) - t_bar(k))**3
                    tut(i,k)=tut(i,k)+(theta(i-1+l,j,k) - t_bar(k))* &
                    (u(i-1+l,j,k) - u_bar(k))
                    tvt(i,k)=tvt(i,k)+(theta(i-1+l,j,k) - t_bar(k))* &
                    (v(i-1+l,j,k) - v_bar(k))

                    if (k==2 .AND. vfact==0) then
                        arg1=0.
                    else
                        arg1=(theta(i-1+l,j,k)+theta(i-1+l,j,k-1) - &
                        t_bar(k) - t_bar(k-1))*0.5d0
                    end if

                    twt(i,k)=twt(i,k)+(w(i-1+l,j,k) - w_bar(k))*arg1

                enddo
            enddo
        enddo
    enddo



    do kk=2,Nzb+1
    ! c NOTE rightnow just matching with the location of jump (no interpolation)
        do i=1,aNx
            if(aNx /= 1)then
                l=ilow(i)
                h=ilow(i)
                if(l == nx)then
                    h=nx
                endif
            else
                l=1
                h=1
            endif

            if(nprocs == 1)then
                k=kk-1
            else
                k=kk
            endif
            	    
            at(i,k)=at(i,k)+ftn*((1-wgx)*tt1(l,k)+wgx*tt1(h,k))
            t2(i,k)=t2(i,k)+ftn*((1-wgx)*tt2(l,k)+wgx*tt2(h,k))
            t3(i,k)=t3(i,k)+ftn*((1-wgx)*tt3(l,k)+wgx*tt3(h,k))
            asgs_t1(i,k)=asgs_t1(i,k)+ftn*((1-wgx)*tsgst1(l,k)+ &
            wgx*tsgst1(h,k))
            asgs_t2(i,k)=asgs_t2(i,k)+ftn*((1-wgx)*tsgst2(l,k)+ &
            wgx*tsgst2(h,k))
            asgs_t3(i,k)=asgs_t3(i,k)+ftn*((1-wgx)*tsgst3(l,k)+ &
            wgx*tsgst3(h,k))
            aut(i,k)=aut(i,k)+ftn*((1-wgx)*tut(l,k)+ &
            wgx*tut(h,k))
            avt(i,k)=avt(i,k)+ftn*((1-wgx)*tvt(l,k)+ &
            wgx*tvt(h,k))
            awt(i,k)=awt(i,k)+ftn*((1-wgx)*twt(l,k)+ &
            wgx*twt(h,k))
            adtdx(i,k)=adtdx(i,k)+ftn*((1-wgx)*tdtdx(l,k)+ &
            wgx*tdtdx(h,k))
            adtdy(i,k)=adtdy(i,k)+ftn*((1-wgx)*tdtdy(l,k)+ &
            wgx*tdtdy(h,k))
            adtdz(i,k)=adtdz(i,k)+ftn*((1-wgx)*tdtdz(l,k)+ &
            wgx*tdtdz(h,k))
            arg2=((1-wgx)*tCs2Pr(l,k)+wgx*tCs2Pr(h,k))
            if(arg2 == 0)then
                aPr(i,k)=aPr(i,k)
            else
                aPr(i,k)=aPr(i,k)+fr*((1-wgx)*tCs2(l,k)+ &
                wgx*tCs2(h,k))/arg2
            endif
            aCs2Pr(i,k)=aCs2Pr(i,k)+ftn*((1-wgx)*tCs2Pr(l,k)+ &
            wgx*tCs2Pr(h,k))
            abeta2(i,k)=abeta2(i,k)+ftn*((1-wgx)*tbeta2(l,k)+ &
            wgx*tbeta2(h,k))
            aET(i,k)=aET(i,k)+ftn*((1-wgx)*tET(l,k)+ &
            wgx*tET(h,k))

        enddo
    enddo

    if(vfact == 0)then
        do j=1,nyb
            do i=1,nx
                aqz_s(i,j)=aqz_s(i,j)+fr*sgs_t3(i,j,2)
            enddo
        enddo

    endif

    return
    end subroutine scalar_slice
