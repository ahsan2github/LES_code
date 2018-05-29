!....subroutine to compute lagrangian averaging for scale-dependent model
!....flag is 0 for 2delta scale and 1 for 4delta scale
    subroutine lagrng_sd_s(a,b,c,d,e,a_old,b_old,c_old,d_old,e_old, &
    beta_old,un,vn,wn,i,j,k,flag,sig_t)

    use globals
    implicit none

    interface
    include './interfaces/interp3D.f90'
    end interface


    real*8, dimension(:,:,:):: a_old,b_old,c_old,d_old,e_old

    real*8 :: a,b,c,d,e,un,vn,wn,eps,Tn,beta_old,TLM,TMM, &
    xi,yi,zi,aint,bint,cint,dint,eint,X(size(a_old,1)+1), &
    Y(size(a_old,2)+1),Z(size(a_old,3)),beta_used,sig_t,iw,jw,kw
          
    integer*4 :: i,j,k,ii,jj,kk,iplus,jplus,flag

    if(flag == 0)then
        beta_used=beta_old
    else
        beta_used=beta_old*beta_old
    endif
             
    j=j+1

    TLM=a_old(i,j,k)+beta_used*b_old(i,j,k)
    TMM=c_old(i,j,k)+beta_used*d_old(i,j,k)+ &
    e_old(i,j,k)*beta_used*beta_used

    if(TLM <= 0.0 .OR. TMM < (1e-10) .OR. sig_t <= 0.0)then
        eps=0.0
    else
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        Tn=1.5*delta*(TLM*TMM)**(-1./4.)*sig_t
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        eps=(dtl/Tn)/(1.+dtl/Tn)
    endif

    xi=-un*dtl
    yi=-vn*dtl
    zi=-wn*dtl

!@@@@@@@@@@@@@@@@@@@@@New Code @@@@@@@@@@@@@@@@@@@@@@@@@
!     *********** x-comp *************
    if(xi < 0.0)then
        iw=(dx-abs(xi))*idx
        if(i == 1)then
            ii=nx
        else
            ii=i-1
        endif
    else
        iw=(xi)*idx
        ii=i
    endif
          
!     *********** y-comp *************
    if(yi < 0.0)then
        jw=(dy-abs(yi))*idy
        jj=j-1
    else
        jw=(yi)*idy
        jj=j
    endif
          
!     *********** z-comp *************
    if(vfact == 0 .AND. k <= 3 .AND. verticalBC == 0)then
        if(zi < 0.0)then
            kk=2
            if(k == 2)then
                kw=0.0
            else
                if(scl_nodes == 0)then
                    kw=(dz/2.-abs(zi))/(dz/2.)
                else
                    kw=(dz-abs(zi))*idz
                endif
            endif
        else
            if(k == 2)then
                if(scl_nodes == 0)then
                    kw=(zi)/(dz/2.)
                else
                    kw=(zi)*idz
                endif
                kk=2
            else
                kw=(zi)*idz
                kk=3
            endif
        endif
    elseif(vfact == vprocs-1 .AND. k == nzb+1 .AND. verticalBC == 0)then
        if(zi < 0.0)then
            kw=(dz-abs(zi))*idz
            kk=nzb
        else
            kw=0
            kk=nzb+1
        endif
    else
        if(zi < 0.0)then
            kw=(dz-abs(zi))*idz
            kk=k-1
        else
            kw=(zi)*idz
            kk=k
        endif
    endif
!     ********************************
          
!      if(t.le.1.and.INITU.eq.0)then
!         kw=0.0
!         kk=k
!      endif
          
    if(t <= 1 .AND. INITS == 0)then
        kw=0.0
        kk=k
    endif
              
    if(ii == nx)then
        iplus=-(nx-1)
    else
        iplus=1
    endif
          
    jplus=1
     
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,a_old,aint)
    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,b_old,bint)
    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,c_old,cint)
    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,d_old,dint)
    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,e_old,eint)

    a=eps*a+(1.-eps)*aint
    b=eps*b+(1.-eps)*bint
    c=eps*c+(1.-eps)*cint
    d=eps*d+(1.-eps)*dint
    e=eps*e+(1.-eps)*eint

    j=j-1

    return
          
    end subroutine lagrng_sd_s
          










