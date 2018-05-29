!....subroutine to compute lagrangian averaging for dynamic model

    subroutine lagrng_dyn(a,b,a_old,b_old,un,vn,wn,i,j,k)
          
    use globals
    implicit none

    interface
    include './interfaces/interp3D.f90'
    end interface

    real*8, dimension(:,:,:):: a_old,b_old

    real*8 :: a,b,un,vn,wn,eps,Tn,TLM,TMM, &
    xi,yi,zi,aint,bint,iw,jw,kw
          
    integer*4 :: i,j,k,ii,jj,kk,iplus,jplus

    j=j+1

    TLM=a_old(i,j,k)
    TMM=b_old(i,j,k)

    if(TLM <= 0.0 .OR. TMM < (1e-10))then
        eps=0.0
    else
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        Tn=1.5*delta*(TLM*TMM)**(-1./8.)
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        eps=(dtl/Tn)/(1.+dtl/Tn)
    endif

    xi=-un*dtl
    yi=-vn*dtl
    zi=-wn*dtl
          
!@@@@@@@@@@@@@@@@@@@@@ New Code @@@@@@@@@@@@@@@@@@@@@@@@@
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
                if(mom_nodes == 0)then
                    kw=(dz/2.-abs(zi))/(dz/2.)
                else
                    kw=(dz-abs(zi))*idz
                endif
            endif
        else
            if(k == 2)then
                if(mom_nodes == 0)then
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
          
    if(t <= 1 .AND. INITU == 0)then
        kw=0.0
        kk=k
    endif
          
!      if(t.le.1.and.INITS.eq.0)then
!         kw=0.0
!         kk=k
!      endif
              
    if(ii == nx)then
        iplus=-(nx-1)
    else
        iplus=1
    endif
          
    jplus=1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,a_old,aint)
    call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,b_old,bint)

          
    a=eps*a+(1.-eps)*aint
    b=eps*b+(1.-eps)*bint

    j=j-1

    return
          
    end subroutine lagrng_dyn
          










