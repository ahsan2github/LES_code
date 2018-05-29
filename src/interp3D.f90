!  recieves an array and a set of three global points and returns the value
!  interpolated using linear interpolation for the new point
!  note that ot=0 for momentum and 1 for scalars...

    subroutine interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,U,ui)

    use globals
    implicit none

    real*8 :: ui,U(:,:,:),iw,jw,kw, &
    ui_hgh,ui_low
    integer*4 :: ii,jj,kk,iplus,jplus
    	
!       use periodic boundary condition in x and y

    if(kw == 0.0)then
        ui =    U(ii,jj,kk)*(1.-iw)*(1.-jw)+ &
        U(ii+iplus,jj,kk)*iw*(1.-jw)+ &
        U(ii+iplus,jj+jplus,kk)*iw*jw+ &
        U(ii,jj+jplus,kk)*(1.-iw)*jw
    else
        ui_low = U(ii,jj,kk)*(1.-iw)*(1.-jw)+ &
        U(ii+iplus,jj,kk)*iw*(1.-jw)+ &
        U(ii+iplus,jj+jplus,kk)*iw*jw+ &
        U(ii,jj+jplus,kk)*(1.-iw)*jw 

        ui_hgh = U(ii,jj,kk+1)*(1.-iw)*(1.-jw)+ &
        U(ii+iplus,jj,kk+1)*iw*(1.-jw)+ &
        U(ii+iplus,jj+jplus,kk+1)*iw*jw+ &
        U(ii,jj+jplus,kk+1)*(1.-iw)*jw

        ui = ui_low*(1.-kw)+ui_hgh*kw
    endif

    return
    	
    end subroutine interp3D



