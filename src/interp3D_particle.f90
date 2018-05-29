    subroutine interp3D_particle(ui,iw,jw,kw,U)

    implicit none
     
    real*8 :: iw,jw,kw,u_low,u_high,ui
    real*8,dimension(:,:,:) :: U

    if(kw == 0.d0)then
        ui = (1-iw)*(1-jw)*U(1,1,1)+iw*(1-jw)*U(2,1,1)+ &
        iw*jw*U(2,2,1)+(1-iw)*jw*U(1,2,1)
    else
        u_low = (1-iw)*(1-jw)*U(1,1,1)+iw*(1-jw)*U(2,1,1)+ &
        iw*jw*U(2,2,1)+(1-iw)*jw*U(1,2,1)
        u_high = (1-iw)*(1-jw)*U(1,1,2)+iw*(1-jw)*U(2,1,2)+ &
        iw*jw*U(2,2,2)+(1-iw)*jw*U(1,2,2);
        ui = (u_high-u_low)*kw+u_low
    endif

!$$$      ui = ((1-iw)*(1-jw)*U(1,1,2)+iw*(1-jw)*U(2,1,2)+iw*jw*U(2,2,2)+
!$$$     +     (1-iw)*jw*U(1,2,2) - (1-iw)*(1-jw)*U(1,1,1)+
!$$$     +     iw*(1-jw)*U(2,1,1)+iw*jw*U(2,2,1)+(1-iw)*jw*U(1,2,1))*kw+
!$$$     +     (1-iw)*(1-jw)*U(1,1,1)+iw*(1-jw)*U(2,1,1)+iw*jw*U(2,2,1)+
!$$$     +     (1-iw)*jw*U(1,2,1)

    return
    end subroutine interp3D_particle
