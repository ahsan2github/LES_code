    subroutine scalar_RHS(s,u_m,v_m,w_m,qx,qy,qz,txz,tyz,dsdx,dsdy, &
    dsdz,RHS,Surf_flux)
    real*8, dimension(:,:,:):: dsdx,dsdy,dsdz,RHS,s,txz,tyz, &
    qx,qy,qz,u_m,v_m,w_m
    real*8, dimension(:,:)::Surf_flux
    end subroutine scalar_RHS
