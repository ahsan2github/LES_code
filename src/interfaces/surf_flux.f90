    subroutine surf_flux(theta,thetaFactor,u,v,t_flux,Psi,Psi0,fi, &
    fiH,zo,t_s,coolrate,scalarIndex)
    integer*4 :: scalarIndex
    real*8,dimension(:,:,:)::theta,thetaFactor,u,v,Psi,Psi0,fi,fiH
    real*8,dimension(:,:)::t_flux,zo,t_s,coolrate
    end subroutine surf_flux
