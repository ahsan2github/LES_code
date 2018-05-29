    subroutine interpLOG(ui,uu,kw,zo)

    use globals
    use wallBoundaryConditions
    implicit none

    real*8 :: kw,ui,uu,zo

    ui = uu*log((kw*dz*z_i)/zo)/log((u_star*dz*z_i)/zo)

    return
    end subroutine interpLOG
