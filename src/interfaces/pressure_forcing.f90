    subroutine pressure_forcing(force_x,force_y,u,v,RHSx,RHSy,RHSx_f, &
    RHSy_f)
    real*8,dimension(:) :: force_x, force_y
    real*8,dimension(:,:,:) :: u,v,RHSx,RHSy,RHSx_f,RHSy_f
    end subroutine pressure_forcing
