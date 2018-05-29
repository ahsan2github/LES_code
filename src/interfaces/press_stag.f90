    subroutine Press_stag(P, RHSx, RHSy, RHSz, RHSx_f, RHSy_f, &
    RHSz_f,u,v,w,DFDX,DFDY,divtz)
    real*8,dimension(:,:,:):: RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f, &
    u,v,w,DFDX,DFDY,divtz,P
    end subroutine Press_stag
