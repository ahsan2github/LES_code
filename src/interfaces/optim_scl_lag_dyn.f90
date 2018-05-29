    subroutine optim_scl_lag_dyn(Pr2,S11,S33,S22,S12,S13,S23,S_hat, &
    u_,v_,w_,L,t_,tx,ty,tz,KX_old,XX_old)
    real*8,dimension(:,:,:):: tz,tx,ty,S11,S22,S33,S12,S13,S23, &
    S_hat,u_,v_,w_,t_,KX_old,XX_old,Pr2
    real*8 :: L(:)
    end subroutine optim_scl_lag_dyn
