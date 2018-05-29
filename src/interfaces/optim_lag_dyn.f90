    subroutine optim_lag_dyn(Cs2,S11,S33,S22,S12,S13,S23,S,S_hat, &
    u_,v_,w_,L,LM_old,MM_old)
    real*8,dimension(:,:,:):: &
    Cs2,S11,S22,S33,S12,S13,S23,S,S_hat,w_,u_,v_,LM_old,MM_old
    real*8 :: L(:)
    end subroutine optim_lag_dyn
