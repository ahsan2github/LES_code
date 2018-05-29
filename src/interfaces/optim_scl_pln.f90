    subroutine optim_scl_pln(Pr2,S11,S33,S22,S12,S13,S23,S_hat, &
    S_hatd,u_,v_,w_,L,t_,tx,ty,tz,betaa1)
    real*8,dimension(:,:,:):: S11,S22,S33,S12,S13,S23,S_hat,S_hatd, &
    u_,v_,w_,t_,tx,ty,tz,betaa1,Pr2
    real*8 :: L(:)
    end subroutine optim_scl_pln
