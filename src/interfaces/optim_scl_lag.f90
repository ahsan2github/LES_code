    subroutine optim_scl_lag(Pr2,S11,S33,S22,S12,S13,S23,S_hat, &
    S_hatd,u_,v_,w_,L,t_,tx,ty,tz,betaa1,a2_old,b2_old, &
    c2_old,d2_old,e2_old,a4_old,b4_old,c4_old,d4_old,e4_old)
    real*8,dimension(:,:,:):: Pr2,tz,tx,ty,S11,S22,S33,S12,S13,S23, &
    S_hat,t_,S_hatd,w_,u_,v_,betaa1, &
    a2_old,b2_old,c2_old,d2_old,e2_old,a4_old,b4_old,c4_old, &
    d4_old,e4_old
    real*8,dimension(:) :: L
    end subroutine optim_scl_lag
