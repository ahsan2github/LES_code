    subroutine optim_lag(Cs2,S11,S33,S22,S12,S13,S23,S,S_hat, &
    S_hatd,u_,v_,w_,L,betaa,a1_old,b1_old,c1_old, &
    d1_old,e1_old,a2_old,b2_old,c2_old,d2_old,e2_old,TL2)
    real*8,dimension(:,:,:):: &
    Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd,u_,v_,w_,betaa, &
    a1_old,b1_old,c1_old,d1_old,e1_old,a2_old,b2_old,c2_old, &
    d2_old,e2_old,TL2
    real*8,dimension(:) :: L
    end subroutine optim_lag
