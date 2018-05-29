    subroutine optim_loc(Cs2,S11,S33,S22,S12,S13,S23,S,S_hat,S_hatd &
    ,u_,v_,w_,L,betaa)   !,bad_cs1)
    real*8,dimension(:,:,:):: &
    S,S11,S22,S33,S12,S13,S23,S_hat,S_hatd,w_,u_,v_,Cs2,betaa
    real*8,dimension(:)::L    !,bad_cs1
    end subroutine optim_loc
