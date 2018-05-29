    subroutine lagrng_sd_s(a,b,c,d,e,a_old,b_old,c_old,d_old,e_old, &
    beta_old,un,vn,wn,i,j,k,flag,sig_t)
    real*8, dimension(:,:,:):: a_old,b_old,c_old,d_old,e_old
    integer*4 :: flag
    real*8 :: a,b,c,d,e,un,vn,wn,sig_t,beta_old
    end subroutine lagrng_sd_s
