    subroutine lagrng_sd(a,b,c,d,e,a_old,b_old,c_old,d_old,e_old, &
    beta_old,Tn,un,vn,wn,i,j,k,flag)
    real*8, dimension(:,:,:):: a_old,b_old,c_old,d_old,e_old
    real*8 :: a,b,c,d,e,un,vn,wn,beta_old,Tn
    integer*4 :: i,j,k,flag
    end subroutine lagrng_sd
