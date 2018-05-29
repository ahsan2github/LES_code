    subroutine lagrng_dyn_s(a,b,a_old,b_old, &
    un,vn,wn,i,j,k,sig_t)
    real*8, dimension(:,:,:):: a_old,b_old
    real*8 :: a,b,un,vn,wn,sig_t
    integer*4 :: i,j,k,ii
    end subroutine lagrng_dyn_s
