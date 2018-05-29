    SUBROUTINE tridag(a,b,c,r,u)
    implicit none
    INTEGER*4 :: j,n
    REAL*8 :: bet
    REAL*8,dimension(:) ::  a,b,c,r,u
    REAL*8,dimension(size(a,1)) :: gam

    n=size(a,1)

    bet=b(1)
    u(1)=r(1)/bet
    do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet == 0.) then
            print *, 'tridag failed at k=',j
            print *,'a, b, c, gam, and bet=',a(j),b(j),c(j),gam(j),bet
        end if
        u(j)=(r(j)-a(j)*u(j-1))/bet
    11 END DO
    do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    12 END DO
    return
    end SUBROUTINE tridag
!  (C) Copr. 1986-92 Numerical Recipes Software ]2#"0>Ya%.
