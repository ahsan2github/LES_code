    SUBROUTINE cyclic(a,b,c,alpha,beta,r,x)
    implicit none
    interface
    include 'interfaces/tridag.f90'
    end interface
    integer*4 :: i,n
    real*8 :: alpha,beta,fact,gamma
    real*8,dimension(:) :: a,b,c,r,x
    real*8,dimension(size(a,1)) :: bb,u,z

    n=size(a,1)

    gamma=-b(1)
    bb(1)=b(1)-gamma
    bb(n)=b(n)-alpha*beta/gamma
    do i=2,n-1
        bb(i)=b(i)
    enddo
    call tridag(a,bb,c,r,x)
    u(1)=gamma
    u(n)=alpha
    do i=2,n-1
        u(i)=0.
    enddo
    call tridag(a,bb,c,u,z)
    fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
    do i=1,n
        x(i)=x(i)-fact*z(i)
    enddo
    return
    end SUBROUTINE cyclic

!  (C) Copr. 1986-92 Numerical Recipes Software
          
