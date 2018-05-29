    Subroutine CONVEC (Cx,Cy,Cz,u1_m,u2_m,u3_m, &
    du1,du2,du3,du4,du5,du6)

    use globals
    implicit none

    interface
    include './interfaces/dealias1.f90'
    include './interfaces/dealias2.f90'
    include './interfaces/update1.f90'
    include './interfaces/update3.f90'
    end interface

    integer*4 :: i,j,k
    Real*8,dimension (:,:,:):: du1,du2,du3,du4,du5,du6, &
    cx,cy,cz
    real*8,dimension (:,:,:):: u1_m,u2_m,u3_m
    real*8,dimension (size(u1_m,1),size(u1_m,2),size(u1_m,3)):: &
    cc,du1_m,du2_m,du3_m,du4_m,du5_m,du6_m
    real*8 :: arg1, arg2a, arg2b, arg2
          
    call dealias1(du1,du1_m)
    call dealias1(du2,du2_m)
    call dealias1(du3,du3_m)
    call dealias1(du4,du4_m)
    call dealias1(du5,du5_m)
    call dealias1(du6,du6_m)

    call update3(du2_m,du4_m,du5_m)
    call update1(du6_m)
!     ************************************************************
    Do k=2,nzb+1
        Do j=1,nyb2
            Do i=1,nx2
                arg1=u2_m(i,j,k)*(du1_m(i,j,k)-du3_m(i,j,k))
                arg2a=u3_m(i,j,k+1)*(du2_m(i,j,k+1)-du5_m(i,j,k+1))
                arg2b=u3_m(i,j,k)*(du2_m(i,j,k)-du5_m(i,j,k))
                cc(i,j,k)=arg1+0.5*(arg2a+arg2b)
            End Do
        End Do
    End Do
          
    if(verticalBC == 0)then
        Do j=1,nyb2
            Do i=1,nx2
                            
                IF (vfact==0) then
                    k=2
                    arg1=u2_m(i,j,k)*(du1_m(i,j,k)-du3_m(i,j,k))
                    arg2a=u3_m(i,j,k+1)*(du2_m(i,j,k+1)-du5_m(i,j,k+1))
                    arg2b=0.
                    cc(i,j,k)=arg1+0.5*(arg2a+arg2b)
                END IF
                               
                IF (vfact==vprocs-1) then
                    k=Nzb+1
                    arg1=u2_m(i,j,k)*(du1_m(i,j,k)-du3_m(i,j,k))
                    arg2a=u3_m(i,j,k)*(du2_m(i,j,k)-du5_m(i,j,k))
                    cc(i,j,k)=arg1+arg2a
                END IF
                               
            End Do
        End Do
    endif

    call dealias2(Cx,cc)

!     **********************************************************
          
    Do k=2,nzb+1
        Do j=1,nyb2
            Do i=1,nx2
                arg1=u1_m(i,j,k)*(du3_m(i,j,k)-du1_m(i,j,k))
                arg2a=u3_m(i,j,k+1)*(du4_m(i,j,k+1)-du6_m(i,j,k+1))
                arg2b=u3_m(i,j,k)*(du4_m(i,j,k)-du6_m(i,j,k))
                cc(i,j,k)=arg1+0.5*(arg2a+arg2b)
            End Do
        End Do
    End Do
          
    if(verticalBC == 0)then
        Do j=1,nyb2
            Do i=1,nx2
                               
                if (vfact==0) then
                    k=2
                    arg1=u1_m(i,j,k)*(du3_m(i,j,k)-du1_m(i,j,k))
                    arg2a=u3_m(i,j,k+1)*(du4_m(i,j,k+1)-du6_m(i,j,k+1))
                    arg2b=0.
                    cc(i,j,k)=arg1+0.5*(arg2a+arg2b)
                end if
                               
                if (vfact==vprocs-1) then
                    k=Nzb+1
                    arg1=u1_m(i,j,k)*(du3_m(i,j,k)-du1_m(i,j,k))
                    arg2a=u3_m(i,j,k)*(du4_m(i,j,k)-du6_m(i,j,k))
                    cc(i,j,k)=arg1+arg2a
                end if
                               
            End Do
        End Do
    endif
          
    call dealias2(Cy,cc)

!     *******************************************************

    Do k=2,nzb+1
        Do j=1,nyb2
            Do i=1,nx2
                arg1=(0.5*(u1_m(i,j,k)+u1_m(i,j,k-1)))* &
                (du5_m(i,j,k)-du2_m(i,j,k))
                arg2=(0.5*(u2_m(i,j,k)+u2_m(i,j,k-1)))* &
                (du6_m(i,j,k)-du4_m(i,j,k))
                cc(i,j,k)=arg1+arg2
            End Do
        End Do
    End Do

    if(verticalBC == 0)then
        do j=1,nyb2
            do i=1,nx2
                               
                if(vfact==0) then
                    k=2
                    cc(i,j,k)=0.
                end if
                               
                if(vfact==vprocs-1)then
                    k=nzb+1
                    cc(i,j,k)=0
                end if
                               
            end do
        end do
    endif
          
    call dealias2(Cz,cc)
          
    Return
    end Subroutine CONVEC
          
