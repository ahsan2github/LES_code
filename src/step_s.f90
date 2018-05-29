    subroutine STEP_S (S, RHS, RHS_f, scalarInd)

    use globals
    use scalars
    implicit none

    integer*4 :: i,j,k,scalarInd
    real*8,dimension(:,:,:):: s, RHS, RHS_f

    do k=2,Nzb+1
        do j=1,Nyb
            do i=1,Nx
                s(i,j,k)= s(i,j,k)+DT*(1.5*RHS(i,j,k)-0.5*RHS_f(i,j,k))
            end do
        end do
    end do

!...  No-stress top
    if (vfact==vprocs-1 .AND. verticalBC == 0) then
        do j=1,Nyb
            do i=1,Nx
                s(i,j,Nzb+1)=s(i,j,Nzb)+inversion(scalarInd)*dz
            end do
        end do
    end if

    return
    end subroutine STEP_S
