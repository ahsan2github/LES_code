    Subroutine CalcBeta (scalar,b)
    	
    use globals
    use scalars
    implicit none

    interface
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: i, j, k, ll
    real*8,dimension(:,:,:) :: b
    real*8,dimension(:,:,:,:) :: scalar
    real*8 :: above, below
    real*8,dimension(size(scalar,3)) :: t_bar,q_bar
    real*8,dimension(size(scalar,1),size(scalar,2),size(scalar,3)) :: &
    theta,q
    	
!...    Note Beta is stored on W nodes, but Scalar is on UVP nodes

    do ll=1,scalarCount
        if(ll == temperatureIndex)then
            theta=scalar(:,:,:,ll)
            call plane_avg(theta,t_bar)
        elseif(ll == moistureIndex)then
            q=scalar(:,:,:,ll)
            call plane_avg(q,q_bar)
        endif
    enddo

    do k=2,nzb+1
        if (t_bar(k) == (0.d0))then
            b(:,:,k) = 0.
        elseif(q_bar(k) == (0.d0) .OR. moistureIndex == 0)then
            do j=1,Nyb
                do i=1,nx
                    above=(theta(i,j,k)-t_bar(k))/(theta_0)
                    below=(theta(i,j,k-1)-t_bar(k-1))/(theta_0)
                    b(i,j,k)=g_hat*(above + below)*0.5
                end do
            end do
        else
            do j=1,Nyb
                do i=1,nx
                    above=(theta(i,j,k)-t_bar(k) + &
                    &    		 0.61d0*theta_0*(q(i,j,k)-q_bar(k)))/(theta_0)
                    below=(theta(i,j,k-1)-t_bar(k-1) + &
                    &    		 0.61d0*theta_0*(q(i,j,k-1)-q_bar(k-1)))/(theta_0)
                    b(i,j,k)=g_hat*(above + below)*0.5
                end do
            end do
        endif
    end do
    	
    if (vfact == 0) then
        	   
        do j=1,Nyb
            do i=1,nx
                b(i,j,2)=0
            enddo
        enddo
        	   
    endif
         
    return
    	
    end Subroutine CalcBeta

    				

