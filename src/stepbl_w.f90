    subroutine STEPBL_w (ui,RHSi,RHSi_f)

    use globals
    implicit none

    interface
    subroutine plane_reduce(XX)
    real*8 :: XX
    end subroutine plane_reduce
    end interface
      
    integer*4 :: i,j,k
    real*8,dimension(:,:,:)::ui,RHSi,RHSi_f
    real*8 :: w_bar

    do k=2,Nzb+1

        do j=1,Nyb
            do i=1,Nx
                ui(i,j,k)=ui(i,j,k)+DT*(1.5d0*RHSi(i,j,k)-0.5d0*RHSi_f(i,j,k))
            end do
        end do

    !	  w_bar=0.0
    !          do j=1,Nyb
    !              do i=1,Nx
    !		w_bar= w_bar+ui(i,j,k)
    !              end do
    !          end do
    !          call plane_reduce(w_bar)
    !          w_bar=w_bar*iNxNy
        	  
    !........And, remove the mean vertical velocity
    !          do j=1,Nyb
    !              do i=1,Nx
    !		  ui(i,j,k)= ui(i,j,k)-w_bar
    !	      end do
    !          end do

    end do
    	
!...No slip wall and no flow through top

    if(verticalBC == 0)then
        if (vfact==0) then
            ui(:,:,2)=0.d0
        end if
                 
        if (vfact==vprocs-1) then
            ui(:,:,nzb+1)=0.d0
        end if
    elseif(verticalBC == 1)then
        if(vfact == 0)then
            call MPI_SEND(ui(1,1,2),size(ui,1)*size(ui,2), &
            MPI_DOUBLE_PRECISION,(nprocs-1)-(hprocs-1-me), &
            (nprocs-1)-(hprocs-1-me),nall,ierr )
        elseif(vfact == vprocs-1)then
            call MPI_RECV(ui(1,1,nzb+1),size(ui,1)*size(ui,2), &
            MPI_DOUBLE_PRECISION,me-hfact,me,nall, &
            MPI_STATUS_IGNORE,ierr )
        endif
    endif
    	
    return
    end subroutine  STEPBL_w
