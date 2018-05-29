function getBLH(auw, avw, atxz, atyz) result(outBLH)
    
    use globals
    use sundry
    implicit none
    
    real(kind=8), dimension(:,:) :: auw, avw, atxz, atyz
    real(kind=8), dimension(nzb) :: muw, mvw, u_shear, zh2 
    real(kind=8), dimension(nzb) :: planeAvg_muw, planeAvg_mvw
    real(kind=8) ::  a1, b1, a2, b2, a, outBLH, maxshear ! output is BLH, BLH is defined in sundry module
    integer(kind = 4) :: i

    muw = (sum(auw(:, 2:nzb+1),1) + sum(atxz(:, 2:nzb+1),1))/aNx
    mvw = (sum(avw(:, 2:nzb+1),1) + sum(atyz(:, 2:nzb+1),1))/aNx

    call MPI_ALLREDUCE(muw,planeAvg_muw,nzb,MPI_DOUBLE_PRECISION,MPI_SUM, &
    &      MPI_COMM_LEVEL,ierr)
    planeAvg_muw = planeAvg_muw / hprocs

    call MPI_ALLREDUCE(mvw,planeAvg_mvw,nzb,MPI_DOUBLE_PRECISION,MPI_SUM, &
    &      MPI_COMM_LEVEL,ierr)
    planeAvg_mvw = planeAvg_mvw / hprocs

    u_shear = (planeAvg_muw**2+planeAvg_mvw**2)**0.25;
    zh2 = (/(i, i=(vfact*nzb),(vfact+1)*nzb-1,1)/)/real(Nz)

    ! process 0 always contains the lowest portions of the boundary layer
    if (me == 0) then 
        groundShear = u_shear(1)        
    end if

    call MPI_BCAST(groundShear, 1, MPI_Double_Precision, 0, nall, ierr)
 
    do i=1, nzb
        ! this if condition will be true for multiple processes
        if (u_shear(i)**2 <= 0.05*groundShear**2) then
            a1 = zH2(i)
            b1 = u_shear(i)**2            
            a2 = zh2(i-1)
            b2 = u_shear(i-1)**2
            a = a2 + (a1-a2)*(0.05*groundShear**2-b2)/(b1-b2)
            outBLH = a /0.95d0 ! non dimensional Boundary layer height
            exit
        end if        
    end do
    call MPI_ALLREDUCE(outBLH, BLH, 1, MPI_Double_Precision, MPI_MIN, nall, ierr)
    outBLH = BLH
    call MPI_BARRIER(nall,ierr)
!    write(*, *) 'groundShear :: BLH  === ', groundShear, '::', BLH * l_z

end function getBLH




! muw      = mauw+mtuw;
! mvw      = mavw+mtvw;
! u_sL     = (muw.^2+mvw.^2).^0.25;
! u_s      = u_sL(1);
! 
! Ind = find(u_sL.^2 <= 0.05*(u_sL(1).^2));
!     a1 = zH_2(Ind(1))*L_z; b1 = u_sL(Ind(1))^2; a2 = zH_2(Ind(1)-1)*L_z; b2 = u_sL(Ind(1)-1)^2;
!     a = a2 + (a1-a2)*(0.05*(u_sL(1)^2)-b2)/(b1-b2); mBLA=a/0.95; clear Ind;
!     display(['boundary layer height=',num2str(mBLA)])
