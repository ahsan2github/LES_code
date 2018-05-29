    subroutine scalar_slice(u,v,w,theta,sgs_t1,sgs_t2,sgs_t3,dtdx, &
    dtdy,dtdz,Pr2,Cs2,beta2,at,t2,t3,asgs_t1,asgs_t2,asgs_t3, &
    aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aqz_s, &
    ET3D,aET,ilow,wgx)
    real*8,dimension(:,:,:):: u,v,w,theta,sgs_t1,sgs_t2,sgs_t3, &
    dtdx,dtdy,dtdz,Pr2,Cs2,beta2,ET3D
    real*8,dimension(:,:) :: at,t2,t3,asgs_t1,asgs_t2,asgs_t3, &
    aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aET
    real*8,dimension(:,:) :: aqz_s
    integer*4,dimension(:) :: ilow
    real*8 :: wgx
    end subroutine scalar_slice
