    subroutine  s_zeroslice(at,t2,t3,asgs_t1,asgs_t2,asgs_t3,aut,avt, &
    awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aqz_s,aET)

    implicit none
          
    real*8, dimension(:,:,:) :: at,t2,t3,asgs_t3,asgs_t1, &
    asgs_t2,aut,avt,awt,adtdx,adtdy,adtdz,aPr,aCs2Pr,abeta2,aET
    real*8, dimension(:,:,:) :: aqz_s

    at = 0.0
    t2 = 0.0
    t3 = 0.0
    asgs_t1 = 0.0
    asgs_t2 = 0.0
    asgs_t3 = 0.0
    aut = 0.0
    avt = 0.0
    awt = 0.0
    adtdx = 0.0
    adtdy = 0.0
    adtdz = 0.0
    aPr = 0.0
    aCs2Pr = 0.0
    abeta2 = 0.0
    aqz_s = 0.0
    aET = 0.0

    return
    end subroutine 
