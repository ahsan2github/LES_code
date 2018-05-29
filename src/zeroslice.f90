    subroutine  zeroslice(au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy, &
    atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz,adwdx, &
    e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL)
          
    implicit none

    real*8,dimension(:,:) :: au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy, &
    atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz,adwdx, &
    e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL

    au = 0.d0
    av = 0.d0
    aw = 0.d0
    ap = 0.d0
    u2 = 0.d0
    v2 = 0.d0
    w2 = 0.d0
    w3 = 0.d0
    atxx = 0.d0
    atxz = 0.d0
    atyy = 0.d0
    atyz = 0.d0
    atzz = 0.d0
    atxy = 0.d0
    p2 = 0.d0
    auw = 0.d0
    avw = 0.d0
    auv = 0.d0
    adudz = 0.d0
    adudx = 0.d0
    advdz = 0.d0
    adwdz = 0.d0
    adwdx = 0.d0
    e = 0.d0
    aCs2 = 0.d0
    aCs = 0.d0
    abeta1 = 0.d0
    atxz_s   = 0.d0
    aESGS = 0.d0
    aDSGS = 0.d0
    aTL = 0.d0

    return
    end subroutine 
