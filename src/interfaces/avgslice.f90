    subroutine avgslice(u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz,dudx, &
    dvdz,dwdz,dwdx,Cs2,beta1,ESGS,DSGS,TL,au,av,aw,ap,u2,v2,w2, &
    p2,w3,atxx,atxz,atyy,atyz,atzz,atxy,auw,avw,auv,adudz,adudx, &
    advdz,adwdz,adwdx,e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL, &
    ilow,wgx)
    real*8,dimension(:,:,:)::u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz, &
    dudx,dvdz,dwdz,dwdx,Cs2,beta1,ESGS,DSGS,TL
    real*8,dimension(:,:)::au,av,aw,ap,u2,v2,w2,p2,w3,atxx,atxz, &
    atyy,atyz,atzz,atxy,auw,avw,auv,adudz,adudx,advdz,adwdz, &
    adwdx,e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL
    integer*4,dimension(:)::ilow
    real*8 :: wgx
    end subroutine avgslice
