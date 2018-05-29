    subroutine sgs_stag(qx,qy,qz,u,v,w,dudx,dudy,dudz,dvdx,dvdy, &
    dvdz,dwdx,dwdy,dwdz,dsdx,dsdy,dsdz,pt,Cs2,Cs2_m,Pr2,Pr2_m, &
    cst,txx,txy,txz,tyy,tyz,tzz,scalar, &
    a1_old,b1_old,c1_old,d1_old,e1_old, &
    a2_old,b2_old,c2_old,d2_old,e2_old, &
    a4_old,b4_old,c4_old,d4_old,e4_old, &
    a8_old,b8_old,c8_old,d8_old,e8_old, &
    beta1,beta2,ESGS3D,DSGS3D,ET3D,TL)
    integer*4 :: pt,cst
    real*8,dimension(:,:,:):: dudx,dudy,dudz,dvdx,dvdy,dvdz, &
    dwdx,dwdy,dwdz,u,v,w,txx,txy,txz,tyy,tyz,tzz, &
    a1_old,b1_old,c1_old,d1_old,e1_old, &
    a2_old,b2_old,c2_old,d2_old,e2_old,Cs2,Cs2_m,beta1,ESGS3D, &
    DSGS3D,TL
    real*8,dimension(:,:,:,:)::qx,qy,qz,dsdx,dsdy,dsdz, &
    scalar,ET3D,Pr2,Pr2_m,beta2, &
    a4_old,b4_old,c4_old,d4_old,e4_old, &
    a8_old,b8_old,c8_old,d8_old,e8_old
    end subroutine sgs_stag
