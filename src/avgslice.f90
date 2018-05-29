      subroutine avgslice(u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz,dudx,    &
     &     dvdz,dwdz,dwdx,Cs2,beta1,ESGS,DSGS,TL,au,av,aw,ap,u2,v2,w2,  &
     &     p2,w3,atxx,atxz,atyy,atyz,atzz,atxy,auw,avw,auv,adudz,adudx, &
     &     advdz,adwdz,adwdx,e,aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL,  &
     &     ilow,wgx)                                                    
                                                                        
      use globals 
      implicit none 
                                                                        
      interface 
         include './interfaces/plane_avg.f90' 
      end interface 
                                                                        
      real*8,dimension(:,:,:)::u,v,w,p,txx,txz,tyy,tyz,tzz,txy,dudz,    &
     &     dudx,dvdz,dwdz,dwdx,Cs2,beta1,ESGS,DSGS,TL                   
                                                                        
      real*8,dimension(:,:)::au,av,aw,ap,u2,v2,w2,p2,w3,atxx,atxz,atyy, &
     &     atyz,atzz,atxy,auw,avw,auv,adudz,adudx,advdz,adwdz,adwdx,e,  &
     &     aCs2,aCs,abeta1,atxz_s,aESGS,aDSGS,aTL                       
                                                                        
      integer*4,dimension(:)::ilow 
                                                                        
      integer*4 i,j,k,l,kk,h 
                                                                        
      real*8,dimension(size(au,1),size(au,2)) :: tu1,tv1,tw1,tp1,tu2,   &
     &     tv2,tw2,tp2,tw3,ttxx,ttxz,ttyy,ttyz,ttzz,ttxy,tuw,tvw,tuv,   &
     &     tdudz,tdudx,te,tCs2,tbeta1,tCs,tdwdz,tdwdx,tdvdz,tESGS,tDSGS,&
     &     tTL                                                          
      real*8,dimension(size(u,3)) :: u_bar,v_bar,w_bar 
      real*8 fr,arg1,arg2,norm,DDD,xpnt,ftn,wgx 
                                                                        
      fr   =(1./p_count)*c_count 
      norm =1.d0/(Ny*(1+Nx-aNx)) 
      ftn=fr*norm 
                                                                        
      DDD = dmod(dble(t)*dt*Ugal,dx) 
      wgx=(DDD)/dx 
                                                                        
!cc compute the plane averages of terms involved in products ccc        
                                                                        
      call plane_avg(u,u_bar) 
      call plane_avg(v,v_bar) 
      call plane_avg(w,w_bar) 
                                                                        
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
                                                                        
      tu1=0. 
      tv1=0. 
      tw1=0. 
      tp1=0. 
      tu2=0. 
      tv2=0. 
      tw2=0. 
      tp2=0. 
      tw3=0. 
      ttxx=0. 
      ttxz=0. 
      ttyy=0. 
      ttyz=0. 
      ttzz=0. 
      ttxy=0. 
      tuw=0. 
      tvw=0. 
      tuv=0. 
      tdudz=0. 
      tdudx=0. 
      tdvdz=0. 
      tdwdz=0. 
      tdwdx=0. 
      te=0. 
      tCs2=0. 
      tCs=0. 
      tbeta1=0. 
      tESGS=0. 
      tDSGS=0. 
      tTL=0. 
                                                                        
      do k=2,Nzb+1 
         do i=1,aNx 
            do j=1,Nyb  	 
               do l=1,1+Nx-aNx 
                  tu1(i,k)=tu1(i,k)+u(i-1+l,j,k)+Ugal 
                  tv1(i,k)=tv1(i,k)+v(i-1+l,j,k)+Vgal 
                  tw1(i,k)=tw1(i,k)+w(i-1+l,j,k) 
                  tp1(i,k)=tp1(i,k)+p(i-1+l,j,k) 
                  tp2(i,k)=tp2(i,k)+p(i-1+l,j,k)*p(i-1+l,j,k) 
                  ttxx(i,k)=ttxx(i,k)+txx(i-1+l,j,k) 
                  ttxz(i,k)=ttxz(i,k)+txz(i-1+l,j,k) 
                  ttyy(i,k)=ttyy(i,k)+tyy(i-1+l,j,k) 
                  ttyz(i,k)=ttyz(i,k)+tyz(i-1+l,j,k) 
                  ttzz(i,k)=ttzz(i,k)+tzz(i-1+l,j,k) 
                  ttxy(i,k)=ttxy(i,k)+txy(i-1+l,j,k) 
                  tdudz(i,k)=tdudz(i,k)+dudz(i-1+l,j,k) 
                  tdudx(i,k)=tdudx(i,k)+dudx(i-1+l,j,k) 
                  tdvdz(i,k)=tdvdz(i,k)+dvdz(i-1+l,j,k) 
                  tdwdz(i,k)=tdwdz(i,k)+dwdz(i-1+l,j,k) 
                  tdwdx(i,k)=tdwdx(i,k)+dwdx(i-1+l,j,k) 
                  tCs2(i,k)=tCs2(i,k)+Cs2(i-1+l,j,k) 
                  tCs(i,k)=tCs(i,k)+Cs2(i-1+l,j,k)**0.5 
                  tbeta1(i,k)=tbeta1(i,k)+beta1(i-1+l,j,k) 
                  tESGS(i,k)=tESGS(i,k)+ESGS(i-1+l,j,k) 
                  tDSGS(i,k)=tDSGS(i,k)+DSGS(i-1+l,j,k) 
                  tTL(i,k)=tTL(i,k)+TL(i-1+l,j,k) 
                                                                        
                  if (k==2.and.vfact==0) then 
                     arg1=0. 
                     arg2=0. 
                  else 
                     arg1=0.5d0*(u(i-1+l,j,k)+u(i-1+l,j,k-1) -          &
     &                           u_bar(k) - u_bar(k-1))                 
                     arg2=0.5d0*(v(i-1+l,j,k)+v(i-1+l,j,k-1) -          &
     &                           v_bar(k) - v_bar(k-1))                 
                  end if 
                                                                        
                  tuw(i,k)=tuw(i,k)+(w(i-1+l,j,k) - w_bar(k))*arg1 
                  tvw(i,k)=tvw(i,k)+(w(i-1+l,j,k) - w_bar(k))*arg2 
                  te(i,k)=te(i,k)+arg1*arg1 + arg2*arg2 +               &
     &                  (w(i-1+l,j,k) - w_bar(k))**2                    
                                                                        
                  tu2(i,k)=tu2(i,k)+(u(i-1+l,j,k) - u_bar(k))**2 
                  tv2(i,k)=tv2(i,k)+(v(i-1+l,j,k) - v_bar(k))**2 
                  tw2(i,k)=tw2(i,k)+(w(i-1+l,j,k) - w_bar(k))**2 
                  tw3(i,k)=tw3(i,k)+(w(i-1+l,j,k) - w_bar(k))**3 
                  tuv(i,k)=tuv(i,k)+(u(i-1+l,j,k) - u_bar(k))*          &
     &                              (v(i-1+l,j,k) - v_bar(k))           
                                                                        
               enddo 
            enddo 
         enddo 
      enddo 
                                                                        
      do kk=2,Nzb+1 
!ccc NOTE WE ARE JUST USING STEPS HERE (NO INTERPOLATION)               
         do i=1,aNx 
            if(aNx.ne.1)then 
               l=ilow(i) 
               h=ilow(i) 
               if(l.eq.nx)then 
                  h=nx 
               endif 
            else 
               l=1 
               h=1 
            endif 
            if(nprocs.eq.1)then 
               k=kk-1 
            else 
               k=kk 
            endif 
                                                                        
            au(i,k)=au(i,k)+ftn*((1-wgx)*tu1(l,k)+wgx*tu1(h,k)) 
            av(i,k)=av(i,k)+ftn*((1-wgx)*tv1(l,k)+wgx*tv1(h,k)) 
            aw(i,k)=aw(i,k)+ftn*((1-wgx)*tw1(l,k)+wgx*tw1(h,k)) 
            ap(i,k)=ap(i,k)+ftn*((1-wgx)*tp1(l,k)+wgx*tp1(h,k)) 
            u2(i,k)=u2(i,k)+ftn*((1-wgx)*tu2(l,k)+wgx*tu2(h,k)) 
            v2(i,k)=v2(i,k)+ftn*((1-wgx)*tv2(l,k)+wgx*tv2(h,k)) 
            w2(i,k)=w2(i,k)+ftn*((1-wgx)*tw2(l,k)+wgx*tw2(h,k)) 
            p2(i,k)=p2(i,k)+ftn*((1-wgx)*tp2(l,k)+wgx*tp2(h,k)) 
            w3(i,k)=w3(i,k)+ftn*((1-wgx)*tw3(l,k)+wgx*tw3(h,k)) 
                                                                        
            atxx(i,k)=atxx(i,k)+ftn*((1-wgx)*ttxx(l,k)+                 &
     &           wgx*ttxx(h,k))                                         
            atxz(i,k)=atxz(i,k)+ftn*((1-wgx)*ttxz(l,k)+                 &
     &           wgx*ttxz(h,k))                                         
            atyy(i,k)=atyy(i,k)+ftn*((1-wgx)*ttyy(l,k)+                 &
     &           wgx*ttyy(h,k))                                         
            atyz(i,k)=atyz(i,k)+ftn*((1-wgx)*ttyz(l,k)+                 &
     &           wgx*ttyz(h,k))                                         
            atzz(i,k)=atzz(i,k)+ftn*((1-wgx)*ttzz(l,k)+                 &
     &           wgx*ttzz(h,k))                                         
            atxy(i,k)=atxy(i,k)+ftn*((1-wgx)*ttxy(l,k)+                 &
     &           wgx*ttxy(h,k))                                         
            auw(i,k)=auw(i,k)+ftn*((1-wgx)*tuw(l,k)+                    &
     &           wgx*tuw(h,k))                                          
            avw(i,k)=avw(i,k)+ftn*((1-wgx)*tvw(l,k)+                    &
     &           wgx*tvw(h,k))                                          
            auv(i,k)=auv(i,k)+ftn*((1-wgx)*tuv(l,k)+                    &
     &           wgx*tuv(h,k))                                          
            adudz(i,k)=adudz(i,k)+ftn*((1-wgx)*tdudz(l,k)+              &
     &           wgx*tdudz(h,k))                                        
            adudx(i,k)=adudx(i,k)+ftn*((1-wgx)*tdudx(l,k)+              &
     &           wgx*tdudx(h,k))                                        
            advdz(i,k)=advdz(i,k)+ftn*((1-wgx)*tdvdz(l,k)+              &
     &           wgx*tdvdz(h,k))                                        
            adwdz(i,k)=adwdz(i,k)+ftn*((1-wgx)*tdwdz(l,k)+              &
     &           wgx*tdwdz(h,k))                                        
            adwdx(i,k)=adwdx(i,k)+ftn*((1-wgx)*tdwdx(l,k)+              &
     &           wgx*tdwdx(h,k))                                        
            aCs2(i,k)=aCs2(i,k)+ftn*((1-wgx)*tCs2(l,k)+                 &
     &           wgx*tCs2(h,k))                                         
            aCs(i,k)=aCs(i,k)+ftn*((1-wgx)*tCs(l,k)+                    &
     &           wgx*tCs(h,k))                                          
            abeta1(i,k)=abeta1(i,k)+ftn*((1-wgx)*tbeta1(l,k)+           &
     &           wgx*tbeta1(h,k))                                       
            aESGS(i,k)=aESGS(i,k)+ftn*((1-wgx)*tESGS(l,k)+              &
     &           wgx*tESGS(h,k))                                        
            aDSGS(i,k)=aDSGS(i,k)+ftn*((1-wgx)*tDSGS(l,k)+              &
     &           wgx*tDSGS(h,k))                                        
            aTL(i,k)=aTL(i,k)+ftn*((1-wgx)*tTL(l,k)+                    &
     &           wgx*tTL(h,k))                                          
            e(i,k)=e(i,k)+ftn*((1-wgx)*te(l,k)+wgx*te(h,k)) 
                                                                        
         enddo 
      enddo 
                                                                        
      if(vfact.eq.0 .and. verticalBC.eq.0)then 
         do j=1,nyb 
            do i=1,nx 
               atxz_s(i,j)=atxz_s(i,j)+fr*                              &
     &              (txz(i,j,2)**2+tyz(i,j,2)**2)**0.25                 
            enddo 
         enddo 
                                                                        
      endif 
                                                                        
      return 
      END                                           
