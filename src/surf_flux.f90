    subroutine surf_flux(theta,thetaFactor,u,v,t_flux,Psi,Psi0,fi, &
    fiH,zo,t_s,coolrate,scalarIndex)

    use globals
    use scalars
    use wallBoundaryConditions
    implicit none

    interface
    include './interfaces/filter_2dsl.f90'
    include './interfaces/plane_reduce.f90'
    include './interfaces/plane_avg.f90'
    end interface

    integer*4 :: scalarIndex
    integer*4 :: i,j,k,kk,ll,pos_cnt,neg_cnt,ddd
    real*8,dimension(:,:,:)::theta,thetaFactor,u,v,Psi,Psi0,fi,fiH
    real*8,dimension(:,:)::t_flux,zo,t_s,coolrate

    real*8, dimension(size(theta,1),size(theta,2)) :: ustar,t_sn,OB_L, &
    u_res,t_res,u_hat,v_hat,t_hat,clrtn,zon
    real*8, dimension(size(theta,1),size(theta,2),2) :: PsiH,PsiH0
    real*8 :: t_res_sum,u_res_sum,denom,t_flux_avg,z,x,y,denomH, &
    OBL_avg,ustar_avg,t_bar

    if(surfaceFlags(scalarIndex) == 0)then

        t_flux=surfaceFluxes(scalarIndex)*10**(-4.)

    else
             
        ddd = nint(10*dmod((dble(t+nrsub)*dt*Ugal),dx)/dx)
                    
        do i=1,nx
            if(ddd == swtx)then
                if(i == nx)then
                    t_sn(i,:)=t_s(1,:)
                    clrtn(i,:)=coolrate(1,:)
                    zon(i,:)=zo(1,:)
                else
                    t_sn(i,:)=t_s(i+1,:)
                    clrtn(i,:)=coolrate(i+1,:)
                    zon(i,:)=zo(i+1,:)
                endif
            else
                t_sn(i,:)=t_s(i,:)
                clrtn(i,:)=coolrate(i,:)
                zon(i,:)=zo(i,:)
            endif
                                   
        end do
                 
        if(t+nrsub > nint(8.*3600.*u_star/z_i/dt))then
            do j=1,nyb
                do i=1,nx
                    t_s(i,j)=t_sn(i,j)-(dsdtHomogeneous(scalarIndex)/ &
                    scalarScales(scalarIndex)/3600)*(z_i/u_star)*dt
                end do
            end do
        else
            do j=1,nyb
                do i=1,nx
                    t_s(i,j)=t_sn(i,j)-(clrtn(i,j)/ &
                    scalarScales(scalarIndex)/3600)*(z_i/u_star)*dt
                end do
            end do
        endif
                 
        coolrate = clrtn
        zo       = zon

    endif
            
!       Option: Homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!         U_res_sum = 0.
!         t_res_sum = 0.
!         do j=1,Nyb
!            do i=1,Nx
!               U_res_sum = U_res_sum + ((u(i,j,2)+Ugal)**2.+
!     +              v(i,j,2)**2.)**0.5
!               t_res_sum = t_res_sum + theta(i,j,2)
!            end do
!         end do
!         call plane_reduce(U_res_sum)
!         call plane_reduce(t_res_sum)
!         endif
!         do j=1,Nyb
!            do i=1,Nx
!               U_res(i,j)=U_res_sum*iNxNy
!               t_res(i,j)=t_res_sum*iNxNy
!            end do
!         end do
!        t_bar = t_res_sum*iNxNy
!       Option: Heterogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call Filter_2dsl(t_res,theta(:,:,2))
    call Filter_2dsl(u_hat,u(:,:,2))
    call Filter_2dsl(v_hat,v(:,:,2))

!     if t > 1 Filter_2dsl does the following  (?? verify this ??)
!        t_hat=theta(:,:,2)
!        u_hat=u(:,:,2)
!        v_hat=v(:,:,2)
!     when t == 1 AND flag == 1 Filter_2dsl also determines what fftw method to use later

! convert to virtual potential temperature (if theta is temperature) by multiplying by temperature*(1 + 0.61*r)
! do for t_s and t_res before computing the flux ( only if scalar is temperature !!! )
! done below at computation of t_flux (????)
! for scalars other than temperature, thetaFactor = 1
            

    t_bar = 0.d0
    do j=1,Nyb
        do i=1,Nx
            U_res(i,j) = ((u_hat(i,j)+Ugal)**2.+ &
            (v_hat(i,j)+Vgal)**2.)**0.5
        !              	    t_bar      = t_bar+theta(i,j,2)
        end do
    end do

!        call plane_reduce(t_bar)
!        t_bar=t_bar*inxny

!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(t == 1)then
        Psi=0.d0
        Psi0=0.d0
        PsiH0=0.d0
        PsiH=0.d0
    endif

    do kk=1,4

        t_flux_avg=0.d0
        ustar_avg=0.d0
        do j=1,Nyb
            do i=1,Nx
                denom = dlog((dz/2.*1.)/zo(i,j))+Psi(i,j,1) &
                -Psi0(i,j,1)
                ustar(i,j) = (u_res(i,j)*vonk/denom)
                if(surfaceFlags(scalarIndex) == 1)then
                    denomH = dlog((dz/2.*1.)/(0.1*zo(i,j)))+ &
                    PsiH(i,j,1)-PsiH0(i,j,1)
                    t_flux(i,j) = (t_s(i,j)-t_res(i,j))*ustar(i,j)* &
                    vonk/denomH
                endif
                t_flux_avg = t_flux_avg+t_flux(i,j)
                ustar_avg=ustar_avg+ustar(i,j)
            end do
        end do

        call plane_reduce(t_flux_avg)
        call plane_reduce(ustar_avg)
        t_flux_avg=t_flux_avg*inxny
        ustar_avg=ustar_avg*inxny

        pos_cnt=0
        neg_cnt=0
                 
        do j=1,nyb
            do i=1,nx
                OB_L(i,j)=-ustar(i,j)**3.*t_res(i,j)* &
                (1.0/(vonk*g_hat*t_flux(i,j)))

                do k=1,2
                !                  if(scalarFlags(scalarIndex).eq.1)then !if temperature (0 is passive, 2 is moisture)
                    if(scalarFlags(scalarIndex) /= 0)then ! if scalar is not passive
                        z=+k*0.5*dz
                        if(z/OB_L(i,j) > 5.)then
                            OB_L(i,j)=z/5.
                        elseif(z/OB_L(i,j) < -5.)then
                            OB_L(i,j)=-z/5.
                        endif
                        if ((t_flux(i,j)) > 0.) then
                            if(k == 1) pos_cnt=pos_cnt+1
                            x=+(1.-(15.*z/OB_L(i,j)))**0.25
                            y=+(1.-(15.*zo(i,j)/OB_L(i,j)))**0.25
                            if(scalarFlags(scalarIndex) == 1)then !if temperature (0 is passive, 2 is moisture)
                                Psi(i,j,k)=-2.*dlog(0.5*(1+x))- &
                                dlog(0.5*(1+x**2.))+2.*atan(x)-Pi/2.
                                Psi0(i,j,k)=-2.*dlog(0.5*(1+y))- &
                                dlog(0.5*(1+y**2.))+2.*atan(y)-Pi/2.
                                fi(i,j,k)=1./x
                            endif
                        !     Equation 11.9 and 11.14 from Arya
                            psiH(i,j,k)=-2*dlog(0.5*(1+x**2.))
                            psiH0(i,j,k)=-2*dlog(0.5*(1+y**2.))
                            fiH(i,j,k)=fi(i,j,k)**2.0
                        else if ((t_flux(i,j)) < (-0.)) then
                            if(k == 1) neg_cnt=neg_cnt+1
                            if(scalarFlags(scalarIndex) == 1)then !if temperature (0 is passive, 2 is moisture)
                                Psi(i,j,k)=+4.8*z/OB_L(i,j)
                                Psi0(i,j,k)=+4.8*zo(i,j)/OB_L(i,j)
                                fi(i,j,k)=1.+psi(i,j,k)
                            endif
                            PsiH(i,j,k)=7.8*z/OB_L(i,j)
                            PsiH0(i,j,k)=7.8*zo(i,j)/OB_L(i,j)
                            fiH(i,j,k)=1.+psiH(i,j,k)
                        else
                            if(scalarFlags(scalarIndex) == 1)then !if temperature (0 is passive, 2 is moisture)
                                Psi(i,j,k)=0.
                                Psi0(i,j,k)=0.
                                fi(i,j,k)=1.
                            endif
                            PsiH(i,j,k) = Psi(i,j,k)
                            PsiH0(i,j,k)= Psi0(i,j,k)
                            fiH(i,j,k) = fi(i,j,k)
                        end if
                                             
                    else ! if scalar is passive

                        Psi(i,j,k)=0.
                        Psi0(i,j,k)=0.
                        fi(i,j,k)=1.
                        PsiH(i,j,k)=0.
                        PsiH0(i,j,k)=0.
                        fiH(i,j,k)=0.74

                    endif

                end do
                              
            enddo
        enddo
    enddo

!      OBL_avg=0.d0
!      do j=1,Nyb
!         do i=1,Nx
!            OBL_avg=OBL_avg+OB_L(i,j)
!         enddo
!      enddo
!      call plane_reduce(OBL_avg)
!      write(1525,*) OBL_avg*z_i*inxny,ustar_avg
!      write (1999,*) -t_flux_avg/ustar_avg,pos_cnt,neg_cnt
!      write(*,*) 'tflux,u*',t_flux_avg,ustar_avg
    	        

    return

    end subroutine surf_flux
