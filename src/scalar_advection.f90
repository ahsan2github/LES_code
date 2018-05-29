    subroutine scalar_advection(scalarRHS)
    use globals
    use scalars
    use wallBoundaryConditions
    implicit none

    real*8,dimension(:,:,:,:) :: scalarRHS

    integer*4 :: k,l
    real*8 :: adv,ztemp
    real*8,dimension(size(scalarRHS,3)) :: scalarAdv

    do l=1,scalarCount
        if(S_advec(l) == 1)then
            UTC_hrs = startUTC+(dble(ttt)*dt*z_i/u_star)/3600.d0
            if( l==moistureIndex)then
                if(UTC_hrs < 2.d0)then
                    adv = 0.0d0 &
                    /(u_star*scalarScales(moistureIndex)/z_i)
                elseif(UTC_hrs < 5.d0)then
                    adv = -0.00000008d0 &
                    /(u_star*scalarScales(moistureIndex)/z_i)
                else
                    adv = 0.00000000d0 &
                    /(u_star*scalarScales(moistureIndex)/z_i)
                endif
            elseif( l == temperatureIndex )then
                if(UTC_hrs < 1.d0)then
                    adv = -0.000025d0 &
                    /(u_star*scalarScales(temperatureIndex)/z_i)
                elseif(UTC_hrs < 6.d0)then
                    adv = 0.000075d0 &
                    /(u_star*scalarScales(temperatureIndex)/z_i)
                else
                    adv = 0.00000000d0 &
                    /(u_star*scalarScales(temperatureIndex)/z_i)
                endif
            endif
                        
            scalarAdv = 0
            do k = 2,nzb+1
                ztemp=(me*nzb+k-1)*dz*z_i-dz*z_i/2.0d0
                if(ztemp >= 200)then
                    scalarAdv(k)=adv
                else
                    scalarAdv(k)=adv*ztemp/200.d0
                endif
                scalarRHS(:,:,k,l)=scalarRHS(:,:,k,l)+scalarAdv(k)
            enddo
        endif

    enddo

    return
    end subroutine scalar_advection
