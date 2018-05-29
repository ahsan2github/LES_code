    subroutine newroots(rtnewt,aa,bb,cc,dd,ee,ff)
    INTEGER :: m
    PARAMETER (m=5)
    REAL*8 :: rtnewt,a(m+1),aa,bb,cc,dd,ee,ff,rtr(m),rti(m)
    	
    a(1) = aa
    a(2) = bb
    a(3) = cc
    a(4) = dd
    a(5) = ee
    a(6) = ff
    if (aa == 0 .AND. bb == 0 .AND. cc == 0 .AND. dd == 0. &
       .and.ee == 0 .AND. ff == 0) then
        rtnewt = -1
    else
        CALL zrhqr(a,m,rtr,rti)
        rtnewt = 0
        do i=1,m
            if (abs(rti(i)) <= 1e-10 .AND. rtr(i) > rtnewt) then
                rtnewt = rtr(i)
            end if
        end do
    end if
    	
    return
    end subroutine newroots
    	
    SUBROUTINE zrhqr(a,m,rtr,rti)
    INTEGER :: m
    PARAMETER (MAXM=10)
    REAL*8 :: a(m+1),rtr(m),rti(m)
    INTEGER :: j,k
    REAL*8 :: hess(MAXM,MAXM),xr,xi
    do k=1,m
        hess(1,k)=-a(m+1-k)/a(m+1)
        do j=2,m
            hess(j,k)=0.
        end do
        if (k /= m) hess(k+1,k)=1.
    end do

    call balanc(hess,m,MAXM)
    call hqr(hess,m,MAXM,rtr,rti)
    do j=2,m
        xr=rtr(j)
        xi=rti(j)
        do k=j-1,1,-1
            if(rtr(k) <= xr)goto 1
            rtr(k+1)=rtr(k)
            rti(k+1)=rti(k)
        end do
        k=0
        1 rtr(k+1)=xr
        rti(k+1)=xi
    end do
    return
    end SUBROUTINE zrhqr
    	
    SUBROUTINE hqr(a,n,np,wr,wi)
    INTEGER :: n,np
    REAL*8 :: a(np,np),wi(np),wr(np)
    INTEGER :: i,its,j,k,l,m,nn
    REAL*8 :: anorm,p,q,r,s,t,u,v,w,x,y,z
    anorm=0.
    do i=1,n
        do j=max(i-1,1),n
            anorm=anorm+abs(a(i,j))
        enddo
    enddo
    nn=n
    t=0.
    1	if(nn >= 1)then
    its=0
    2 do l=nn,2,-1
        s=abs(a(l-1,l-1))+abs(a(l,l))
        if(s == 0.)s=anorm
        if(abs(a(l,l-1))+s == s) then
            a(l,l-1)=0.
            goto 3
        end if
    end do
    l=1
    3 x=a(nn,nn)
    if(l == nn)then
        wr(nn)=x+t
        wi(nn)=0.
        nn=nn-1
        	      
    else
        y=a(nn-1,nn-1)
        w=a(nn,nn-1)*a(nn-1,nn)
        if(l == nn-1)then
            p=0.5*(y-x)
            q=p**2+w
            z=sqrt(abs(q))
            x=x+t
            if(q >= 0.)then
                z=p+sign(z,p)
                wr(nn)=x+z
                wr(nn-1)=wr(nn)
                if(z /= 0.)wr(nn)=x-w/z
                wi(nn)=0.
                wi(nn-1)=0.
            else
                wr(nn)=x+p
                wr(nn-1)=wr(nn)
                wi(nn)=z
                wi(nn-1)=-z
            endif
            nn=nn-2
        else
        !		 if(its.eq.30) pause 'too many iterations in hqr'
            if(its == 30)then
                do ii=1,np
                    wr(ii)=0
                    wi(ii)=0
                enddo
                return
            endif
            if(its == 10 .OR. its == 20)then
                t=t+x
                do i=1,nn
                    a(i,i)=a(i,i)-x
                enddo
                s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
                x=0.75*s
                y=x
                w=-0.4375*s**2
            endif
            its=its+1
            do m=nn-2,l,-1
                z=a(m,m)
                r=x-z
                s=y-z
                p=(r*s-w)/a(m+1,m)+a(m,m+1)
                q=a(m+1,m+1)-z-r-s
                r=a(m+2,m+1)
                s=abs(p)+abs(q)+abs(r)
                p=p/s
                q=q/s
                r=r/s
                if(m == l)goto 4
                u=abs(a(m,m-1))*(abs(q)+abs(r))
                v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
                if(u+v == v) goto 4
            enddo
            4 do i=m+2,nn
                a(i,i-2)=0.
                if (i /= m+2) a(i,i-3)=0.
            enddo
            do k=m,nn-1
                if(k /= m)then
                    p=a(k,k-1)
                    q=a(k+1,k-1)
                    r=0.
                    if(k /= nn-1)r=a(k+2,k-1)
                    x=abs(p)+abs(q)+abs(r)
                    if(x /= 0.)then
                        p=p/x
                        q=q/x
                        r=r/x
                    endif
                endif
                s=sign(sqrt(p**2+q**2+r**2),p)
                if(s /= 0.)then
                    if(k == m)then
                        if(l /= m)a(k,k-1)=-a(k,k-1)
                    else
                        a(k,k-1)=-s*x
                    endif
                    p=p+s
                    x=p/s
                    y=q/s
                    z=r/s
                    q=q/p
                    r=r/p
                    do j=k,nn
                        p=a(k,j)+q*a(k+1,j)
                        if(k /= nn-1)then
                            p=p+r*a(k+2,j)
                            a(k+2,j)=a(k+2,j)-p*z
                        endif
                        a(k+1,j)=a(k+1,j)-p*y
                        a(k,j)=a(k,j)-p*x
                    enddo
                    do i=l,min(nn,k+3)
                        p=x*a(i,k)+y*a(i,k+1)
                        if(k /= nn-1)then
                            p=p+z*a(i,k+2)
                            a(i,k+2)=a(i,k+2)-p*r
                        endif
                        a(i,k+1)=a(i,k+1)-p*q
                        a(i,k)=a(i,k)-p
                    enddo
                endif
            enddo
            goto 2
        endif
    endif
    goto 1
endif
    return
    end SUBROUTINE hqr

    SUBROUTINE balanc(a,n,np)
    INTEGER :: n,np
    REAL*8 :: a(np,np),RADIX,SQRDX
    PARAMETER (RADIX=2.,SQRDX=RADIX**2)
    INTEGER :: i,j,last
    REAL*8 :: c,f,g,r,s
    1	continue
    last=1
    do i=1,n
        c=0.
        r=0.
        do j=1,n
            if(j /= i)then
                c=c+abs(a(j,i))
                r=r+abs(a(i,j))
            endif
        enddo
        if(c /= 0. .AND. r /= 0.)then
            g=r/RADIX
            f=1.
            s=c+r
            2 if(c < g)then
                f=f*RADIX
                c=c*SQRDX
                goto 2
            endif
            g=r*RADIX
            3 if(c > g)then
                f=f/RADIX
                c=c/SQRDX
                goto 3
            endif
            if((c+r)/f < 0.95*s)then
                last=0
                g=1./f
                do j=1,n
                    a(i,j)=a(i,j)*g
                enddo
                do j=1,n
                    a(j,i)=a(j,i)*f
                enddo
            endif
        endif
    enddo
    if(last == 0)goto 1
    return
    end SUBROUTINE balanc
