    subroutine root8(rtnewt,aa,bb,cc,dd,ee,ff)

    INTEGER :: JMAX,j
    REAL*8 :: rtnewt,xinit,xacc,b,aa,bb,cc,dd,ee,ff,df,dx,f
    REAL*8 :: ccp,ddp,eep,ffp
    PARAMETER (JMAX=100)
    PARAMETER(xinit=100.0,xacc=0.0001)

    rtnewt=xinit

    ccp = cc*2.
    ddp = dd*3.
    eep = ee*4
    ffp = ff*5.

    do 11 j=1,JMAX
        b=rtnewt

        f=aa+bb*b+cc*b**2.+dd*b**3.+ee*b**4.+ff*b**5.
    !        f=aa+bb*b+cc*b*b+dd*b*b*b+ee*b*b*b*b+ff*b*b*b*b*b

        df=bb+ccp*b+ddp*b**2.+eep*b**3.+ffp*b**4.
    !        df=bb+ccp*b+ddp*b*b+eep*b*b*b+ffp*b*b*b*b

        dx=f/df

        rtnewt=rtnewt-dx

        if(abs(dx) < xacc) return
    11 END DO

    222	rtnewt=0.

    end subroutine root8
