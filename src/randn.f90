    FUNCTION randn(idum)
    use globals
    implicit none
    integer :: idum
    real :: randn

    integer :: iset
    real :: ran1,fac,gset,rsq,v1,v2,rand
    save iset,gset
    data iset/0/

    if(idum < 0) iset=0
    if(iset == 0) then
        1 v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq >= 1. .OR. rsq == 0)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        randn=v2*fac
        iset=1
    else
        randn=gset
        iset=0
    endif
    end FUNCTION randn


    FUNCTION ran1(idum)
    implicit none
    integer :: idum,IA,IM,IQ,IR,NTAB,NDIV
    real :: ran1,AM,EPS,RNMX
    parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
    NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    integer :: j,k,iv(NTAB),iy
    save iv,iy
    data iv /NTAB*0/, iy /0/

    if (idum <= 0 .OR. iy == 0)then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if(idum < 0) idum=idum+IM
            if(j <= NTAB) iv(j)=idum
        enddo
        iy=iv(1)
    endif
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if(idum < 0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
    ran1=min(AM*iy,RNMX)
    return
    end FUNCTION ran1



