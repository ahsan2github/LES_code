    Subroutine openfiles (dsdxUnits,dsdyUnits,dsdzUnits, &
    scalarMeanUnits,scalar2Units,scalar3Units,xScalarFluxUnits, &
    yScalarFluxUnits,zScalarFluxUnits,awsUnits,ausUnits,avsUnits, &
    scalarSpectraUnits,scalarPrandtlUnits,cs2prUnits,beta2Units, &
    ETUnits,obukovUnits,scalarStarUnits,surfaceFluxUnits)

    use globals
    use scalars
    use mainmodule

    integer*4,dimension(:):: dsdxUnits, &
    dsdyUnits, dsdzUnits,scalarMeanUnits,scalar2Units,scalar3Units, &
    xScalarFluxUnits,yScalarFluxUnits,zScalarFluxUnits,awsUnits, &
    ausUnits,avsUnits,scalarSpectraUnits,scalarPrandtlUnits, &
    cs2prUnits,beta2Units,ETUnits,obukovUnits,scalarStarUnits, &
    surfaceFluxUnits
    integer*4 :: scalarIter
    character(11) :: ff
    character(6) :: pp
    character(7) :: ss
    character(100) :: fileStr
    character(15) :: scalarName
    integer*4 :: nameLength,l
! ...  the following is an example of what fileStr is for scalarN for each file
!      fileStr = 'output/scalar3/prandtl.bin'
    ff='unformatted'

    if(initu == 0 .AND. inits == 0)then
        pp='asis'
        ss='replace'
    else
        pp='append'
        ss='unknown'
    end if

    scalarIter = 0
!...  Open output files
    Open (unit=11,file='output/ke.out',status=ss,position=pp)
    Open (unit=39,file='output/tke.bin',status=ss,form=ff,position=pp)
    Open (unit=71,file='output/au.bin',status=ss,form=ff,position=pp)
    Open (unit=72,file='output/av.bin',status=ss,form=ff,position=pp)
    Open (unit=73,file='output/aw.bin',status=ss,form=ff,position=pp)
    Open (unit=74,file='output/u2.bin',status=ss,form=ff,position=pp)
    Open (unit=75,file='output/v2.bin',status=ss,form=ff,position=pp)
    Open (unit=76,file='output/w2.bin',status=ss,form=ff,position=pp)
    Open (unit=40,file='output/w3.bin',status=ss,form=ff,position=pp)
    Open (unit=57,file='output/atxx.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=78,file='output/atxz.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=79,file='output/atyy.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=80,file='output/atyz.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=87,file='output/atzz.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=90,file='output/atxy.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=81,file='output/p2.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=82,file='output/auw.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=83,file='output/avw.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=77,file='output/auv.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=85,file='output/ap.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=88,file='output/dudz.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=91,file='output/dudx.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=95,file='output/dvdz.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=93,file='output/dwdz.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=92,file='output/dwdx.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=299,file='output/spectru.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=298,file='output/spectrv.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=297,file='output/spectrw.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=295,file='output/spectrp.bin',status=ss,form=ff, &
    position=pp)
    Open(unit=398,file='output/Cs2_ALL.bin',status=ss,form=ff, &
    position=pp)
    Open(unit=399,file='output/Cs_ALL.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=978,file='output/beta1.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=888,file='output/ESGS.bin',status=ss,form=ff, &
    position=pp)

    open(unit=132,file='output/ustar.bin',status=ss,form=ff, &
    position=pp)

    Open (unit=320,file='output/TL.bin',status=ss,form=ff, &
    position=pp)
    Open (unit=321,file='output/dsgs.bin',status=ss,form=ff, &
    position=pp)

    If (scalarCount >= 1) then
        do l=1,scalarCount
            if (scalarFlags(l) == 1)then
                write(scalarName,'(A)') 'temperature'
            elseif(scalarFlags(l) == 2)then
                write(scalarName,'(A)') 'moisture'
            else
                scalarIter = scalarIter + 1
                write(scalarName,'(A,I1)') 'scalar',scalarIter
            endif
        ! system() may be compiler dependent and may not work on other systems
        !  however it may not be needed for some compilers
        !  most will create the directory if it doesn't exist when open is called
            call system("mkdir -p ./output/"//TRIM(scalarName))
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/dsdz.bin'
            Open (unit=dsdzUnits(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! old unit = 89
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/dsdx.bin'
            Open (unit=dsdxUnits(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! 114
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/dsdy.bin'
            Open (unit=dsdyUnits(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! 115
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/scalarMean.bin'
            Open (unit=scalarMeanUnits(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! 100
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/t2.bin'
            Open (unit=scalar2Units(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! 102
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/t3.bin'
            Open (unit=scalar3Units(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! 41
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/flux_t1.bin'
            Open (unit=xScalarFluxUnits(l),file=TRIM(fileStr), &
            status=ss,form=ff,position=pp)!112
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/flux_t2.bin'
            Open (unit=yScalarFluxUnits(l),file=TRIM(fileStr), &
            status=ss,form=ff,position=pp)!116
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/flux_t3.bin'
            Open (unit=zScalarFluxUnits(l),file=TRIM(fileStr), &
            status=ss,form=ff,position=pp)!104
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/awt.bin'
            Open(unit=awsUnits(l),file=TRIM(fileStr),status=ss,form=ff, &
            position=pp)! 106
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/aut.bin'
            Open(unit=ausUnits(l),file=TRIM(fileStr),status=ss,form=ff, &
            position=pp)! 113
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/avt.bin'
            Open(unit=avsUnits(l),file=TRIM(fileStr),status=ss,form=ff, &
            position=pp)! 117
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/spectrt.bin'
            open (unit=scalarSpectraUnits(l),file=TRIM(fileStr), &
            status=ss,form=ff,position=pp)! 300
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/prandtl.bin'
            open(unit=scalarPrandtlUnits(l),file=TRIM(fileStr), &
            status=ss,form=ff,position=pp)! 1526
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/cs2pr.bin'
            open(unit=cs2prUnits(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)! 1527
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/beta2.bin'
            Open (unit=beta2Units(l),file=TRIM(fileStr),status=ss, &
            form=ff,position=pp)!979
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/ET.bin'
            Open (unit=ETUnits(l),file=TRIM(fileStr),status=ss,form=ff, &
            position=pp)!999

            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/obukov.out'
            open(unit=obukovUnits(l),file=TRIM(fileStr), &
            status=ss,position=pp)!1525
            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/t_star.out'
            open(unit=scalarStarUnits(l),file=TRIM(fileStr), &
            status=ss,position=pp)!1999

            write(fileStr, '(A,A,A)' ) 'output/',TRIM(scalarName), &
            '/qz_surf.bin'
            open(unit=surfaceFluxUnits(l),file=TRIM(fileStr), &
            status=ss,form=ff,position=pp)!133
        enddo
    end if

    end subroutine openfiles

    Subroutine closefiles (dsdxUnits,dsdyUnits,dsdzUnits, &
    scalarMeanUnits,scalar2Units,scalar3Units,xScalarFluxUnits, &
    yScalarFluxUnits,zScalarFluxUnits,awsUnits,ausUnits,avsUnits, &
    scalarSpectraUnits,scalarPrandtlUnits,cs2prUnits,beta2Units, &
    ETUnits,obukovUnits,scalarStarUnits,surfaceFluxUnits)
    use globals
    use scalars
    integer*4,dimension(:):: dsdxUnits, &
    dsdyUnits,dsdzUnits,scalarMeanUnits,scalar2Units,scalar3Units, &
    xScalarFluxUnits,yScalarFluxUnits,zScalarFluxUnits,awsUnits, &
    ausUnits,avsUnits,scalarSpectraUnits,scalarPrandtlUnits, &
    cs2prUnits,beta2Units,ETUnits,obukovUnits,scalarStarUnits, &
    surfaceFluxUnits

    Close(unit=11)
    Close(unit=39)
    Close(unit=71)
    Close(unit=72)
    Close(unit=73)
    Close(unit=74)
    Close(unit=77)
    Close(unit=75)
    Close(unit=76)
    Close(unit=40)
    Close(unit=57)
    Close(unit=78)
    Close(unit=79)
    Close(unit=80)
    Close(unit=87)
    Close(unit=81)
    Close(unit=82)
    Close(unit=83)
    Close(unit=85)
    Close(unit=88)
    Close(unit=90)
    Close(unit=91)
    Close(unit=95)
    Close(unit=92)
    Close(unit=93)
    Close(unit=295)
    Close(unit=297)
    Close(unit=298)
    Close(unit=299)
    Close(unit=398)
    Close(unit=399)
    Close(unit=978)
    Close(unit=888)
    Close(unit=132)
    Close(unit=320)
    Close(unit=321)

    If (scalarCount >= 1) then
        do l=1,scalarCount
            Close (unit=dsdzUnits(l))
            Close (unit=scalarMeanUnits(l))
            Close (unit=scalar2Units(l))
            Close (unit=scalar3Units(l))
            Close (unit=zScalarFluxUnits(l))
            Close (unit=awsUnits(l))
            Close (unit=xScalarFluxUnits(l))
            Close (unit=ausUnits(l))
            Close (unit=dsdxUnits(l))
            Close (unit=dsdyUnits(l))
            Close (unit=yScalarFluxUnits(l))
            Close (unit=avsUnits(l))
            Close (unit=scalarPrandtlUnits(l))
            Close (unit=cs2prUnits(l))
            Close (unit=obukovUnits(l))
            Close (unit=scalarStarUnits(l))
            Close (unit=scalarSpectraUnits(l))
            Close (unit=beta2Units(l))
            Close (unit=ETUnits(l))
            Close (unit=surfaceFluxUnits(l))
        enddo
    end if

    end subroutine closefiles
