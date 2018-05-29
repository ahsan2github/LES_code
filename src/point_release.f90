    subroutine point_release(particle,release_pos,zo,idum)

    use globals
    use particleModule
    use canopyModule
    use wallBoundaryConditions
    implicit none

    logical :: dir_ext
    integer :: idum
    integer*4 :: i,j,f,ipart0
    real :: randn,ran1
    real*8 :: ranx,rany,ranz
    real*8,dimension(:,:) :: release_pos,zo
    real*8,dimension(:,:,:) :: particle
          
    if(me == 0)then
        if(ttt == start_release)then
            inquire(file="./output/part_frame/.",exist=dir_ext)
            if( .NOT. dir_ext)then
                write(*,*)'making output/part_frame'
                call system('mkdir output/part_frame')
            endif
            inquire(file="./output/part_frame/attribs_nomodel.bin", &
            exist=dir_ext)
            if(dir_ext)then
                write(*,*)'removing old attribs.bin files'
                call system('rm output/part_frame/attribs_nomodel.bin')
                if(part_model >= 2)then
                    call &
                    system('rm output/part_frame/attribs_weiliso.bin')
                endif
                if(part_model == 3)then
                    call &
                    system('rm output/part_frame/attribs_weilaniso.bin')
                endif
            endif
        endif

        open(unit=895,file='output/part_frame/attribs_nomodel.bin' &
        ,form='unformatted',access='stream',position='append')
        open(unit=896,file='output/part_frame/attribs_weiliso.bin' &
        ,form='unformatted',access='stream',position='append')
        open(unit=897,file='output/part_frame/attribs_weilaniso.bin' &
        ,form='unformatted',access='stream',position='append')
                 
    end if

    ipart0=ipart+1

    do i=1,freq
        do j=1,nr
            ipart = ipart + 1

            if(me == 0)then
                ranx=randn(idum)*1.e-3
                rany=randn(idum)*1.e-3
                ranz=randn(idum)*1.e-3
            end if
            call MPI_BARRIER(nall,ierr)
            call MPI_BCAST(ranx,1,MPI_DOUBLE_PRECISION,0,nall,ierr)
            call MPI_BCAST(rany,1,MPI_DOUBLE_PRECISION,0,nall,ierr)
            call MPI_BCAST(ranz,1,MPI_DOUBLE_PRECISION,0,nall,ierr)

            particle(ipart,1:3,1) = release_pos(j,:)+(/ranx,rany,ranz/)
            particle(ipart,4,1) = 0
            particle(ipart,5,1) = j
            if(deposition == 1)then
                particle(ipart,6,1) = 1
                write(895)j,1
            else
                particle(ipart,6,1) = 0
                write(895)j,0
            endif
            if(particle(ipart,3,1) <= zo(1,1))then
                particle(ipart,3,1)=zo(1,1)+abs(ranz)
            elseif(particle(ipart,3,1) >= h_canopy .AND. c_flag == 1)then
                particle(ipart,3,1)=h_canopy-abs(ranz)
            endif

            if(part_model >= 2)then
                particle(ipart,1:3,2) = release_pos(j,:)+ &
                (/ranx,rany,ranz/)
                particle(ipart,4,2) = 0
                particle(ipart,5,2) = j
                if(deposition == 1)then
                    particle(ipart,6,2) = 1
                    write(896)j,1
                else
                    particle(ipart,6,2) = 0
                    write(896)j,0
                endif
                if(particle(ipart,3,2) <= zo(1,1))then
                    particle(ipart,3,2)=zo(1,1)+abs(ranz)
                elseif(particle(ipart,3,2) >= h_canopy .AND. c_flag == 1) &
                    then
                    particle(ipart,3,2)=h_canopy-abs(ranz)
                endif

            endif

            if(part_model == 3)then
                particle(ipart,1:3,3) = release_pos(j,:)+ &
                (/ranx,rany,ranz/)
                particle(ipart,4,3) = 0
                particle(ipart,5,3) = j
                if(deposition == 1)then
                    particle(ipart,6,3) = 1
                    write(897)j,1
                else
                    particle(ipart,6,3) = 0
                    write(897)j,0
                endif
                if(particle(ipart,3,3) <= zo(1,1))then
                    particle(ipart,3,3)=zo(1,1)+abs(ranz)
                elseif(particle(ipart,3,3) >= h_canopy .AND. c_flag == 1) &
                    then
                    particle(ipart,3,3)=h_canopy-abs(ranz)
                endif

            endif

            if(ipart == npart)then
                return
            endif

        enddo
    enddo

    if(nodep_copy == 1)then
        particle(ipart+1:2*ipart-ipart0+1,:,:)= &
        particle(ipart0:ipart,:,:)
        particle(ipart+1:2*ipart-ipart0+1,6,:)=0
        ipart=2*ipart-ipart0+1
    endif

    if(me == 0)then
        close(895)
        close(896)
        close(897)
    end if

    return
    end subroutine point_release
