    subroutine write_frames(u,v,w,p,esgs,dsgs,TL,lambda2,particle, &
    scalar1)
! ossible frames are:
!1. u-velocity (u_frame)
!2. v-velocity (v_frame)
!3. w-velocity (w_frame)
!4. pressure (p_frame)
!5. SGS tke (esgs_frame)
!6. SGS dissipation (dsgs_frame)
!7. SGS Lagrangian timescale (TL_frame)
!8. lambda2 (lambda2_frame)
!9. particle positions (part_frame)
! 0. scalar1 (scalar1_frame)

    use globals
    use mainModule
    use particleModule
    use frameModule
    implicit none

    real*8,dimension(:,:,:) :: u,v,w,p,esgs,dsgs,TL,lambda2,particle, &
    scalar1

    integer*4 :: j, NNz, fh
    real*8,dimension(size(p,1),size(p,2),size(p,3)) :: p_dyn
    character(50) :: frameStr
    integer (kind=MPI_OFFSET_KIND) :: view_disp
    logical :: dir_ext, writeBovHeader

    if( t == 1 )then

        output_frames(:)= .FALSE. 

    ! check to see which frames we'll be outputting
                 
        do j=1,framefiles
                     
        ! _frame
            if(framenames(j)=='u_frame')then
                output_frames(1)= .TRUE. 
            ! _frame
            elseif(framenames(j)=='v_frame')then
                output_frames(2)= .TRUE. 
            ! _frame
            elseif(framenames(j)=='w_frame')then
                output_frames(3)= .TRUE. 
            ! _frame
            elseif(framenames(j)=='p_frame')then
                output_frames(4)= .TRUE. 
            ! sgs_frame
            elseif(framenames(j)=='esgs_frame')then
                output_frames(5)= .TRUE. 
            ! sgs_frame
            elseif(framenames(j)=='dsgs_frame')then
                output_frames(6)= .TRUE. 
            ! L_frame
            elseif(framenames(j)=='TL_frame')then
                output_frames(7)= .TRUE. 
            ! ambda2_frame
            elseif(framenames(j)=='lambda2_frame')then
                output_frames(8)= .TRUE. 
            ! art_frame
            elseif(framenames(j)=='part_frame')then
                output_frames(9)= .TRUE. 
            ! calar1_frame
            elseif(framenames(j)=='scalar1_frame')then
                output_frames(10)= .TRUE. 
            end if

        end do

    ! ake directory structure

        if(me == 0)then

        ! _frame
            if(output_frames(1))then
                inquire(file="./output/u_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/u_frame'
                    call system('mkdir output/u_frame')
                endif
            end if
        ! _frame
            if(output_frames(2))then
                inquire(file="./output/v_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/v_frame'
                    call system('mkdir output/v_frame')
                endif
            end if
        ! _frame
            if(output_frames(3))then
                inquire(file="./output/w_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/w_frame'
                    call system('mkdir output/w_frame')
                endif
            end if
        ! _frame
            if(output_frames(4))then
                inquire(file="./output/p_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/p_frame'
                    call system('mkdir output/p_frame')
                endif
            end if
        ! sgs_frame
            if(output_frames(5))then
                inquire(file="./output/esgs_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/esgs_frame'
                    call system('mkdir output/esgs_frame')
                endif
            end if
        ! sgs_frame
            if(output_frames(6))then
                inquire(file="./output/dsgs_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/dsgs_frame'
                    call system('mkdir output/dsgs_frame')
                endif
            end if
        ! L_frame
            if(output_frames(7))then
                inquire(file="./output/TL_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/TL_frame'
                    call system('mkdir output/TL_frame')
                endif
            end if
        ! ambda2_frame
            if(output_frames(8))then
                inquire(file="./output/lambda2_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/lambda2_frame'
                    call system('mkdir output/lambda2_frame')
                endif
            end if
        ! art_frame
            if(output_frames(9))then
                inquire(file="./output/part_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/part_frame'
                    call system('mkdir output/part_frame')
                endif
            end if
        ! calar1_frame
            if(output_frames(10))then
                inquire(file="./output/scalar1_frame/.",exist=dir_ext)
                if( .NOT. dir_ext)then
                    write(*,*)'making output/scalar1_frame'
                    call system('mkdir output/scalar1_frame')
                endif
            end if

        end if

    end if

    call MPI_BARRIER(nall,ierr)

! write frames

    if(ttt >= sframe .AND. ttt <= sframe+(nframe-1)*framestep)then

        if(ttt == sframe .OR. mod(ttt-sframe,framestep) == 0)then

            if((vfact+1)*nzb <= frameh)then

                if(frameh > vfact*nzb .AND. frameh <= (vfact+1)*nzb)then
                    NNz=frameh-vfact*nzb
                else
                    NNz=nzb
                endif

                view_disp=0

            ! u_frame
                if(output_frames(1))then
                    write(frameStr, '(A,I4.4,A)')'output/u_frame/u_frame', &
                    frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,u(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! v_frame
                if(output_frames(2))then
                    write(frameStr, '(A,I4.4,A)')'output/v_frame/v_frame', &
                    frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,v(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! w_frame
                if(output_frames(3))then
                    write(frameStr, '(A,I4.4,A)')'output/w_frame/w_frame', &
                    frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,w(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! p_frame
                if(output_frames(4))then
                    p_dyn=p-0.5d0*(u**2+v**2+w**2)
                    write(frameStr, '(A,I4.4,A)')'output/p_frame/p_frame', &
                    frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,p(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! esgs_frame
                if(output_frames(5))then
                    write(frameStr, '(A,I4.4,A)') &
                    'output/esgs_frame/esgs_frame',frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,esgs(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! dsgs_frame
                if(output_frames(6))then
                    write(frameStr, '(A,I4.4,A)') &
                    'output/dsgs_frame/dsgs_frame',frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,dsgs(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! TL_frame
                if(output_frames(7))then
                    write(frameStr, '(A,I4.4,A)') &
                    'output/TL_frame/TL_frame',frame_cnt+1,'.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,TL(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! lambda2_frame
                if(output_frames(8))then
                    write(frameStr, '(A,I4.4,A)') &
                    'output/lambda2_frame/lambda2_frame',frame_cnt+1, &
                    '.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,lambda2(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if
            ! scalar1_frame
                if(output_frames(10))then
                    write(frameStr, '(A,I4.4,A)') &
                    'output/scalar1_frame/scalar1_frame',frame_cnt+1, &
                    '.bin'
                    call MPI_FILE_OPEN(MPI_COMM_FRAMES,frameStr, &
                    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL, &
                    fh,ierr)
                    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
                    fileview_3D,"native",MPI_INFO_NULL,ierr)
                    call MPI_FILE_WRITE_ALL(fh,scalar1(1,1,2), &
                    Nx*Nyb*NNz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, &
                    ierr)
                    call MPI_FILE_CLOSE(fh,ierr)
                end if

            end if

        ! part_frame
            if(me == 0)then
                if(output_frames(9))then
                    write(frameStr, '(A,I4.4,A)') &
                    'output/part_frame/part_frame_nomodel', &
                    frame_cnt+1,'.bin'

                    open(unit=890,file=frameStr,status='replace', &
                    form='unformatted')
                    write(890)particle(1:ipart,1:3,1)
                    close(890)
                                   
                    if(part_model >= 2)then
                        write(frameStr, '(A,I4.4,A)') &
                        'output/part_frame/part_frame_weiliso', &
                        frame_cnt+1,'.bin'
                        open(unit=890,file=frameStr,status='replace', &
                        form='unformatted')
                        write(890)particle(1:ipart,1:3,2)
                        close(890)
                    endif
                end if
            end if

            frame_cnt=frame_cnt+1
!             dir_ext = writeBovHeader()

        end if

    end if

    end subroutine write_frames
