!! This subroutine writes the .bov header files at each time step that enables VisIt and Paraview read
!! the binary frames written by MPI rourtines. For details about BOV header files 
!! refer to the VisIt help file "Getting data into VisIt"
!! call this function right after  write_frames() or within write_frames() after the line frame_cnt = frame_cnt + 1
!! the function returns a logical value, if succesful retruns .TRUE. otherwise .FALSE.
!! written by Ahsan , Oct 16, 2014

function writeBovHeader() result(retVal)
    use globals
    use mainModule
    use particleModule
    use frameModule
    use wallBoundaryConditions
    use visualizationModule
    implicit none
    character(100) :: frameStr, uStr, vStr, wStr, pStr, lambda2Str, scalarStr, endian, &
                     lambda2FiltStr
    logical :: retVal
    integer :: noCurFrame, looper
    real :: dimTime
    dimTime = (ttt * dt) / (u_star/z_i) ! time in SEC
    retVal = .FALSE.
    endian = 'LITTLE'
    if(ttt >= sframe .AND. ttt <= sframe+(nframe-1)*framestep)then
        noCurFrame = int((ttt-sframe)/framestep)+1
        if(me == 0)then

        ! open the header file on output diretory to write velocity vector
            if(output_frames(1) .and. output_frames(2) .and.  output_frames(3)) then
                write(frameStr, '(A,I4.4,A)')'output/u_', &
                    noCurFrame,'.bov'
                write(uStr, '(A,I4.4,A)')'output/u_frame/u_frame', &
                    frame_cnt-1,'.bin'
                INQUIRE( FILE=uStr, EXIST=retVal)
                write(*, *) uStr, retVal
                open(unit=101, file=frameStr, status='unknown')
            
                write(101, '(A,F16.4)') "TIME: ", dimTime
                write(101, '(A,A)') "DATA_FILE: ", uStr
                write(101, '(A,I5.5,I5.5,I5.5)') "DATA_SIZE: ", Nx, Ny, Nz
                !! Allowable values for DATA_FORMAT are: BYTE, INT, FLOAT, DOUBLE
                write(101, '(A,A)') "DATA_FORMAT: ", 'DOUBLE'
                write(101, '(A,A)') "VARIABLE: ", 'u'
                !! Endian representation of the computer that created the data.
                !! Intel is LITTLE, many other processors are BIG.
                write(101,'(A,A)') "DATA_ENDIAN: ", endian
                !! Centering refers to how the data is distributed in a cell. If you
                !! give “zonal” then it’s 1 data value per zone. Otherwise the data
                !! will be centered at the nodes.
                !! write(101, *) 'CENTERING:', 'zonal'
                !! BRICK_ORIGIN lets you specify a new coordinate system origin for
                !! the mesh that will be created to suit your data.
                write(101, '(A,F10.5,F10.5,F10.5)') "BRICK_ORIGIN: ", 0.0, 0.0, (dz/2.0)*z_i
                !! BRICK_SIZE lets you specify the size of the brick.
                write(101, '(A,F16.4, F16.4, F16.4)') "BRICK_SIZE: ", 2.0*pi*z_i, 2.0*pi*z_i, L_z
                !! Additional BOV options:
                !! BYTE_OFFSET: is optional and lets you specify some number of
                !! bytes to skip at the front of the file. This can be useful for
                !! skipping the 4-byte header that Fortran tends to write to files.
                !! If your file does not have a header then DO NOT USE BYTE_OFFSET.

                !! write(101, *) "BYTE_OFFSET:  "

                !! DIVIDE_BRICK: is optional and can be set to “true” or “false”.
                !! When DIVIDE_BRICK is true, the BOV reader uses the values stored
                !! in DATA_BRICKLETS to divide the data into chunks that can be
                !! processed in parallel.

                !! write(101, *) "DIVIDE_BRICK: ", "false"

                !! DATA_BRICKLETS: is optional and requires you to specify 3 integers
                !! that indicate the size of the bricklets to create when you have
                !! also specified the DIVIDE_BRICK option. The values chosen for
                !! DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE.

                !! write(101, *) "DATA_BRICKLETS: " 

                !! DATA_COMPONENTS: is optional and tells the BOV reader how many
                !! components your data has. 1=scalar, 2=complex number, 3=vector,
                !! 4 and beyond indicate an array variable. You can use “COMPLEX”
                !! instead of “2” for complex numbers.
                write(101, '(A,I2.2)') "DATA_COMPONENTS: ", 1
                close(101)
            end if
            INQUIRE( FILE=frameStr, EXIST=retVal ) 



        ! open the header file on output diretory to write velocity vector
            if(output_frames(1) .and. output_frames(2) .and.  output_frames(3)) then
                write(frameStr, '(A,I4.4,A)')'output/v_', &
                    noCurFrame,'.bov'
                write(vStr, '(A,I4.4,A)')'output/v_frame/v_frame', &
                    frame_cnt-1,'.bin'
                INQUIRE( FILE=vStr, EXIST=retVal)
                write(*, *) vStr, retVal
                if(retVal) then 
                    open(unit=101, file=frameStr, status='unknown')
                
                    write(101, '(A,F16.4)') "TIME: ", dimTime
                    write(101, '(A,A)') "DATA_FILE: ", vStr
                    write(101, '(A,I5.5,I5.5,I5.5)') "DATA_SIZE: ", Nx, Ny, Nz
                    !! Allowable values for DATA_FORMAT are: BYTE, INT, FLOAT, DOUBLE
                    write(101, '(A,A)') "DATA_FORMAT: ", 'DOUBLE'
                    write(101, '(A,A)') "VARIABLE: ", 'v'
                    !! Endian representation of the computer that created the data.
                    !! Intel is LITTLE, many other processors are BIG.
                    write(101,'(A,A)') "DATA_ENDIAN: ", endian
                    !! Centering refers to how the data is distributed in a cell. If you
                    !! give “zonal” then it’s 1 data value per zone. Otherwise the data
                    !! will be centered at the nodes.
                    !! write(101, *) 'CENTERING:', 'zonal'
                    !! BRICK_ORIGIN lets you specify a new coordinate system origin for
                    !! the mesh that will be created to suit your data.
                    write(101, '(A,F10.5,F10.5,F10.5)') "BRICK_ORIGIN: ", 0.0, 0.0, (dz/2.0)*z_i
                    !! BRICK_SIZE lets you specify the size of the brick.
                    write(101, '(A, F16.4, F16.4, F16.4)') "BRICK_SIZE: ", 2.0*pi*z_i, 2.0*pi*z_i, L_z
                    !! Additional BOV options:
                    !! BYTE_OFFSET: is optional and lets you specify some number of
                    !! bytes to skip at the front of the file. This can be useful for
                    !! skipping the 4-byte header that Fortran tends to write to files.
                    !! If your file does not have a header then DO NOT USE BYTE_OFFSET.

                    !! write(101, *) "BYTE_OFFSET:  "

                    !! DIVIDE_BRICK: is optional and can be set to “true” or “false”.
                    !! When DIVIDE_BRICK is true, the BOV reader uses the values stored
                    !! in DATA_BRICKLETS to divide the data into chunks that can be
                    !! processed in parallel.

                    !! write(101, *) "DIVIDE_BRICK: ", "false"

                    !! DATA_BRICKLETS: is optional and requires you to specify 3 integers
                    !! that indicate the size of the bricklets to create when you have
                    !! also specified the DIVIDE_BRICK option. The values chosen for
                    !! DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE.

                    !! write(101, *) "DATA_BRICKLETS: " 

                    !! DATA_COMPONENTS: is optional and tells the BOV reader how many
                    !! components your data has. 1=scalar, 2=complex number, 3=vector,
                    !! 4 and beyond indicate an array variable. You can use “COMPLEX”
                    !! instead of “2” for complex numbers.
                    write(101, '(A,I1.1)') "DATA_COMPONENTS: ", 1
                    close(101)
                end if ! retVal
            end if
            INQUIRE( FILE=frameStr, EXIST=retVal ) 


        ! open the header file on output diretory to write velocity vector
            if(output_frames(1) .and. output_frames(2) .and.  output_frames(3)) then
                write(frameStr, '(A,I4.4,A)')'output/w_', &
                    noCurFrame,'.bov'
                write(wStr, '(A,I4.4,A)')'output/w_frame/w_frame', &
                    frame_cnt-1,'.bin'
                INQUIRE( FILE=wStr, EXIST=retVal)
                write(*, *) wStr, retVal
                if(retVal) then 
                    open(unit=101, file=frameStr, status='unknown')
                
                    write(101, '(A,F16.4)') "TIME: ", dimTime
                    write(101, '(A,A)') "DATA_FILE: ", wStr
                    write(101, '(A,I5.5,I5.5,I5.5)') "DATA_SIZE: ", Nx, Ny, Nz
                    !! Allowable values for DATA_FORMAT are: BYTE, INT, FLOAT, DOUBLE
                    write(101, '(A,A)') "DATA_FORMAT: ", 'DOUBLE'
                    write(101, '(A,A)') "VARIABLE: ", 'w'
                    !! Endian representation of the computer that created the data.
                    !! Intel is LITTLE, many other processors are BIG.
                    write(101,'(A,A)') "DATA_ENDIAN: ", endian
                    !! Centering refers to how the data is distributed in a cell. If you
                    !! give “zonal” then it’s 1 data value per zone. Otherwise the data
                    !! will be centered at the nodes.
                    !! write(101, *) 'CENTERING:', 'zonal'
                    !! BRICK_ORIGIN lets you specify a new coordinate system origin for
                    !! the mesh that will be created to suit your data.
                    write(101, '(A,F10.5,F10.5,F10.5)') "BRICK_ORIGIN: ", 0.0, 0.0, 0.0
                    !! BRICK_SIZE lets you specify the size of the brick.
                    write(101, '(A, F16.4, F16.4, F16.4)') "BRICK_SIZE: ", 2.0*pi*z_i, 2.0*pi*z_i, L_z
                    !! Additional BOV options:
                    !! BYTE_OFFSET: is optional and lets you specify some number of
                    !! bytes to skip at the front of the file. This can be useful for
                    !! skipping the 4-byte header that Fortran tends to write to files.
                    !! If your file does not have a header then DO NOT USE BYTE_OFFSET.

                    !! write(101, *) "BYTE_OFFSET:  "

                    !! DIVIDE_BRICK: is optional and can be set to “true” or “false”.
                    !! When DIVIDE_BRICK is true, the BOV reader uses the values stored
                    !! in DATA_BRICKLETS to divide the data into chunks that can be
                    !! processed in parallel.

                    !! write(101, *) "DIVIDE_BRICK: ", "false"

                    !! DATA_BRICKLETS: is optional and requires you to specify 3 integers
                    !! that indicate the size of the bricklets to create when you have
                    !! also specified the DIVIDE_BRICK option. The values chosen for
                    !! DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE.

                    !! write(101, *) "DATA_BRICKLETS: " 

                    !! DATA_COMPONENTS: is optional and tells the BOV reader how many
                    !! components your data has. 1=scalar, 2=complex number, 3=vector,
                    !! 4 and beyond indicate an array variable. You can use “COMPLEX”
                    !! instead of “2” for complex numbers.
                    write(101, '(A,I1.1)') "DATA_COMPONENTS: ", 1
                    close(101)
                end if ! retVal
            end if
            INQUIRE( FILE=frameStr, EXIST=retVal ) 

     !!! write pressure frame bov header
        ! open the header file on output diretory to write velocity vector
            if(output_frames(4)) then
                write(frameStr, '(A,I4.4,A)')'output/p_', &
                    noCurFrame,'.bov'
                write(pStr, '(A,I4.4,A)')'output/p_frame/p_frame', &
                    frame_cnt-1,'.bin'
                INQUIRE( FILE=pStr, EXIST=retVal)
                write(*, *) pStr, retVal
                if(retVal) then
                    open(unit=101, file=frameStr, status='unknown')                
                    write(101, '(A,F16.4)') "TIME: ", dimTime
                    write(101, '(A,A)') "DATA_FILE: ", pStr
                    write(101, '(A,I5.5,I5.5,I5.5)') "DATA_SIZE: ", Nx, Ny, Nz
                    !! Allowable values for DATA_FORMAT are: BYTE, INT, FLOAT, DOUBLE
                    write(101, '(A,A)') "DATA_FORMAT: ", 'DOUBLE'
                    write(101, '(A,A)') "VARIABLE: ", 'p'
                    !! Endian representation of the computer that created the data.
                    !! Intel is LITTLE, many other processors are BIG.
                    write(101,'(A,A)') "DATA_ENDIAN: ", endian
                    !! Centering refers to how the data is distributed in a cell. If you
                    !! give “zonal” then it’s 1 data value per zone. Otherwise the data
                    !! will be centered at the nodes.
                    !! write(101, *) 'CENTERING:', 'zonal'
                    !! BRICK_ORIGIN lets you specify a new coordinate system origin for
                    !! the mesh that will be created to suit your data.
                    write(101, '(A,I5.5,I5.5,I5.5)') "BRICK_ORIGIN: ", 0.0, 0.0, dz/2.0*z_i
                    !! BRICK_SIZE lets you specify the size of the brick.
                    write(101, '(A, F16.4, F16.4, F16.4)') "BRICK_SIZE: ", 2.0*pi*z_i, 2.0*pi*z_i, L_z
                    !! Additional BOV options:
                    !! BYTE_OFFSET: is optional and lets you specify some number of
                    !! bytes to skip at the front of the file. This can be useful for
                    !! skipping the 4-byte header that Fortran tends to write to files.
                    !! If your file does not have a header then DO NOT USE BYTE_OFFSET.

                    !! write(101, *) "BYTE_OFFSET:  "

                    !! DIVIDE_BRICK: is optional and can be set to “true” or “false”.
                    !! When DIVIDE_BRICK is true, the BOV reader uses the values stored
                    !! in DATA_BRICKLETS to divide the data into chunks that can be
                    !! processed in parallel.

                    !! write(101, *) "DIVIDE_BRICK: ", "false"

                    !! DATA_BRICKLETS: is optional and requires you to specify 3 integers
                    !! that indicate the size of the bricklets to create when you have
                    !! also specified the DIVIDE_BRICK option. The values chosen for
                    !! DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE.

                    !! write(101, *) "DATA_BRICKLETS: " 

                    !! DATA_COMPONENTS: is optional and tells the BOV reader how many
                    !! components your data has. 1=scalar, 2=complex number, 3=vector,
                    !! 4 and beyond indicate an array variable. You can use “COMPLEX”
                    !! instead of “2” for complex numbers.
                    write(101, '(A,I1.1)') "DATA_COMPONENTS: ", 1
                    close(101)
                end if ! retVal
            end if
            INQUIRE( FILE=frameStr, EXIST=retVal ) 

     !!! write lambda2 frame bov header
        ! open the header file on output diretory to write velocity vector
            if(output_frames(8)) then
                write(frameStr, '(A,I4.4,A)')'output/lambda2_', &
                    noCurFrame,'.bov'
                write(lambda2Str, '(A,I4.4,A)')'output/lambda2_frame/lambda2_frame', &
                    frame_cnt-1,'.bin'
                INQUIRE( FILE=lambda2Str, EXIST=retVal)
                write(*, *) lambda2Str, retVal
                if(retVal) then 
                    open(unit=101, file=frameStr, status='unknown')                
                    write(101, '(A,F16.4)') "TIME: ", dimTime
                    write(101, '(A,A)') "DATA_FILE: ", lambda2Str
                    write(101, '(A,I5.5,I5.5,I5.5)') "DATA_SIZE: ", Nx, Ny, Nz
                    !! Allowable values for DATA_FORMAT are: BYTE, INT, FLOAT, DOUBLE
                    write(101, '(A,A)') "DATA_FORMAT: ", 'DOUBLE'
                    write(101, '(A,A)') "VARIABLE: ", 'lambda2'
                    !! Endian representation of the computer that created the data.
                    !! Intel is LITTLE, many other processors are BIG.
                    write(101,'(A,A)') "DATA_ENDIAN: ", endian
                    !! Centering refers to how the data is distributed in a cell. If you
                    !! give “zonal” then it’s 1 data value per zone. Otherwise the data
                    !! will be centered at the nodes.
                    !! write(101, *) 'CENTERING:', 'zonal'
                    !! BRICK_ORIGIN lets you specify a new coordinate system origin for
                    !! the mesh that will be created to suit your data.
                    write(101, '(A,F10.5,F10.5,F10.5)') "BRICK_ORIGIN: ", 0.0, 0.0, (dz/2.0)*z_i
                    !! BRICK_SIZE lets you specify the size of the brick.
                    write(101, '(A, F16.4, F16.4, F16.4)') "BRICK_SIZE: ", 2.0*pi*z_i, 2.0*pi*z_i, L_z
                    !! Additional BOV options:
                    !! BYTE_OFFSET: is optional and lets you specify some number of
                    !! bytes to skip at the front of the file. This can be useful for
                    !! skipping the 4-byte header that Fortran tends to write to files.
                    !! If your file does not have a header then DO NOT USE BYTE_OFFSET.

                    !! write(101, *) "BYTE_OFFSET:  "

                    !! DIVIDE_BRICK: is optional and can be set to “true” or “false”.
                    !! When DIVIDE_BRICK is true, the BOV reader uses the values stored
                    !! in DATA_BRICKLETS to divide the data into chunks that can be
                    !! processed in parallel.

                    !! write(101, *) "DIVIDE_BRICK: ", "false"

                    !! DATA_BRICKLETS: is optional and requires you to specify 3 integers
                    !! that indicate the size of the bricklets to create when you have
                    !! also specified the DIVIDE_BRICK option. The values chosen for
                    !! DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE.

                    !! write(101, *) "DATA_BRICKLETS: " 

                    !! DATA_COMPONENTS: is optional and tells the BOV reader how many
                    !! components your data has. 1=scalar, 2=complex number, 3=vector,
                    !! 4 and beyond indicate an array variable. You can use “COMPLEX”
                    !! instead of “2” for complex numbers.
                    write(101, '(A,I1.1)') "DATA_COMPONENTS: ", 1
                    close(101)
                end if !! retVal
            end if
            INQUIRE( FILE=frameStr, EXIST=retVal ) 


     !!! write filtered lambda2 frame bov header
        ! open the header file on output diretory to write velocity vector
            if(filtLevelsLambda > 0) then
                do looper = 1,filtLevelsLambda                     
                    write(lambda2FiltStr, '(A,I4.4,A,I4.4,A)')'output/lambda2_frame/lambdaFrame_cutoff_', &
                        nint(filtsLambda(looper)),'_frame',(frame_cnt-1),'.bin'
                    INQUIRE( FILE=lambda2FiltStr, EXIST=retVal)
                    write(*, *) lambda2FiltStr, retVal
                    write(frameStr, '(A,I4.4,A,I4.4,A)')'output/lambdaFrame_cutoff_', &
                        nint(filtsLambda(looper)),'_frame',(frame_cnt-1),'.bov'
                    if(retVal) then  
                        open(unit=101, file=frameStr, status='unknown')                    
                        write(101, '(A,F16.4)') "TIME: ", dimTime
                        write(101, '(A,A)') "DATA_FILE: ", lambda2FiltStr
                        write(101, '(A,I5.5,I5.5,I5.5)') "DATA_SIZE: ", Nx, Ny, Nz
                        !! Allowable values for DATA_FORMAT are: BYTE, INT, FLOAT, DOUBLE
                        write(101, '(A,A)') "DATA_FORMAT: ", 'DOUBLE'
                        write(101, '(A,A)') "VARIABLE: ", 'lambda2'
                        !! Endian representation of the computer that created the data.
                        !! Intel is LITTLE, many other processors are BIG.
                        write(101,'(A,A)') "DATA_ENDIAN: ", endian
                        !! Centering refers to how the data is distributed in a cell. If you
                        !! give “zonal” then it’s 1 data value per zone. Otherwise the data
                        !! will be centered at the nodes.
                        !! write(101, *) 'CENTERING:', 'zonal'
                        !! BRICK_ORIGIN lets you specify a new coordinate system origin for
                        !! the mesh that will be created to suit your data.
                        write(101, '(A,F10.5,F10.5,F10.5)') "BRICK_ORIGIN: ", 0.0, 0.0, (dz/2.0)*z_i
                        !! BRICK_SIZE lets you specify the size of the brick.
                        write(101, '(A, F16.4, F16.4, F16.4)') "BRICK_SIZE: ", 2.0*pi*z_i, 2.0*pi*z_i, L_z
                        !! Additional BOV options:
                        !! BYTE_OFFSET: is optional and lets you specify some number of
                        !! bytes to skip at the front of the file. This can be useful for
                        !! skipping the 4-byte header that Fortran tends to write to files.
                        !! If your file does not have a header then DO NOT USE BYTE_OFFSET.

                        !! write(101, *) "BYTE_OFFSET:  "

                        !! DIVIDE_BRICK: is optional and can be set to “true” or “false”.
                        !! When DIVIDE_BRICK is true, the BOV reader uses the values stored
                        !! in DATA_BRICKLETS to divide the data into chunks that can be
                        !! processed in parallel.

                        !! write(101, *) "DIVIDE_BRICK: ", "false"

                        !! DATA_BRICKLETS: is optional and requires you to specify 3 integers
                        !! that indicate the size of the bricklets to create when you have
                        !! also specified the DIVIDE_BRICK option. The values chosen for
                        !! DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE.

                        !! write(101, *) "DATA_BRICKLETS: " 

                        !! DATA_COMPONENTS: is optional and tells the BOV reader how many
                        !! components your data has. 1=scalar, 2=complex number, 3=vector,
                        !! 4 and beyond indicate an array variable. You can use “COMPLEX”
                        !! instead of “2” for complex numbers.
                        write(101, '(A,I1.1)') "DATA_COMPONENTS: ", 1
                        close(101)
                    end if !! retVal
                end do !! i = 1,filtLevelsLambda  
            end if
            INQUIRE( FILE=frameStr, EXIST=retVal ) 
        end if !! me == 0
    end if

end function writeBovHeader



