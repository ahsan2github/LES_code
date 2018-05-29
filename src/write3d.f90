!      file created by mohammad, Feb 12, 2013
    subroutine write3D(fh,u,i,j,k,filetyp,view_disp)
        use globals
        implicit none
        real(kind=8), dimension(:,:,:)  :: u
        integer(kind=4)                 :: i,j,k
        integer(kind=MPI_OFFSET_KIND)   :: view_disp
        integer(kind=4)                 :: fh
        integer(kind=4)                 ::filetyp
        call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
        filetyp,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_WRITE_ALL(fh,u(i,j,k),Nx*Nyb*Nzb, &
        MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    end subroutine write3D
