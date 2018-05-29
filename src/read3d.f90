!     Read one 3D variable from a binary file using MPI I/O

    subroutine Read3D(u,i,j,k,filename,fileview_3D)

    use globals
    implicit none

    real(kind=8), dimension(:,:,:)  :: u
    integer(kind=4)                 ::fileview_3D,i,j,k
    character(len=* )               ::filename

    integer(kind=MPI_OFFSET_KIND)   :: view_disp
    integer(kind=4)                 ::fh

    call MPI_FILE_OPEN(nall,filename, &
    MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    view_disp=0
    call MPI_FILE_SET_VIEW(fh,view_disp,MPI_DOUBLE_PRECISION, &
    fileview_3D,"native",MPI_INFO_NULL,ierr)
    call MPI_FILE_READ_ALL(fh,u(i,j,k), &
    Nx*Nyb*Nzb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FILE_CLOSE(fh,ierr)

    return

    end subroutine Read3D
