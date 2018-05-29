!     ... Read one 3D variable from a binary file using MPI I/O

    subroutine read3D(u,i,j,k,filename,fileview_3D)

    use globals
    implicit none

    real*8, dimension(:,:,:):: u
    integer*4 :: fileview_3D,i,j,k
    character(len=*) :: filename

    integer(kind=MPI_OFFSET_KIND) :: view_disp
    integer*4 :: fh

    end subroutine read3D
