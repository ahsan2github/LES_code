!      file created by mohammad, Feb 12, 2013
    subroutine write3D(fh,u,i,j,k,fileview,view_disp)
    use globals
    implicit none
    real(kind=8), dimension(:,:,:)  :: u
    integer(kind = 4)                         :: i,j,k
    integer(kind=MPI_OFFSET_KIND)   ::view_disp
    integer(kind = 4)                          :: fh
    integer(kind = 4)                         ::fileview
    end subroutine write3D
