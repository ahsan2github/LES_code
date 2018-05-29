    subroutine sponge_layer(scalar,u,v,w,scalarRHS,RHSx,RHSy,RHSz, &
    rdmp,flag)
    integer*4 :: flag
    real*8,dimension(:,:,:) :: rdmp,u,v,w,RHSx,RHSy,RHSz
    real*8,dimension(:,:,:,:) :: scalar,scalarRHS
    end subroutine sponge_layer
