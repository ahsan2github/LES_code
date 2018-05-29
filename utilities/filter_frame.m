function uf=filter_frame(u,fgr1,fgr2,Nx,Ny,Nz)


%filter
uf=zeros(Nx,Ny,Nz);
for k=1:Nz
    %wavenumbers
    
    kx=[(0:Nx/2),(-Nx/2+1:-1)];ky=[(0:Ny/2),(-Ny/2+1:-1)];
    f_kernel=zeros(Nx,Ny);
    for j=1:Ny
        for i=1:Nx
            if(abs(kx(i))>Nx/2/fgr1 || abs(ky(j))>Ny/2/fgr1)
                f_kernel(i,j)=0;
            else
                f_kernel(i,j)=1;
            end
%             if(abs(kx(i))<Nx/2/fgr2 || abs(ky(j))<Ny/2/fgr2)
%                 f_kernel(i,j)=0;
%             else
%                 f_kernel(i,j)=1;
%             end
        end
    end
    
    u_hat=fftn(u(:,:,k));
    u_hat=u_hat.*f_kernel;
    uf(:,:,k)=ifftn(u_hat,'symmetric');
end

return
end