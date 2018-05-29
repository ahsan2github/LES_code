function [dudx,dudy,dudz]=ddx(u,var,dz,z_i,l_r,NX)

Nx=size(u,1);
Ny=size(u,2);
Nz=size(u,3);

l_x=Nx/NX*z_i;
l_y=z_i/l_r;

%---------DUDX, DUDY------------%

%wavenumbers
Kx = zeros(Nx,Ny);
for j=1:Ny
    Kx(:,j)=[(0:Nx/2),(-Nx/2+1:-1)]'/l_x;
end
Ky = zeros(Nx,Ny);
for i=1:Nx
    Ky(i,:)=[(0:Ny/2),(-Ny/2+1:-1)]'/l_y;
end

Kx(Nx/2+1,:) = 0.; %set Nyquist to zero
Ky(:,Ny/2+1) = 0.;

%2D Fourier Transform
u_hat=zeros(Nx,Ny,Nz);dudx_hat=u_hat;dudy_hat=u_hat;
for k=1:Nz
    u_hat(:,:,k)=fftn(u(:,:,k));
    %horizontal derivatives
    dudx_hat(:,:,k) = Kx.*u_hat(:,:,k)*1i;
    dudy_hat(:,:,k) = Ky.*u_hat(:,:,k)*1i;
end

dudx=zeros(Nx,Ny,Nz);dudy=zeros(Nx,Ny,Nz);
for k=1:Nz
    dudx(:,:,k)=ifftn(dudx_hat(:,:,k),'symmetric');
    dudy(:,:,k)=ifftn(dudy_hat(:,:,k),'symmetric');
end

%vertical derivatives
dudz=zeros(Nx,Ny,Nz);
for k=1:Nz-1
   if((strcmp(var,'u')||strcmp(var,'v')) && k==1)

   else
       dudz(:,:,k)=(u(:,:,k+1)-u(:,:,k))/dz;
   end
end

return
end