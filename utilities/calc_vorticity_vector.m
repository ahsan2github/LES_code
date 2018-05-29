function omega=calc_vorticity_vector(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
%NOTE: output is on w-nodes

omega=cell(3,1);

omega{1}=dwdy-dvdz;
omega{2}=dudz-dwdx;

omega3=zeros(size(dudy,1),size(dudy,2),size(dudy,3));
for k=size(dudy,3):-1:2
    omega3(:,:,k)=(dvdx(:,:,k)+dvdx(:,:,k-1))/2-(dudy(:,:,k)+dudy(:,:,k-1))/2;
end
omega{3}=omega3;


return
end