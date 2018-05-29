function u=loadframe(SA,Nx,Ny,Nz,endian)

%check that file exists
if(~exist(SA))
    error(['File ',SA,' does not exist.'])
end

fid=fopen(SA,'r',endian);
u=reshape(fread(fid,Nx*Ny*Nz,'double'),Nx,Ny,Nz);
fclose(fid);

return
end