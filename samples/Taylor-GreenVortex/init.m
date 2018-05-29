clear,clc
input_path='../../input/';

%input parameters

dimension=2;    %dimension of problem (2=2D, 3=3D)
dir_flag=1;     %direction of symmetry if dimension=2

nsteps=500;
dt=0.001;

Nx=32;
Ny=32;
Nz=32;

l_r=1;
z_i=1;
u_star=1;
l_z=2*pi;
press_cor=0;
f_p=0;
verticalBC=1;
nu=0.001;
theta_mean_wind=0;

cs_count=1;
nframe=50;
sframe=1;
framestep=10;
frameh=Nz;

npart=0;
c_flag=0;
scalar_count=0;

theta=pi/2;

%calculated parameters
dz=l_z/(Nz-1);

%add utilities to the path temorarily
addpath ../../utilities

%check if input_path has a trailing slash
if( ~strcmp(input_path(end),'/') )
    input_path(end+1)='/';
end

%write inputs to LESinputs.txt
if( ~exist(input_path,'dir') )
    error('Input directory does not exist.')
end
system(['rm ',input_path,'LESinputs.txt']);
system(['cp ',input_path,'LESinputs.txt.master ',input_path,'LESinputs.txt']);
writeinput('nsteps',nsteps,input_path);
writeinput('dt',dt,input_path);
writeinput('Nx',Nx,input_path);
writeinput('Ny',Ny,input_path);
writeinput('Nz',Nz,input_path);
writeinput('l_r',l_r,input_path);
writeinput('z_i',z_i,input_path);
writeinput('u_star',u_star,input_path);
writeinput('l_z',l_z,input_path);
writeinput('press_cor',press_cor,input_path);
writeinput('f_p',f_p,input_path);
writeinput('verticalBC',verticalBC,input_path);
writeinput('nu',nu,input_path);
writeinput('theta_mean_wind',theta_mean_wind,input_path);
writeinput('cs_count',cs_count,input_path);
writeinput('nframe',nframe,input_path);
writeinput('sframe',sframe,input_path);
writeinput('framestep',framestep,input_path);
writeinput('frameh',frameh,input_path);
writeinput('npart',npart,input_path);
writeinput('c_flag',c_flag,input_path);
writeinput('scalar_count',scalar_count,input_path);

%check that write was successful
if( ~exist([input_path,'LESinputs.txt'],'file') )
    error('Creation of LESinputs.txt failed.')
end

%build u,v,w fields
u=zeros(Nx,Ny,Nz);v=u;w=u;
if(dimension==2)
    if(dir_flag == 1)
        [X,Y]=meshgrid([0:Nx-1]*2*pi/Nx,[0:Ny-1]*2*pi/Ny);
        for k=1:Nz
            u(:,:,k) = sin(X).*cos(Y);
            v(:,:,k) = -cos(X).*sin(Y);
        end
    elseif(dir_flag == 2) 
        [X_u,Z_u]=meshgrid([0:Nx-1]*2*pi/Nx,[0.5:Nz-0.5]*l_z/(Nz-1));
        [X_w,Z_w]=meshgrid([0:Nx-1]*2*pi/Nx,[0:Nz-1]*l_z/(Nz-1));
        for j=1:Ny
            u(:,j,:) = sin(X_u).*cos(Z_u);
            w(:,j,:) = -cos(X_w).*sin(Z_w);
        end
    elseif(dir_flag == 3)
        [Y,Z]=meshgrid([0:Ny-1]*2*pi/Ny,[0:Nz-1]*2*pi/Nz);
        for i=1:Nx
            v(i,:,:) = sin(Y).*cos(Z);
            w(i,:,:) = -cos(Y).*sin(Z);
        end
    end
elseif(dimension==3)
    z_uv=linspace(dz/2,l_z+dz/2,Nz);
    z_w=linspace(0,l_z,Nz);
    y=linspace(0,2*pi*z_i/l_r,Ny);
    x=linspace(0,2*pi*z_i,Nx);
    for k=1:Nz
        for j=1:Ny
            for i=1:Nx
                %             u(i,j,k)=2/sqrt(3)*sin(theta+2*pi/3)*sin(x(i))*cos(y(j))*cos(z_uv(k));
                %             v(i,j,k)=2/sqrt(3)*sin(theta-2*pi/3)*sin(x(i))*cos(y(j))*cos(z_uv(k));
                %             w(i,j,k)=2/sqrt(3)*sin(theta)*sin(x(i))*cos(y(j))*cos(z_w(k));
                u(i,j,k)=sin(x(i))*cos(y(j))*cos(z_uv(k));
                v(i,j,k)=-cos(x(i))*sin(y(j))*cos(z_uv(k));
                w(i,j,k)=0;
            end
        end
    end
end

[~,input_contents]=system(['find ',input_path,' -name \*.ini']);
if( ~isempty(input_contents) )
    system(['rm ',input_path,'*.ini']);
end

%write to file
fw = fopen([input_path,'vel.ini'],'w','l');

fwrite(fw,u(:,:,:),'double');
fwrite(fw,v(:,:,:),'double');
fwrite(fw,w(:,:,:),'double');

fclose(fw);
disp('done creating initial files')