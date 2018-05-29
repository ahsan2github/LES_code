clear,clc
input_path='../../input/';

%input parameters

z_o=0.16;

nsteps=4000;
dt=3.25;

Nx=40;
Ny=40;
Nz=40;
hprocs=1;

l_r=1;
z_i=1100.0789667;
u_star=0.45;
l_z=2400;

press_cor=0;
f_p=0;
Ug=0.0;
Vg=0.0;
f_c=0;

verticalBC=0;
nu=0;
theta_mean_wind=0;

cs_count=1;
nframe=0;

c_flag=0;

npart=0;

scalarcount=1;
scalarFlags=1;
surfaceFlags=0;
surfaceFluxes=0.06;
scalarScales=1;
S_advec=0;
inversion=0.003;
theta_0=300;
Ri_flag=1;
sponge=1;
z_d=1800;
rlx_time=300;

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
system(['cp ',input_path,'LESinputs.txt.master ',input_path,'LESinputs.txt']);
writeinput('nsteps',nsteps,input_path);
writeinput('dt',dt,input_path);
writeinput('Nx',Nx,input_path);
writeinput('Ny',Ny,input_path);
writeinput('Nz',Nz,input_path);
writeinput('hprocs',hprocs,input_path);
writeinput('l_r',l_r,input_path);
writeinput('z_i',z_i,input_path);
writeinput('u_star',u_star,input_path);
writeinput('l_z',l_z,input_path);
writeinput('verticalBC',verticalBC,input_path);
writeinput('Ug',Ug,input_path);
writeinput('Vg',Vg,input_path);
writeinput('press_cor',press_cor,input_path);
writeinput('f_p',f_p,input_path);
writeinput('theta_mean_wind',theta_mean_wind,input_path);
writeinput('f_c',f_c,input_path);
writeinput('nu',nu,input_path);
writeinput('cs_count',cs_count,input_path);
writeinput('nframe',nframe,input_path);
writeinput('c_flag',c_flag,input_path);
writeinput('npart',npart,input_path);

writeinput('scalarcount',scalarcount,input_path);
writeinput('scalarFlags',scalarFlags,input_path);
writeinput('surfaceFlags',surfaceFlags,input_path);
writeinput('surfaceFluxes',surfaceFluxes,input_path);
writeinput('scalarScales',scalarScales,input_path);
writeinput('S_advec',S_advec,input_path);
writeinput('inversion',inversion,input_path);
writeinput('theta_0',theta_0,input_path);
writeinput('sponge',sponge,input_path);
writeinput('z_d',z_d,input_path);
writeinput('rlx_time',rlx_time,input_path);
writeinput('Ri_flag',Ri_flag,input_path);

system(['cp ',input_path,'LESinputs.txt  LESinputs.txt']);

%check that write was successful
if( ~exist([input_path,'LESinputs.txt'],'file') )
    error('Creation of LESinputs.txt failed.')
end

%create temperature and velocity fields
u=zeros(Nx,Ny,Nz);
v=zeros(Nx,Ny,Nz);
w=zeros(Nx,Ny,Nz);
T=zeros(Nx,Ny,Nz);

z  = [1:Nz]*l_z/(Nz-1)-l_z/(2*(Nz-1));
zw = [0:Nz-1]*l_z/(Nz-1);

gamma = 0.003;
g=9.81;
zio = 1600;
zi1 = 0.844*zio;
To = 300;
Os = 0.06;
wso = ((g/To)*Os*zio)^(1/3);
Tso = Os/wso;

for k=1:Nz
    if(z(k)<zi1)
        T(:,:,k) = To+0.1*(rand(Nx,Ny)-0.5)*(1-z(k)/zi1)*Tso;
    else
        T(:,:,k) = To+(z(k)-zi1)*gamma;
    end
    if(k == 1)
        w(:,:,k) = zeros(Nx,Ny);
    elseif(zw(k)<zi1)
        w(:,:,k) = 0.1 *(rand(Nx,Ny)-0.5)*(1-zw(k)/zi1)*wso;
    else
        w(:,:,k) = zeros(Nx,Ny);
    end
    if(k == 1)
        u(:,:,k) = Ug + 0.1*(rand(Nx,Ny)-0.5)*wso;
        v(:,:,k) = 0.1*(rand(Nx,Ny)-0.5)*wso;
    else
        u(:,:,k) = Ug + zeros(Nx,Ny);
        v(:,:,k) = zeros(Nx,Ny);
    end
end

zo=ones(Nx,Ny)*z_o;

figure;
for k=1:Nz
    MaT(k)=sum(sum(T(:,:,k)))/(Nx*Ny);
end
plot(MaT,z)

%write to file
fw = fopen([input_path,'vel.ini'],'w','l');
fwrite(fw,u,'double');
fwrite(fw,v,'double');
fwrite(fw,w,'double');
fclose(fw);

fw = fopen([input_path,'zo.ini'],'w','l');
fwrite(fw,zo,'double');
fclose(fw);

fw = fopen([input_path,'temperature.ini'],'w','l');
fwrite(fw,T,'double');
fclose(fw);

disp('done creating initial files')
