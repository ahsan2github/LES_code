clear,clc
input_path='../../input/';

%input parameters

z_o=0.02;

nsteps=30000;
dt=0.05;

Nx=96;
Ny=48;
Nz=30;
hprocs=1;

l_r=2;
z_i=20;
u_star=0.5;
l_z=60;

press_cor=3;
f_p=0;
press_step=5;
u_avg=2;

verticalBC=0;
nu=0;
theta_mean_wind=0;

cs_count=2;
nframe=0;

c_flag=1;
Cd=0.15;
h_canopy=20;

npart=0;
scalarcount=0;

SS_LAD=[2.1236 2.3199 2.7751 3.7024 5.5269 6.8554 7.5469 7.4832 6.8768 5.1139 0]/h_canopy;
SS_U=[0.2081 0.2003 0.1743 0.1455 0.1436 0.1657 0.215 0.2972 0.4454 0.7467 ...
    0.9639 1.0582 1.1524 1.2106 1.2626 1.3028 1.34 1.3742 1.4113 1.4364 1.4705...
    1.4926 1.5238 1.5457 1.562 1.5721 1.5793 1.5805 1.5875 1.6006];
    
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
writeinput('press_cor',press_cor,input_path);
writeinput('press_step',press_step,input_path);
writeinput('u_avg',u_avg,input_path);
writeinput('f_p',f_p,input_path);
writeinput('verticalBC',verticalBC,input_path);
writeinput('nu',nu,input_path);
writeinput('theta_mean_wind',theta_mean_wind,input_path);
writeinput('cs_count',cs_count,input_path);
writeinput('nframe',nframe,input_path);
writeinput('c_flag',c_flag,input_path);
writeinput('Cd',Cd,input_path);
writeinput('h_canopy',h_canopy,input_path);
writeinput('npart',npart,input_path);
writeinput('scalarcount',scalarcount,input_path);

system(['cp ',input_path,'LESinputs.txt  LESinputs.txt']);

%check that write was successful
if( ~exist([input_path,'LESinputs.txt'],'file') )
    error('Creation of LESinputs.txt failed.')
end

%create velocity fields
u=zeros(Nx,Ny,Nz);
for k=1:Nz
    u(:,:,k)=SS_U(k)*u_avg+randn(Nx,Ny)*u_star;
end
v=randn(Nx,Ny,Nz)*u_star;
w=randn(Nx,Ny,Nz)*u_star;
zo=ones(Nx,Ny)*z_o;

%create LAD profiles (u-nodes)
LAD=zeros(Nx,Ny,Nz);
LAD(:,:,1)=SS_LAD(1);
for k=2:ceil(h_canopy/dz)
    LAD(:,:,k)=(SS_LAD(k)+SS_LAD(k+1))/2;
end

%write to file
fw = fopen([input_path,'u.ini'],'w','l');
fwrite(fw,u,'double');
fclose(fw);
fw = fopen([input_path,'v.ini'],'w','l');
fwrite(fw,v,'double');
fclose(fw);
fw = fopen([input_path,'w.ini'],'w','l');
fwrite(fw,w,'double');
fclose(fw);

fw = fopen([input_path,'zo.ini'],'w','l');
fwrite(fw,zo,'double');
fclose(fw);

fw = fopen([input_path,'PlantDensity.ini'],'w','l');
fwrite(fw,LAD,'double');
fclose(fw);

disp('done creating initial files')
