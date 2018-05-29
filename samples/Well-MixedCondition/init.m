clear,clc
input_path='../../input/';

%input parameters

z_o=0.1;

N_turnovers=25;

Nx=36;
Ny=36;
Nz=36;
hprocs=2;

l_r=1;
z_i=1000;
u_star=0.45;
l_z=1000;

press_cor=0;
f_p=u_star^2/l_z;
mom_nodes=0;

verticalBC=0;
nu=0;
theta_mean_wind=0;

cs_count=2;
nframe=0;

c_flag=0;
scalarcount=0;

npart=100000;
part_model=2;
release_type=1;
start_release=20000;
partstep=1;
skip_step=1;
freq=100000;
nr=100000;

%calculated parameters
dz=l_z/(Nz-1);

Umax=u_star/0.4*log(l_z/z_o);
dt=dz/Umax/10;

nsteps=round(N_turnovers/(dt*u_star/l_z));

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
writeinput('f_p',f_p,input_path);
writeinput('mom_nodes',mom_nodes,input_path);
writeinput('verticalBC',verticalBC,input_path);
writeinput('nu',nu,input_path);
writeinput('theta_mean_wind',theta_mean_wind,input_path);
writeinput('cs_count',cs_count,input_path);
writeinput('nframe',nframe,input_path);
writeinput('c_flag',c_flag,input_path);
writeinput('scalarcount',scalarcount,input_path);
writeinput('npart',npart,input_path);
writeinput('part_model',part_model,input_path);
writeinput('release_type',release_type,input_path);
writeinput('start_release',start_release,input_path);
writeinput('partstep',partstep,input_path);
writeinput('skip_step',skip_step,input_path);
writeinput('freq',freq,input_path);
writeinput('nr',nr,input_path);

system(['cp ',input_path,'LESinputs.txt  LESinputs.txt']);

%check that write was successful
if( ~exist([input_path,'LESinputs.txt'],'file') )
    error('Creation of LESinputs.txt failed.')
end

z_u=linspace(dz/2,l_z+dz/2,Nz);
z_w=linspace(0,l_z,Nz);

%create velocity fields
u=zeros(Nx,Ny,Nz);
v=zeros(Nx,Ny,Nz);
w=zeros(Nx,Ny,Nz);

for k=1:Nz
    
    ur=randn(Nx,Ny);
    vr=randn(Nx,Ny);
    wr=randn(Nx,Ny);
    
   u(:,:,k)=u_star/0.4*log(z_u(k)/z_o)+ur/u_star^2;
   v(:,:,k)=vr/u_star^2;
   w(:,:,k)=wr/u_star^2;
    
end
zo=ones(Nx,Ny)*z_o;

%write to file
fw = fopen([input_path,'vel.ini'],'w','l');
fwrite(fw,u,'double');
fwrite(fw,v,'double');
fwrite(fw,w,'double');
fclose(fw);

fw = fopen([input_path,'zo.ini'],'w','l');
fwrite(fw,zo,'double');
fclose(fw);

disp('done creating initial files')
