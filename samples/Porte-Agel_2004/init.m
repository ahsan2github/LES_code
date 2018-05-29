clear,clc
input_path='../../input/';

%input parameters

z_o=0.1;

N_turnovers=15;

Nx=32;
Ny=32;
Nz=32;
hprocs=4;

l_r=1;
z_i=1000;
u_star=0.45;
l_z=1000;

mom_nodes=0;
press_cor=0;
f_p=u_star^2/l_z;

verticalBC=0;
nu=0;
theta_mean_wind=0;

cs_count=2;
nframe=0;

c_flag=0;

scalarcount=1;
scalarFlags=0;
surfaceFlags=0;
surfaceFluxes=1;
scalarScales=1;
S_advec=0;
Ri_flag=0;

npart=0;

%calculated parameters
dz=l_z/(Nz-1);

Umax=u_star/0.4*log(l_z/z_o);
dt=dz/Umax/5;

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
writeinput('mom_nodes',mom_nodes,input_path);
writeinput('verticalBC',verticalBC,input_path);
writeinput('press_cor',press_cor,input_path);
writeinput('f_p',f_p,input_path);
writeinput('verticalBC',verticalBC,input_path);
writeinput('nu',nu,input_path);
writeinput('theta_mean_wind',theta_mean_wind,input_path);
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
writeinput('Ri_flag',Ri_flag,input_path);

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
    uw=-(l_z-z_w(k))/l_z*u_star^2;

    if(k<Nz-6)
        varu=0.25*u_star^2;
        varv=0.25*u_star^2;
        varw=0.25*u_star^2;
    else
        varu=0;
        varv=0;
        varw=0;
    end
    
    ur=randn(Nx,Ny);
    vr=randn(Nx,Ny);
    wr=randn(Nx,Ny);

    u(:,:,k)=u_star/0.4*log(z_u(k)/z_o)+ur*varu;
    v(:,:,k)=vr*varv;
    w(:,:,k)=wr*varw;
    
end

for k=1:Nz
    U(k)=sum(sum(u(:,:,k)))/(Nx*Ny);
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

fw = fopen([input_path,'scalar1.ini'],'w','l');
fwrite(fw,ones(Nx,Ny,Nz),'double');
fclose(fw);

disp('done creating initial files')
