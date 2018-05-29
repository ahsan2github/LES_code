clear,clc
input_path='../../input/';

%input parameters

z_o=0.16;

nsteps=35000;
dt=0.5;

Nx=48;
Ny=48;
Nz=48;
hprocs=1;

l_r=1;
z_i=795.775;
u_star=0.5;
l_z=2000;

press_cor=1;
Ug=3.6;
Vg=0.0;
f_c=0.5;

verticalBC=0;
nu=0;
theta_mean_wind=0;

cs_count=1;
nframe=500;
sframe=25000;
framestep=20;

c_flag=0;

npart=72900;
part_model=2;
start_release=25000;
partstep=225;
skip_step=3;
freq=60;
nr=81;
deposition=0;
wd=0.0;

scalarcount=1;
scalarFlags=1;
surfaceFlags=0;
surfaceFluxes=0.24;
scalarScales=1;
S_advec=0;
inversion=0.003;
theta_0=294;
Ri_flag=0;
sponge=2;
z_d=1500;
rlx_time=60;

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
writeinput('theta_mean_wind',theta_mean_wind,input_path);
writeinput('f_c',f_c,input_path);
writeinput('nu',nu,input_path);
writeinput('cs_count',cs_count,input_path);
writeinput('nframe',nframe,input_path);
writeinput('sframe',sframe,input_path);
writeinput('framestep',framestep,input_path);
writeinput('c_flag',c_flag,input_path);
writeinput('npart',npart,input_path);
writeinput('part_model',part_model,input_path);
writeinput('start_release',start_release,input_path);
writeinput('partstep',partstep,input_path);
writeinput('skip_step',skip_step,input_path);
writeinput('freq',freq,input_path);
writeinput('nr',nr,input_path);
writeinput('deposition',deposition,input_path);
writeinput('wd',wd,input_path);
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

for k=1:Nz
    u(:,:,k)=u_star/1.3*log(z(k)/z_o)+0.25*randn(Nx,Ny)*u_star;
end
v=0.25*randn(Nx,Ny,Nz)*u_star;

zi=900;
gamma1 = 0.06;
gamma2 = 0.003;
To    = 300.0; 
wso   = ((9.81/To)*surfaceFluxes*zi)^(1/3);
Tso   = surfaceFluxes/wso;
H = round(zi/dz);

Tprof(1:H) = To+0.1*(rand-0.5)*(1-z(1:H)/zi)*Tso;
Tprof(H+1:H+6) = To+(z(H+1:H+6)-zi)*gamma1;
Tprof(H+7:Nz) = Tprof(H+6)+(z(H+7:Nz)-zi)*gamma2;

for k=1:Nz
    T(:,:,k) = Tprof(k);
    
    if(k == 1)
        w(:,:,k) = zeros(Nx,Ny);
    elseif(zw(k)<zi)
        w(:,:,k) = 0.1 *(rand(Nx,Ny)-0.5)*(1-zw(k)/zi)*wso;
    else
        w(:,:,k) = zeros(Nx,Ny);
    end
    if(k == Nz)
        u(:,:,k) = Ug + 0.1*(rand(Nx,Ny)-0.5)*wso;
        v(:,:,k) = 0.1*(rand(Nx,Ny)-0.5)*wso;
    end
end

zo=ones(Nx,Ny)*z_o;

%particle release positions
zos=0.07*zi;
DX=2*pi*z_i/9;
DY=2*pi*z_i/9;
[X,Y]=meshgrid(DX/2:DX:2*pi*z_i-DX/2,DY/2:DY:2*pi*z_i-DY/2);
for i=1:length(X(:))
    release_pos(i,:)=[X(i),Y(i),zos];
end

disp(release_pos)

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

fw = fopen([input_path,'release_pos.ini'],'w','l');
fwrite(fw,release_pos,'double');
fclose(fw);

disp('done creating initial files')
