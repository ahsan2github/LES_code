clear,clc
input_path='.';

%input parameters

continuous=0;  %0=homogeneous, 1=row-resolved
LAI=0.077;
row_space=6.0;

z_o=0.005;
row_width=0.5;
h_canopy=2.0;
n_width=2;

nsteps=250000;
dt=0.004;

Nx=192;
Ny=192;
Nz=160;
hprocs=12;

l_r=1;
u_star=0.5;
z_i=row_width/n_width*Nx/(2*pi);
l_z=8*h_canopy;

press_cor=0;
f_p=0.025;
press_step=0;
u_avg=0.0;

verticalBC=0;
nu=0;
theta_mean_wind=0;

c_count=100;
p_count=500;

cs_count=2;
nframe=2000;
sframe=50000;
framestep=175;
framefiles={'part_frame'};

c_flag=1;
Cd=0.5;

npart=1000000;
part_model=2;
start_release=50000;
skip_step=1;
freq=1;
deposition=0;

scalarcount=0;

LAD_prof=[2.0 2.0 2.0 2.0 2.0 2.0 2.0 3.0 3.8 6.0 12.2 12.9 12.4 11.6 10.3 9.2 8.5 7.4 5.9 2.0];
    
%calculated parameters
dz=l_z/(Nz-1);

%add utilities to the path temorarily
addpath ../utilities

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
writeinput('c_count',c_count,input_path);
writeinput('p_count',p_count,input_path);
writeinput('cs_count',cs_count,input_path);
writeinput('nframe',nframe,input_path);
writeinput('sframe',sframe,input_path);
writeinput('framestep',framestep,input_path);
writeinput('framefiles',framefiles,input_path);
writeinput('c_flag',c_flag,input_path);
writeinput('Cd',Cd,input_path);
writeinput('h_canopy',h_canopy,input_path);
writeinput('npart',npart,input_path);
writeinput('part_model',part_model,input_path);
writeinput('start_release',start_release,input_path);
writeinput('skip_step',skip_step,input_path);
writeinput('freq',freq,input_path);
writeinput('deposition',deposition,input_path);
writeinput('scalarcount',scalarcount,input_path);

%system(['cp ',input_path,'LESinputs.txt  LESinputs.txt']);

%check that write was successful
if( ~exist([input_path,'LESinputs.txt'],'file') )
    error('Creation of LESinputs.txt failed.')
end

%create velocity fields
u=randn(Nx,Ny,Nz)*u_star;
v=randn(Nx,Ny,Nz)*u_star;
w=randn(Nx,Ny,Nz)*u_star;
zo=ones(Nx,Ny)*z_o;

dx=2*pi*z_i/Nx;

%scale LAD profile to have correct LAI

LAD_top=0.25;
Z = linspace(0,h_canopy,length(LAD_prof));

%LAD_prof(end)=LAD_top;

g = 0:0.001:1;
for i=1:length(g)
    scl=LAD_prof*g(i);
    %scl(end)=LAD_top;
    res(i)=abs(LAI-sum(scl*dz));
end
[err,ind] = min(res);

prof = LAD_prof*g(ind);
%prof(end)=LAD_top;

%create LAD field
LAD=zeros(Nx,Ny,Nz);

n_space = round(row_space/dx);
n_rows = ceil(Nx/(n_space+n_width));
        
start = round((Nx-(n_rows*n_width+(n_rows-1)*n_space))/2);

if(start<n_space/4)
    n_rows=n_rows-1;
    start=start+round((n_space+n_width)/2);
end

if(continuous==1)
    for i=1:n_rows
        for j=1:Ny
            for ii=start+(i-1)*(n_space+n_width):start+(i-1)*(n_space+n_width)+n_width-1;
                LAD(ii,j,1:length(LAD_prof)) = prof(:);
            end
        end
    end
else
    for j=1:Ny
        for i=1:Nx
             LAD(i,j,1:length(LAD_prof))=prof(:);
        end
    end
end
    
%Particle release positions
h=linspace(0.4,h_canopy,5);
y = linspace(0,2*pi*z_i/l_r,7);y=y(2:end-1);
    
dr=(row_space+row_width)/4;

count=1;release_pos=[0,0,0];
for k=1:length(h)
    for i=1:n_rows
        for ii=1:4
            for j=1:length(y)
                x_pos=(start-1)*dx-dx/2+row_width/2+(row_width+row_space)*(i-1)+dr*(ii-1);
                release_pos(count,:) = [x_pos,y(j),h(k)];
                count=count+1;
            end
        end
    end
end
nr=count-1;
writeinput('nr',nr,input_path);
disp(['nr=',num2str(count-1)])
partstep=round(nsteps*(count-1)/npart);
writeinput('partstep',partstep,input_path);
disp(['partstep=',num2str(round(nsteps*(count-1)/npart))])

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

fw = fopen([input_path,'release_pos.ini'],'w','l');
fwrite(fw,release_pos,'double');
fclose(fw);

disp('done creating initial files')
