function readinputs(SA)
% READINPUTS reads inputs from LES code
%   READINPUTS('DIRECTORYNAME') loads input parameters from LESinputs.txt 
%   file into the MATLAB workspace, where DIRECTORYNAME is a string 
%   containing the full system path to the directory containing the 
%   LESinputs.txt file. 

%check that directory was input with trailing slash
if(strcmp(SA(end),'/')~=1)
    SA=[SA,'/']; %if not, fix it
end

%check that directory exists
if(exist(SA,'dir')~=7)
    error('%s%s%s','Error: Directory ',SA,' does not exist.')
end

%check that LESinputs.txt file exists
if(exist([SA,'input/LESinputs.txt'],'file')~=2)
    error('%s%s','Error: input/LESinputs.txt file does not exist in directory ',...
        SA)
end

%read nsteps
name='nsteps';
readFile(name,SA)

%read dt
name='dt';
readFile(name,SA)

%read Nx
name='Nx';
readFile(name,SA)

%read Ny
name='Ny';
readFile(name,SA)

%read Nz
name='Nz';
readFile(name,SA)

%read l_r
name='l_r';
readFile(name,SA)

%read z_i
name='z_i';
readFile(name,SA)

%read u_star
name='u_star';
readFile(name,SA)

%read l_z
name='l_z';
readFile(name,SA)

%read nu
name='nu';
readFile(name,SA)

%read h_canopy
name='h_canopy';
readFile(name,SA)

%read p_count
name='p_count';
readFile(name,SA)

%read nframe
name='nframe';
readFile(name,SA)
% 
% if(nframe>0)
% 
    %read sframe
    name='sframe';
    readFile(name,SA)
     
    %read framestep
    name='framestep';
    readFile(name,SA)
    
    %read frameh
    name='frameh';
    readFile(name,SA)
% 
% end

%read scalarcount
name='scalarcount';
readFile(name,SA)

if(scalarcount>0)
    
    %read surfaceFluxes
    name='surfaceFluxes';
    readFile(name,SA)
        
    %read scalarScales
    name='scalarScales';
    readFile(name,SA)
    
    %read theta_0
    name='theta_0';
    readFile(name,SA)
    
end

%read npart
name='npart';
readFile(name,SA)

if(npart>0)
   
    %read start_release
    name='start_release';
    readFile(name,SA)
    
    %read partstep
    name='partstep';
    readFile(name,SA)
    
    %read skip_step
    name='skip_step';
    readFile(name,SA)
    
    %read freq
    name='freq';
    readFile(name,SA)
    
    %read nr
    name='nr';
    readFile(name,SA)
    
end

%read canopy_h
name='h_canopy';
readFile(name,SA)

%read zo
if(exist([SA,'input/zo.ini'],'file'))
    fid=fopen([SA,'input/zo.ini'],'r','l');
    zo=fread(fid,1,'double');
    fclose(fid);
    assignin('base','zo',zo)
    assignin('caller','zo',zo)
else
    disp('Warning: zo.ini file does not exist.')
end

%calculated values
if(exist('Nx','var')==1 && exist('z_i','var')==1)
    dx=2*pi*z_i/Nx;
    assignin('base','dx',dx)
else
    disp('Warning: dx could not be calculated')
end
if(exist('Ny','var')==1 && exist('z_i','var')==1)
    dy=2*pi*z_i/Ny;
    assignin('base','dy',dy)
else
    disp('Warning: dy could not be calculated')
end
if(exist('Nz','var')==1 && exist('l_z','var')==1)
    dz=l_z/(Nz-1);
    assignin('base','dz',dz)
else
    disp('Warning: dz could not be calculated')
end

return
end

function readFile(name,SA)
str=['grep ''= ',name,''' ',SA,'input/LESinputs.txt | cut -f1 -d''='''];
[status,result]=system(str);
if(status==0)
    %check if result is a double (i.e., contains 'd0')
    if(isempty(strfind(result,'d'))==0)
        result(strfind(result,'d'))='0';
    end
    %assign value to base and caller workspaces
    assignin('base',name,str2double(result))
    assignin('caller',name,str2double(result))
else
    sprintf('%s%s%s','Warning: Variable "',name,'" was not found')
end

return
end
