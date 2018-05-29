function writeinput(var_str,var_val,input_path)
%WRITEINPUT This is a utility for replacing input values in a current
%LESinputs.txt file. 
%
%WRITEINPUT(string,value,path) replaces the input designated by
%'string' with 'value'. 'input_path' specifies the absolute or relative path to the
%LESinputs.txt file. By default, WRITEINPUT will place a copy of the
%modified LESinputs.txt file in the current directory as well as in
%'input_path'.
%
%Examples:
%   writeinput('Nx',96../../input/')

%check if input_path has a trailing slash
if( ~strcmp(input_path(end),'/') )
    input_path(end+1)='/';
end

if( strcmp(var_str,'framefiles') )
    write_str=var_val;
else
    write_str=num2str(var_val);
end
    
%MODEL PARAMETERS AND FLAGS
if( strcmp(var_str,'mom_nodes') )
    search_str='.*= mom_nodes,.*';
    write_str=[write_str,' \t\t\t= mom_nodes, vertical level for computation of dynamic momentum optimization (0=w-nodes, 1=uvp-nodes)'];
elseif( strcmp(var_str,'nsteps') )
    search_str='.*= nsteps,.*';
    write_str=[write_str,' \t\t\t= nsteps, number of timesteps'];
elseif( strcmp(var_str,'dt') )
    write_str=add_double(write_str);
    search_str='.*= dt,.*';
    write_str=[write_str,' \t\t= dt, physical timestep (sec)'];
elseif( strcmp(var_str,'Nx') )
    search_str='.*= Nx,.*';
    write_str=[write_str,' \t\t\t= Nx, x-direction # of grid points, streamwise direction'];
elseif( strcmp(var_str,'Ny') )
    search_str='.*= Ny,.*';
    write_str=[write_str,' \t\t\t= Ny, y-direction # of grid points, spanwise direction'];
elseif( strcmp(var_str,'Nz') )
    search_str='.*= Nz,.*';
    write_str=[write_str,' \t\t\t= Nz, z-direction # of grid points, wall-normal direction'];
elseif( strcmp(var_str,'hprocs') )
    search_str='.*= hprocs,.*';
    write_str=[write_str,' \t\t\t= hprocs, number of processors, horizontal direction'];
elseif( strcmp(var_str,'l_r') )
    search_str='.*= l_r,.*';
    write_str=[write_str,' \t\t\t= l_r, ratio of streamwise to spanwise (an integer), l_r = Lx/Ly'];
elseif( strcmp(var_str,'z_i') )
    write_str=add_double(write_str);
    search_str='.*= z_i,.*';
    write_str=[write_str,' \t\t\t= z_i, normalization height, dimensionless depth of domain=L_z/z_i (recall: Lx/z_i=2*pi, Ly/z_i=2*pi/l_r)'];
elseif( strcmp(var_str,'u_star') )
    write_str=add_double(write_str);
    search_str='.*= u_star,.*';
    write_str=[write_str,' \t\t\t= u_star, normalization velocity, dimensionless velocity=u/u_star'];
elseif( strcmp(var_str,'l_z') )
    write_str=add_double(write_str);
    search_str='.*= l_z,.*';
    write_str=[write_str,' \t\t\t= l_z, dimensioned depth of domain'];
elseif( strcmp(var_str,'verticalBC') )
    search_str='.*= verticalBC,.*';
    write_str=[write_str,' \t\t\t= verticalBC, type of vertical boundary condition used, 0=rigid lid/wall, 1=periodic'];
elseif( strcmp(var_str,'Ug') )
    write_str=add_double(write_str);
    search_str='.*= Ug,.*';
    write_str=[write_str,' \t\t\t= Ug, geostrophic wind x-component [m/s]'];
elseif( strcmp(var_str,'Vg') )
    write_str=add_double(write_str);
    search_str='.*= Vg,.*';
    write_str=[write_str,' \t\t\t= Vg, geostrophic wind y-component [m/s]'];
elseif( strcmp(var_str,'press_cor') )
    search_str='.*= press_cor,.*';
    write_str=[write_str,' \t\t\t= press_cor, 0=constant pressure forcing 1=const geo wind 2=variable geo wind 3=const Uavg forcing'];
elseif( strcmp(var_str,'f_p') )
    write_str=add_double(write_str);
    search_str='.*= f_p,.*';
        write_str=[write_str,' \t\t\t= f_p, pressure gradient (dynamic) if press_cor=0 (m/s^2)'];
elseif( strcmp(var_str,'press_step') )
    search_str='.*= press_step,.*';
    write_str=[write_str,' \t\t\t= press_step, frequency of pressure forcing calculations if press_cor=3']; 
elseif( strcmp(var_str,'u_avg') )
    write_str=add_double(write_str);
    search_str='.*= u_avg,.*';
    write_str=[write_str,' \t\t\t= u_avg, desired vertical average of u-profile (m/s) if press_cor=3'];     
elseif( strcmp(var_str,'theta_mean_wind') )
    write_str=add_double(write_str);
    search_str='.*= theta_mean_wind,.*';
    write_str=[write_str,' \t\t\t= theta_mean_wind, meand wind direction (degrees), takes values between +/- 90 degrees']; 
elseif( strcmp(var_str,'f_c') )
    write_str=add_double(write_str);
    search_str='.*= f_c,.*';
    write_str=[write_str,' \t\t\t= f_c, coriolis parameter, f_c=1.45*sin(latitude)']; 
elseif( strcmp(var_str,'nu') )
    write_str=add_double(write_str);
    search_str='.*= nu,.*';
    write_str=[write_str,' \t\t\t= nu, kinematic viscosity (m^2/s^2) (=0 for inviscid flow)'];
elseif( strcmp(var_str,'c_count') )
    search_str='.*= c_count,.*';
    write_str=[write_str,' \t\t\t= c_count, frequency of computation of statistics'];
elseif( strcmp(var_str,'p_count') )
    search_str='.*= p_count,.*';
    write_str=[write_str,' \t\t\t= p_count, frequency of printing statistics'];
elseif( strcmp(var_str,'cs_count') )
    search_str='.*= cs_count,.*';
    write_str=[write_str,' \t\t\t= cs_count, frequency of computing sub-grid scale coefficients'];
elseif( strcmp(var_str,'ruler') )
    search_str='.*= ruler,.*';
    write_str=[write_str,' \t\t\t= ruler, frequency of updating checkpoint files'];


%SCALAR PARAMETERS
elseif( strcmp(var_str,'scalarcount') )
    search_str='.*= scalarcount,.*';
    write_str=[write_str,' \t\t\t= scalarcount, number of scalars to compute (ie  0...N)'];
elseif( strcmp(var_str,'scalarFlags') )
    search_str='.*= scalarFlags,.*';
    write_str=[write_str,' \t\t\t= scalarFlags, type of scalar, 1=temperature, 2=moisture, 0=passive (one flag for each scalar to be computed)'];
elseif( strcmp(var_str,'surfaceFlags') )
    search_str='.*= surfaceFlags,.*';
    write_str=[write_str,' \t\t\t= surfaceFlags, boundary condition of scalar, 0=constant flux, 1=cooling rate, 2=SEB'];
elseif( strcmp(var_str,'surfaceFluxes') )
    write_str=add_double(write_str);
    search_str='.*= surfaceFluxes,.*';
    write_str=[write_str,' \t\t\t= surfaceFluxes, surface flux of scalars when surfaceFlag==0 [K-m/s]'];
elseif( strcmp(var_str,'scalarScales') )
     write_str=add_double(write_str);
    search_str='.*= scalarScales,.*';
    write_str=[write_str,' \t\t\t= scalarScales, normalization scale of scalar - used for non-dimensional computation'];
elseif( strcmp(var_str,'S_advec') )
    search_str='.*= S_advec,.*';
    write_str=[write_str,' \t\t\t= S_advec,'];
elseif( strcmp(var_str,'inversion') )
    write_str=add_double(write_str);
    search_str='.*= inversion,.*';
    write_str=[write_str,' \t\t\t= inversion, inversion strength at domain top (lapse rate) [K/m]'];
elseif( strcmp(var_str,'theta_0') )
     write_str=add_double(write_str);
    search_str='.*= theta_0,.*';
    write_str=[write_str,' \t\t\t= theta_0, buoyancy flux reference temperature [K]'];
elseif( strcmp(var_str,'sponge') )
    search_str='.*= sponge,.*';
    write_str=[write_str,' \t\t\t= sponge, 1=dampening layer starting at z_d, 0=none'];
elseif( strcmp(var_str,'z_d') )
    write_str=add_double(write_str);
    search_str='.*= z_d,.*';
    write_str=[write_str,' \t\t\t= z_d, height to start dampening layer if sponge==1 [m]'];
elseif( strcmp(var_str,'rlx_time') )
    write_str=add_double(write_str);
    search_str='.*= rlx_time,.*';
    write_str=[write_str,' \t\t\t= rlx_time, relaxation time for dampening layer [sec]'];
elseif( strcmp(var_str,'Ri_flag') )
    search_str='.*= Ri_flag,.*';
    write_str=[write_str,' \t\t\t= Ri_flag, Richardson # criteria, if pointwise gradient Ri # exceeds pointwise SGS prandtl #, then SGS terms switched OFF, 1=criteria on, 0=off'];


%PARTICLE DISPERSION PARAMETERS
elseif( strcmp(var_str,'npart') )
    search_str='.*= npart,.*';
    write_str=[write_str,' \t\t\t= npart, total number of particles to be released'];
elseif( strcmp(var_str,'part_model') )
    search_str='.*= part_model,.*';
    write_str=[write_str,' \t\t\t= part_model, 1=no model, 2=Weil SFS model'];
elseif( strcmp(var_str,'start_release') )
    search_str='.*= start_release,.*';
    write_str=[write_str,' \t\t\t= start_release, timestep to begin releases'];
elseif( strcmp(var_str,'partstep') )
    search_str='.*= partstep,.*';
    write_str=[write_str,' \t\t\t= partstep, # of timesteps between releases'];
elseif( strcmp(var_str,'skip_step') )
    search_str='.*= skip_step,.*';
    write_str=[write_str,' \t\t\t= skip_step, # of timesteps between trajectory calcs'];
elseif( strcmp(var_str,'freq') )
    search_str='.*= freq,.*';
    write_str=[write_str,' \t\t\t= freq, # of particles per release (per release position)'];
elseif( strcmp(var_str,'nr') )
    search_str='.*= nr,.*';
    write_str=[write_str,' \t\t\t= nr, number of release points']; 
elseif( strcmp(var_str,'deposition') )
    search_str='.*= deposition,.*';
    write_str=[write_str,' \t\t\t= deposition, (0=no deposition model, 1=deposition model)']; 
elseif( strcmp(var_str,'wd') )
    search_str='.*= wd,.*';
    write_str=[write_str,' \t\t\t= wd, particle drift velocity (m/s)']; 



%CANOPY PARAMETERS
elseif( strcmp(var_str,'c_flag') )
    search_str='.*= c_flag,.*';
    write_str=[write_str,' \t\t\t= c_flag, (= 1 if canopy is present)'];
elseif( strcmp(var_str,'Cd') )
    write_str=add_double(write_str);
    search_str='.*= Cd,.*';
    write_str=[write_str,' \t\t\t= Cd, canopy drag coefficient'];
elseif( strcmp(var_str,'h_canopy') )
    write_str=add_double(write_str);
    search_str='.*= h_canopy,.*';
    write_str=[write_str,' \t\t\t= h_canopy, height of canopy (m)'];



%WRITING AND STATISTICS
elseif( strcmp(var_str,'nframe') )
    search_str='.*= nframe,.*';
    write_str=[write_str,' \t\t\t= nframe, number of instantaneous frames to output (0 will output nothing, skip frame outputting)'];
elseif( strcmp(var_str,'sframe') )
    search_str='.*= sframe,.*';
    write_str=[write_str,' \t\t\t= sframe, timestep at which to start saving frames'];
elseif( strcmp(var_str,'framestep') )
    search_str='.*= framestep,.*';
    write_str=[write_str,' \t\t\t= framestep, number of timesteps between frames'];
elseif( strcmp(var_str,'frameh') )
    search_str='.*= frameh,.*';
    write_str=[write_str,' \t\t\t= frameh, height of frames in grid points'];
elseif( strcmp(var_str,'framefiles') )
    search_str='.*= framefiles,.*';
    write_str_t=[num2str(length(write_str)),' \t\t\t= framefiles, number of frame variables to output (=0 outputs nothing)'];
    for i=1:length(write_str)
        write_str_t=[write_str_t,'\n',char(write_str{i})];
    end
    write_str=write_str_t;


%ERROR
else
    error(['Invalid input variable: ',var_str])
end

sed_str=['sed ''s;',search_str,';',write_str,';'' ',input_path, ...
         'LESinputs.txt > tempfile; mv tempfile ',input_path,'LESinputs.txt'];

system(sed_str);

[~,result]=system('cat LESinputs.txt');
if(isempty(result))
    error('%s%s%s','Write for ',var_str,' did not succeed.\n')
end

return
end

function write_str=add_double(var_str)

    if(isempty(strfind(var_str,'.')))
        write_str=[var_str,'.d0'];
    else
         write_str=[var_str,'d0'];
    end

end
