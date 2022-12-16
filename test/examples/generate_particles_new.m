clc; clear;
addpath('/scratch365/cswarm/kramos/interface/generate_pgfem3d_input_stack_with_interfaces_local_new/generate_layer_spherical_pack')
% initialization of library path
c_dir = pwd;
path(path, '..');
setenv('LD_LIBRARY_PATH', '');
generate_new_particles = 1;
%% create output directory
fid_wdir = fopen('working_dir.txt', 'r');
wdir = fgetl(fid_wdir);
fclose(fid_wdir);
in = [wdir,'/input/Rocpack'];
check_dir(in);
out = [wdir,'/input/t3d'];
check_dir(out);
%%%%%%%%%%%%%%% MODIFY PARAMETERS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sampling parameters
file_base = 'particles';     % file_base name 
layer_thickness = 0.01;      % in [mm] 
particle_gap = 0.001;
vf0 = 0.40;                % particle fraction
r= [0.0148];  % particle radius size 
box = [1000e-3, 1000e-3, 1000e-3];           % cell size Lx, Ly, Lz
N = [20];            % number of particles to create 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_base = 'particles';
ext_in = '.in';
ext_rp = '.rp';
ext_rp2stat3D = '.rp2stat3D';
ext_rp2t3d = '.rp2t3d';
ext_t3d = '.t3d';
ext_regions = '.out.regions';

%%generate particles
rng('shuffle'); % seeds the random number generator 

if(generate_new_particles)
  ia=1; 
  fn_rp      = [in, '/', file_base, ext_rp];
  fn_in      = [in, '/', file_base, ext_in];
  fn_rp2stat = [in, '/', file_base, ext_rp2stat3D]; 
  fn_rp2t3d  = [in, '/', file_base, ext_rp2t3d]
  fn_t3d     = [out, '/', file_base, ext_t3d]
  fn_regions = [out, '/', file_base, ext_regions]
  seed_number = round( 1 + 1000*rand(1));
  
  fid = fopen(fn_in, 'w');
  fprintf(fid, 'set seed = %d;\n\n', seed_number);
  fprintf(fid, 'set packing_fraction = %e;\n\n', vf0);
  fprintf(fid, 'boundary { box %e %e %e periodic}\n\n', box(1), box(2), box(3));
  for i=1:length(r)
  fprintf(fid, 'create %d sphere size %e ;\n', N(i), r(i));
  end  
  fclose(fid);
  
  cmd = ['module load gcc; /cswarm/tools/bin/pack-4-8-0-SHAKE -o', fn_rp, ' ', fn_in];
  fprintf('run: %s\n', cmd);
  
  if(system(cmd)==0)
    fprintf('run: [%s] was successful.\n', cmd);
  else
    fprintf('run: [%s] failed.\n', cmd);
    return;
  end
  
  %Generate the rp2stat3D file which decreases the particle sizes
  resize_particles(fn_rp,fn_rp2stat, particle_gap); 
  %Genereate the rp2t3d file 
  get_particle_positions(fn_rp2stat, fn_rp2t3d);
  %Generate the t3d file 
  gen_T3D_layer_spherical_pack(fn_rp2t3d, fn_t3d, layer_thickness)
  %Generate the regions file from t3d file 
  write_regions(fn_t3d, fn_regions)
  
end 
fprintf('\n');  
  
