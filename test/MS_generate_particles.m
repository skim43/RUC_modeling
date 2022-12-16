clear all; close all; clc;

% initialization of library path
c_dir = pwd;
path(path, '..');
setenv('LD_LIBRARY_PATH', '');

%% create output directory
fid_wdir = fopen('working_dir.txt', 'r');
wdir = fgetl(fid_wdir);
fclose(fid_wdir);
in = [wdir,'/input/Rocpack'];
check_dir(in);
out = [wdir,'/input/MS_PGFem3D'];
check_dir(out);

%% file extension
file_base = 'pack';
ext_in = '.in';
ext_rp = '.rp';
ext_rp2stat3D = '.rp2stat3D';
ext_rp2t3d = '.rp2t3d';
ext_t3d = '.t3d';
ext_regions = '.out.regions';

fn_rp2t3d  = [in, '/', file_base, ext_rp2t3d];
fn_t3d     = [out, '/', file_base];
% fn_regions = [out, '/', file_base, ext_regions];

%% Generate particles
S = 0.0182;
lc=1.0e-1;
box  = [1.0e-1 1.0e-1 lc];
cut_box = [-box(1)/2, box(1)/2
           -box(2)/2, box(2)/2
           -box(3)/2, box(3)/2];

p_top = cut_box(3,2)/2;
p_btm = cut_box(3,1)/2;
p_left = cut_box(2,1);
p_right = cut_box(2,2);
p_front = cut_box(1,2);
p_rear = cut_box(1,1);

tol = 1.0e-8;
P = [
    % tol tol tol;
     0.0 p_right-tol 0.0; % right
     p_front-tol 0.0 0.0; % front
     p_front-tol p_right-tol p_top; % top
     p_front-tol p_right-tol p_btm; % bottom
     ];% R-RB

 vf0 = box(1)*box(2)*box(3);
  vf = S^3*4/3*pi * size(P,1);
  vf = vf/vf0;
  fprintf('INFO: volume fraction %4.8f\n',vf);   
  fprintf('INFO: reference volume fraction %4.8f\n',(S^3*4/3*pi)/(0.1^3));  
S = S*ones(size(P,1),1);

%% print Rockpack2T3D file, enclosure
fid = fopen(fn_rp2t3d, 'w');
fprintf(fid, 'MinDist 1.0e-6\n');
fprintf(fid, 'RefinementDist 0.1\n');
fprintf(fid, 'Cohesive 0\n');
fprintf(fid, 'Holes 0\n');
fprintf(fid, '\n');

box_encls = (cut_box(:, 2) - cut_box(:, 1))/2;
box_cent  = (cut_box(:, 2) + cut_box(:, 1))/2;
P(:, 1) = P(:, 1) - box_cent(1);
P(:, 2) = P(:, 2) - box_cent(2);
P(:, 3) = P(:, 3) - box_cent(3);

fprintf(fid, 'Enclosure %d %e %e %e %e %e %e\n', 1, 0.0, 0.0, 0.0, ...
	  box_encls(1), box_encls(2), box_encls(3));
fprintf(fid, '\n');

for b=1:size(P,1)
  fprintf(fid, 'Particle %d %e %e %e %e %d %d\n', b, P(b, 1),...
          P(b, 2),...
          P(b, 3),...
          S(b, 1),...
          0,0); % material IDs are all 0
end
fprintf(fid, '\n');
fprintf(fid, 'End\n');
fclose(fid);

%Generate the t3d file 
system(sprintf('/cswarm/tools/bin/RoP2T3d -i %s -o %s', fn_rp2t3d, fn_t3d));
fprintf('\n');



