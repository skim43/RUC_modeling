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

fn_rp      = [in, '/', file_base, ext_rp];
fn_in      = [in, '/', file_base, ext_in];
fn_rp2t3d  = [in, '/', file_base, ext_rp2t3d];
fn_t3d     = [out, '/', file_base];


%% Generate particles - RocPack
r = 0.010;
lc = 1.0e-1;
box = [1.0e-1 1.0e-1 lc-2*r];
%box = [2.0e-1+2*r 2.0e-1+2*r lc-2*r];
%seed_number = round( 1 + 1000*rand(1)); 816 - 318 - 951 - 35 - 440
seed_number = 440;
N = 4;

vf = box(1)*box(2)*box(3);
vf0 = r^3*4/3*pi * N;
vf0 = vf0/vf;

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
return
end

%% Read particles -
fid = fopen(fn_rp, 'r');
pno = 0;
while(~feof(fid))
    line = fgetl(fid);
    count = strfind(line,'sphere');
    if(numel(count)>=1)
        if(numel(count)==1)
            pno = pno + 1;
        else
            fprintf('Rocpack file format is weird: multiple spheres are read in a line\n');
        end
    end
end
fclose(fid);

fid = fopen(fn_rp, 'r');
P = zeros(pno, 3); % translate "position"
count = 0; % number of particle

while (~feof(fid))
  n = 0;
  while(n==0)
    if(feof(fid))
      break;
    end
    line = fgetl(fid);
    [~, n] = sscanf(line, 'sphere %s\n');
  end
  if(n~=0)
    count = count + 1;
    if(count>pno)
      fprintf('Something is wrong on reading particles: number of particle (%d) is incorrect\n', pno);
    end

    line = fgetl(fid);
    [p, ~] = sscanf(line, '     translate <%f, %f, %f>\n');
    P(count, :) = p;
    line = fgetl(fid);
   end
end
fclose(fid);

%% Generate particles - RoP2T3d
S = r;
% S = 0.0182;
tol = 1.0e-5;
% lc = 1.0e-1;

small_box = [1.0e-1 1.0e-1 lc];
box  = [2.0e-1 2.0e-1 lc];
nbox = 2;
cut_box = [-box(1)/2, box(1)/2
           -box(2)/2, box(2)/2
           -box(3)/2, box(3)/2];
P

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





