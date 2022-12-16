%% matlab initialization
c_dir = pwd;
path(path, '..');

% directory setup
filebase = 'pack';

initial_disp = 1.0e-4;
density = [8.88e-9, 2.71e-9]; % particle matrix

% mat = 0- matrix, 1-interface (not here), 2-particle properties 
% homogeneous material properties 
mat = [2 2.4e+3 18.254e+2 0.0 8.9552e+2 0.0 0.0 0.34 0.0 0.0 1 2
       0 8.0e+2 6.5299e+2 0.0 2.9851e+2 0.0 0.0 0.34 0.0 0.0 1 2];

BC = [5,1,1,1 % Z-bottom
      6,1,1,1]; % Z-top

Mat.mat = mat;

% find regions
fn_region = sprintf('%s.out.regions', filebase);
fid_region = fopen(['./input/MS_PGFem3D/',fn_region], 'r');
tline = fgetl(fid_region);
region_no = sscanf(tline, '%d');
region = zeros(region_no, 4);
for ib=1: region_no
    tline = fgetl(fid_region);
    region(ib, :) = sscanf(tline, '%d %d %d %d');
end
fclose(fid_region);

Mat.mat_ids = [region(:, 2), zeros(region_no, 1)];
Mat.mat_ids(end, 2) = 1;    

MS_write_PGFem3D_mat_bc_json_files(['./input/MS_PGFem3D/',filebase], BC, Mat, initial_disp, density);

%% Update periodic.json for interface multiscale simulation
% Make sure the directory is right
%
% PERIODICITY
filename_ss = ['./input/MS_PGFem3D/',filebase,'.out.periodic'];
filename_to = './input/MS_PGFem3D/periodic.csv'; %temporary file
filename_ms = './input/MS_PGFem3D/periodic.json';

n_data = zeros(6,1);

% STEP 1 : READ SS OUTPUS
fid = fopen(filename_ss);
data = fscanf(fid,'%ld');
fclose(fid);

% STEP 1.5 : CONVERT TO CSV FILE
fid= fopen(filename_to,'w');
for ii = 1:6
    count = 0;
    data_start = data(1);
    n_data(ii) = data_start;
    while count < data_start*2-1
        data(1) = [];
        fprintf(fid,'%ld, ',data(1));
        count = count +1;
    end
    data(1) = [];
    fprintf(fid,'%ld\n',data(1));
    data(1) = [];
end
fclose(fid);
fprintf('... csv file is written\n');
fprintf('... size of data = %d, %d, %d, %d, %d, %d\n',n_data);
clear data

% STEP 2 : WRITE MS INPUTS
data = csvread('./input/MS_PGFem3D/periodic.csv');
fid = fopen(filename_ms,'w');
fprintf(fid,'{\n');
fprintf(fid,'    "pairs":\n');
fprintf(fid,'    [\n');
for ii=1:6
    fprintf(fid,'        [ ');
    for jj=1:n_data(ii)
        fprintf(fid,'{"type":%ld, "ID":%ld}',data(ii,2*jj-1),data(ii,2*jj));
        if jj<n_data(ii)
            fprintf(fid,',');
        end
    end
    fprintf(fid,']');
    
    if ii==6
        fprintf(fid,'\n');
    else
        fprintf(fid,',\n');
    end
end
fprintf(fid,'    ]\n');
fprintf(fid,'}\n');
fclose(fid);

system(['rm ',filename_to]);
fprintf('... complete!\n');

