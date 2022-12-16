function MS_write_PGFem3D_mat_bc_json_files(filename_base, BC, Mat, initial_disp, density)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write interface multiscale PGFem3D input file for material and boundary condition
% input:
%   filename_base : file name base. {filename_base}_mat.json and {filename_base}_bc.json files will be created
%   Mat           : material object
%   BC            : boundary conditions
%   initial_disp  : prescribed initial displacement
%   denstiy       : material density  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODIFIED 01-Apr-2022
% by Sion Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 21-Aug-2018
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% matlab initialization
  c_dir = pwd;
  path(path, '..');

  dist = 10;
  
  if ~isfile('t3d')
      fprintf("ERROR:: T3D EXCUTABLE REQUIRED! \n\n");
      system('cp /cswarm/tools/bin/t3d ./');
  end
  % generate T3D outputs to read mesh
  fn_in = 'temp.out';
  cmd = sprintf('rm -rf *.out; t3d -i %s.t3d -o %s -d %d', filename_base, fn_in, dist);
  if(system(cmd)~=0)
    fprintf('T3D is failed\n');
    return;
  end
  % read number of region to set material
  fn_region = sprintf('%s.out.regions', filename_base);
  fid_region = fopen(fn_region, 'r');
  tline = fgetl(fid_region);
  region_no = sscanf(tline, '%d');
  region = zeros(region_no, 4);
  for ib=1: region_no
    tline = fgetl(fid_region);
    region(ib, :) = sscanf(tline, '%d %d %d %d');
  end
  fclose(fid_region);

  Mat.solidno = region_no;
  if(~isfield(Mat,'mat_ids'))
      Mat.mat_ids = region(:, 2:3);
  end
  %% generate header files
  [grid, node_info, ~] = load_from_t3d_out(fn_in, 0);

  xmin = zeros(grid.nsd,1);
  xmax = zeros(grid.nsd,1);
  for a=1: grid.nsd
    xmin(a) = min(grid.NODE(:, a));
    xmax(a) = max(grid.NODE(:, a));
  end
 
  % compute smallest distance of nodes
  nx = (grid.nodeno)^(1/3);
  min_r = min((xmax - xmin)./nx);
  tol = min_r*1.0e-3;
  Lx_ids = cell(grid.nsd, 1);
  Hx_ids = cell(grid.nsd, 1);
  for a = 1: grid.nsd
    m = node_info(abs(grid.NODE(:, a) - xmin(a)) < tol, :);
    Lx_ids{a} = unique(m(:, [1,2]), 'rows');
    clear m;
    m = node_info(abs(grid.NODE(:, a) - xmax(a)) < tol, :);
    Hx_ids{a} = unique(m(:, [1,2]), 'rows');
  end

  ids = cell(3,6);
  flag = [1 2 5];

  for flagid = 1: numel(flag)
    for xid = 1: grid.nsd
      m = Lx_ids{xid}(:, 1) == flag(flagid);
      ids{flagid, 2*xid-1} = Lx_ids{xid}(m, 2);
      clear m;
      m = Hx_ids{xid}(:, 1) == flag(flagid);
      ids{flagid, 2*xid} = Hx_ids{xid}(m, 2);
    end
  end

  vertex.LX = ids{1,1};
  vertex.HX = ids{1,2};
  vertex.LY = ids{1,3};
  vertex.HY = ids{1,4};
  vertex.LZ = ids{1,5};
  vertex.HZ = ids{1,6};

  edge.LX = ids{2,1};
  edge.HX = ids{2,2};
  edge.LY = ids{2,3};
  edge.HY = ids{2,4};
  edge.LZ = ids{2,5};
  edge.HZ = ids{2,6};

  face.LX = ids{3,1};
  face.HX = ids{3,2};
  face.LY = ids{3,3};
  face.HY = ids{3,4};
  face.LZ = ids{3,5};
  face.HZ = ids{3,6};

  Bnd.face = face;
  Bnd.edge = edge;
  Bnd.vertex = vertex;
  delete(fn_in);
%  geom = GEOM_GEN;
write_mat_bc_files(filename_base, Mat, Bnd, BC, initial_disp, density);
end



function write_mat_bc_files(str, Mat, Bnd, BC, disp, density)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write_mat_bc_files(str, Mat, Bnd, BC, disp, density)
%
% write material and boundary condition json files
% input:
%   str     : file name base. str_mat.json and str_bc.json files will be created
%   Mat     : material object
%   Bnd     : bunndary info
%   BC      : boundary conditions
%   disp    : initial displacement
%   density : material density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 21-Aug-2018
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  BC_ids = get_bounary_ids4BC(Bnd, BC);

  E = 100.0e+9;
  nu = 0.25;
  mu = E/2/(1+nu);
  
  if(nargin<6)
    density = [];
  end
  if(nargin<5)
    density = [];
    disp = [];
  end
  if(nargin<4)
    density = [];
    disp = [];
    BC = [];
  end
  if(nargin<3)
    density = [];
    disp = [];
    BC = [];
    Bnd = [];
    mat = [E mu/2 0 mu 0 0 nu 0 0 1.0 1.0 1.0 1.7e+03 1 2];
    Mat.mat = mat;
    Mat.solidno = 1;
    Mat.mat_ids = 0;
  end
  
  mat = Mat.mat;
  solidno = Mat.solidno;
  mat_ids = Mat.mat_ids;
  
  if(numel(mat)==0)
    mat = [E mu/2 0 mu 0 0 nu 0 0 1.0 1.0 1.0 1.7e+03 1 2];
  end
  if(numel(density)==0)
    density = zeros(size(mat, 1), 1);
  end
  
  % write material file
  matno = size(mat, 1);
  fid = fopen([str,'_mat.json'], 'w');
  
  fprintf(fid, '{\n');
  fprintf(fid, '"number_of_materials": %d,\n', size(mat, 1));
  fprintf(fid, '"number_of_volume_fraction": 1,\n');
  fprintf(fid, '"volume_fraction": [1],\n');
  fprintf(fid, '"number_of_bases": 1,\n');
  fprintf(fid, '"basis_vectors": [\n');
  fprintf(fid, '{\n');
  fprintf(fid, '  "e1": [1, 0, 0],\n');
  fprintf(fid, '  "e2": [0, 1, 0],\n');
  fprintf(fid, '  "e3": [0, 0, 1]\n');
  fprintf(fid, '}\n');
  fprintf(fid, '],\n');
  fprintf(fid, '"density": [');
  if(numel(density)>1)
    fprintf(fid, '%e, ', density(1:end-1));
  end
  fprintf(fid, '%e],\n', density(end));
  
  fprintf(fid, '"number_of_regions_to_set_material": %d,\n', solidno);
  fprintf(fid, '"material_regions": [\n');
  n = size(mat_ids,2);
  if(n==2)
    for a=1:solidno-1
      fprintf(fid, '\t["region", %d, %d, %d],\n', mat_ids(a, 1), mat_ids(a, 2), mat_ids(a, 2));
    end
    fprintf(fid, '\t["region", %d, %d, %d]\n],\n', mat_ids(solidno, 1), mat_ids(solidno, 2), mat_ids(solidno, 2));
  else
    for a=1:solidno
      fprintf(fid, '\t["region", %d, %d, %d],\n', a, mat_ids(a), mat_ids(a));
    end
    fprintf(fid, '\t["region", %d, %d, %d]\n],\n', solidno, mat_ids(solidno), mat_ids(solidno));
  end
  fprintf(fid, '"materials": [\n');
  
  for a=1: size(mat, 1)
    fprintf(fid, '{\n');
    fprintf(fid, '"name": "A",\n');
    fprintf(fid, '"ID": %d,\n', a-1);
    fprintf(fid, '"youngs_modulus": [%e, %e, %e],\n', mat(a, 2), mat(a, 3), mat(a, 4));
    fprintf(fid, '"shear_modulus": [%e, %e, %e],\n', mat(a, 5), mat(a, 6), mat(a, 7));
    fprintf(fid, '"poissons_ratio": [%e, %e, %e],\n', mat(a, 8), mat(a, 9), mat(a, 10));
    fprintf(fid, '"coefficient_of_thermal_expansion": [0.0, 0.0, 0.0],\n');
    fprintf(fid, '"sig": 1.0e+7,\n');
    fprintf(fid, '"strain_energy_function_dev": %d,\n', mat(a, 11));
    fprintf(fid, '"strain_energy_function_vol": %d,\n', mat(a, 12));
    fprintf(fid, '"fraction_of_heat_from_mechanical": 0.0\n');
    if(a<size(mat, 1))
      fprintf(fid, '},\n');
    else
      fprintf(fid, '}\n');
    end
  end
  fprintf(fid, ']\n');
  fprintf(fid, '}\n');
  fclose(fid);
  
  % write boundary condition file
  t3d_gid = {'"vertex"', '"curve"', '"surface"', '"region"', '"patch"', '"shell"', '"interface"'};
  fid = fopen([str,'_bc.json'], 'w');
  
  fprintf(fid, '{\n');
  if(~(size(BC, 1) ==0 || isnumeric(Bnd)))
    bc_cnt = 0;
    fprintf(fid, '"bc_data": [\n');
    for a=1: size(BC_ids.vertex, 1)
      fprintf(fid, '[%s, %d, %d, %d, %d]',t3d_gid{1}, ...
        BC_ids.vertex(a, 1), BC_ids.vertex(a, 2), BC_ids.vertex(a, 3), BC_ids.vertex(a, 4));
      
      bc_cnt = bc_cnt + 1;      
      if(bc_cnt<BC_ids.BCno)
        fprintf(fid, ',\n');
      else
        fprintf(fid, '\n');
      end
    end
    for a=1: size(BC_ids.edge, 1)
      fprintf(fid, '[%s, %d, %d, %d, %d]',t3d_gid{2}, ...
        BC_ids.edge(a, 1), BC_ids.edge(a, 2), BC_ids.edge(a, 3), BC_ids.edge(a, 4));

      bc_cnt = bc_cnt + 1;        
      if(bc_cnt<BC_ids.BCno)
        fprintf(fid, ',\n');
      else
        fprintf(fid, '\n');
      end
    end
    
    for a=1: size(BC_ids.face, 1)
      fprintf(fid, '[%s, %d, %d, %d, %d]',t3d_gid{5}, ...
        BC_ids.face(a, 1), BC_ids.face(a, 2), BC_ids.face(a, 3), BC_ids.face(a, 4));

      bc_cnt = bc_cnt + 1;        
      if(bc_cnt<BC_ids.BCno)
        fprintf(fid, ',\n');
      else
        fprintf(fid, '\n');
      end
    end
  end

  if(numel(disp)>0)
      fprintf(fid, '],\n');
      fprintf(fid, '"replacements": [ 0.0, 0.0, ');
      if(numel(disp)>1)
          fprintf(fid, '%e, ', disp(1:end-1));
      end
      fprintf(fid, '%e, 0.0, 0.0, 0.0]\n', disp(end));
  else
      fprintf(fid, ']\n');
  end

  fprintf(fid, '}\n');
  fclose(fid);
end

function BC_ids = get_bounary_ids4BC(Bnd, BC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BC_ids = get_bounary_ids4BC(Bnd, BC)
%
% compute boundary condition ids
% input:
%   Bnd     : bunndary info
%   BC      : boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 21-Aug-2018
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FF = [];
  EE = [];
  VV = [];
  F = {Bnd.face.LX,   Bnd.face.HX, ...
    Bnd.face.LY,   Bnd.face.HY, ...
    Bnd.face.LZ,   Bnd.face.HZ};
  E = {Bnd.edge.LX,   Bnd.edge.HX, ...
    Bnd.edge.LY,   Bnd.edge.HY, ...
    Bnd.edge.LZ,   Bnd.edge.HZ};
  V = {Bnd.vertex.LX, Bnd.vertex.HX, ...
    Bnd.vertex.LY, Bnd.vertex.HY, ...
    Bnd.vertex.LZ, Bnd.vertex.HZ};
  
  for a=1: size(BC, 1)
    id = BC(a, 1);
    no = numel(F{id});
    FF = [FF; F{id}(:), BC(a, 2)*ones(no, 1), BC(a, 3)*ones(no, 1), BC(a, 4)*ones(no, 1)];
    no = numel(E{id});
    EE = [EE; E{id}(:), a*(BC(a, 2)~=0)*ones(no, 1), a*(BC(a, 3)~=0)*ones(no, 1), a*(BC(a, 4)~=0)*ones(no, 1)];
    no = numel(V{id});
    VV = [VV; V{id}(:), a*(BC(a, 2)~=0)*ones(no, 1), a*(BC(a, 3)~=0)*ones(no, 1), a*(BC(a, 4)~=0)*ones(no, 1)];
  end
  
  list = VV(:, 1);
  [~, Io, Ic] = unique(list);
  nV = VV(Io, :);
  
  BCid_2 = abs(VV(:, 2));
  BCid_3 = abs(VV(:, 3));
  BCid_4 = abs(VV(:, 4));
  
  J = BCid_2~=0;
  nV(Ic(J), 2) = BC(BCid_2(J), 2);
  J = BCid_3~=0;
  nV(Ic(J), 3) = BC(BCid_3(J), 3);
  J = BCid_4~=0;
  nV(Ic(J), 4) = BC(BCid_4(J), 4);
  
  list = EE(:, 1);
  [~, Io, Ic] = unique(list);
  
  nE = EE(Io, :);
  BCid_2 = abs(EE(:, 2));
  BCid_3 = abs(EE(:, 3));
  BCid_4 = abs(EE(:, 4));
  
  J = BCid_2~=0;
  nE(Ic(J), 2) = BC(BCid_2(J), 2);
  J = BCid_3~=0;
  nE(Ic(J), 3) = BC(BCid_3(J), 3);
  J = BCid_4~=0;
  nE(Ic(J), 4) = BC(BCid_4(J), 4);
  
  BC_ids.BC = BC;
  BC_ids.vertex = nV;
  BC_ids.edge = nE;
  BC_ids.face = FF;
  BC_ids.BCno = size(nV, 1) + size(nE, 1) + size(FF, 1);
end

