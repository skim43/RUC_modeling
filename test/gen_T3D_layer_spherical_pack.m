function gen_T3D_layer_spherical_pack(fn_rp2t3d, fn_t3d, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gen_T3D_layer_spherical_pack(fn_rp2t3d, fn_t3d, t)
% 
% This function will generate T3D file (spherical particle pack with layer.
% fn_rp2t3d format (RoP2T3d input file format) is needed.
% T3D and T3DAssemble classes are used in the folder to read T3D files and
% assemble the geometries.
%
% e.g.
% gen_T3D_layer_spherical_pack('example.pack/pack.rp2t3d', 'example.pack/pack.t3d', 0.01); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \param[in] fn_rp2t3d Rocpack to T3D file
% \param[in] fn        T3D file (this file will be created)
% \param[in] t         layer thickness
%
% 16-Feb-2020
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  fn = write_inner_particle_Ropack(fn_rp2t3d, t);
  system(sprintf('/cswarm/tools/bin/RoP2T3d -i %s -o %s', fn_rp2t3d, fn_rp2t3d));
  system(sprintf('/cswarm/tools/bin/RoP2T3d -i %s -o %s', fn,        fn));
  t3d_o = T3D;
  t3d_i = T3D;
  t3d   = T3DAssemble;
  t3d_o.read_t3d_file(sprintf('%s.t3d', fn_rp2t3d));
  t3d_i.read_t3d_file(sprintf('%s.t3d', fn));
  t3d.assemble_two_packs(t3d_o, t3d_i);
  t3d.write_t3d_file(fn_t3d);
end

function fn = write_inner_particle_Ropack(fn_rp2t3d, t)
  fn = sprintf('%s.tmp', fn_rp2t3d);
  rp = read_rp2t3d(fn_rp2t3d);
  fid = fopen(fn, 'w');
  if(~fid)
    fprintf('Cannot create [%s]\n', fn_rp2t3d);
    return;
  end
  
  for ia = 1: numel(rp.head)
    fprintf(fid, '%s\n', rp.head{ia});
  end

  for ia = 1: size(rp.P, 1)
    fprintf(fid, 'Particle %d ', rp.P(ia, 1));
    fprintf(fid, '%g ', rp.P(ia, 2:4));
    fprintf(fid, '%g ', rp.P(ia, 5)-t);
    fprintf(fid, '%d ', rp.P(ia, 6:7));
    fprintf(fid, '\n');
  end

  for ia = 1: numel(rp.tail)
    fprintf(fid, '%s\n', rp.tail{ia});
  end

  fclose(fid);
end

function rp = read_rp2t3d(fn_rp2t3d)
  fid = fopen(fn_rp2t3d, 'r');

  if(~fid)
    fprintf('Cannot open [%s]\n', fn_rp2t3d);
    return;
  end

  pno     = 0;
  headno  = 0;
  lineno  = 0;
  while(~feof(fid))
    l = fgetl(fid);
    lineno = lineno + 1;
    if(numel(l)>0)
      s = split(l);
      if(strcmp(s{1}, 'Particle'))
        pno = pno + 1;
      end
    end

    if(pno==0)
      headno = headno + 1;
    end
  end
  
  frewind(fid);
  rp.head = cell(headno, 1);

  for ia = 1: headno
    rp.head{ia} = fgetl(fid);
  end

  rp.P = zeros(pno, 7);
  for ia = 1: pno
    l = fgetl(fid);
    rp.P(ia, :) = str2num(l(9:end));
  end

  rp.tail = cell(lineno - headno - pno, 1);
  cnt = 0;
  while(~feof(fid))
    cnt = cnt + 1;
    rp.tail{cnt} = fgetl(fid);
  end

  fclose(fid);
end

