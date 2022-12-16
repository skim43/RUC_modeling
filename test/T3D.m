classdef T3D < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T3D class
% Provide data structure for T3D entities and means of reading and writing
% T3D files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  properties
    E;
    map;
    t3d;
    tol;
  end
  methods 
    function obj = T3D
      obj.t3d = interpret;
      obj.tol = 1.0e-12;
    end
    function read_t3d_file(obj, fn_t3d)
    
      fid = fopen(fn_t3d, 'r');
      if(~fid)
        fprintf('Cannot open [%s]\n', fn_t3d);
        return;
      end
      
      % set size of t3d entities
      [no, limit] = count_number_of_entities(fid);
      obj.map = cell(obj.t3d.ids(end), 1);
      obj.E   = cell(obj.t3d.ids(end), 1);
      for ia = 1: obj.t3d.ids(end)
        if(limit(ia)>0)
          obj.map{ia} = sparse(limit(ia), 1);
        end
        if(no(ia)>0)
          obj.E{ia} = cell(no(ia), 1);
        end
      end
      
      % read t3d file and set entities
      cnt = zeros(obj.t3d.ids(end), 1);
      while(~feof(fid))
        l = fgetl(fid);
        o = interpret(l, fid);
        if(o.out.valid)
          cnt(o.out.t3d) = cnt(o.out.t3d) + 1;           % count each entity to save in order
          obj.E{o.out.t3d}{cnt(o.out.t3d), 1} = o.out;   % save read entity in order
          obj.map{o.out.t3d}(o.out.id) = cnt(o.out.t3d); % set map id. T3D has arbtrary ID numbers
                                                         % indexing using T3D id is enabled by this map
        end
      end
      fclose(fid);
      for ia = 1: numel(obj.E{4})-1
        obj.E{4}{ia}.property = obj.E{4}{ia}.property + 1;
      end
    end

    function write_t3d_file(obj, fn_t3d)
      fid = fopen(fn_t3d, 'w');
      write_T3D_vertex( fid, obj.E{1});
      write_T3D_curve(  fid, obj.E{2});
      write_T3D_surface(fid, obj.E{3});
      write_T3D_patch(  fid, obj.E{5});
      write_T3D_shell(  fid, obj.E{6});
      write_T3D_region( fid, obj.E{4});
      fclose(fid);
    end
    function min_d = compute_smallest_edge(obj)
      min_d = 1.0e+15;
      vno = numel(obj.E{1});
      for ia = 1: vno-1
        if(~obj.E{1}{ia}.valid); continue; end
        for ib = ia+1: vno
          if(~obj.E{1}{ib}.valid); continue; end
          d = norm(obj.E{1}{ia}.e - obj.E{1}{ib}.e);
          if(d<min_d)
            min_d = d;
          end
        end
      end
    end    
  end
end

function [no, limit] = count_number_of_entities(fid)
  % count number of T3D entities
  e = interpret;
  no = zeros(e.ids(end), 1);
  limit = zeros(e.ids(end), 1);

  while(~feof(fid))
    l = fgetl(fid);
    if(numel(l)<1); continue; end
    if(l(1)=='#');  continue; end

    s = split(l);

    for ia = 1: e.ids(end)
      if(strcmp(s{1}, e.entities{ia}))
        no(ia) = no(ia) + 1;
        limit(ia) = max(limit(ia), str2num(s{2}));
        break;
       end
    end
  end
  frewind(fid);
end
function out = interpret(l, fid)
  % set basic data structure for T3D entities 
  out.out.valid = false;
  out.ids = 1:7;
  out.entities = {'vertex' 'curve' 'surface' 'region' 'patch' 'shell' 'interface'};
  
  out.out.t3d    = [];
  out.out.id     = [];
  out.out.name   = '';
  out.out.option = '';
  out.out.order  = [];
  out.out.e      = [];
  out.out.sub    = [];
  out.out.n      = [];
  out.out.mirror = [];
  out.out.property = [];  

  if(nargin<1);   return; end
  if(numel(l)<1); return; end
  if(l(1)=='#');  return; end
  
  
  tmp = split(l);
  s = {};

  newline = true;
  while(newline)
    newline = false;
    end_id = [];
    for ia = 1: numel(tmp)
      if(strcmp(tmp{ia}, '\'))
        newline = true;
        tmp(ia:end) = [];
        break;
      end
    end
    if(newline)
      s = {s{:}, tmp{:}};
      l = fgetl(fid);
      tmp = split(l);      
    end
  end

  s = {s{:}, tmp{:}};
  
  for ia = 1: out.ids(end)
    % if line conttains T3D entity name
    % read all information
    if(strcmp(s{1}, out.entities{ia}))
      % set basic data structure for T3D entity
      out.out.valid = true;
      out.out.t3d = ia;
      out.out.id  = str2double(s{2});
      out.out.name = s{1};

      switch(out.entities{ia})
        case 'vertex'
          [out.out.e, out.out.option] = read_T3D_vertex(s);
        case 'curve'
          [out.out.e, out.out.option, out.out.order] = read_T3D_curve(s);
        case 'surface'
          [out.out.e, out.out.option, out.out.order] = read_T3D_surface(s);
        case 'region'
          [out.out.e, out.out.property] = read_T3D_region(s);      
        case 'patch'
          [out.out.e, out.out.n, out.out.mirror] = read_T3D_patch(s);
        case 'shell'
          out.out.e = read_T3D_shell(s);
      end
      break 
    end
  end
  % read if curve or surface have control points
  if(numel(out.out.order)>0)
    out.out.sub = cell(numel(out.out.order), 2);
    for ia = 1: numel(out.out.order)
      out.out.sub{ia, 1} = get_ctr_pt(split(fgetl(fid)), out.out.name);
      if(out.out.order(ia)==4)
        out.out.sub{ia, 2} = get_ctr_pt(split(fgetl(fid)), out.out.name);
      end
    end
  end
end
function [e, o] = read_T3D_vertex(s)
  o = '';
  e = [str2double(s{4}), str2double(s{5}), str2double(s{6})];
  if(numel(s)>6)
    o = sprintf('%s ', s{7:end});
  end
end
function write_T3D_vertex(fid, E)
  fprintf(fid, '# Vertex\n');
  for ia = 1: numel(E)
    if(~E{ia}.valid); continue; end
    fprintf(fid, 'vertex %d xyz ', E{ia}.id);
    fprintf(fid, '%g ',            E{ia}.e);
    fprintf(fid, '%s\n',           E{ia}.option);
  end
end
function [e, o, order] = read_T3D_curve(s)
  o = '';
  order = [];
  have_ov = 0;
  for ib = 3: numel(s)
    if(strcmp(s{ib}, 'order'))
      have_ov = have_ov + 2;
      order = str2double(s{ib+1});
    end
    if(strcmp(s{ib}, 'vertex'))
      have_ov = have_ov + 3;
      e = [str2double(s{ib+1}), str2double(s{ib+2})];
    end
  end
  if(numel(s)-have_ov > 2)
    o = sprintf('%s ', s{3+have_ov:end});
  end
end
function write_T3D_curve(fid, E)
  fprintf(fid, '# Curve\n');
  for ia = 1: numel(E)
    if(~E{ia}.valid); continue; end
    fprintf(fid, 'curve %d ', E{ia}.id);
    if(numel(E{ia}.order)>0)
      fprintf(fid, 'order '); fprintf(fid, '%d ', E{ia}.order);
    end
    fprintf(fid, 'vertex ');
    fprintf(fid, '%d ',              E{ia}.e);
    fprintf(fid, '%s\n',             E{ia}.option);
    print_cp(fid, E{ia});
  end
end
function [e, o, order] = read_T3D_surface(s)
  o = '';
  order = [];
  have_ov = 0;
  for ib = 3: numel(s)
    if(strcmp(s{ib}, 'order'))
      have_ov = have_ov + 3;
      order = [str2double(s{ib+1}), str2double(s{ib+2})];
    end
    if(strcmp(s{ib}, 'Curve') || strcmp(s{ib}, 'curve'))
      have_ov = have_ov + 5;
      e = [str2double(s{ib+1}), str2double(s{ib+2}), str2double(s{ib+3}), str2double(s{ib+4})];
    end
  end
  if(numel(s)-have_ov > 2)
    o = sprintf('%s ', s{3+have_ov:end});
  end
end
function write_T3D_surface(fid, E)
  fprintf(fid, '# Surface\n');
  for ia = 1: numel(E)
    if(~E{ia}.valid); continue; end
    fprintf(fid, 'surface %d ', E{ia}.id);
    if(numel(E{ia}.order)>0)
      fprintf(fid, 'order '); fprintf(fid, '%d ', E{ia}.order);
    end
    fprintf(fid, 'curve ');
    fprintf(fid, '%d ',  E{ia}.e);
    fprintf(fid, '%s\n', E{ia}.option);
    print_cp(fid, E{ia});
  end
end
function [e, p] = read_T3D_region(s)
  e = cell(3,1);
  have_sp = zeros(3, 2);
  ic = 0; id = 0;
  for ib = 3: numel(s)
    if(strcmp(s{ib}, 'size'))
      have_sp(ic, 2) = ib-1;
    end
    if(strcmp(s{ib}, 'boundary'))
      if(strcmp(s{ib+1}, 'surface'))
        id = 1;
      end
      if(strcmp(s{ib+1}, 'shell'))
        id = 2;
      end
      if(strcmp(s{ib+1}, 'patch'))
        id = 3;
      end
    else
      continue;
    end
    have_sp(id, 1) = ib+2;
    if(ic>0)
      have_sp(ic, 2) = ib-1;
    end
    ic = id;
  end
  for ic = 1: size(have_sp, 1)
    if(have_sp(ic, 1) > 0)
      e{ic} = str2num(sprintf('%s ', s{have_sp(ic, 1):have_sp(ic, 2)}));
    end
  end
  p = str2double(s{end});
end
function write_T3D_region(fid, E)
  fprintf(fid, '# Region\n');
  surf_names = {'surface', 'shell', 'patch'};

  for ia = 1: numel(E)
    cnt = 0;
    if(~E{ia}.valid); continue; end
    fprintf(fid, 'region %d ', E{ia}.id);
    for ib = 1: numel(E{ia}.e)
      fno = numel(E{ia}.e{ib});
      if(fno==0); continue; end
      fprintf(fid, 'boundary %s ', surf_names{ib});
      for ic = 1: fno
        cnt = cnt + 1;
        fprintf(fid, '%d ', E{ia}.e{ib}(ic));
        if(mod(cnt, 40)==0 && ic ~= fno); fprintf(fid, '\\\n'); end
      end
    end
    fprintf(fid, 'size def property %d\n', E{ia}.property);
  end
end
function [e, n, mirror] = read_T3D_patch(s)
  mirror = [];
  have_mirror = 0;
  for ib = 3: numel(s)
    n = [str2double(s{4}), str2double(s{5}), str2double(s{6})];
    if(strcmp(s{ib}, 'mirror'))
      have_mirror = have_mirror + 2;
      mirror = str2double(s{ib+1});
    end
  end
  e = str2num(sprintf('%s ', s{9:(end-2-have_mirror)}));
end
function write_T3D_patch(fid, E)
  fprintf(fid, '# Patch\n');
  for ia = 1: numel(E)
    if(~E{ia}.valid); continue; end
    fprintf(fid, 'patch %d normal %g %g %g boundary curve ', E{ia}.id, E{ia}.n);
    
    cno = numel(E{ia}.e);
    for ib = 1: cno
      fprintf(fid, '%d ', E{ia}.e(ib));
      if(mod(ib, 40)==0 && ib ~= cno); fprintf(fid, '\\\n'); end
    end    
    
    if(numel(E{ia}.mirror))
      fprintf(fid, 'mirror %d ', E{ia}.mirror);
    end
    fprintf(fid, 'size def\n');
  end
end
function e = read_T3D_shell(s)
  e = str2num(sprintf('%s ', s{4}, s{7:end-2}));
end
function write_T3D_shell(fid, E)
  fprintf(fid, '# Shell\n');
  for ia = 1: numel(E)
    if(~E{ia}.valid); continue; end
    fprintf(fid, 'shell %d bgsurface %d boundary curve ', E{ia}.id, E{ia}.e(1));
    fprintf(fid, '%d ', E{ia}.e(2:end));
    fprintf(fid, 'size def\n');
  end
end
function cp = get_ctr_pt(s, name)
  if(strcmp(name, 'surface'))
    cp = {[str2double(s{2}), str2double(s{3})]
          [str2double(s{5}), str2double(s{6}), str2double(s{7})]
           str2double(s{9})};
  else
    cp = { str2double(s{2})
          [str2double(s{4}), str2double(s{5}), str2double(s{6})]
           str2double(s{8})}; 
  end
end
function print_cp(fid, E)
  if(numel(E.order)>0)
    for ia = 1: numel(E.order)
      fprintf(fid, 'polygon '); fprintf(fid, '%d ', E.sub{ia, 1}{1});
      fprintf(fid, 'xyz ');      fprintf(fid, '%g ', E.sub{ia, 1}{2});
      fprintf(fid, 'weight %g\n',                   E.sub{ia, 1}{3});
      if(E.order(ia)==4)
        fprintf(fid, 'polygon '); fprintf(fid, '%d ', E.sub{ia, 2}{1});
        fprintf(fid, 'xyz ');      fprintf(fid, '%g ', E.sub{ia, 2}{2});
        fprintf(fid, 'weight %g\n', E.sub{ia, 2}{3});
      end
    end
  end
end
