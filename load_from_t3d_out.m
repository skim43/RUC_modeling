function [grid, node_info, elem_info] = load_from_t3d_out(fn_in, load_elem)
  if(nargin<2)
    load_elem = 1;
  end

  fprintf('read [%s] ...', fn_in);
  fid = fopen(fn_in, 'r');
  if(fid==-1)
    fprintf('%s is not found\n', str);
    return;
  end
  fgets(fid);
  line = fgets(fid);
  num = str2num(line);
  
  nodeno = num(1);
  elemno = num(5);
  if(elemno==0)
    elemno = num(end);
  end
  
  fgets(fid);
  temp = fscanf(fid, '%f', [7, nodeno]);
  temp = temp';
  
  node = temp(:, 2:4);
  node_info = temp(:, 5:7);
  
  temp = [];
  elem = [1 2 3 4];
  elem_info = [];
 
  if(load_elem)
    % determine number of ids per element
    fgets(fid);
    fgets(fid);
    line = fgets(fid);
    num = str2num(line);
    numno = numel(num);
    nne = 4; % for tetrahedron. Other types are not supported.
    temp = fscanf(fid, '%d', [numno, elemno-1]);
    temp = temp';
    elem = [num(2: nne+1); temp(:, 2:nne+1)];
    elem_info = [num(nne+2:end); temp(:, nne+2:end)];
  end
  
  grid.NODE = node;
  grid.ELM  = elem;
  grid.nodeno = nodeno;
  grid.elmno  = elemno;
  grid.nsd = size(node, 2);

  fclose(fid);
  fprintf('done.\n');
end
