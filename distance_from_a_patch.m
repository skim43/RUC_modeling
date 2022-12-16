function dist = distance_from_a_patch(p, x_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = distance_from_a_patch(p, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12-Jun-2016
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: p : list of points constructing a pahch
%             e.g) p = [0 1 0; 1 1 0; 1 0 0; 0 0 0];
% input: x_in : position to be checked
% output: dist : shortest distance form a patch p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how to use:
% p = [0 1 0; 1 1 0; 1 0 0; 0 0 0];
% x = [0.5, 0.5, 0.5];
% dist = distance_from_a_patch(p, x);
% dist =
%    0.5000
  x = [x_in(1); x_in(2); x_in(3)];
  vtxno = size(p, 1);
  %% compute projection of x on the patch: y 
  v0 = x - p(1, :)';
  v1 = p(2, :)' - p(1, :)';
  v2 = p(3, :)' - p(2, :)';
  n = cross(v1, v2);
  n = n./sqrt(n'*n);
  d = v0'*n;
  y = x - d*n;
  
  %% deside the y on the patch by computing area 
  % compute original patch area
  q = [p; p(1, :)];
  v = q(2:end, :) - q(1:end-1, :);
  A0 = 0;
  for a = 1: vtxno-1
    v1 = p(a+1, :) - p(1, :);
    v2 = v(a+1, :);
    A0 = A0 + 0.5*(n'*cross(v1', v2'));
  end
  
  % compute patch area with y
  %%
  z = [y(1) - q(:, 1), y(2) - q(:, 2), y(3) - q(:, 3)];
  A = 0;
  for a = 1: vtxno
    v1 = v(a, :);
    v2 = z(a+1, :);
    A = A + abs(0.5*(n'*cross(v1', v2')));
  end  
  
  if(A>A0) % out of the patch
    % compute distance on the line
    dist = 1.0e+15;
    for ia=1:vtxno
      line = [q(ia+1, :); q(ia, :)];
      dist = min([dist, distance_from_a_line(line, x)]);
    end
  else % on the patch
    dist = abs(d);
  end
end