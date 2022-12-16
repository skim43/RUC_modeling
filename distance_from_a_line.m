function dist = distance_from_a_line(L, x_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = distance_from_a_line(L, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 12-Jun-2016
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: L : list of points constructing a line
%             e.g) p = [0 0 0; 1,0,0];
% input: x : position to be checked
% output: dist : shortest distance form a patch p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how to use:
% L = [0 0 0; 1 0 0];
% x = [0.5, 0.5, 0];
% dist = distance_from_a_line(L, x);
% dist =
%    0.5000

  dist = [];
  [p, q] = size(L);
  s = numel(x_in);

  if(q~=s)
    fprintf('dimension is not consistance: L(%dx%d) and x_in(%d)\n', p, q, s);
    return;
  end
  
  dim = q;
  if(0 > dim || dim > 3)
    fprintf('dimension (%d) must be in [1,3]\n', dim);
    return;
  end
  
  line = L;
  if(dim == 1)
    line = [L(1), 0, 0
            L(2), 0, 0];
    x = [x_in(1), 0, 0];
  end
  
  if(dim == 2)
    line = [L(1, :), 0
            L(2, :), 0];
    x = [x_in(1), x_in(2), 0];
  end

  if(dim==3)
     x = [x_in(1), x_in(2), x_in(3)];
  end

  a = line(2, :) - line(1, :);
  v = x - line(1, :);
  ma = sqrt(a*a');
  mm = (v*a')/ma;
  m = a*mm/ma;
  
  if(m*a'<0)
    dist = norm(v);
  else
    if(mm<norm(a))
      dist = norm(v - m);
    else
      dist = norm(x-line(2, :));
    end
  end
end