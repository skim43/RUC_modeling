function q = sphere_box_overlap_test(p, out_box, tol)
%% sphere-box overlap test
% input: p = [cx, cy, cz, r
%              :
%              :
%            ];
%        out_box = [xmin xmax
%                   ymin ymax
%                   zmin zmax];
%        tol = tolarence when the decision is made, default = 1.0e-12
%
% output: q = [ 0 or 1
%                  :
%             ];
% example:
% box = [0 1; 0 1; 0 1];
% p = [-0.51  1.01, 0, 0.5
%      -0.49 -0.49, 0, 0.5];  
% q = sphere_box_overlap_test(p, box)
% q =
%     1
%     0
% %plot results
% ang=0:0.01:2*pi; 
% xp=0.5*cos(ang);
% yp=0.5*sin(ang);
% plot(xp-0.51, yp+1.01), hold on;
% plot(xp-0.49, yp-0.49), hold on;
% plot([0, 0, 1, 1, 0], [0, 1, 1, 0, 0], 'k'); hold off; axis equal;

  if(nargin<3)
    tol = 1.0e-12;
  end

  [pno, test_size] = size(p);
  if(test_size~=4)
    fprintf('Wrong input on particle positions and radii, p = %d x %d\n', pno, test_size);
  end
  q = zeros(pno, 1);
  
  box = [out_box(1,1), out_box(2,1), out_box(3,1)
         out_box(1,2), out_box(2,1), out_box(3,1) 
         out_box(1,2), out_box(2,2), out_box(3,1)
         out_box(1,1), out_box(2,2), out_box(3,1)
         out_box(1,1), out_box(2,1), out_box(3,2)
         out_box(1,2), out_box(2,1), out_box(3,2) 
         out_box(1,2), out_box(2,2), out_box(3,2)
         out_box(1,1), out_box(2,2), out_box(3,2)];

  patch = {box((4:-1:1),:), ...
           box((5:8),:), ...
          [box((1:2),:); box((5:6),:)], ...
          [box((2:3),:); box((6:7),:)], ...
          [box((3:4),:); box((7:8),:)], ...
          [box([4,1],:); box([8,5],:)]};

  patchno = numel(patch);        
  for ip = 1: pno        
    c = p(ip, 1:3);
    r = p(ip, 4);

    is_bounded = true;
    for ia=1: 3
      if((c(ia)<out_box(ia, 1)) || (c(ia)>out_box(ia, 2)))
        is_bounded = false;
      end
    end
  
    if(is_bounded)
      q(ip) = 1;
      continue;
    end    
    
    dist = 1.0e+10;
    for ia = 1: patchno
      dist = min([dist, distance_from_a_patch(patch{ia}, c)]);
    end % for ia

    if(dist < r-tol)
      q(ip) = 1;
    end  
  end
end