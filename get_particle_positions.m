function get_particle_positions(fn_in, fn_out, cut_box)
%% This script read packing data from Rocpack output
 % It supports only sphere particles and box boundary

plot_particles = 0;

%% count number of particles
fid = fopen(fn_in, 'r');
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

%%
fid = fopen(fn_in, 'r');
n=0; % temporal counting for reading variables
%% read scaling factor
scailing = 0.0;
while(n==0)
  if(feof(fid))
    break;
  end
  line = fgetl(fid);
  [scailing, n] = sscanf(line, 'set scaling = %f\n');
end

% check scale needs to be 1 otherwise there is a problem
if(abs(1.0-scailing)>1.0e-3)
  fprintf('Boundary sized is modified, scaling = %e\n', scailing);
end
%% read packing fraction
n=0;
packing_fraction = 0.0;
while(n==0)
  if(feof(fid))
    break;
  end
  line = fgetl(fid);
  [packing_fraction, n] = sscanf(line, 'set packing_fraction =  %f\n');
end
%% read boundary

n=0;
box = [1.0, 1.0, 1.0];
while(n==0)
  if(feof(fid))
    break;
  end
  line = fgetl(fid);
  [box, n] = sscanf(line, 'boundary { box %f %f %f }\n');
end

%% read sphere particles
P = zeros(pno, 3); % translate
R = zeros(pno, 4); % rotate
S = zeros(pno, 1); % scale
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
    [p, ~] = sscanf(line, '	translate <%f, %f, %f>\n');
    P(count, :) = p;
    line = fgetl(fid);
    [r, ~] = sscanf(line, '	rotate < %f, %f, %f> %f\n');
    R(count, :) = r;
    line = fgetl(fid);
    [s, ~] = sscanf(line, '	scale %f color %f %f %f tag %f\n');
    S(count) = s(1);    
   end
end
fclose(fid);
%% rmove particles outside of the boundary

%% expand z direction
maxz = max(P(:, 3));
minz = min(P(:, 3));
P(:, 3) = (box(3)-S*2.05)./(maxz - minz).*P(:, 3);

if(nargin<3)
  cut_box = [0, box(1)/2
             0, box(2)/2
             -box(3)/3 box(3)/3];
end

list = sphere_box_overlap_test([P, S], cut_box);

P(list==0, :) = [];
S(list==0, :) = [];
R(list==0, :) = [];

pno = size(P, 1);
%% print T3d2Rpcpack file
fid = fopen(fn_out, 'w');
%fprintf(fid, '# file name : %s\n', fn_out);
%fprintf(fid, '# number of particles: %d\n', pno);
%fprintf(fid, '# packing fraction: %e\n', packing_fraction);
%fprintf(fid, '# scaling: %e\n', scailing);
%fprintf(fid, '# First pariticle''s radius: %e\n', S(1));

fprintf(fid, 'MinDist 1.0e-6\n');
fprintf(fid, 'RefinementDist 1.0e-6\n');
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

for b=1: pno
  fprintf(fid, 'Particle %d %e %e %e %e %d %d\n', b, P(b, 1),...
                                                     P(b, 2),...
                                                     P(b, 3),...
                                                     S(b, 1),...
                                                     0,0); % material IDs are all 0
end
fprintf(fid, '\n');
fprintf(fid, 'End\n');

fclose(fid);
%% temporal plot particles
if(plot_particles)
  %%
  box = cut_box(:,2) - cut_box(:, 1);
    
  bnd_box1 = [-box(1)/2, -box(2)/2, -box(3)/2;  box(1)/2, -box(2)/2, -box(3)/2;
               box(1)/2,  box(2)/2, -box(3)/2; -box(1)/2,  box(2)/2, -box(3)/2;
              -box(1)/2, -box(2)/2, -box(3)/2; -box(1)/2, -box(2)/2,  box(3)/2;
               box(1)/2, -box(2)/2,  box(3)/2;  box(1)/2,  box(2)/2,  box(3)/2;
              -box(1)/2,  box(2)/2,  box(3)/2; -box(1)/2, -box(2)/2,  box(3)/2];
            
  bnd_box2 = [-box(1)/2,  box(2)/2, -box(3)/2; -box(1)/2,  box(2)/2,  box(3)/2];
  bnd_box3 = [ box(1)/2,  box(2)/2, -box(3)/2;  box(1)/2,  box(2)/2,  box(3)/2];
  bnd_box4 = [-box(1)/2,  box(2)/2, -box(3)/2; -box(1)/2,  box(2)/2,  box(3)/2];
  bnd_box5 = [ box(1)/2, -box(2)/2, -box(3)/2;  box(1)/2, -box(2)/2,  box(3)/2];
  
  figure; hold on;
  view(42,14);
  set(gcf,'color','w')
  
  plot3(bnd_box1(:, 1), bnd_box1(:, 2), bnd_box1(:, 3), '-k');
  plot3(bnd_box2(:, 1), bnd_box2(:, 2), bnd_box2(:, 3), '-k');
  plot3(bnd_box3(:, 1), bnd_box3(:, 2), bnd_box3(:, 3), '-k');
  plot3(bnd_box4(:, 1), bnd_box4(:, 2), bnd_box4(:, 3), '-k');
  plot3(bnd_box5(:, 1), bnd_box5(:, 2), bnd_box5(:, 3), '-k');
  
  [x,y,z] = sphere(20);
  for k= 1: size(P, 1)
    d = P(k,:);
    x2 = x*S(k)+d(1);
    y2 = y*S(k)+d(2);
    z2 = z*S(k)+d(3);
    figs = surf(x2,y2,z2,ones(size(z2))*S(k));
  end
  
  set(figs,'FaceAlpha',1);
  set(figs,'DiffuseStrength', 0.9);
  set(figs,'SpecularStrength',0.1);
  set(gca,'dataaspectratio',[1 1 1])
  hold off;
end
end



