function resize_particles(fn_in, fn_out, d_increase) 
%% This script reads packing data from Rocpack output, resizes "deflates" ...
%  particles to account for enforced inter-particle distances
%
%% Count number of particles

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

%Section 1 : read scaling & packing
fid = fopen(fn_in, 'r');
n=0; % temporal counting for reading variables
scaling = 0.0;
while(n==0)
  if(feof(fid))
    break;
  end
  line = fgetl(fid);
  [scaling, n] = sscanf(line, 'set scaling = %f\n');
end

%%Section 2 : read seed
seed = 0.0;
while(n==0)
  if(feof(fid))
    break;
  end
  line = fgetl(fid);
  [seed, n] = sscanf(line, 'set seed  = %f\n');
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

%%Section 3 : read boundary information
n=0;
box = [1.0, 1.0, 1.0];
while(n==0)
  if(feof(fid))
    break;
  end
  line = fgetl(fid);
  [box, n] = sscanf(line, 'boundary { box %f %f %f periodic }\n');
end

%%Section 4 : details from particles

P = zeros(pno, 3); % translate "position"
R = zeros(pno, 4); % rotate
S = zeros(pno, 1); % scale
C = zeros(pno, 3); % color
T = zeros(pno, 1); % tag
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
    C(count,:) = [s(2),s(3),s(4)];
    T(count) = s(5); 
   end
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --   print Rocpack2Stat3D file   -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Resize particles to account for enforced inter-particle spacing

S = S./(1.0+d_increase); 

%% print Rocpack2Stat3D file
fid = fopen(fn_out, 'w');

fprintf(fid, 'set scaling = %f;\n', scaling);
fprintf(fid, '\n');
fprintf(fid, 'set packing_fraction = %f;\n', packing_fraction);
fprintf(fid, '\n');
fprintf(fid, 'boundary { box %f %f %f periodic}\n\n', box(1), box(2), box(3));
fprintf(fid, '\n');
fprintf(fid, '\n');
for b=1:pno
    fprintf(fid, 'sphere{\n');
    fprintf(fid, '	translate <%f, %f, %f>\n', P(b,1),P(b,2),P(b,3));
    fprintf(fid, '	rotate < %f, %f, %f> %f\n',R(b,1),R(b,2),R(b,3), R(b,4));
    fprintf(fid, '	scale %f color %f %f %f tag %.1f\n',S(b,1),C(b,1), C(b,2), C(b,3),T(b,1));
    fprintf(fid, '}\n\n');
end
fclose(fid);

