function write_regions(fn_in, fn_out) 
%fn_in='/scratch365/cswarm/kramos/interface/generate_pgfem3d_input_stack_with_interfaces_local/input/t3d/particles.t3d';
%fn_out='/scratch365/cswarm/kramos/interface/generate_pgfem3d_input_stack_with_interfaces_local/input/t3d/test.out.regions';
%num_prop=3; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write regions input file for material and boundary condition
% input:
%   fn_in : file name base. {filename_base}_mat.json and {filename_base}_bc.json files will be created
% output:
%   fn_out: file containing regions information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13-May-2020
% by Katherine Ramos
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fn_in, 'r');
rno = 0;
while(~feof(fid))
  line = fgetl(fid);
  count = strfind(line,'region');
  if(numel(count)>=1)
    if(numel(count)==1)
      rno = rno + 1;
    end
  end 
end
fclose(fid);
fprintf('The number of regions is = %i\n',rno); 


%%Section 2 : read shape ID
fid = fopen(fn_in, 'r');
region_id = zeros(rno, 1);
count = 1;
while(count<=rno)
    line = fgetl(fid);
    [r, ~] = sscanf(line, 'region %d \n');
    if ~isempty(r)
        region_id(count, 1) = r;
        count=count+1;
    end
end
fclose(fid);
fprintf('The region id = %i\n',region_id); 

%assign property to vector counter to define regions
fid=fopen(fn_in,'r');
text=textscan(fid,'%s','Delimiter','','endofline','');
text=text{1}{1};
fid=fclose(fid);

property_id=regexp(text,'property[\s]+(\d+)','tokens');
property_id=str2double([property_id{:}]);

%pno=num_prop; 
%prop_id=zeros(pno,1);
%start=0; 
%for i=1:pno
%    prop_id(i)=start;
%    start=start+1;
%end 
% diff_regions=(rno-1)/2;
% property_id=zeros(rno,1);
% property_id(1:diff_regions)=prop_id(2);
% property_id(diff_regions+1)=prop_id(1);
% property_id(diff_regions+2:rno)=prop_id(3);

% fprintf('The properties id = %i\n',property_ids); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --     print out.region   file   -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% write the regions file 
region_type=4; 
fid = fopen(fn_out, 'w');
fprintf(fid, '%d\n', rno);
for i=1:rno
fprintf(fid, '%d %d %d %d \n',region_type,region_id(i),property_id(i),property_id(i));
end
fclose(fid); 
exit 
