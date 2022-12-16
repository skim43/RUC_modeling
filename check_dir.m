function status = check_dir(dir_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% status = check_dir(dir_name)
% dir_name: target dir
% if this directory exists, status will be 1 otherwise 0
% if status is 0, a folder which has name as dir_name will be created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11-Oct-2011
% by Sangmin Lee, Ph.D.
% Post-doctoral Research Fellow
% Advanced Materials Systems Laboratory
% 1000 Victors Way, Ann Arbor, MI 48108
% Department of Mechanical Engineering
% University of Michigan, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirinf = dir(dir_name);

status = 1;
if(numel(dirinf)==0)
  mkdir(dir_name);
  status = 0;
end
end