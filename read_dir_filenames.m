function filenames = read_dir_filenames(target_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenames = read_dir_filenames(target_dir)
% read list of file names in a folder whose address is given in target_dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 27-Mar-2012
% by Sangmin Lee, Ph.D.
% Post-doctoral Research Fellow
% Advanced Materials Systems Laboratory
% 1000 Victors Way, Ann Arbor, MI 48108
% Department of Mechanical Engineering
% University of Michigan, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    crt_dir = pwd;
    cd(target_dir);
    path = what;
    dirinf = dir(path.path);
    dirnames = {dirinf.name};
    fnames = {};
    cnt = 0;
    for a=1: numel(dirnames)
        fname = dirnames{a};
        if(fname(1)~='.')
            cnt = cnt + 1;
            fnames{cnt, 1} = fname;
        end
    end
    filenames = fnames;
    cd(crt_dir);
end