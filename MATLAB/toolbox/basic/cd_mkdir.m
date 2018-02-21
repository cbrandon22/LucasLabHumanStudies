function [] = cd_mkdir(pathname)
if ~exist(pathname,'dir')
    mkdir(pathname)
end
cd(pathname)
    