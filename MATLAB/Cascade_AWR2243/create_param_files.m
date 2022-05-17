clc; clear; close all;

paramfile = 'test.mmwave.json';
src = ['/mnt/HDD04/Projection_data/Cascade_AWR2243/cascadeParams/' paramfile];

fol_path = '/mnt/HDD04/Projection_data/Cascade_AWR2243';
folders = dir(fol_path);
folders(1:3) = [];


for i = 1:length(folders)
        disp([int2str(i) '/' int2str(length(folders))]);
        subfol = folders(i).name;
        dst = [folders(i).folder '/' folders(i).name '/' folders(i).name '.mmwave.json'];
        copyfile(src, dst);
end






