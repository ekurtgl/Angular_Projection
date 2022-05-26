clc; clear; close all;

datapath = '/mnt/HDD04/Projection_data/Cascade_AWR2243/*';
% outpath = '/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/RDCs/';
outpath = '/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/microDoppler/';
subfolds = dir(datapath);
subfolds(1:3) = [];

for i = 1:length(subfolds)
        tic
        disp([int2str(i) '/' int2str(length(subfolds))]);
        fname = [subfolds(i).folder '/' subfolds(i).name];
%         fOut = [outpath subfolds(i).name '.mat'];
        fOut = [outpath subfolds(i).name '.png'];
        RDC = RDC_extract_cascade_AWR2243(fname);
        RDC_to_microDopp_cascade( RDC, fOut, 1)
        RDC_to_microDopp_cascade( RDC, fOut, 0)
%         RDC_to_microDopp_cascade_2Dfft( RDC, fOut, 0)
        toc
%         tic
%         disp('Saving...');
%         save(fOut, 'RDC', '-v7.3');
%         toc

end




