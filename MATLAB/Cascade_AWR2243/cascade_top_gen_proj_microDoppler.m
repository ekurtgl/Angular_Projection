clc; clear; close all;

datapath = '/mnt/HDD04/Projection_data/Cascade_AWR2243/*';
% outpath = '/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/RDCs/';
outpath = '/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/projected_microDoppler/';
subfolds = dir(datapath);
subfolds(1:3) = [];
NPpF = 60;
antenna_azimuth_only = [1;2;3;4;17;18;19;20;33;34;35;5;6;7;8;21;22;23;24;37;38;39;40;53;54;55;56;69;70;71;72;85;86;87;88;101;
        102;103;104;117;118;119;120;133;134;135;9;10;11;12;13;14;15;16;29;30;31;32;45;46;47;48;61;62;63;64;77;78;79;80;93;94;
        95;96;109;110;111;112;125;126;127;128;141;142;143;144]; % retrieved from caliboj of TI's script, for 12 TXs

for i = 1:length(subfolds)
        tic
        disp([int2str(i) '/' int2str(length(subfolds))]);
        fname = [subfolds(i).folder '/' subfolds(i).name];
%         fOut = [outpath subfolds(i).name '.mat'];

        underScores = strfind(subfolds(i).name, '_');
        angStart = strfind(subfolds(i).name, 'ang');
        angles = subfolds(i).name(angStart+3:underScores(end)-1);
        plusses = strfind(angles, '+');
        minuses = strfind(angles, '-');
        all_ang_ids = sort([plusses minuses length(angles)+1]);
        if all_ang_ids(1) == 1
                angs = zeros(1, length(all_ang_ids) - 1);
                for id = 1:length(all_ang_ids) - 1
                        angs(id) = str2num(angles(all_ang_ids(id): all_ang_ids(id+1)-1)); 
                end
        else
                angs = zeros(1, length(all_ang_ids));
                angs(1) = str2num(angles(1: all_ang_ids(1)-1)); 
                for id = 1:length(all_ang_ids) - 1
                        angs(id+1) = str2num(angles(all_ang_ids(id): all_ang_ids(id+1)-1)); 
                end
        end
        
        if length(angs) == 1
                continue
        end
        RDC = RDC_extract_cascade_AWR2243(fname);
        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
        
        fOut = [outpath subfolds(i).name '.png'];
%         RDC_to_microDoppler_projection(RDC, NPpF, fOut)
        RDC_to_microDoppler_projection_fun(RDC, NPpF, angs, fOut)
        toc

end




