clc; clear; close all


datapath = '/mnt/HDD04/Projection_data/Cascade_AWR2243/*';
subfolds = dir(datapath);
subfolds(1:3) = [];

formatSpec = '%s\n';
train_files = fopen('train_files.txt','w');
test_files = fopen('test_files.txt','w');
for i = 1:length(subfolds)
        disp([int2str(i) '/' int2str(length(subfolds))]);
        
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
        
        if length(angs) > 1
                fprintf(test_files, formatSpec, subfolds(i).name);
        else
                r = rand(1);
                if r > 0.3
                        fprintf(train_files, formatSpec, subfolds(i).name);
                else
                        fprintf(test_files, formatSpec, subfolds(i).name);
                end
        end
        
end
fclose(train_files);
fclose(test_files);





