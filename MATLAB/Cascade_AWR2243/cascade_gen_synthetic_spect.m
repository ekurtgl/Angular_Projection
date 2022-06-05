clc; clear; close all; warning off

mode = {'train', 'test'};
% num_twoAct_samples = [300 150];
% num_threeAct_samples = [50 50];
num_twoAct_samples = [200 100];
num_threeAct_samples = [50 50];
datapath = '/mnt/HDD04/Projection_data/Cascade_AWR2243/';
classes = {'1', '2', '3', '4', '5', '6', '7', '13', '14'};
angles = [-45, 0, 45];
NPpF = 60;
antenna_azimuth_only = [1;2;3;4;17;18;19;20;33;34;35;5;6;7;8;21;22;23;24;37;38;39;40;53;54;55;56;69;70;71;72;85;86;87;88;101;
        102;103;104;117;118;119;120;133;134;135;9;10;11;12;13;14;15;16;29;30;31;32;45;46;47;48;61;62;63;64;77;78;79;80;93;94;
        95;96;109;110;111;112;125;126;127;128;141;142;143;144]; % retrieved from caliboj of TI's script, for 12 TXs
% classes_int = cellfun(@str2num, classes);

for m = 1:length(mode) % for train or test
        
        outpath = ['/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/' mode{m} '_twoAct/'];
        fileID = fopen([mode{m} '_files.txt'],'r');
        
        file = cell(1,1);
        cnt = 1;
        while ~feof(fileID)
                file{cnt} = fgetl(fileID);
                cnt = cnt + 1;
        end
        
        for i = 1:num_twoAct_samples(m)
                tic
                disp([mode{m} ', two act, ' num2str(i) '/' num2str(num_twoAct_samples(m))]);
                firstAct = classes{randi(length(classes))};
                firstAng = angles(randi(length(angles)));
                secondAct = classes{randi(length(classes))};
                secondAng = firstAng;
                
                while secondAng == firstAng
                        secondAng = angles(randi(length(angles)));
                end
                
                if firstAng == 45
                        firstAng = '+45';
                else
                        firstAng = num2str(firstAng);
                end
                
                if secondAng == 45
                        secondAng = '+45';
                else
                        secondAng = num2str(secondAng);
                end
                
                fnames = cell(1,2);
                while true
                        pos = randi(length(file));
                        fnames{1} = file{pos};
                        if contains(fnames{1}, ['class' firstAct '_ang' firstAng])
                                break
                        end
                end
                
                while true
                        pos = randi(length(file));
                        fnames{2} = file{pos};
                        if contains(fnames{2}, ['class' secondAct '_ang' secondAng])
                                break
                        end
                end
                
                RDCs = cell(1, length(fnames));
                for f = 1:length(fnames)
                        RDC = RDC_extract_cascade_AWR2243([datapath fnames{f}]);
                        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
                        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
                        RDCs{f} = RDC; % normalize
                end
                
                if size(RDCs{1},2) ~= size(RDCs{2},2)
                        num_zeros = abs(size(RDCs{1},2) - size(RDCs{2},2));
                        if size(RDCs{1},2) > size(RDCs{2},2)
                              RDCs{2} = cat(2, RDCs{2}, zeros(size(RDCs{2},1), num_zeros, size(RDCs{2},3)));
                        else
                              RDCs{1} = cat(2, RDCs{1}, zeros(size(RDCs{1},1), num_zeros, size(RDCs{1},3)));
                        end
                end
                
                com_RDC = RDCs{1} + RDCs{2};
                fNameOut = [outpath fnames{1} '_' fnames{2} '.png'];
                RDC_to_microDoppler_projection_fun(com_RDC, NPpF, [str2num(firstAng) str2num(secondAng)], fNameOut);
                toc
        end
        
        outpath = ['/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/' mode{m} '_threeAct/'];
        %% three acts
        for i = 1:num_threeAct_samples(m)
                tic
                disp([mode{m} ', three act, ' num2str(i) '/' num2str(num_threeAct_samples(m))]);
                firstAct = classes{randi(length(classes))};
                firstAng = angles(randi(length(angles)));
                secondAct = classes{randi(length(classes))};
                secondAng = firstAng;
                thirdAct = classes{randi(length(classes))};
                thirdAng = firstAng;
                
                while secondAng == firstAng
                        secondAng = angles(randi(length(angles)));
                end
                
                while thirdAng == firstAng || thirdAng == secondAng 
                        thirdAng = angles(randi(length(angles)));
                end
                
                if firstAng == 45
                        firstAng = '+45';
                else
                        firstAng = num2str(firstAng);
                end
                
                if secondAng == 45
                        secondAng = '+45';
                else
                        secondAng = num2str(secondAng);
                end
                
                if thirdAng == 45
                        thirdAng = '+45';
                else
                        thirdAng = num2str(thirdAng);
                end
                
                fnames = cell(1,3);
                while true
                        pos = randi(length(file));
                        fnames{1} = file{pos};
                        if contains(fnames{1}, ['class' firstAct '_ang' firstAng])
                                break
                        end
                end
                
                while true
                        pos = randi(length(file));
                        fnames{2} = file{pos};
                        if contains(fnames{2}, ['class' secondAct '_ang' secondAng])
                                break
                        end
                end
                
                while true
                        pos = randi(length(file));
                        fnames{3} = file{pos};
                        if contains(fnames{3}, ['class' thirdAct '_ang' thirdAng])
                                break
                        end
                end
                
                RDCs = cell(1, 3);
                for f = 1:length(fnames)
                        RDC = RDC_extract_cascade_AWR2243([datapath fnames{f}]);
                        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
                        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
                        RDCs{f} = RDC;
%                         RDCs{f} = RDC / max(max(abs(RDC(:))));
                end
                
                max_len = max([size(RDCs{1},2), size(RDCs{2},2) , size(RDCs{3},2)]);
                min_len = min([size(RDCs{1},2), size(RDCs{2},2) , size(RDCs{3},2)]);
                
                if max_len ~= min_len
                        num_zeros = [max_len - size(RDCs{1},2) max_len - size(RDCs{2},2) max_len - size(RDCs{3},2)];
                        for z = 1:length(num_zeros)
                                if num_zeros(z) ~= 0
                                        RDCs{z} = cat(2, RDCs{z}, zeros(size(RDCs{z},1), num_zeros(z), size(RDCs{z},3)));
                                end
                        end
                end
                
                com_RDC = RDCs{1} + RDCs{2} + RDCs{3};
                fNameOut = [outpath fnames{1} '_' fnames{2} '_' fnames{3} '.png'];
                RDC_to_microDoppler_projection_fun(com_RDC, NPpF, [str2num(firstAng) str2num(secondAng) str2num(thirdAng)], fNameOut);
                toc
        end
        
end









