clc; clear; close all;

datapath = '/mnt/HDD04/Projection_data/Cascade_AWR2243/';
subpath = '/mnt/HDD04/Projection_data/OUTPUTS/Cascade_AWR2243/';
subfolds = {'test_threeAct/', 'test_twoAct/', 'train_threeAct/' , 'train_twoAct/'};
NPpF = 60;
antenna_azimuth_only = [1;2;3;4;17;18;19;20;33;34;35;5;6;7;8;21;22;23;24;37;38;39;40;53;54;55;56;69;70;71;72;85;86;87;88;101;
        102;103;104;117;118;119;120;133;134;135;9;10;11;12;13;14;15;16;29;30;31;32;45;46;47;48;61;62;63;64;77;78;79;80;93;94;
        95;96;109;110;111;112;125;126;127;128;141;142;143;144]; % retrieved from caliboj of TI's script, for 12 TXs

for s = 1:length(subfolds)
        
        orig_files = dir([subpath subfolds{s} '*_orig.png']);
        outpath = [subpath 'est_' subfolds{s}];
                
        for o = 1:length(orig_files)
                tic
                disp([subfolds{s} ', ' int2str(o) '/' int2str(length(orig_files))]);
                cur_org_file = [orig_files(o).folder '/' orig_files(o).name];
                slashes = find(cur_org_file == '/');
                underscores = find(cur_org_file(slashes(end)+1:end) == '_');
                
                if length(underscores) == 12 % if three act
                        fnames = cell(1,3);
                        fnames{1} = cur_org_file(slashes(end)+1:slashes(end)+underscores(4)-1);
                        fnames{2} = cur_org_file(slashes(end)+underscores(4)+1:slashes(end)+underscores(8)-1);
                        fnames{3} = cur_org_file(slashes(end)+underscores(8)+1:slashes(end)+underscores(end)-1);
                else
                        fnames = cell(1,2);
                        fnames{1} = cur_org_file(slashes(end)+1:slashes(end)+underscores(4)-1);
                        fnames{2} = cur_org_file(slashes(end)+underscores(4)+1:slashes(end)+underscores(8)-1);
                end
                
                RDCs = cell(1, length(fnames));
                lens = [];
                for f = 1:length(fnames)
                        RDC = RDC_extract_cascade_AWR2243([datapath fnames{f}]);
                        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
                        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
                        RDCs{f} = RDC;
                        lens = [lens size(RDC,2)];
                end
                
                if max(lens) ~= min(lens)
                        for z = 1:length(lens)
                                num_zeros = max(lens) - size(RDCs{z},2);
                                if num_zeros ~= 0
                                        RDCs{z} = cat(2, RDCs{z}, zeros(size(RDCs{z},1), num_zeros, size(RDCs{z},3)));
                                end
                        end
                end
                
                com_RDC = RDCs{1};
                
                for rdc = 2:length(RDCs)
                        com_RDC = com_RDC + RDCs{rdc};
                end
                
                proj_angles = RDC_to_rangeDopp_with_angle_est_v2(RDC);
                
                fNameOut = [outpath fnames{1} '_' fnames{2}];
                if length(fnames) == 3
                        fNameOut = [fNameOut '_' fnames{3}];
                end
                fNameOut = [fNameOut  '.png'];
                RDC_to_microDoppler_projection_fun(com_RDC, NPpF, proj_angles, fNameOut);
                toc
        end
end













