clc; clear; close all
num_elem = [4 16 32 72 86];
main = '/mnt/HDD04/Projection_data/Cascade_Alperen/';
all_files = dir(main);
all_files(1:2) = [];
fnames = cell(1,2);
f1 = all_files(11).name;
fnames{1} = [main f1];
all_files(11:12) = [];
NPpF = 60;
antenna_azimuth_only = [1;2;3;4;17;18;19;20;33;34;35;5;6;7;8;21;22;23;24;37;38;39;40;53;54;55;56;69;70;71;72;85;86;87;88;101;
        102;103;104;117;118;119;120;133;134;135;9;10;11;12;13;14;15;16;29;30;31;32;45;46;47;48;61;62;63;64;77;78;79;80;93;94;
        95;96;109;110;111;112;125;126;127;128;141;142;143;144]; % retrieved from caliboj of TI's script, for 12 TXs

RDC1 = RDC_extract_cascade_AWR2243(fnames{1});
RDC1 = reshape(RDC1,size(RDC1,1), size(RDC1,2), size(RDC1,3)*size(RDC1,4));
RDC1 = RDC1(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
RDCs = cell(1, length(fnames));

ang1 = -45;

for f = 3:length(all_files)
        f2 = all_files(f).name; 
        fnames{2} = [main f2];
        RDC = RDC_extract_cascade_AWR2243(fnames{2});
        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
        
        underScores = strfind(f2, '_');
        angStart = strfind(f2, 'ang');
        ang2 = str2num(f2(angStart+3:underScores(end)-1));
        
        for n = 1:length(num_elem)
                disp(['File ' int2str(f) '/' int2str(length(all_files)) ', num_elem: ' int2str(num_elem(n))]);
                fNameOut = ['/mnt/HDD04/Projection_data/OUTPUTS/Cascade_Alperen/projected/' f1 '_' f2 '_numelem' int2str(num_elem(n)) '.png'];
                
                RDCs{1} = RDC1(:,:,1:num_elem(n));
                RDCs{2} = RDC(:,:,1:num_elem(n));
                
                if size(RDCs{1},2) ~= size(RDCs{2},2)
                        num_zeros = abs(size(RDCs{1},2) - size(RDCs{2},2));
                        if size(RDCs{1},2) > size(RDCs{2},2)
                              RDCs{2} = cat(2, RDCs{2}, zeros(size(RDCs{2},1), num_zeros, size(RDCs{2},3)));
                        else
                              RDCs{1} = cat(2, RDCs{1}, zeros(size(RDCs{1},1), num_zeros, size(RDCs{1},3)));
                        end
                end
                
                com_RDC = RDCs{1};
                
                if length(RDCs) > 1
                        for i = 2:length(RDCs)
                                com_RDC = com_RDC + RDCs{i};
                        end
                end
                
                
                
                RDC_to_microDoppler_projection_fun(com_RDC, NPpF, [ang1 ang2], fNameOut);

        end
end
