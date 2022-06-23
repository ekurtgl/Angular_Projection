clc; clear; close all;

main = '/mnt/HDD04/Projection_data/Cascade_AWR2243/';
fnames = cell(1,3);
fnames{1} = 'pRony_class1_ang0_iter3';
fnames{2} = 'pRony_class6_ang-45_iter1';
fnames{3} = 'pAhmed_class7_ang+45_iter4';

num_elem = [4 8 16 32 64 72 86];
ang = [-45 0 45];
NPpF = 60;
antenna_azimuth_only = [1;2;3;4;17;18;19;20;33;34;35;5;6;7;8;21;22;23;24;37;38;39;40;53;54;55;56;69;70;71;72;85;86;87;88;101;
        102;103;104;117;118;119;120;133;134;135;9;10;11;12;13;14;15;16;29;30;31;32;45;46;47;48;61;62;63;64;77;78;79;80;93;94;
        95;96;109;110;111;112;125;126;127;128;141;142;143;144]; % retrieved from caliboj of TI's script, for 12 TXs

RDCs = cell(1, 3);
lens =[];
for i = 1:length(fnames)
        disp(['RDC ' int2str(i) '/' int2str(length(fnames))]);
        RDC = RDC_extract_cascade_AWR2243([main fnames{i}]);
        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
        RDCs{i} = RDC;
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

if length(RDCs) > 1
        for i = 2:length(RDCs)
                com_RDC = com_RDC + RDCs{i};
        end
end

for i = 1:length(num_elem)
        disp([int2str(i) '/' int2str(length(num_elem)) ', num elem: ' int2str(num_elem(i))]);
        fNameOut = ['num_elem_sim/' fnames{1} '_' fnames{2} '_' fnames{3} 'numelem_' int2str(num_elem(i)) '.png'];
        RDC_to_microDoppler_projection_fun(com_RDC(:,:,1:num_elem(i)), NPpF, ang, fNameOut);
end








