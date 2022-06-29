clc; clear; close all

main = '/mnt/HDD04/Projection_data/Cascade_AWR2243/';
fnames = cell(1,2);
f1 = 'pAhmed_class1_ang+45_iter1';
f2 = 'pAhmed_class2_ang-45_iter3';
fnames{1} = [main f1];
fnames{2} = [main f2];
fNameOut = ['results/rd ra md sample/' f1 '_' f2 '.png'];
ang = [-45 45];
NPpF = 60;
antenna_azimuth_only = [1;2;3;4;17;18;19;20;33;34;35;5;6;7;8;21;22;23;24;37;38;39;40;53;54;55;56;69;70;71;72;85;86;87;88;101;
        102;103;104;117;118;119;120;133;134;135;9;10;11;12;13;14;15;16;29;30;31;32;45;46;47;48;61;62;63;64;77;78;79;80;93;94;
        95;96;109;110;111;112;125;126;127;128;141;142;143;144]; % retrieved from caliboj of TI's script, for 12 TXs

RDCs = cell(1, length(fnames));

for i = 1:length(RDCs)
        disp(['RDC ' int2str(i) '/' int2str(length(RDCs))]);
        RDC = RDC_extract_cascade_AWR2243(fnames{i});
        RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));
        RDC = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt
        RDCs{i} = RDC;
end

com_RDC = RDCs{1};

if length(RDCs) > 1
        for i = 2:length(RDCs)
                com_RDC = com_RDC + RDCs{i};
        end
end

RDC_to_RDM_RAM_microDoppler_projection_fun(com_RDC, NPpF, ang, fNameOut)


