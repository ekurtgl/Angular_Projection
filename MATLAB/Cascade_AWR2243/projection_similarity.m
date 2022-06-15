clc; clear; close all;

num_channels = 86;
main = 'C:\Users\emrek\Desktop\Publications\Projection Journal\proj\data\';
origpath = [main 'microDoppler\'];
twoact = [main '*_twoAct\*png'];
threeact = [main '*_threeAct\*png'];

two_files = dir(twoact);
idx = [];
for i = 1:length(two_files)
   if contains(two_files(i).folder, '\est_')
       idx = [idx i];
   end
end
two_files(idx) = []; % remove estimated files
idx_orig = [];
for i = 1:length(two_files)
   if contains(two_files(i).name, '_orig')
       idx_orig = [idx_orig i];
   end
end
two_origs = two_files(idx_orig);

three_files = dir(threeact); % remove estimated files
idx = [];
for i = 1:length(three_files)
   if contains(three_files(i).folder, '\est_')
       idx = [idx i];
   end
end
three_files(idx) = [];
idx_orig = [];
for i = 1:length(three_files)
   if contains(three_files(i).name, '_orig')
       idx_orig = [idx_orig i];
   end
end
three_origs = three_files(idx_orig);

%% Two act

twoact_sim_img = zeros(4, 3, length(two_origs));
for i = 1:length(two_origs)
    disp(['Two act: ' num2str(i) '/' num2str(length(two_origs))]);
    underScores = strfind(two_origs(i).name, '_');
    angs = strfind(two_origs(i).name, 'ang');
    ang1 = int2str(str2num(two_origs(i).name(angs(1)+3:underScores(3)-1)));
    ang2 = int2str(str2num(two_origs(i).name(angs(2)+3:underScores(end-1)-1)));
    orig1 = imread([origpath two_origs(i).name(1:underScores(4)-1) '.png']);
    orig2 = imread([origpath two_origs(i).name(underScores(4)+1:underScores(8)-1) '.png']);
    
    im1 = imread([two_origs(i).folder '\' two_origs(i).name(1:end-9) '_proj' ang1 '.png']);
    im2 = imread([two_origs(i).folder '\' two_origs(i).name(1:end-9) '_proj' ang2 '.png']);
    
    [twoact_sim_img(1, 1, i),~] = ssim(rgb2gray(im1),rgb2gray(orig1));
    [twoact_sim_img(2, 1, i),~] = ssim(rgb2gray(im2),rgb2gray(orig1));
    [twoact_sim_img(3, 1, i),~] = ssim(rgb2gray(im1),rgb2gray(orig2));
    [twoact_sim_img(4, 1, i),~] = ssim(rgb2gray(im2),rgb2gray(orig2));

    twoact_sim_img(1, 2, i) = immse(rgb2gray(im1),rgb2gray(orig1));
    twoact_sim_img(2, 2, i) = immse(rgb2gray(im2),rgb2gray(orig1));
    twoact_sim_img(3, 2, i) = immse(rgb2gray(im1),rgb2gray(orig2));
    twoact_sim_img(4, 2, i) = immse(rgb2gray(im2),rgb2gray(orig2));
    
    twoact_sim_img(1, 3, i) = psnr(rgb2gray(im1),rgb2gray(orig1));
    twoact_sim_img(2, 3, i) = psnr(rgb2gray(im2),rgb2gray(orig1));
    twoact_sim_img(3, 3, i) = psnr(rgb2gray(im1),rgb2gray(orig2));
    twoact_sim_img(4, 3, i) = psnr(rgb2gray(im2),rgb2gray(orig2));
end

%% Three act

threeact_sim_img = zeros(9, 3, length(three_origs));
for i = 1:length(three_origs)
    disp(['Three act: ' num2str(i) '/' num2str(length(three_origs))]);
    underScores = strfind(three_origs(i).name, '_');
    angs = strfind(three_origs(i).name, 'ang');
    ang1 = int2str(str2num(three_origs(i).name(angs(1)+3:underScores(3)-1)));
    ang2 = int2str(str2num(three_origs(i).name(angs(2)+3:underScores(7)-1)));
    ang3 = int2str(str2num(three_origs(i).name(angs(3)+3:underScores(11)-1)));
    orig1 = imread([origpath three_origs(i).name(1:underScores(4)-1) '.png']);
    orig2 = imread([origpath three_origs(i).name(underScores(4)+1:underScores(8)-1) '.png']);
    orig3 = imread([origpath three_origs(i).name(underScores(8)+1:underScores(12)-1) '.png']);
    
    im1 = imread([three_origs(i).folder '\' three_origs(i).name(1:end-9) '_proj' ang1 '.png']);
    im2 = imread([three_origs(i).folder '\' three_origs(i).name(1:end-9) '_proj' ang2 '.png']);
    im3 = imread([three_origs(i).folder '\' three_origs(i).name(1:end-9) '_proj' ang3 '.png']);
    
    [threeact_sim_img(1, 1, i),~] = ssim(rgb2gray(im1),rgb2gray(orig1));
    [threeact_sim_img(2, 1, i),~] = ssim(rgb2gray(im2),rgb2gray(orig1));
    [threeact_sim_img(3, 1, i),~] = ssim(rgb2gray(im3),rgb2gray(orig1));
    [threeact_sim_img(4, 1, i),~] = ssim(rgb2gray(im1),rgb2gray(orig2));
    [threeact_sim_img(5, 1, i),~] = ssim(rgb2gray(im2),rgb2gray(orig2));
    [threeact_sim_img(6, 1, i),~] = ssim(rgb2gray(im3),rgb2gray(orig2));
    [threeact_sim_img(7, 1, i),~] = ssim(rgb2gray(im1),rgb2gray(orig3));
    [threeact_sim_img(8, 1, i),~] = ssim(rgb2gray(im2),rgb2gray(orig3));
    [threeact_sim_img(9, 1, i),~] = ssim(rgb2gray(im3),rgb2gray(orig3));

    threeact_sim_img(1, 2, i) = immse(rgb2gray(im1),rgb2gray(orig1));
    threeact_sim_img(2, 2, i) = immse(rgb2gray(im2),rgb2gray(orig1));
    threeact_sim_img(3, 2, i) = immse(rgb2gray(im3),rgb2gray(orig1));
    threeact_sim_img(4, 2, i) = immse(rgb2gray(im1),rgb2gray(orig2));
    threeact_sim_img(5, 2, i) = immse(rgb2gray(im2),rgb2gray(orig2));
    threeact_sim_img(6, 2, i) = immse(rgb2gray(im3),rgb2gray(orig2));
    threeact_sim_img(7, 2, i) = immse(rgb2gray(im1),rgb2gray(orig3));
    threeact_sim_img(8, 2, i) = immse(rgb2gray(im2),rgb2gray(orig3));
    threeact_sim_img(9, 2, i) = immse(rgb2gray(im3),rgb2gray(orig3));
    
    threeact_sim_img(1, 3, i) = psnr(rgb2gray(im1),rgb2gray(orig1));
    threeact_sim_img(2, 3, i) = psnr(rgb2gray(im2),rgb2gray(orig1));
    threeact_sim_img(3, 3, i) = psnr(rgb2gray(im3),rgb2gray(orig1));
    threeact_sim_img(4, 3, i) = psnr(rgb2gray(im1),rgb2gray(orig2));
    threeact_sim_img(5, 3, i) = psnr(rgb2gray(im2),rgb2gray(orig2));
    threeact_sim_img(6, 3, i) = psnr(rgb2gray(im3),rgb2gray(orig2));
    threeact_sim_img(7, 3, i) = psnr(rgb2gray(im1),rgb2gray(orig3));
    threeact_sim_img(8, 3, i) = psnr(rgb2gray(im2),rgb2gray(orig3));
    threeact_sim_img(9, 3, i) = psnr(rgb2gray(im3),rgb2gray(orig3));
end

save([num2str(num_channels) '_elem_similarity_results.mat'], ...
    'twoact_sim_img', 'threeact_sim_img');
disp(['File ' num2str(num_channels) '_elem_similarity_results.mat is created!'])


