clc; clear; close all;

num_channels = [4 16 32 72 86];
angs = {'-30' '-15' '0' '+15' '+30' '+45'};
ang_ax = [-30 -15 0 15 30 45];
main = 'C:\Users\emrek\Desktop\Publications\Projection Journal\proj\matlab\Cascade_Alperen\';
main2 = 'C:\Users\emrek\Desktop\Publications\Projection Journal\proj\data\';
origpath = [main 'microDoppler\'];
origpath2 = [main2 'microDoppler\'];
twoact = [main 'projected3\'];

f1 = 'pMehedi_class5_ang-45_iter1.png';

sim_img = zeros(4, 3, length(num_channels), length(angs));
for i = 1:length(angs)
    disp(['File: ' num2str(i) '/' num2str(length(angs))]);
    f2 = ['pAlperen_class2_ang' angs{i} '_iter2.png'];
    orig1 = imread([origpath2 f1]);
    orig2 = imread([origpath f2]);
    
    for n = 1:length(num_channels)

        im1 = imread([twoact f1(1:end-4) '_' f2(1:end-4) '_numelem' int2str(num_channels(n)) '_proj-45.png']);
        im2 = imread([twoact f1(1:end-4) '_' f2(1:end-4) '_numelem' int2str(num_channels(n)) '_proj' int2str(str2num(angs{i})) '.png']);
        
        [sim_img(1, 1, n, i),~] = ssim(rgb2gray(im1),rgb2gray(orig1));
        [sim_img(2, 1, n, i),~] = ssim(rgb2gray(im2),rgb2gray(orig1));
        [sim_img(3, 1, n, i),~] = ssim(rgb2gray(im1),rgb2gray(orig2));
        [sim_img(4, 1, n, i),~] = ssim(rgb2gray(im2),rgb2gray(orig2));

        sim_img(1, 2, n, i) = immse(rgb2gray(im1),rgb2gray(orig1));
        sim_img(2, 2, n, i) = immse(rgb2gray(im2),rgb2gray(orig1));
        sim_img(3, 2, n, i) = immse(rgb2gray(im1),rgb2gray(orig2));
        sim_img(4, 2, n, i) = immse(rgb2gray(im2),rgb2gray(orig2));

        sim_img(1, 3, n, i) = psnr(rgb2gray(im1),rgb2gray(orig1));
        sim_img(2, 3, n, i) = psnr(rgb2gray(im2),rgb2gray(orig1));
        sim_img(3, 3, n, i) = psnr(rgb2gray(im1),rgb2gray(orig2));
        sim_img(4, 3, n, i) = psnr(rgb2gray(im2),rgb2gray(orig2));
    end
end

fname = 'alperen_elem_similarity_results.mat';
save(fname, 'sim_img');
disp(['File ' fname '_elem_similarity_results.mat is created!'])

mse_mat = squeeze(sim_img(1, 2, :, :));
sssim_mat = squeeze(sim_img(4, 1, :, :));
% sssim_mat = squeeze(sim_img(4, 1, :, :) ./ sim_img(3, 1, :, :));

figure
fontsize = 15;
hold on; grid on
p1 = plot(ang_ax, mse_mat(1,:));
p2 = plot(ang_ax, mse_mat(2,:));
p3 = plot(ang_ax, mse_mat(3,:));
p4 = plot(ang_ax, mse_mat(4,:));
p5 = plot(ang_ax, mse_mat(5,:));
% p1 = plot(ang_ax, sssim_mat(1,:));
% p2 = plot(ang_ax, sssim_mat(2,:));
% p3 = plot(ang_ax, sssim_mat(3,:));
% p4 = plot(ang_ax, sssim_mat(4,:));
% p5 = plot(ang_ax, sssim_mat(5,:));
xticks(ang_ax);
xticklabels(angs);
legend([p1 p2 p3 p4, p5], {'4 elements', '16 elements', '32 elements', '72 elements', '86 elements'}, 'fontname', 'times new roman', 'location', 'best')
set(gca, 'fontsize', fontsize)
set(gca, 'TickLabelInterpreter', 'none')
set(gcf, 'color', 'white')
set(gca, 'fontname', 'times new roman')
axis tight

xlim([-32 47])
ylim([0.4 0.9])

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print(gcf,'-dpdf','C:\Users\emrek\Desktop\Publications\Projection Journal\proj\results\figs\aspect_angle_vs_num_elem.pdf');


