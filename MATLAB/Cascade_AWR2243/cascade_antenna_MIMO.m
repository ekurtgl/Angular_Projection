clc; clear; close all
% This programme is written based on
% C:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\MatlabExamples\signal_processing_4chip_cascade.pdf
% Figures 1 and 2
fontsize = 15;
c = ['k', 'y', 'm', 'b'];
num_dev = 4;
num_tx = 3;
num_rx = 4;

antenna_dist = 0.5; % in lambda
tot_rx_len = 26.5 / antenna_dist;
tot_tx_len = 8*2 / antenna_dist;
rx_x_offset = (tot_tx_len - tot_rx_len) / 2 + 0.5; % 0.5 is arbitrary to get rid of 0.5 indices
rx_y_offset = 30; % arbitrary

%% device 4 (slave)
dev4_tx = [1 1; 5 1; 9 1];
dev4_rx = [rx_x_offset + 1 rx_y_offset; rx_x_offset + 2 rx_y_offset; rx_x_offset + 3 rx_y_offset; rx_x_offset + 4 rx_y_offset];
figure;
hold on
grid on
% title('TX - RX Positions')
% xlabel('Azimuth ($\frac{\lambda}{2}$)', 'fontname', 'times new roman', 'fontsize', fontsize)
% ylabel('Elevation ($\frac{\lambda}{2}$)', 'fontname', 'times new roman', 'fontsize', fontsize)
% xlabel('Azimuth ({\lambda}/2)', 'fontname', 'times new roman', 'fontsize', fontsize)
% ylabel('Elevation ({\lambda}/2)', 'fontname', 'times new roman', 'fontsize', fontsize)
xlabel('Azimuth (\lambda/2)', 'fontname', 'times new roman', 'fontsize', fontsize)
ylabel('Elevation (\lambda/2)', 'fontname', 'times new roman', 'fontsize', fontsize)
for i = 1:size(dev4_tx,1)
   ptx4 = scatter(dev4_tx(i,1), dev4_tx(i,2), 100, c(4)); 
end
for i = 1:size(dev4_rx,1)
   prx = scatter(dev4_rx(i,1), dev4_rx(i,2), 100, 'r'); 
end

%% device 1 (master)
dev1_tx = [10 2; 11 5; 12 7];
dev1_rx = [rx_x_offset + 12 rx_y_offset; rx_x_offset + 13 rx_y_offset; rx_x_offset + 14 rx_y_offset; rx_x_offset + 15 rx_y_offset];
for i = 1:size(dev1_tx,1)
   ptx1 = scatter(dev1_tx(i,1), dev1_tx(i,2), 100, c(1)); 
end
for i = 1:size(dev1_rx,1)
   scatter(dev1_rx(i,1), dev1_rx(i,2), 100, 'r'); 
end

%% device 3 (slave)
dev3_tx = [13 1; 17 1; 21 1];
dev3_rx = [rx_x_offset + 47 rx_y_offset; rx_x_offset + 48 rx_y_offset; rx_x_offset + 49 rx_y_offset; rx_x_offset + 50 rx_y_offset];
for i = 1:size(dev3_tx,1)
   ptx3 = scatter(dev3_tx(i,1), dev3_tx(i,2), 100, c(3)); 
end
for i = 1:size(dev3_rx,1)
   scatter(dev3_rx(i,1), dev3_rx(i,2), 100, 'r'); 
end

%% device 2 (slave)
dev2_tx = [25 1; 29 1; 33 1];
dev2_rx = [rx_x_offset + 51 rx_y_offset; rx_x_offset + 52 rx_y_offset; rx_x_offset + 53 rx_y_offset; rx_x_offset + 54 rx_y_offset];
for i = 1:size(dev2_tx,1)
   ptx2 = scatter(dev2_tx(i,1), dev2_tx(i,2), 100, c(2)); 
end
for i = 1:size(dev2_rx,1)
   scatter(dev2_rx(i,1), dev2_rx(i,2), 100, 'r'); 
end
legend([ptx1 ptx2 ptx3 ptx4, prx], {'Dev-1 TX', 'Dev-2 TX', 'Dev-3 TX', 'Dev-4 TX', 'RX'}, 'fontname', 'times new roman')
set(gca, 'fontsize', fontsize)
set(gca, 'TickLabelInterpreter', 'none')
set(gcf, 'color', 'white')
set(gca, 'fontname', 'times new roman')
axis tight

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print(gcf,'-dpdf','C:\Users\emrek\Desktop\Publications\Projection Journal\proj\results\figs\mimo array_v2.pdf');

%% Combine all devices

com_tx = cell(num_dev, 1);
com_tx{4,1} = dev4_tx;
com_tx{1,1} = dev1_tx;
com_tx{3,1} = dev3_tx;
com_tx{2,1} = dev2_tx;

com_rx = cell(num_dev, 1);
com_rx{4,1} = dev4_rx;
com_rx{1,1} = dev1_rx;
com_rx{3,1} = dev3_rx;
com_rx{2,1} = dev2_rx;

%% MIMO Channels

virt_array = cell(num_dev, num_tx, num_rx, num_dev);
figure;
hold on
grid on
% title('Virtual Array')
xlabel('Azimuth (\lambda/2)', 'fontname', 'times new roman', 'fontsize', fontsize)
ylabel('Elevation (\lambda/2)', 'fontname', 'times new roman', 'fontsize', fontsize)

for d = 1:num_dev
    for tx = 1:num_tx
        cur_tx = com_tx{d,1}(tx,:);
        for rx = 1:num_rx
            for dev = 1:num_dev
                cur_rx = com_rx{dev,1}(rx,:);
                x = cur_tx(1) - cur_rx(1);
                y = cur_tx(2) - cur_rx(2);
                virt_array{d, tx, rx, dev} = [x y];
                s(d) = scatter(x, y, 100, c(d)); 
            end
        end
    end
end
legend([s(1) s(2) s(3) s(4)], {'Dev-1', 'Dev-2', 'Dev-3', 'Dev-4'}, 'fontname', 'times new roman')
set(gca, 'fontsize', fontsize)
set(gca, 'fontname', 'times new roman')
set(gca, 'TickLabelInterpreter', 'none')
set(gcf, 'color', 'white')
axis tight

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print(gcf,'-dpdf','C:\Users\emrek\Desktop\Publications\Projection Journal\proj\results\figs\mimo array_v22.pdf');










