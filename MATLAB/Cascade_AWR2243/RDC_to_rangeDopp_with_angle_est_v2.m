function [proj_angles] = RDC_to_rangeDopp_with_angle_est_v2(RDC)
%         fNameOut = 'test.avi';
%         numTX = 1;
%         numRX = 4;
%         NTS = size(RDC,1); %64 Number of time samples per sweep
        NoC = 60; % Number of chirp loops
%         NPpF = numTX*NoC; % Number of pulses per frame
%         fstart = 77e9; % Start Frequency
%         fstop = fstart+4e9;%1.79892e9;%   Stop Frequency
%         sampleFreq = 6.25e6; % 2e6 ADC Sampling frequency
%         slope = 66.578e12; %29.982e12; % Mhz / us = e6/e-6 = e12
        %     numADCBits = 16; % number of ADC bits per sample
        
%         fc = (fstart+fstop)/2; % Center Frequency
%         c = physconst('LightSpeed'); % Speed of light
%         lambda = c/fc; % Lambda
%         d = lambda/2; % element spacing (in wavelengths)
%         SweepTime = 40e-3; % Time for 1 frame=sweep
%         numChirps = size(RDC,2);
%         NoF = round(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
%         Bw = fstop - fstart; % Bandwidth
        
%         dT = SweepTime/NPpF; %
%         prf = 1/dT;
%         timeAxis = linspace(0,SweepTime*NoF,numChirps);%[1:NPpF*NoF]*SweepTime/NPpF ; % Time
%         duration = max(timeAxis);
%         
%         idletime = 100e-6;
%         adcStartTime = 6e-6;
%         rampEndTime = 60e-6;
        
        %% Range-Velocity Map
        
%         Rmax = sampleFreq*c/(2*slope);
%         Tc = idletime+adcStartTime+rampEndTime;
%         Tf = SweepTime;
%         velmax = lambda/(Tc*4); % Unambiguous max velocity
%         DFmax = velmax/(c/fc/2);
%         rResol = c/(2*Bw);
%         vResol = lambda/(2*Tf);
        % define frame size
%         PN = NTS; %10 equally time spaced matricex: 10X500=5000
%         RANGE_FFT_SIZE = NTS;
%         DOPPLER_FFT_SIZE = PN*2; %*2
        
        
%         RNGD2_GRID = linspace(0, Rmax, RANGE_FFT_SIZE);
%         DOPP_GRID = linspace(DFmax, -DFmax, DOPPLER_FFT_SIZE);
        
%         V_GRID = (c/fc/2)*DOPP_GRID;
        
%         RCData = RDC(:,:,1);
%         fps = 1/SweepTime;
        n_frames = floor(size(RDC,2)/NoC);
        shft = NoC;
        
      %%  CA-CFAR params
        numGuard = 4;
        numTrain = numGuard*2;
        P_fa = 1e-5; % Prob of false alarm
        SNR_OFFSET = -10; % -5
        
%         figure('Visible','off')%,
        % set(gcf,  'units', 'normalized','position', [0.2 0.2 0.4 0.6])
        ang_list = cell(n_frames,1);
%         targets_list = cell(n_frames,1);
%         num_targets_list = zeros(n_frames,1);
%         modes_list = zeros(n_frames,1);
        target_det_threshold = 87; % half of num frames
        max_ang_sep = 15; % degrees
%         ang_separation = 10; 
%         range_separation = 10;
        tot_ang_list = [];
        tot_cnt_list = [];
        for k = 1:n_frames
%                 if mod(k,50) == 1 
%                         disp(['Frame ' int2str(k) '/' int2str(n_frames)]);
%                 end
                RData_frame = RDC(:, 1+(k-1)*shft:k*shft,:);
%                 RData_frame = RCData(:, 1+(k-1)*shft:k*shft);
%                 RData_frame = bsxfun(@minus, RData_frame, mean(RData_frame,2));   % subtract stationary objects
                G_frame = zeros(size(RData_frame));
                for i = 1:size(G_frame, 3)
                        G_frame(:,:,i) = fft(RData_frame(:,:,i));
                        G_frame(:,:,i) = G_frame(:,:,i) - repmat(mean(G_frame(:,:,i), 2), [1, size(G_frame(:,:,i), 2)]); % DC removal
                        G_frame(:,:,i) = fftshift(fft(G_frame(:,:,i), [], 2), 2);
                end
%                 G_frame = fftshift(fft2(RData_frame, RANGE_FFT_SIZE,DOPPLER_FFT_SIZE),2); % 2 adjust  shift
                RDM_dB = 10*log10(abs(G_frame(:,:,1)) / max(max(abs(G_frame(:,:,1)))));
                %         time_Counter = (k/n_frames)*duration;
                [RDM_mask, cfar_ranges, cfar_dopps] = ca_cfar(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);
                
                if ~isempty(cfar_ranges)
%                         cfar_bins(1,k) = min(cfar_ranges);
%                         cfar_bins(2,k) = max(cfar_ranges);
                        %% angle estimation for the detected RD points
                        angles = zeros(1, length(cfar_ranges));
                        for i = 1:length(cfar_ranges)
                                angles(i) = music_ang_est(squeeze(G_frame(cfar_ranges(i), cfar_dopps(i), :)));
                        end
                        ang_list{k} = angles;
                        
                        for a = 1:length(angles)
                                if any(tot_ang_list == angles(a))
%                                         idx = find(tot_ang_list == angles(a));
                                        tot_cnt_list(tot_ang_list == angles(a)) = tot_cnt_list(tot_ang_list == angles(a)) + 1;
                                else
                                        tot_ang_list(end + 1) = angles(a);
                                        tot_cnt_list(end + 1) = 1;
                                end
                        end
                        
%                         modes_list(k) = mode(angles);
%                         targets = cell(1,1);
%                         mid_trgt = mode(angles); % initialize with most possible target
%                         mid_idx = find(angles==mid_trgt);
%                         if length(mid_idx) > 1
%                                 mid_idx = mid_idx(1);
%                         end
%                         targets{1} = [cfar_ranges(mid_idx) angles(mid_idx)];
%                         if length(angles) > 1
%                                 for t = 2:length(angles)
%                                         add_target = 1;
%                                         for tt = 1:length(targets)
%                                                 if abs(targets{tt}(1) - cfar_ranges(t)) <  range_separation && ...
%                                                                 abs(targets{tt}(2) - angles(t)) <  ang_separation
%                                                         add_target = 0; % don't add new target
%                                                         break
%                                                 end
%                                         end
%                                         
%                                         if add_target == 1
%                                                 targets{length(targets) + 1} = [cfar_ranges(t) angles(t)]; % add new target
%                                         end
%                                 end
%                         end
%                         targets_list{k} = targets;
%                         num_targets_list(k) = length(targets);
                end
%                 imagesc(V_GRID,RNGD2_GRID,RDM_dB);
%                 %         xlabel('Radial Velocity (m/s)','FontSize',13, 'FontName','Times')
%                 %         ylabel('Range (meter)','FontSize',13, 'FontName','Times')
%                 %         title({'Range-Velocity Map';num2str(time_Counter,'%.2f')},'FontSize',13, 'FontName','Times')
%                 %         colorbar
% %                 set(gca, 'CLim',[-10,0]); % [-35,0],
%                 set(gca, 'CLim',[-35,0]); % [-35,0],
%                 colormap(jet) % jet
%                 %         caxis([90 130]) % 90 130
%                 %         axis xy;
%                 axis([-velmax/numTX velmax/numTX 0 6])
%                 %         set(gcf, 'Position',  [100, 100, size(G_frame,1), size(G_frame,2)])
%                 
%                 drawnow
%                 F(k) = getframe(gca); % gcf returns the current figure handle
%                 
%                 %       colormap(gray)
%                 %       F2(k) =  getframe(gca);
%                 %% figure cfar map
%                 
%                 imagesc(V_GRID,RNGD2_GRID,RDM_mask);
%                 axis([-velmax/numTX velmax/numTX 0 6])
%                 set(gca, 'Visible', 'off')
%                 drawnow
%                 F2(k) = getframe(gca); % gcf returns the current figure handle
% %                 F3(:,:,k) = getimage(gca);
                
        end
%         num_targets = mode(num_targets_list);
%         proj_angles = zeros(1, length(num_targets));
%         for i = 1:num_targets
%               proj_angles(i) = mode(modes_list);
%               modes_list(modes_list == proj_angles(i)) = [];
%         end
        [sorted_cnt, idx] = sort(tot_cnt_list, 'descend');
        sorted_ang = tot_ang_list(idx);
        sorted_ang(sorted_cnt < target_det_threshold) = [];
        sorted_ang = -sorted_ang;
%         sorted_cnt(sorted_cnt < target_det_threshold) = [];
        if ~isempty(sorted_ang)
                proj_angles(1) = sorted_ang(1);
                for i = 1:length(sorted_ang)
                        if any(abs(proj_angles - sorted_ang(i)) < max_ang_sep)
                                continue
                        else
                                proj_angles(end + 1) = sorted_ang(i);
                        end
                end
        else
                proj_angles = 'no_target';
        end
        
%         writerObj = VideoWriter(fNameOut);
%         writerObj.FrameRate = fps;
%         open(writerObj);
%         
%         writerObj2 = VideoWriter([fNameOut(1:end-4) '_cfar.avi']);
%         writerObj2.FrameRate = fps;
%         open(writerObj2);
%         
%         for i=1:length(F)
%                 % convert the image to a frame
%                 frame = F(i) ;
%                 writeVideo(writerObj, frame);
%                 frame2 = F2(i);
%                 writeVideo(writerObj2, frame2);
%         end
%         close(writerObj);
%         close(writerObj2);
        close all
end