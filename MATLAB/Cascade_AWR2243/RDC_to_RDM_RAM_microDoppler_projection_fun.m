function [] = RDC_to_RDM_RAM_microDoppler_projection_fun(RDC, NPpF, ang, fNameOut)
        SweepTime = 40e-3;
        prf = 1/(SweepTime / NPpF);
        %% MTI
%         h = [1 -2 1]; % [1 -2 1]
%         if MTI
%                 RDC = bsxfun(@minus, RDC, mean(RDC,2));   % subtract stationary objects
%                 RDC  = filter(h,1,RDC,[],2);
%         end
        RDC_raw = RDC;
         %% Range profile
 
         % take FFT on each channel separately, otherwise range profile gets corrupted
         % limit range bins after plotting
%          num_el = 86;
%          RDC = RDC(:,:,1:num_el);
         for i = 1:size(RDC,3)
                 RDC(:,:, i) = fft(RDC(:,:, i));
         end
         
        RDC = RDC - repmat(mean(RDC, 2), [1, size(RDC, 2)]); % DC removal
         
        timeAxis = 0:SweepTime*size(RDC,2)/NPpF;
        %% original spect
        rBin = 1:256;
        nfft = 2^12;window = 256;noverlap = 200;shift = window - noverlap;
        sx = myspecgramnew(sum(RDC(rBin,:,1)),window,nfft,shift); 
        sx2 = abs(flipud(fftshift(sx,1)));
        figure('visible','off');
        colormap(jet(256));
        imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2 / max(sx2(:))));
%         colorbar;
%         set(gcf,'units','normalized','outerposition',[0,0,1,1]);
         caxis([-40 0]) % 40
%         clim = get(gca,'CLim');
%         set(gca, 'YDir','normal','clim',[clim(1)+90 clim(2)])
%         axis([0 timeAxis(end) -prf/2 prf/2])
        
        set(gca,'xtick',[],'ytick',[])
        set(gca, 'YDir','normal')
        frame = frame2im(getframe(gca));
        %         imwrite(frame,[fNameOut(1:end-4) '_right_' int2str(i) '.png']);
        imwrite(frame,[fNameOut(1:end-4) '_orig.png']);
        
        %% orig RDM
        
        numTX = 2;
        numRX = 4;
        NTS = size(RDC_raw,1); %64 Number of time samples per sweep
        NoC = 60; % Number of chirp loops
        NPpF = NoC; % Number of pulses per frame
        fstart = 77e9; % Start Frequency
        fstop = fstart+4e9;%1.79892e9;%   Stop Frequency
        sampleFreq = 6.25e6; % 2e6 ADC Sampling frequency
        slope = 66.578e12; %29.982e12; % Mhz / us = e6/e-6 = e12
        %     numADCBits = 16; % number of ADC bits per sample
        
        fc = (fstart+fstop)/2; % Center Frequency
        c = physconst('LightSpeed'); % Speed of light
        lambda = c/fc; % Lambda
        d = lambda/2; % element spacing (in wavelengths)
        SweepTime = 40e-3; % Time for 1 frame=sweep
        numChirps = size(RDC_raw,2);
        NoF = round(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
        Bw = fstop - fstart; % Bandwidth
        
        dT = SweepTime/NPpF; %
        prf = 1/dT;
        timeAxis = linspace(0,SweepTime*NoF,numChirps);%[1:NPpF*NoF]*SweepTime/NPpF ; % Time
        duration = max(timeAxis);
        
        idletime = 100e-6;
        adcStartTime = 6e-6;
        rampEndTime = 60e-6;
        
        %% Range-Velocity Map
        
        Rmax = sampleFreq*c/(2*slope);
        Tc = idletime+adcStartTime+rampEndTime;
        Tf = SweepTime;
        velmax = lambda/(Tc*4); % Unambiguous max velocity
        DFmax = velmax/(c/fc/2);
        rResol = c/(2*Bw);
        vResol = lambda/(2*Tf);
        % define frame size
        PN = NTS; %10 equally time spaced matricex: 10X500=5000
        RANGE_FFT_SIZE = NTS;
        DOPPLER_FFT_SIZE = NoC*2; %*2
        
        
        RNGD2_GRID = linspace(0, Rmax, RANGE_FFT_SIZE);
        DOPP_GRID = linspace(DFmax, -DFmax, DOPPLER_FFT_SIZE);
        
        V_GRID = (c/fc/2)*DOPP_GRID;
        
        RCData = RDC_raw(:,:,1);
        fps = 25;%1/SweepTime;
        n_frames = floor(size(RDC_raw,2)/NPpF);
        shft = NPpF;
        
      %%  CA-CFAR params
        numGuard = 4;
        numTrain = numGuard*2;
        P_fa = 1e-5; % Prob of false alarm
        SNR_OFFSET = -5; % -10
        %     cfar_bins = ones(2,n_frames);
        figure('Visible','off')%,
        % set(gcf,  'units', 'normalized','position', [0.2 0.2 0.4 0.6])
%         F3 = zeros(RANGE_FFT_SIZE,DOPPLER_FFT_SIZE,n_frames);
        for k = 1:n_frames
                
                RData_frame = RCData(:, 1+(k-1)*shft:k*shft);
                RData_frame = bsxfun(@minus, RData_frame, mean(RData_frame,2));   % subtract stationary objects
                G_frame = fftshift(fft2(RData_frame, RANGE_FFT_SIZE,DOPPLER_FFT_SIZE),2); % 2 adjust  shift
                RDM_dB = 10*log10(abs(G_frame)./max(abs(G_frame(:))));
                %         time_Counter = (k/n_frames)*duration;
                [RDM_mask, cfar_ranges, ~] = ca_cfar(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);
                if ~isempty(cfar_ranges)
                        cfar_bins(1,k) = min(cfar_ranges);
                        cfar_bins(2,k) = max(cfar_ranges);
                else
                        cfar_bins(1,k) = 17;
                        cfar_bins(2,k) = 35;
                end
                
                imagesc([],[],RDM_dB);
                %         xlabel('Radial Velocity (m/s)','FontSize',13, 'FontName','Times')
                %         ylabel('Range (meter)','FontSize',13, 'FontName','Times')
                %         title({'Range-Velocity Map';num2str(time_Counter,'%.2f')},'FontSize',13, 'FontName','Times')
                %         colorbar
                set(gca, 'CLim',[-10,0]); % [-35,0],
                set(gca, 'Visible', 'off')
                colormap(jet) % jet
                %         caxis([90 130]) % 90 130
                %         axis xy;
%                 axis([-velmax/numTX velmax/numTX 0 6])
                ylim([0 90])
                %         set(gcf, 'Position',  [100, 100, size(G_frame,1), size(G_frame,2)])
                
                drawnow
                F(k) = getframe(gca); % gcf returns the current figure handle
%                 F(k) = getimage(gca); % gcf returns the current figure handle
                
                %% figure cfar map
                imagesc([],[],RDM_mask);
%                 axis([-velmax/numTX velmax/numTX 0 6])
                ylim([0 90])
                colormap(parula)
                set(gca, 'Visible', 'off')
                drawnow
                F2(k) = getframe(gca); % gcf returns the current figure handle
%                 F3(:,:,k) = getimage(gca);
                
                %       colormap(gray)
                %       F2(k) =  getframe(gca);
                
        end
%         save([fcfar(1:end-3),'mat'],'F3');
        
        %     fGray = [fNameOut(1:end-4) '_gray.avi'];
        
        writerObj = VideoWriter([fNameOut(1:end-3) 'avi']);
        writerObj.FrameRate = fps;
        open(writerObj);
        
        writerObj2 = VideoWriter([fNameOut(1:end-4) '_cfar.avi']);
        writerObj2.FrameRate = fps;
        open(writerObj2);
        
        for i=1:length(F)
                % convert the image to a frame
                frame = F(i) ;
                writeVideo(writerObj, frame);
                frame2 = F2(i);
                writeVideo(writerObj2, frame2);
        end
        close(writerObj);
        close(writerObj2);
        close all
        %% Steering Matrix
        
        ang_ax = -90:90;
        d = 0.5;
        
        for k=1:length(ang_ax)
                a1(:,k)=exp(-1i*2*pi*(d*(0:size(RDC,3)-1)'*sin(ang_ax(k).'*pi/180)));
        end
        
        %% orig RAM
        MTI = 0;
        M_pulse = 60; % num. pulses to be used in cov matrix
        K = 1; % num. of targets
        OF = 0; % optical flow
        RDC_az = RDC;
        numGuard = 4;
        numTrain = numGuard*2;
        P_fa = 1e-5; % Prob of false alarm
        SNR_OFFSET = -1; % -5
        
        figure('visible','off')
        colormap(jet)
        for j = 1:n_frames
                disp(['Frame ' int2str(j) '/' int2str(n_frames)]);
                for i = 1:90 % size(RDC_az,1)
                        
                        Rxx = zeros(size(a1,1),size(a1,1));
                        
                        for mp = 1:M_pulse
                                p_idx = (j-1)*NPpF+mp;
                                if j == n_frames
                                        p_idx = (j-2)*NPpF+mp;
                                end
                                A = squeeze(RDC_az(i,p_idx,:));
                                Rxx = Rxx + 1/(M_pulse) * (A*A');
                        end
                        
                        [Q,D] = eig(Rxx); % Q: eigenvectors (columns), D: eigenvalues
                        [D2, I] = sort(diag(D),'descend');
                        Q = Q(:,I); % Sort the eigenvectors to put signal eigenvectors first
                        Qs = Q(:,1:K); % Get the signal eigenvectors
                        Qn = Q(:,K+1:end); % Get the noise eigenvectors
                        
                        for k=1:length(ang_ax)
                                music_spectrum2(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
                        end
                        
                        range_az_music(i,:) = music_spectrum2;
                end
                
                if OF
                        f2 = abs(range_az_music); %imresize(frame2im(F2), 0.5);
                        flow = estimateFlow(opticFlow,f2);
                        magnitudes = flow.Magnitude;
                        imagesc(ang_ax,RNGD2_GRID(1:rangelimMatrix),20*log10(abs(range_az_music(1:rangelimMatrix,:)) .* magnitudes(1:rangelimMatrix,:) ./ max(max(abs(range_az_music((1:rangelimMatrix),:))))));
                        set(gca, 'CLim',[-25,0]);
                else
                        %             imagesc(ang_ax,[],10*log10(abs(range_az_music(1:rangelimMatrix,:))./max(abs(range_az_music(:)))));
                        RAM_dB = 10*log10(abs(range_az_music)./max(abs(range_az_music(:))));
                        [RAM_mask, cfar_ranges, cfar_dopps] = ca_cfar(RAM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);
                        imagesc(ang_ax,[],RAM_dB);
                end
                
                xlabel('Azimuth')
                colormap(jet)
                ylabel('Range (m)')
                %                 set(gca, 'CLim',[-35,0]); % [-35,0]
                axis([-60 60 0 90])
                %                 title('MUSIC Range-Angle Map')
                %                 clim = get(gca,'clim');
                drawnow
                F(j) = getframe(gca); % gcf returns the current figure handle
                
                imagesc(ang_ax,[],RAM_mask);
                axis([-60 60 0 90])
                colormap(parula)
                drawnow
                F2(j) = getframe(gca); % gcf returns the current figure handle
        end
        
        
        % fname = [fNameOut(1:end-4) '_K' int2str(K) '_Mpulse' ...
        % int2str(M_pulse) '_MTI' int2str(MTI) '_OF' int2str(OF) '_azimuth.avi'];
        fname = [fNameOut(1:end-4) '_DoA.avi'];
        writerObj = VideoWriter(fname);
        writerObj.FrameRate = fps;
        open(writerObj);
        
        writerObj2 = VideoWriter([fname(1:end-4) '_cfar.avi']);
        writerObj2.FrameRate = fps;
        open(writerObj2);
        
        
        for i=1:length(F)
                frame = F(i);
                writeVideo(writerObj, frame);
                frame2 = F2(i);
                writeVideo(writerObj2, frame2);
        end
        close(writerObj);
        close(writerObj2);
        close all
        %% Projection Matrix
        
            B_proj = cell(1, length(ang));
            B_herm = cell(1, length(ang));
            RDCs = cell(1, length(ang));
            for a = 1:length(ang)
                  B_proj{a} = a1(:,90-ang(a): 92-ang(a));
                  B_herm{a} = B_proj{a}*inv(B_proj{a}'*B_proj{a})*B_proj{a}';
                  RDCs{a} = zeros(size(RDC,1), size(RDC,2), size(RDC,3)); % use only 1 channel for microDoppler
            end
            
        weight = 2;
        
        for r = 1:size(RDC,1)
%                 disp([int2str(r) '/' int2str(size(RDC,1))]);
                x = squeeze(RDC_raw(r,:,:));
                y = cell(1, length(ang));
                sums = zeros(1, length(ang));
                for proj = 1:length(ang)
                        y{proj} = B_herm{proj}*x.';
                        sums(proj) = sum(abs(y{proj}(:)));
                end
                max_sum = max(sums);
                
                for proj = 1:length(ang)
                        RDCs{proj}(r,:,:) = (y{proj}*((sums(proj)/max_sum)^weight)).';
                end
        end
        
        %% Projection Results
        for a = 1:length(ang)
                raw_RDC = RDCs{a};
                %% RAM
                figure('visible','off')
                colormap(jet)
                for j = 1:n_frames
                        disp(['Frame ' int2str(j) '/' int2str(n_frames)]);
                        for i = 1:90 % size(RDC_az,1)
                                
                                Rxx = zeros(size(a1,1),size(a1,1));
                                
                                for mp = 1:M_pulse
                                        p_idx = (j-1)*NPpF+mp;
                                        if j == n_frames
                                                p_idx = (j-2)*NPpF+mp;
                                        end
                                        A = squeeze(raw_RDC(i,p_idx,:));
                                        Rxx = Rxx + 1/(M_pulse) * (A*A');
                                end
                                
                                [Q,D] = eig(Rxx); % Q: eigenvectors (columns), D: eigenvalues
                                [D2, I] = sort(diag(D),'descend');
                                Q = Q(:,I); % Sort the eigenvectors to put signal eigenvectors first
                                Qs = Q(:,1:K); % Get the signal eigenvectors
                                Qn = Q(:,K+1:end); % Get the noise eigenvectors
                                
                                for k=1:length(ang_ax)
                                        music_spectrum2(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
                                end
                                
                                range_az_music(i,:) = music_spectrum2;
                        end
                        
                        if OF
                                f2 = abs(range_az_music); %imresize(frame2im(F2), 0.5);
                                flow = estimateFlow(opticFlow,f2);
                                magnitudes = flow.Magnitude;
                                imagesc(ang_ax,RNGD2_GRID(1:rangelimMatrix),20*log10(abs(range_az_music(1:rangelimMatrix,:)) .* magnitudes(1:rangelimMatrix,:) ./ max(max(abs(range_az_music((1:rangelimMatrix),:))))));
                                set(gca, 'CLim',[-25,0]);
                        else
                                %             imagesc(ang_ax,[],10*log10(abs(range_az_music(1:rangelimMatrix,:))./max(abs(range_az_music(:)))));
                                RAM_dB = 10*log10(abs(range_az_music)./max(abs(range_az_music(:))));
                                [RAM_mask, cfar_ranges, cfar_dopps] = ca_cfar(RAM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);
                                imagesc(ang_ax,[],RAM_dB);
                        end
                        
                        xlabel('Azimuth')
                        colormap(jet)
                        ylabel('Range (m)')
                        %                 set(gca, 'CLim',[-35,0]); % [-35,0]
                        axis([-60 60 0 90])
                        %                 title('MUSIC Range-Angle Map')
                        %                 clim = get(gca,'clim');
                        drawnow
                        F(j) = getframe(gca); % gcf returns the current figure handle
                        
                        imagesc(ang_ax,[],RAM_mask);
                        axis([-60 60 0 90])
                        colormap(parula)
                        drawnow
                        F2(j) = getframe(gca); % gcf returns the current figure handle
                end
                
                
                % fname = [fNameOut(1:end-4) '_K' int2str(K) '_Mpulse' ...
                % int2str(M_pulse) '_MTI' int2str(MTI) '_OF' int2str(OF) '_azimuth.avi'];
                fname = [fNameOut(1:end-4) '_DoA_proj' int2str(ang(a)) '.avi'];
                writerObj = VideoWriter(fname);
                writerObj.FrameRate = fps;
                open(writerObj);
                
                writerObj2 = VideoWriter([fname(1:end-4) '_cfar.avi']);
                writerObj2.FrameRate = fps;
                open(writerObj2);
                
                
                for i=1:length(F)
                        frame = F(i);
                        writeVideo(writerObj, frame);
                        frame2 = F2(i);
                        writeVideo(writerObj2, frame2);
                end
                close(writerObj);
                close(writerObj2);
                close all
                
                %% RDM
                numGuard = 4;
                numTrain = numGuard*2;
                P_fa = 1e-5; % Prob of false alarm
                SNR_OFFSET = -5; % -10
                %     cfar_bins = ones(2,n_frames);
                figure('Visible','off')%,
                RCData = raw_RDC(:,:,1);
                % set(gcf,  'units', 'normalized','position', [0.2 0.2 0.4 0.6])
                %         F3 = zeros(RANGE_FFT_SIZE,DOPPLER_FFT_SIZE,n_frames);
                for k = 1:n_frames
                        
                        RData_frame = RCData(:, 1+(k-1)*shft:k*shft);
                        RData_frame = bsxfun(@minus, RData_frame, mean(RData_frame,2));   % subtract stationary objects
                        G_frame = fftshift(fft2(RData_frame, RANGE_FFT_SIZE,DOPPLER_FFT_SIZE),2); % 2 adjust  shift
                        RDM_dB = 10*log10(abs(G_frame)./max(abs(G_frame(:))));
                        %         time_Counter = (k/n_frames)*duration;
                        [RDM_mask, cfar_ranges, ~] = ca_cfar(RDM_dB, numGuard, numTrain, P_fa, SNR_OFFSET);
                        if ~isempty(cfar_ranges)
                                cfar_bins(1,k) = min(cfar_ranges);
                                cfar_bins(2,k) = max(cfar_ranges);
                        else
                                cfar_bins(1,k) = 17;
                                cfar_bins(2,k) = 35;
                        end
                        
                        imagesc([],[],RDM_dB);
                        %         xlabel('Radial Velocity (m/s)','FontSize',13, 'FontName','Times')
                        %         ylabel('Range (meter)','FontSize',13, 'FontName','Times')
                        %         title({'Range-Velocity Map';num2str(time_Counter,'%.2f')},'FontSize',13, 'FontName','Times')
                        %         colorbar
                        set(gca, 'CLim',[-10,0]); % [-35,0],
                        set(gca, 'Visible', 'off')
                        colormap(jet) % jet
                        %         caxis([90 130]) % 90 130
                        %         axis xy;
                        %                 axis([-velmax/numTX velmax/numTX 0 6])
                        ylim([0 90])
                        %         set(gcf, 'Position',  [100, 100, size(G_frame,1), size(G_frame,2)])
                        
                        drawnow
                        F(k) = getframe(gca); % gcf returns the current figure handle
                        %                 F(k) = getimage(gca); % gcf returns the current figure handle
                        
                        %% figure cfar map
                        imagesc([],[],RDM_mask);
                        %                 axis([-velmax/numTX velmax/numTX 0 6])
                        ylim([0 90])
                        colormap(parula)
                        set(gca, 'Visible', 'off')
                        drawnow
                        F2(k) = getframe(gca); % gcf returns the current figure handle
                        %                 F3(:,:,k) = getimage(gca);
                        
                        %       colormap(gray)
                        %       F2(k) =  getframe(gca);
                        
                end
                %         save([fcfar(1:end-3),'mat'],'F3');
                
                %     fGray = [fNameOut(1:end-4) '_gray.avi'];
                fname = [fNameOut(1:end-4) '_RDM_ang' int2str(ang(a)) '.avi'];
                writerObj = VideoWriter(fname);
                writerObj.FrameRate = fps;
                open(writerObj);
                
                writerObj2 = VideoWriter([fname(1:end-4) '_cfar.avi']);
                writerObj2.FrameRate = fps;
                open(writerObj2);
                
                for i=1:length(F)
                        % convert the image to a frame
                        frame = F(i) ;
                        writeVideo(writerObj, frame);
                        frame2 = F2(i);
                        writeVideo(writerObj2, frame2);
                end
                close(writerObj);
                close(writerObj2);
                close all
                %% mD spectrogram
                for i = 1:size(RDC,3)
                        raw_RDC(:,:, i) = fft(raw_RDC(:,:, i));
                end
                
                raw_RDC = raw_RDC - repmat(mean(raw_RDC, 2), [1, size(raw_RDC, 2)]); % DC removal
                
                sx = myspecgramnew(sum(raw_RDC(rBin,:,1)),window,nfft,shift); 
                sx2 = abs(flipud(fftshift(sx,1)));
                figure('visible','off');
                colormap(jet(256));
                imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2 / max(sx2(:))));
%                 colorbar;
                %         set(gcf,'units','normalized','outerposition',[0,0,1,1]);
                caxis([-40 0]) % 40
                %         clim = get(gca,'CLim');
                %         set(gca, 'YDir','normal','clim',[clim(1)+90 clim(2)])
%                 axis([0 timeAxis(end) -prf/2 prf/2])
%                 ylabel('Frequency (Hz)');
%                 xlabel('Times (s)');
                set(gca,'xtick',[],'ytick',[])
                set(gca, 'YDir','normal')
                frame = frame2im(getframe(gca));
                %         imwrite(frame,[fNameOut(1:end-4) '_right_' int2str(i) '.png']);
                imwrite(frame,[fNameOut(1:end-4) '_proj' int2str(ang(a)) '.png']);
                
                
        end
        close all
end