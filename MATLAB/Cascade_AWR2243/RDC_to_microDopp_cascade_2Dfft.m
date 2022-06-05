function [] = RDC_to_microDopp_cascade_2Dfft( RDC, fOut, is_spike_filter)
        
        SweepTime = 40e-3; % Time for 1 frame
        NTS = size(RDC,1); %256 Number of time samples per sweep
        numTX = 2; % '1' for 1 TX, '2' for BPM
        NPpF = 60; % Number of pulses per frame
        numRX = 4;
        
        % NoF = fileSize/2/NPpF/numRX/NTS; % Number of frames
        numChirps = size(RDC,2);
        NoF = round(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
        dT = SweepTime/NPpF; %
        prf = 1/dT; %
        %% 2d fft params
        Win2D = hanning(size(RDC,1)); 
        Win2D = repmat(Win2D,1,size(RDC,2)); %for windowing the raw signal
        ScaWin = sum(Win2D(:,1)); %this is to scale down the raw signal after windowing
        NFFT = 2^10; %number of FFT points for fast time
        NFFTVel = 2^10; %number of FFT points for slow timerd_window_size = 100;
        rd_window_size = 30; %excluding first 10 bins
        WinVel = hanning(rd_window_size); %for windowing the range_profile
        
        %% Range FFT
        rp = fft(RDC(:,:,1,1));
        
        %% DC subtraction
        rp = rp - repmat(mean(rp, 2), [1, size(rp, 2)]);
        %% IQ balancing
         miux = mean(real(rp));
         miuy = mean(imag(rp));
         I2_bar = mean((real(rp)-miux).^2);
         Q2_bar = mean((imag(rp)-miuy).^2);
         IQ_bar = mean((real(rp)-miux).*(imag(rp)-miuy));
         D_bar = IQ_bar./I2_bar;
         C_bar = sqrt(Q2_bar./I2_bar-D_bar.^2);
         d_ampImb = sqrt(C_bar.^2+D_bar.^2)-1;
         phi = atan(D_bar./C_bar);
         I_rawdata = real(rp) - miux;
         Q_rawdata = ((imag(rp) - miuy)./(1+d_ampImb)...
                 - real(rp).*sin(phi))./cos(phi);
         rp = I_rawdata + 1i*Q_rawdata;
    
        figure;
        imagesc(10*log10(abs(rp) / max(abs(rp(:)))));
        
        n_frames = floor(size(rp,2) / NPpF);
        shift =10;
        Cfg.NLoop = NPpF;
        N_meas = floor((Cfg.NLoop-rd_window_size)/shift);
        
        for loop = 1:n_frames
                for MeasIdx = 1:N_meas
                        i = ((MeasIdx-1)*shift+1)+(loop-1)*Cfg.NLoop;
                        j = ((MeasIdx-1)*shift+rd_window_size)+(loop-1)*Cfg.NLoop;
                        rp2 = rp(:,i:j); %reshaping data for range doppler processing
                        rp3 = fft(rp2, NFFTVel, 2);
                        rp3 = fftshift(rp3, 2); %the range doppler
                        MD(MeasIdx+(loop-1)*N_meas,:) = sum(abs(rp3)); %concatenate range-doppler data to create micro-doppler image
                end
        end
        MD = MD.';
        
        figure;
        colormap(jet)
        imagesc(10*log10(abs(MD) / max(abs(MD(:)))));
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %% if MTI or not
       mti = 0;
       underScores = strfind(fOut, '_');
       lastUnderScore = underScores(end-3);
       classStart = strfind(fOut, 'class');
       mainclass = str2num(fOut(classStart+5:lastUnderScore-1));
       %% MTI v2
       if any([1 2 3 4 13 14]) == mainclass
               mti = 1;
               [b,a]=butter(1, 0.05, 'high'); %  4th order is 24dB/octave slope, 6dB/octave per order of n
               %                                      [B,A] = butter(N,Wn, 'high') where N filter order, b (numerator), a (denominator), ...
               %                                      highpass, Wn is cutoff freq (half the sample rate)
               [m,n]=size(rp(:,:,1));
               rngpro=zeros(m,n);
               for k=1:size(rp,1)
                       rngpro(k,:)=filter(b,a,rp(k,:,1));
               end
               rBin = 20:109;
       else
               rBin = 20:21;
       end
      %% STFT
%             rBin = min(cfar_bins(1)):1+max(cfar_bins(2)); %covid 18:30, front ignore= 7:nts/2, %lab 15:31 for front
        %     for i = 1:size(cfar_bins,2) % fill zeros
        %             if cfar_bins(1,i) == 0 || cfar_bins(2,i) == 0
        %                     cfar_bins(:,i) = cfar_bins(:,i-1);
        %             end
        %     end
        
       % rBin = min(cfar_bins(cfar_bins>0)):median(cfar_bins(2,:));
%         rBin=22:25;  % observed from the range profile
%         weight_factor = 2;
%         range_weights = rBin.^weight_factor;
%         range_weights = repmat(range_weights.', [1 size(rngpro,2)]);
%         rngpro = rngpro(rBin,:) .* range_weights;
        nfft = 2^12;window = 256;noverlap = 200;shift = window - noverlap;
%         sx = myspecgramnew(sum(rngpro(rBin,:)),window,nfft,shift); % mti filter and IQ correction
        sx = myspecgramnew(sum(rngpro),window,nfft,shift); % mti filter and IQ correction
%         sx = myspecgramnew(sum(rp),window,nfft,shift); % mti filter and IQ correction
%         [sx, diffs] = spike_filter(sx);
      %% cfar bins
        %
        %     numrep = floor(ns/size(cfar_bins,2));
        %     b = ones(1,numrep);
        %     extended_bins = kron(cfar_bins,b);
        %     extended_bins(:,end+1:ns) = repmat(extended_bins(:,end),1,ns-size(extended_bins,2))+1;
        %     mask = zeros(size(rngpro));
        %     for i = 1:ns
        %             mask(extended_bins(1,i):extended_bins(2,i),i) = 1;
        %     end
        %     rngpro2 = rngpro.*mask;
        %     num_used = sum(mask);
        %     sx = myspecgramnew(sum(rngpro2)./num_used,window,nfft,shift); % mti filter and IQ correction
        
        sx2 = abs(flipud(fftshift(sx,1)));
        if is_spike_filter == 1
                [sx2, diffs] = spike_filter(sx2);
                savename = [fOut(1:end-4) '.png'];
                disp('spike filtered');
        else
                savename = [fOut(1:end-4) '_spikeFiltered.png'];
        end
        %% Spectrogram
        timeAxis = [1:NPpF*NoF]*SweepTime/NPpF*numTX ; % Time
        figure('visible','off');
        colormap(jet(256));
        imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2./max(sx2(:))));
        
        if mti == 1
                caxis([-35 0]) % 40
        else
                caxis([-45 0]) % 50
        end

        set(gca, 'YDir','normal')
      
%         axis([0 timeAxis(end) -prf/3 prf/3])
   
        set(gca,'xtick',[],'ytick',[])
        frame = frame2im(getframe(gca));
        imwrite(frame,savename);
        close all
        
end