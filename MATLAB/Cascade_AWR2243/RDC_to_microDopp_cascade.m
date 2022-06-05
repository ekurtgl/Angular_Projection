function [] = RDC_to_microDopp_cascade( RDC, fOut, is_spike_filter)
        
%         numADCBits = 16; % number of ADC bits per sample
        SweepTime = 40e-3; % Time for 1 frame
        NTS = size(RDC,1); %256 Number of time samples per sweep
%         numADCSamples = NTS;
        numTX = 2; % '1' for 1 TX, '2' for BPM
%         NoC = 60;%128; % Number of chirp loops
        NPpF = 60; % Number of pulses per frame
%         numRX = 4;
        
%         numLanes = 2; % do not change. number of lanes is always 4 even if only 1 lane is used. unused lanes
        % NoF = fileSize/2/NPpF/numRX/NTS; % Number of frames
        numChirps = size(RDC,2);
        NoF = floor(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
        dT = SweepTime/NPpF; %
        prf = 1/dT; %
        
        %% Range Profile
        rp = fft(RDC(:,:,1,1));
%         rp = rp - repmat(mean(rp, 2), [1, size(rp, 2)]);
        rngpro = rp;
%         figure;
%         imagesc(10*log10(abs(rngpro) / max(abs(rngpro(:)))));
      %% MTI Filter (not working)
      
%         [m,n]=size(rp);
%         %     ns = size(rp,2)+4;
%         h=[1 -2 3 -2 1]';
%         ns = size(rp,2)+length(h)-1;
%         rngpro=zeros(m,ns);
%         for k=1:m
%                 rngpro(k,:)=conv(h,rp(k,:,1));
%         end
        
        %% if MTI or not
%        mti = 0;
%        underScores = strfind(fOut, '_');
%        lastUnderScore = underScores(end-3);
%        classStart = strfind(fOut, 'class');
%        mainclass = str2num(fOut(classStart+5:lastUnderScore-1));
       %% MTI v2
%        if any([1 2 3 4 13 14]) == mainclass
%                mti = 1;
%                [b,a]=butter(1, 0.05, 'high'); %  4th order is 24dB/octave slope, 6dB/octave per order of n
%                %                                      [B,A] = butter(N,Wn, 'high') where N filter order, b (numerator), a (denominator), ...
%                %                                      highpass, Wn is cutoff freq (half the sample rate)
%                [m,n]=size(rp(:,:,1));
%                rngpro=zeros(m,n);
%                for k=1:size(rp,1)
%                        rngpro(k,:)=filter(b,a,rp(k,:,1));
%                end
%                rBin = 20:109;
%        else
%                rBin = 20:21;
%        end
        rngpro = rngpro - repmat(mean(rngpro, 2), [1, size(rngpro, 2)]); % DC removal
      %% STFT
%             rBin = min(cfar_bins(1)):1+max(cfar_bins(2)); %covid 18:30, front ignore= 7:nts/2, %lab 15:31 for front
        %     for i = 1:size(cfar_bins,2) % fill zeros
        %             if cfar_bins(1,i) == 0 || cfar_bins(2,i) == 0
        %                     cfar_bins(:,i) = cfar_bins(:,i-1);
        %             end
        %     end
        
       % rBin = min(cfar_bins(cfar_bins>0)):median(cfar_bins(2,:));
        rBin=1:NTS;  % observed from the range profile
%         weight_factor = 2;
%         range_weights = rBin.^weight_factor;
%         range_weights = repmat(range_weights.', [1 size(rngpro,2)]);
%         rngpro = rngpro(rBin,:) .* range_weights;
        nfft = 2^12;window = 256;noverlap = 200;shift = window - noverlap;
%         nfft = 2^10; window = 56; noverlap = 50; shift = window - noverlap;
        sx = myspecgramnew(sum(rngpro(rBin,:)),window,nfft,shift); % mti filter and IQ correction
%         sx = myspecgramnew(sum(rngpro),window,nfft,shift); % mti filter and IQ correction
%         sx = [];
%         for i = 1:NoF
%                 sx_temp = myspecgramnew(sum(rngpro(:, (i-1) * NPpF + 1: i*NPpF)),window,nfft,shift);
%                 sx = [sx sx_temp];
%         end
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
                savename = [fOut(1:end-4) '_spikeFiltered.png'];
                disp('spike filtered');
        else
                savename = [fOut(1:end-4) '.png'];
        end
        %% Spectrogram
        timeAxis = [1:NPpF*NoF]*SweepTime/NPpF*numTX ; % Time
        figure('visible','off');
        colormap(jet(256));
        imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2./max(sx2(:))));
        caxis([-40 0]) % 40
        
%         if mti == 1
%                 caxis([-35 0]) % 40
%         else
%                 caxis([-45 0]) % 50
%         end

        set(gca, 'YDir','normal')
      
%         axis([0 timeAxis(end) -prf/3 prf/3])
   
        set(gca,'xtick',[],'ytick',[])
        frame = frame2im(getframe(gca));
        imwrite(frame,savename);
        close all
        
end