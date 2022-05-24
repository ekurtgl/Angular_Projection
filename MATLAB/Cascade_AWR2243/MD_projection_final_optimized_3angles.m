clc
close all
clear all

%loading previously saved data file

raw_data_left = zeros(512,20000,16);
raw_data_zero = zeros(512,20000,16);
raw_data_right = zeros(512,20000,16);

%for z = 1:67
    z = 1;
    a = num2str(z);

    address1 = 'Data_UA\spot\';
    address = strcat(address1,a);
    load (address);

    % load Angle_data\left-right_left.mat;
    % raw_data_left = raw_data_2;
    % load Angle_data\person_walking45_right.mat;
    % raw_data_zero = raw_data_2;
    % load Angle_data\person_zero_onehand.mat;
    % raw_data_right = raw_data_2;

    % raw_data_2 = zeros(512,80000,16);
    % raw_data_2(:,1:40000,:) = raw_data_left;
    % raw_data_2(:,40001:80000,:) = raw_data_right;

    % raw_data_2 = zeros(512,40000,16);
    % raw_data_2 = raw_data_left + raw_data_right + raw_data_zero;

    % Configure script
    Disp_FrmNr = 1;
    Disp_TimSig = 0;      % display time signals
    Disp_RP = 0;      % display range profile
    Disp_RD = 1;      % display range-Doppler map
    c0 = 299792458;     %speed of light

    %configurations
    Cfg.fStrt = 76e9; %chirpup_start_freq
    Cfg.fStop = 78e9; %chirpdown_stop_freq
    Cfg.TRampUp = 260e-6; %chirpup_time
    Cfg.TRampDo = 10e-6; %chirpdown_time
    Cfg.TInt = 240e-3; %time for Cfg.NLoop to finish. 800 chirps=800 frames. 
    Cfg.N = 512; %total range bin in one chirp up
    Cfg.Chirps = length(raw_data_2(1,:,1)); %total number of chirps to be collected
    Cfg.CfgTim = 30e-6; %time between one chirp end to another chirp start
    Cfg.IniTim = 10e-3; %initial time for the device before start recording data
    Cfg.IniEve = 0;
    Cfg.NLoop = 800; %number of frames the device will store
    numTX = 1;
    numRX = 16;

    fs_adc = 40*10^6; %calculated from the board the sampling frequency set by the device
    fs_actual = Cfg.N/Cfg.TRampUp; %actual sampling frequency
    R = ((fs_adc)/(fs_actual)); %the device only accepts 
    R_floor = floor((fs_adc)/(fs_actual)); %from the documents the RBK2 only supports sampling frequency which is an integer multiple of fs_adc(line 28)
    fs = fs_adc/R_floor; %sampling frequency that will be used

    % Processing.................................................................
    Win2D = hanning(Cfg.N-1); 
    Win2D = repmat(Win2D,1,Cfg.Chirps); %for windowing the raw signal
    ScaWin = sum(Win2D(:,1)); %this is to scale down the raw signal after windowing
    NFFT = 2^10; %number of FFT points for fast time
    NFFTVel = 2^10; %number of FFT points for slow time
    kf = (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp; %the slope of the chirp
    B = (Cfg.fStop - Cfg.fStrt); %Bandwidth of the FMCW radar
    vRange = [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf); %calculating the maximum unambiguous range the radar can cover w.r.t. range bin
    fc = (Cfg.fStop + Cfg.fStrt)/2; %centre freuency
    Tp = Cfg.TRampUp + Cfg.TRampDo; %time duration of one chirp
    PRI = Cfg.TRampUp + Cfg.TRampDo + Cfg.CfgTim; %Pulse repetation interval

    RMin = 0; %minimum range
    RMax = 5.5; %maximum range that user wants to claculate

    [Val RMinIdx] = min(abs(vRange - RMin)); 
    [Val RMaxIdx] = min(abs(vRange - RMax));
    vRangeExt = vRange(RMinIdx:RMaxIdx-1); %converts range bins to meters

    rd_window_size = 100;
    exclude_bins = 20; %excluding first 10 bins
    WinVel = hanning(rd_window_size); %for windowing the range_profile
    ScaWinVel = sum(WinVel); %this is to scale down the range profile after windowing
    WinVel2D = repmat(WinVel.',numel(vRangeExt)-exclude_bins+1,1);

    vFreqVel = [-NFFTVel./2:NFFTVel./2-1].'./NFFTVel.*(1/Tp); %calculating the maximum unambiguous frequency w.r.t RBK2 configuration
    [Val vFreqVelmin] = min(abs(vFreqVel - (min(vFreqVel)*(Tp/PRI))));
    [Val vFreqVelmax] = min(abs(vFreqVel - (max(vFreqVel)*(Tp/PRI)))); %calculating the maximum unambiguous frequency w.r.t users configuration

    %for line 60-62 the program of the radarbook2 to doesnot calculate config
    %time in the pulse repetation interval sequence. because of that
    %calculating it with Tp gives a flat line both up and down portions of the
    %micro doppler figure. So user needed to correct the unambiguous frequency
    %and velocity measurements.

    freq_min = vFreqVel(vFreqVelmin); 
    freq_max = vFreqVel(vFreqVelmax); %converts FFT points into frequency for plotting
    vVel = vFreqVel*c0/(2.*fc);  %calculating unambiguous velocity
    velocity_min = freq_min*c0/(2.*fc); 
    velocity_max = freq_max*c0/(2.*fc); %converts FFT points into velocity for plotting

    %--------------------------------------------------------------------------
    % Select channel to be
    ChnSel = 16; %only 1 channel data out of 16 possible channels were collected
    RPExt = zeros(RMaxIdx,Cfg.Chirps,ChnSel);
    %%FFT in fast time..................................................
    for l = 1:16
        RP = fft(raw_data_2(2:end,:,l).*Win2D,NFFT,1)/ScaWin;
        RPExt(:,:,l) = RP(RMinIdx:RMaxIdx,:); %taking only the required portion
    end

    %%DC subtraction
    RPExt = RPExt - repmat(mean(RPExt, 2), [1, size(RPExt, 2)]); 

    %IQ balancing
    miux = mean(real(RPExt));
    miuy = mean(imag(RPExt));
    I2_bar = mean((real(RPExt)-miux).^2);
    Q2_bar = mean((imag(RPExt)-miuy).^2);
    IQ_bar = mean((real(RPExt)-miux).*(imag(RPExt)-miuy));
    D_bar = IQ_bar./I2_bar;
    C_bar = sqrt(Q2_bar./I2_bar-D_bar.^2);
    d_ampImb = sqrt(C_bar.^2+D_bar.^2)-1;
    phi = atan(D_bar./C_bar);
    I_rawdata = real(RPExt) - miux;
    Q_rawdata = ((imag(RPExt) - miuy)./(1+d_ampImb)...
        - real(RPExt).*sin(phi))./cos(phi);
    RPExt_2 = I_rawdata + 1i*Q_rawdata;
    figure(1)
    imagesc(20.*log10(abs(RPExt_2(exclude_bins:end,:,1))));
    ylabel('Range in meters');
    xlabel('no of Chirps');
    time = Cfg.Chirps*PRI;

    RPcal = RPExt_2(exclude_bins+1:end,:,:); %reshaping data for DOA

    %% angle

    ang_ax = -90:90;
    d = 0.5;

    for k=1:length(ang_ax)
            a1(:,k)=exp(-1i*2*pi*(d*(0:numTX*numRX-1)'*sin(ang_ax(k).'*pi/180)));
    end
  
    B_right = a1(:,round(end/6)+1:round(2*end/6)-1); % -60 to -30
    B_zero = a1(:,round(end/2)-15:round(end/2)+13); % -15 to 15
    B_left = a1(:,round(4*end/6)+1:round(5*end/6)); % 30 to 60

    B_herm_left = B_left*inv(B_left'*B_left)*B_left';
    B_herm_zero = B_zero*inv(B_zero'*B_zero)*B_zero';
    B_herm_right = B_right*inv(B_right'*B_right)*B_right';

    RDC_left = zeros(size(RPcal));
    RDC_zero = zeros(size(RPcal));
    RDC_right = zeros(size(RPcal));
    RDC = RPcal;

    %%weighting factor
    for c = 1:size(RPcal,2)
        for r = 1:size(RPcal,1)
            x = squeeze(RPcal(r,c,:));

            y_left = B_herm_left*x;
            y_zero = B_herm_zero*x;
            y_right = B_herm_right*x;
            sum_left = sum(sum(abs(y_left(:,:))));
            sum_zero = sum(sum(abs(y_zero(:,:))));
            sum_right = sum(sum(abs(y_right(:,:))));

            maximum = [sum_left sum_zero sum_right];
            sum_max = max(maximum);
    %         if sum_max == sum_left
    %             RDC_right(r,c,:) = y_right*(sum_right/(sum_max+sum_right));
    %             RDC_left(r,c,:) = y_left*(sum_left/sum_max);
    %         else
    %             RDC_right(r,c,:) = y_right*(sum_right/sum_max);
    %             RDC_left(r,c,:) = y_left*(sum_left/(sum_max+sum_left));
    %         end
            RDC_left(r,c,:) = y_left*((sum_left/sum_max)^2);
            RDC_zero(r,c,:) = y_zero*((sum_zero/sum_max)^2);
            RDC_right(r,c,:) = y_right*((sum_right/sum_max)^2);

        end
    end

    %% creating range_doppler..........................................................

    %there is two loops in this measurement because 800 frames of data was
    %captured and stored at once. and before the start of next sequence there
    %was an initializing time. to avoid the stripes because of windowing
    %operation over that region the user used shifted windows in 800 frames
    %then started using it in the next 800 frames without overlapping them
    shift = 20; %the window shifting size
    for loop = 1:Cfg.Chirps/Cfg.NLoop
        N_meas = floor(Cfg.NLoop-rd_window_size)/shift;
        for MeasIdx = 1:N_meas
            % Reshape measurement data for range doppler processing
            % 2D array with [N_meas, NLoop]
            i = ((MeasIdx-1)*shift+1)+(loop-1)*Cfg.NLoop;
            j = ((MeasIdx-1)*shift+rd_window_size)+(loop-1)*Cfg.NLoop;
            RDC_left_2 = RDC_left(:,i:j,1); %reshaping data for range doppler processing
            RD_left = fft(RDC_left_2.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RD_left = fftshift(RD_left, 2); %the range doppler
            MD_left(MeasIdx+(loop-1)*N_meas,:) = sum(abs(RD_left)); %concatenate range-doppler data to create micro-doppler image

            RDC_zero_2 = RDC_zero(:,i:j,1); %reshaping data for range doppler processing
            RD_zero = fft(RDC_zero_2.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RD_zero = fftshift(RD_zero, 2); %the range doppler
            MD_zero(MeasIdx+(loop-1)*N_meas,:) = sum(abs(RD_zero)); %concatenate range-doppler data to create micro-doppler image

            RDC_right_2 = RDC_right(:,i:j,1); %reshaping data for range doppler processing
            RD_right = fft(RDC_right_2.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RD_right = fftshift(RD_right, 2); %the range doppler
            MD_right(MeasIdx+(loop-1)*N_meas,:) = sum(abs(RD_right)); %concatenate range-doppler data to create micro-doppler image

            RDC_2 = RDC(:,i:j,1); %reshaping data for range doppler processing
            RD = fft(RDC_2.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RD = fftshift(RD, 2); %the range doppler
            MD(MeasIdx+(loop-1)*N_meas,:) = sum(abs(RD)); %concatenate range-doppler data to create micro-doppler image
            %the time signal plot. to show the plot change the value of
            %Disp_TimSig to 1 at the top
            if Disp_TimSig > 0
                % Display time signals
                figure(1)
                plot(MeasChn(2:end,:));
                gridy on;
                xlabel('n ( )');
                xlabel('u (LSB)');     
            end

            %the range profile plot. to show the plot change the value of
            %Disp_RP to 1 at the top
            if Disp_RP > 0
                % Display range profile
                figure(2)
                plot(vRangeExt, 20.*log10(abs(RPExt)));
                grid on;
                xlabel('R (m)');
                ylabel('X (dBV)');
                axis([vRangeExt(1) vRangeExt(end) -120 -40])
            end

            %Displaying range profile. Not display: change value of Disp_RD at the start
            if Disp_RD > 0
                % Display range doppler map
                figure(3)
                subplot(2,2,1)
                grid on
                imagesc(vVel, vRangeExt(exclude_bins:end), 20.*log10(abs(RD_left)));
                colormap('jet')
                caxis([-15 35]);
                xlabel('v (m/s) RD left Projection');
                ylabel('R (m)');
                xlim([velocity_min velocity_max]);

                subplot(2,2,2)
                grid on
                imagesc(vVel, vRangeExt(exclude_bins:end), 20.*log10(abs(RD_zero)));
                colormap('jet')
                caxis([-15 35]);
                xlabel('v (m/s) RD zero Projection');
                ylabel('R (m)');
                xlim([velocity_min velocity_max]);

                subplot(2,2,3)
                grid on
                imagesc(vVel, vRangeExt(exclude_bins:end), 20.*log10(abs(RD_right)));
                colormap('jet')
                caxis([-15 35]);
                xlabel('v (m/s) RD right Projection');
                ylabel('R (m)');
                xlim([velocity_min velocity_max]);

                subplot(2,2,4)
                imagesc(vVel, vRangeExt(exclude_bins:end), 20.*log10(abs(RD)));
                grid on
                xlabel('v (m/s) Merged RD');
                ylabel('R (m)');
                xlim([velocity_min velocity_max]);
                colormap('jet')
                caxis([-30 5]);
            end
        end
    end

    %% Creating Micro-Doppler Image
    MD_final_left = MD_left.';
    MD_final_zero = MD_zero.';
    MD_final_right = MD_right.';
    MD_final = MD.';
    time_MD = 1/size(MD_left,2):12/size(MD_final_left,2):12;

    figure(4)
    imagesc(time_MD, vFreqVel, 20*log10(abs(MD_final_right))); %final Spectrogram image
    caxis([25 60]);
    xlabel('time (sec) (Right side_80_100 projection)');
    ylabel('frequency (Hz)');
    axis off
    ylim([freq_min/1.2 freq_max/1.2]);
    colormap jet
    %exportgraphics(gcf,'Data_UA\multiple\.png');
    chr = int2str(z);
    address2 = 'figures\';
    address = strcat(address1,address2,chr,'45.png');
    exportgraphics(gcf,address);

    figure(5)
    imagesc(time_MD, vFreqVel, 20*log10(abs(MD_final_zero))); %final Spectrogram image
    caxis([25 60]);
    xlabel('time (sec) (straight zero_80_100 projection)');
    ylabel('frequency (Hz)');
    axis off
    ylim([freq_min/1.2 freq_max/1.2]);
    colormap jet
    %exportgraphics(gcf,'Data_UA\multiple\class18.png');
    chr = int2str(z);
    address2 = 'figures\';
    address = strcat(address1,address2,chr,'0.png');
    exportgraphics(gcf,address);
    
    figure(6)
    imagesc(time_MD, vFreqVel, 20*log10(abs(MD_final_left))); %final Spectrogram image
    caxis([25 60]);
    xlabel('time (sec) (left side_130_150 projection)');
    ylabel('frequency (Hz)');
    axis off
    ylim([freq_min/1.2 freq_max/1.2]);
    colormap jet
    %exportgraphics(gcf,'Data_UA\multiple\class113.png');
    chr = int2str(z);
    address2 = 'figures\';
    address = strcat(address1,address2,chr,'-45.png');
    exportgraphics(gcf,address);
    
    figure(7)
    imagesc(time_MD, vFreqVel, 20*log10(abs(MD_final))); %final Spectrogram image
    caxis([10 45]);
    xlabel('time (sec) Merged MD');
    ylabel('frequency (Hz)');
    ylim([freq_min/1.2 freq_max/1.2]);
    colormap jet
%end