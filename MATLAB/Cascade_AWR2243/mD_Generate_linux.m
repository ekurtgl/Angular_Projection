%% get the input path and testList
add_paths_fun();
addpath('/mnt/HDD04/Projection_data/Scripts/Cascade_AWR2243/MatlabExamples/4chip_cascade_MIMO_example/main/cascade');
% pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');  % for windows
pro_path = '/mnt/HDD04/Projection_data/Scripts/Cascade_AWR2243/MatlabExamples/4chip_cascade_MIMO_example';
input_path = strcat(pro_path, '/main/cascade/input/');
testList = strcat(input_path, 'testList.txt');
%path for input folder
fidList = fopen(testList,'r');
testID = 1;
PARAM_FILE_GEN_ON = 1;
dataPlatform = 'TDA2';
while ~feof(fidList)
    
    %% get each test vectors within the test list
    % test data file name
    dataFolder_test = fgetl(fidList);    
   
    %calibration file name
    dataFolder_calib = fgetl(fidList);
    
    %module_param_file defines parameters to init each signal processing
    %module
    module_param_file = fgetl(fidList);
    
     %parameter file name for the test
    pathGenParaFile = [input_path,'test',num2str(testID), '_param.m'];
    %important to clear the same.m file, since Matlab does not clear cache
    %automatically
    clear(pathGenParaFile);
    
    %generate parameter file for the test to run
    if PARAM_FILE_GEN_ON == 1     
        parameter_file_gen_json(dataFolder_test, dataFolder_calib, module_param_file, pathGenParaFile, dataPlatform);
    end
    
    %load calibration parameters
    load(dataFolder_calib)
    
    % simTopObj is used for top level parameter parsing and data loading and saving
    simTopObj           = simTopCascade('pfile', pathGenParaFile);
    calibrationObj      = calibrationCascade('pfile', pathGenParaFile, 'calibrationfilePath', dataFolder_calib);
    rangeFFTObj         = rangeProcCascade('pfile', pathGenParaFile);
    DopplerFFTObj       = DopplerProcClutterRemove('pfile', pathGenParaFile);
    detectionObj        = CFAR_CASO('pfile', pathGenParaFile);
    DOAObj              = DOACascade('pfile', pathGenParaFile);
    
    % get system level variables
    platform            = simTopObj.platform;
    numValidFrames      = simTopObj.totNumFrames;
    cnt = 1;
    frameCountGlobal = 0;
    
    
   % Get Unique File Idxs in the "dataFolder_test"   
   [fileIdx_unique] = getUniqueFileIdx(dataFolder_test);
     all_frame=cell(1,148);%should be number of valid frame. but doesn't have to be accurate
     tot=cell(1,148);
    for i_file = 1:(length(fileIdx_unique))
        
       % Get File Names for the Master, Slave1, Slave2, Slave3   
       [fileNameStruct]= getBinFileNames_withIdx(dataFolder_test, fileIdx_unique{i_file});        
       
      %pass the Data File to the calibration Object
      calibrationObj.binfilePath = fileNameStruct;
        
      detection_results = [];  
        
       % Get Valid Number of Frames 
       [numValidFrames dataFileSize] = getValidNumFrames(fullfile(dataFolder_test, fileNameStruct.masterIdxFile));
        %intentionally skip the first frame due to TDA2 
       
        for frameIdx = 2:1:numValidFrames;%numFrames_toRun
            tic
            %read and calibrate raw ADC data            
            calibrationObj.frameIdx = frameIdx;
            frameCountGlobal = frameCountGlobal+1
            adcData = datapath(calibrationObj);
            
            % RX Channel re-ordering
            adcData = adcData(:,:,calibrationObj.RxForMIMOProcess,:);            
            
            %only take TX and RXs required for MIMO data analysis
            % adcData = adcData
            
            if mod(frameIdx, 10)==1
                fprintf('Processing %3d frame...\n', frameIdx);
            end
            
            all_frame{frameIdx}=adcData;% Get frame by frame ADC data
            
            % For micro-Doppler we need one channel, Hence collapsing 3rd
            % and 4th dimension.
            reshape_adc=reshape(adcData,size(adcData,1), size(adcData,2), size(adcData,3)*size(adcData,4));
            tot{frameIdx}=reshape_adc(:,:,1);
        end
    end
end


 %Concate across cloumn direction to stitch all the frames
rdc=horzcat(tot{:});

rp=fft(rdc);
figure; 
imagesc(20*log10(abs(rp)))
 %% MTI v2
[b,a]=butter(1, 0.01, 'high'); %  4th order is 24dB/octave slope, 6dB/octave per order of n
%                                      [B,A] = butter(N,Wn, 'high') where N filter order, b (numerator), a (denominator), ...
%                                      highpass, Wn is cutoff freq (half the sample rate)
[m,n]=size(rp);
rngpro=zeros(m,n);
for k=1:size(rp,1)
    rngpro(k,:)=filter(b,a,rp(k,:));
end

rBin=20:35;

nfft = 2^12;window = 256;noverlap = 200;shift = window - noverlap;
sx = myspecgramnew(sum(rngpro(rBin,:)),window,nfft,shift); % mti filter and IQ correction

sx2 = abs(flipud(fftshift(sx,1)));


  %% Spectrogram

NPpF=60;
NoF=numValidFrames-2;
SweepTime=100e-3;
numTX=12;
timeAxis =  [1:NPpF*NoF]*SweepTime/NPpF; 
prf=128/0.1; % chirp per frame/ sweep time  

figure('visible','on');
colormap(jet(256));
imagesc(timeAxis,((3e8)*[-prf/2 prf/2])/(2*77e9),20*log10(sx2./max(sx2(:))));
     
axis xy
set(gca,'FontSize',10)
xlabel('Time (sec)');
ylabel('Velocity (M/S)');
caxis([-45 0]) % 40

temp = 0;

