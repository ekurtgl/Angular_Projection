function [] = RDC_to_microDoppler_projection(RDC, NPpF, fNameOut)
        SweepTime = 40e-3;
        prf = 1/(SweepTime / NPpF);
        %% MTI
%         h = [1 -2 1]; % [1 -2 1]
%         if MTI
%                 RDC = bsxfun(@minus, RDC, mean(RDC,2));   % subtract stationary objects
%                 RDC  = filter(h,1,RDC,[],2);
%         end
        
         %% Range profile
 
         % take FFT on each channel separately, otherwise range profile gets corrupted
         % limit range bins after plotting
         
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
        colorbar;
%         set(gcf,'units','normalized','outerposition',[0,0,1,1]);
         caxis([-40 0]) % 40
%         clim = get(gca,'CLim');
%         set(gca, 'YDir','normal','clim',[clim(1)+90 clim(2)])
        axis([0 timeAxis(end) -prf/2 prf/2])
        
        ylabel('Frequency (Hz)');
        xlabel('Times (s)');
        frame = frame2im(getframe(gcf));
        imwrite(frame,[fNameOut(1:end-4) '_orig.png']);
        %% Steering Matrix
        
        ang_ax = -90:90;
        d = 0.5;
        
        for k=1:length(ang_ax)
                a1(:,k)=exp(-1i*2*pi*(d*(0:size(RDC,3)-1)'*sin(ang_ax(k).'*pi/180)));
        end
        
        %% Projection Matrix
        
        for i = 1:1%9 % divide into 10 degrees, set limit to 9
            
%         B_zero = a1(:,round(end/2)-15:round(end/2)+13); % -15 to 15
        B_zero = a1(:,90:92); % -1 to 1
        B_herm_zero = B_zero*inv(B_zero'*B_zero)*B_zero';
%         B_left = a1(:,(i-1)*10+1:i*10); % 10 degree intervals
%         B_right = a1(:,(i+8)*10+1:(i+9)*10); % 10 degree intervals
%         B_left = a1(:,1:floor(end/2)); % -90 to 0
%         B_right = a1(:,ceil(end/2)+1:end); % 0 to 90
%         B_left = a1(:,round(end/6)+1:round(end/2)-1); % -60 to 0
%         B_right = a1(:,round(end/2)+1:round(5*end/6)); % 0 to 60
%         B_left = a1(:,round(end/6)+1:round(2*end/6)-1); % -60 to -30
%         B_right = a1(:,round(4*end/6)+1:round(5*end/6)); % 30 to 60
        B_left = a1(:,44:46); % -45
        B_right = a1(:,134:136); % 45
%         B_left = a1(:,1:round(end/3)-1); % -90 to -30
%         B_right = a1(:,round(2*end/3)+1:end); % 30 to 90 --> this worked fine for ahmed's walking iter 2
        
        B_herm_left = B_left*inv(B_left'*B_left)*B_left';
        B_herm_right = B_right*inv(B_right'*B_right)*B_right';
        
        RDC_left = zeros(size(RDC));
        RDC_zero = zeros(size(RDC));
        RDC_right = zeros(size(RDC));
        weight = 2;
        
        for c = 1:size(RDC,2)
                disp([int2str(c) '/' int2str(size(RDC,2))]);
            for r = 1:size(RDC,1)
                x = squeeze(RDC(r,c,:));
                
                y_left = B_herm_left*x;
                y_zero = B_herm_zero*x;
                y_right = B_herm_right*x;
                
                sum_left = sum(sum(abs(y_left(:,:))));
                sum_zero = sum(sum(abs(y_zero(:,:))));
                sum_right = sum(sum(abs(y_right(:,:))));
                
                maximum = [sum_left sum_zero sum_right];
                sum_max = max(maximum);
                
%                 RDC_left(r,c,:) = y_left;
%                 RDC_zero(r,c,:) = y_zero;
%                 RDC_right(r,c,:) = y_right;
                
                RDC_left(r,c,:) = y_left*((sum_left/sum_max)^weight);
                RDC_zero(r,c,:) = y_zero*((sum_zero/sum_max)^weight);
                RDC_right(r,c,:) = y_right*((sum_right/sum_max)^weight);
                
            end
        end
        
        %% right
        sx = myspecgramnew(sum(RDC_right(rBin,:,1)),window,nfft,shift); 
        sx2 = abs(flipud(fftshift(sx,1)));
        figure('visible','off');
        colormap(jet(256));
        imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2 / max(sx2(:))));
        colorbar;
%         set(gcf,'units','normalized','outerposition',[0,0,1,1]);
         caxis([-40 0]) % 40
%         clim = get(gca,'CLim');
%         set(gca, 'YDir','normal','clim',[clim(1)+90 clim(2)])
        axis([0 timeAxis(end) -prf/2 prf/2])
        ylabel('Frequency (Hz)');
        xlabel('Times (s)');
        frame = frame2im(getframe(gcf));
%         imwrite(frame,[fNameOut(1:end-4) '_right_' int2str(i) '.png']);
        imwrite(frame,[fNameOut(1:end-4) '_weight' int2str(weight) '_right_45.png']);
        
        %% left
        sx = myspecgramnew(sum(RDC_left(rBin,:,1)),window,nfft,shift); 
        sx2 = abs(flipud(fftshift(sx,1)));
        figure('visible','off');
        colormap(jet(256));
        imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2 / max(sx2(:))));
        colorbar;
%         set(gcf,'units','normalized','outerposition',[0,0,1,1]);
         caxis([-40 0]) % 40
%         clim = get(gca,'CLim');
%         set(gca, 'YDir','normal','clim',[clim(1)+90 clim(2)])
        axis([0 timeAxis(end) -prf/2 prf/2])
        ylabel('Frequency (Hz)');
        xlabel('Times (s)');
        frame = frame2im(getframe(gcf));
%         imwrite(frame,[fNameOut(1:end-4) '_left_' int2str(i) '.png']);
        imwrite(frame,[fNameOut(1:end-4) '_weight' int2str(weight) '_left_-45.png']);
        
         %% zero
        sx = myspecgramnew(sum(RDC_zero(rBin,:,1)),window,nfft,shift); 
        sx2 = abs(flipud(fftshift(sx,1)));
        figure('visible','off');
        colormap(jet(256));
        imagesc(timeAxis,[-prf/2 prf/2],20*log10(sx2 / max(sx2(:))));
        colorbar;
%         set(gcf,'units','normalized','outerposition',[0,0,1,1]);
         caxis([-40 0]) % 40
%         clim = get(gca,'CLim');
%         set(gca, 'YDir','normal','clim',[clim(1)+90 clim(2)])
        axis([0 timeAxis(end) -prf/2 prf/2])
        ylabel('Frequency (Hz)');
        xlabel('Times (s)');
        frame = frame2im(getframe(gcf));
%         imwrite(frame,[fNameOut(1:end-4) '_left_' int2str(i) '.png']);
        imwrite(frame,[fNameOut(1:end-4) '_weight' int2str(weight) '_zero_-1to+1.png']);
        end
end