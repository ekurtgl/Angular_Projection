clc; clear; close all;

fname = '/mnt/HDD04/Projection_data/Cascade_AWR2243/pAhmed_class1_ang-45_iter1';
MTI = 0;
M_pulse = 60;
K = 1;
OF = 0;

SweepTime = 40e-3;
fps = 1/SweepTime;
NPpF = 60;

RDC = RDC_extract_cascade_AWR2243(fname);

RDC = reshape(RDC,size(RDC,1), size(RDC,2), size(RDC,3)*size(RDC,4));

RDC_az = RDC(:,NPpF:end, antenna_azimuth_only);  % first frame is generally corrupt

  %% MTI
h = [1 -2 1]; % [1 -2 1]
if MTI
    RDC2 = bsxfun(@minus, RDC_az, mean(RDC_az,2));   % subtract stationary objects
    RDC_az  = filter(h,1,RDC2,[],2);
end
       
 %% Range profile
 
 % take FFT on rach channel separately, otherwise range profile gets corrupted
% limit range bins after plotting

for i = 1:size(RDC_az,3)
        RDC_az(:,:, i) = fft(RDC_az(:,:, i));
end

% figure;
% colormap(jet)
% imagesc(10*log10(abs(RDC_az(:,:,1)) / max(max(abs(RDC_az(:))))));


n_frames = floor(size(RDC_az,2) / NPpF);
duration = n_frames / fps;


%% Steering matrix
ang_ax = -90:90;
d = 0.5;

for k=1:length(ang_ax)
        a1(:,k)=exp(-1i*2*pi*(d*(0:size(RDC_az,3)-1)'*sin(ang_ax(k).'*pi/180)));
end

%% Range-Az map

figure('visible','off')
colormap(jet)
for j = 1:n_frames
        disp(['Frame ' int2str(j) '/' int2str(n_frames)]);
        for i = 1:70 % size(RDC_az,1)
                
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
            imagesc(ang_ax,[],10*log10(abs(range_az_music)./max(abs(range_az_music(:)))));
        end

        xlabel('Azimuth')   
        ylabel('Range (m)')
%                 set(gca, 'CLim',[-35,0]); % [-35,0]
        axis([-60 60 0 70])
%                 title('MUSIC Range-Angle Map')
%                 clim = get(gca,'clim');
        drawnow
        F(j) = getframe(gcf); % gcf returns the current figure handle          
end


% fname = [fNameOut(1:end-4) '_K' int2str(K) '_Mpulse' ...
% int2str(M_pulse) '_MTI' int2str(MTI) '_OF' int2str(OF) '_azimuth.avi'];
fname = 'test_noMTI_ang-45.avi';
writerObj = VideoWriter(fname);
writerObj.FrameRate = fps;
open(writerObj);

for i=1:length(F)
        frame = F(i) ;
        writeVideo(writerObj, frame);
        %         frame2 = F2(i);
        %         writeVideo(writerObj2, frame2);
end
close(writerObj);
close all









