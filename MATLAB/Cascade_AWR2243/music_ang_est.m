function angle = music_ang_est(ch_vector)
        
        ang_ax = -90:90;
        d = 0.5;
        
        for k=1:length(ang_ax)
                a1(:,k)=exp(-1i*2*pi*(d*(0:length(ch_vector)-1)'*sin(ang_ax(k).'*pi/180)));
        end
        
        Rxx = ch_vector*ch_vector';
        [Q,D] = eig(Rxx); % Q: eigenvectors (columns), D: eigenvalues
        [D, I] = sort(diag(D),'descend');
        Q = Q(:,I); % Sort the eigenvectors to put signal eigenvectors first
        Qs = Q(:,1); % Get the signal eigenvectors
        Qn = Q(:,2:end); % Get the noise eigenvectors
        
        music_spectrum = zeros(1, length(ang_ax));
        for k=1:length(ang_ax)
                music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
        end
        [~, max_ind] = max(abs(music_spectrum));
        angle = ang_ax(max_ind);
        
%         figure(3)
%         plot(ang_ax, abs(music_spectrum))
end