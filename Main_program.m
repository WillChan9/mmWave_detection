 clear
% close all
% clc
%% parameter setting
tic
c=3.0e8;  
B=900e6;       %bandwidth
K=29.982e12;       %slope(Hz/s)
T=B/K;         %sampling time
Tc=2*57e-6;     %chirp period(s),multiply 2 in MIMO
fs=6000e3;       %sampling rate
f0=77e9;       %initial freq
lambda=c/f0;   %wave length
d=lambda/2;    %array spacing
n_samples=256; %sampling number
N=256;         %range FFT
n_chirps=128;  %chirps in one frame
M=128;         %doppler FFT
n_RX=8;        %Number of receive antenna channels
Q = 180;       %angle FFT
xx = 1;        % xx-th frame
frame_num=40;   %total frames 
frame_period=250e-3;%frame interval time
delt_r=c*fs/(2*n_samples*K); %range resolution(use effective bandwidth)
r_max=fs*c/2/K;%maximum detected range
delt_v=lambda/(2*n_chirps*Tc);%velocity resolution（Use effective time width）
v_max=lambda/4/Tc;%maximum detected velocity
delt_theta0=2/n_RX*180/pi;%angle resolution when angle=0
delt_theta30=2/n_RX/cosd(30)*180/pi;
delt_theta60=2/n_RX/cosd(60)*180/pi;

fprintf('range resolution： %2f m\n',delt_r);
fprintf('max detect range： %2f m\n',r_max);
fprintf('velocity resolution： %2f m/s\n',delt_v);
fprintf('max detect velocity： %2f m/s (%2f km/h)\n',v_max,v_max*3.6);
fprintf('azimuth resolution:\n  0°： %2f°\n  30°：%2f°\n  60°：%2f°\n',delt_theta0,delt_theta30,delt_theta60);

%% data reader
data_radar=zeros(n_samples,n_chirps,n_RX,frame_num);%4D range velocity angle time
for xx=1:frame_num %read frames
    fname='test.bin';
    fid = fopen(fname,'rb');
    sdata = fread(fid,xx*n_samples*n_chirps*4*2*2,'int16');    
    sdata = sdata((xx-1)*n_samples*n_chirps*4*2*2+1:xx*n_samples*n_chirps*4*2*2);
    fileSize = size(sdata,1);
    lvds_data = zeros(1,fileSize/2);
    count = 1;
    for i=1:4:fileSize-5
       lvds_data(1,count) = sdata(i) + 1i*sdata(i+2); %1i为复数
       lvds_data(1,count+1) = sdata(i+1)+1i*sdata(i+3); 
       count = count + 2;
    end
    lvds_data = reshape(lvds_data, n_samples*8, n_chirps);
    temp_mimo=zeros(size(lvds_data));
    temp_mimo(1:n_samples,:)=lvds_data(1:n_samples,:);%1w
    temp_mimo(1+n_samples:2*n_samples,:)=lvds_data(2*n_samples+1:3*n_samples,:); %2w
    temp_mimo(1+2*n_samples:3*n_samples,:)=lvds_data(4*n_samples+1:5*n_samples,:); %3w
    temp_mimo(1+3*n_samples:4*n_samples,:)=lvds_data(6*n_samples+1:7*n_samples,:); %4w
    temp_mimo(1+4*n_samples:5*n_samples,:)=lvds_data(1*n_samples+1:2*n_samples,:); %5w
    temp_mimo(1+5*n_samples:6*n_samples,:)=lvds_data(3*n_samples+1:4*n_samples,:); %6w
    temp_mimo(1+6*n_samples:7*n_samples,:)=lvds_data(5*n_samples+1:6*n_samples,:); %7w
    temp_mimo(1+7*n_samples:8*n_samples,:)=lvds_data(7*n_samples+1:8*n_samples,:); %8w
    lvds_data = lvds_data.'; 
    cdata = zeros(n_RX,n_chirps*n_samples);
    for row = 1:n_RX
      for i = 1: n_chirps
          cdata(row,(i-1)*n_samples+1:i*n_samples) = lvds_data(i,(row-1)*n_samples+1:row*n_samples);
      end
    end
    for i=1:n_RX
        
    data_radar(:,:,i,xx)=reshape(cdata(i,:),n_samples,n_chirps);     %data cube
    end
    fclose(fid);
end

%% 3DFFT
% stop_flag=1;
% while(stop_flag)
%     str='input frame number:';%choose the frame to observe
%     xx=input(str)
frame_count=1;
for xx=1:frame_num
    %range FFT
    range_profile=zeros(N,n_chirps,n_RX);
    range_win = hamming(n_samples);   %Hamming window
    for k=1:n_RX
       for m=1:n_chirps
          temp=data_radar(:,m,k,xx).*range_win;    %windows add on frist dimention of n_samples
          temp_fft=fft(temp,N);    % N point FFT to every chirp 
          range_profile(:,m,k)=temp_fft; %temp_fft is 256*1 column vector
       end
    end
    %doppler FFT
    speed_profile=zeros(N,M,n_RX);
    doppler_win = hamming(n_chirps);
    for k=1:n_RX
        for n=1:N
          temp=range_profile(n,:,k).*(doppler_win)';    
          temp_fft=fftshift(fft(temp,M));    % M point FFT to rangeFFT result
          speed_profile(n,:,k)=temp_fft;  %temp_fft is a 1*128 row vector
        end
    end
    %Doppler compensation, MIMO mode
    comp_profile=zeros(size(speed_profile));
    comp_profile(:,:,1:4)=speed_profile(:,:,1:4);
    for k=5:n_RX
        for t=1:M
          dop_comp=exp(-1*sqrt(-1)*pi*(t-M/2)/(M-1));
          comp_profile(:,t,k)=speed_profile(:,t,k).*dop_comp;
        end
    end

    %range-doppler map
    RVmap=zeros(N,M);
    for i=1:n_RX
        RVmap=RVmap+abs(comp_profile(:,:,i));
    end

                
    %% MVDR DoA
    RAmap_MVDR=zeros(N,181);
    Rxx=zeros(n_RX,n_RX);
    for i=1:N
        for j=1:M
            x=reshape(comp_profile(i,j,:),[n_RX,1]);
            R=x*x'/N; %data covariance matrix
            Rxx=Rxx+R;
        end
        Rxx=Rxx/M;
        iRxx = inv(Rxx); % notice the Singular value matrix bug
        ThetaSet = -90:1:90;
        theta = length(ThetaSet);   %scan angles
        A = zeros(n_RX,theta);
        phi = 2*pi*d*sind(ThetaSet)/lambda;
        for k = 0:n_RX-1 
           A(k+1,:) = exp(1j*k*phi); % positive, indicating positive angle at rightside
        end
        Amp = zeros(1,theta);
        for mm = 1: theta
            TempValue = A(:,mm)' *iRxx * A(:,mm);
            Amp(1,mm) = 1 / abs(TempValue);
        end
        RAmap_MVDR(i,:)=Amp(1,:);
    end
    
    %% OS_CFAR
    %CFAR to angle-range map's range 
    N_det=12;
    N_pro=2;
    k=N_det/3;
    Pfa=0.2;
    cfar_RAmap=zeros(size(RAmap_MVDR));
    for i=1:N-3  
        xc=(RAmap_MVDR(i,:));
        if i<30
            [~,XT1]=cfar_os(xc,N_det,k,N_pro,0.1);
        else
            [~,XT1]=cfar_os(xc,N_det,k,N_pro,Pfa);
        end
        for j=1:181
            if(xc(1,j)<XT1(1,j))
                xc(1,j)=0;
            end
        end
        cfar_RAmap(i,:)=xc;
    end
    for i=N-3:N
        cfar_RAmap(i,:)=0;
    end
%     
%     for j=1:181
%         xc=reshape(RAmap_MVDR(:,j),[1,N]);
%         [~,XT2]=cfar_os(xc,N_det,k,N_pro,0.1);
%         for i=1:N
%             if(xc(1,i)<XT2(1,i))
%                 cfar_RAmap(i,j)=0;
%             end
%         end
%     end
    %CFAR to range-doppler map,first range then doppler
    N_det=12;
    N_pro=2;
    k=N_det/2;
    Pfa=0.1;
    cfar_RVmap=zeros(size(RVmap));
    for i=1:N  
        xc=(RVmap(i,:));
        [~,XT3]=cfar_os(xc,N_det,k,N_pro,Pfa);
        for j=1:M
            if(xc(1,j)<XT3(1,j))
                xc(1,j)=0;
            end
        end
        cfar_RVmap(i,:)=xc;
    end
    N_det=12;
    N_pro=2;
    k=N_det/2;
    Pfa=0.2;
    for j=1:M  
        xc=reshape(RVmap(:,j),[1,N]);
        [~,XT4]=cfar_os(xc,N_det,k,N_pro,Pfa);
        for i=1:N
            if(xc(1,i)<XT4(1,i))
                cfar_RVmap(i,j)=0;
            end
        end
    end
    %% peak extraction
    peak_judge=zeros(3,3);
    for i=3:N-2
        for j=3:179
            if(cfar_RAmap(i,j)>0)
                peak_judge=(cfar_RAmap(i-2:i+2,j-2:j+2));
                if(cfar_RAmap(i,j)~=max(max(peak_judge)))
                    cfar_RAmap(i,j)=0;
                end
            end
        end
    end
    %% match doppler
    %limited
    RAmap_v=zeros(size(RAmap_MVDR));%save radial velocity
    for i=1:N
        for j=1:181
            if(cfar_RAmap(i,j)>0)
                [peak,v_temp]=max(cfar_RVmap(i,:));
                if(peak>0)
                    RAmap_v(i,j)=(v_temp-M/2)*lambda/Tc/M/2;
                else
                    RAmap_v(i,j)=0;
                end
            end
        end
    end       

    %% clustering
    cluster1=zeros(size(RAmap_MVDR));
    num=0;%point count
    k=1;
    epsilon=1;%radius
    MinPts=3;%min points
    for i=3:N %exclude too close data 
        for j=31:151% save +-60 angle data
            if(cfar_RAmap(i,j)>0)
                cluster1(i,j)=1;
                num=num+1;
            end
        end
    end
    %transform to rectangular coordinate system 
    cluster2=zeros(num,4);%save point set
    for i=1:N  
        for j=1:181 
            if(cluster1(i,j)==1)
                R_temp=i/N*fs*c/2/K;
                theta_temp=j-90;
                cluster2(k,1)=R_temp*sind(theta_temp);
                cluster2(k,2)=R_temp*cosd(theta_temp);
                cluster2(k,3)=RAmap_MVDR(i,j);%save amplitude
                cluster2(k,4)=RAmap_v(i,j);%save velocity
                k=k+1;
            end
        end
    end
    
    IDX1=Clustering(cluster2,epsilon,MinPts);
    
    %find xx frame:feature vector of cluster_IDX=i
   for i=0:max(IDX1)
        Xi=cluster2(IDX1==i,:);
        if i~=0
            %retangular range
            [Xmax,~]=max(Xi(:,1));
            [Xmin,~]=min(Xi(:,1));
            [Ymax,~]=max(Xi(:,2));
            [Ymin,~]=min(Xi(:,2));
            L_cluster=Xmax-Xmin;
            W_cluster=Ymax-Ymin;
            %mean
            x_mean=mean(Xi(:,1));
            y_mean=mean(Xi(:,2));
            %standard deviation
            %x_std=std(Xi(:,1));
            %y_std=std(Xi(:,2));
            %density
            %density=numel(Xi(:,1))/(L_cluster*W_cluster);
            %mean amplitude
            amp_mean=mean(Xi(:,3));
            %mean velocity
            v_mean=mean(Xi(:,4));
            if(v_mean==0)
                v_mean=0.01;
            end
            %save for matching
            cluster_feature(frame_count,i,:)=[x_mean,y_mean,v_mean,L_cluster*W_cluster,Xmax,Xmin,Ymax,Ymin,amp_mean].';
        end
    end
   
%% Figure showing
    %range-doppler map
    figure
    colormap(jet(256))
    Y=((0:N-1)/N*fs)*c/2/K;
    X=(-M/2:M/2-1)*lambda/Tc/M/2;
    imagesc(X,Y,RVmap);
    colorbar;
    xlabel('velocity[m/s]')
    ylabel('range[m]')
    axis xy
    title(['range-doppler map in Frame: ' num2str(xx)])
    %saveas(gcf,strcat('D:\SCI_video\RD_map\',int2str(xx),'.jpg'));
    
    % range_angle map
    figure
    X=((0:N-1)/N*fs)*c/2/K;
    Y=(0:180)-90;
    colormap(jet(256))
    imagesc(Y,X,RAmap_MVDR);
    colorbar;
    xlabel('angle [degree]')
    ylabel('range [m]')
    axis xy
    title(['Range-Angle Rectangle Heatmap in Frame: ' num2str(xx)])
    %saveas(gcf,strcat('D:\SCI_video\RA_map\',int2str(xx),'.jpg'));
    
%     transform polar coordinates to Cartesian coordinates 
%     Display clustering results
    figure
    X=((0:N-1)/N*fs)*c/2/K;
    Y=(0:180)-90;
    YY = X'*cosd(ThetaSet);
    XX = X'*sind(ThetaSet);
    colormap(jet(256))
    h = pcolor(XX, YY, RAmap_MVDR);
    colorbar;
    shading flat
    axis equal;
    xlabel('Cross-Range [m]');
    ylabel('Range [m]');
    yLim = [0,X(end)];
    xLim = yLim(2)*sin(max(abs(ThetaSet))) * [-1,1];
    ylim(yLim);
    xlim(xLim);
    set(gca, 'Xtick', [-50:1:50]);
    set(gca, 'Ytick', [0:1:100]);
    title(['Range-Angle Polar Heatmap in Frame:' num2str(xx)])
    %saveas(gcf,strcat('D:\SCI_video\polar\',int2str(xx),'.jpg'));
    %grid on;
    
    figure
    X=((0:N-1)/N*fs)*c/2/K;
    Y=(0:180)-90;
    YY = X'*cosd(ThetaSet);
    XX = X'*sind(ThetaSet);
    colormap(jet(256))
    h = pcolor(XX, YY, RAmap_MVDR);
    colorbar;
    shading flat
    axis equal;
    xlabel('Cross-Range [m]');
    ylabel('Range [m]');
    yLim = [0,X(end)];
    xLim = yLim(2)*sin(max(abs(ThetaSet))) * [-1,1];
    ylim(yLim);
    xlim(xLim);
    set(gca, 'Xtick', [-50:1:50]);
    set(gca, 'Ytick', [0:1:100]);
    grid on;
    alpha(0.1)
    hold on
    PlotClusterinResult(cluster2,squeeze(cluster_feature(frame_count,:,:)),IDX1,xx);%cluster2是点集，cluster_feature为聚类集
    %saveas(gcf,strcat('D:\SCI_video\cluster\',int2str(xx),'.jpg'));
    
    % After CFAR
    figure
    X=((0:N-1)/N*fs)*c/2/K;
    Y=(0:180)-90;
    colormap(jet(256))
    imagesc(Y,X,cfar_RAmap);
    colorbar;
    xlabel('angle [degree]')
    ylabel('range [m]')
    axis xy
    title(['CFAR RA Rectangle Heatmap in Frame: ' num2str(xx)])
    %saveas(gcf,strcat('D:\SCI_video\CFAR\',int2str(xx),'.jpg'));
    
    % Doppler velocity
    figure
    X=((0:N-1)/N*fs)*c/2/K;
    Y=(0:180)-90;
    imagesc(Y,X,RAmap_v);
    colorbar;
    xlabel('angle [degree]')
    ylabel('range [m]')
    axis xy
str='continue:1 stop:0 \n';
stop_flag=input(str);
close all
frame_count=frame_count+1;
end
toc
tic