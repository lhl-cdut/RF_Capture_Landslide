%%
clc
clear all
close all
%%
% Reading .sac files for single-component seismic data
filename ='HHZ.sac';
tic
% Reading a .sac file
sacdata = rdsac(filename);
% Getting data
data = sacdata.d;
% data=data(1:15000);
Fs = 1./sacdata.HEADER.DELTA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data preprocessing
fl=1;
fh=14;
wp=[fl./(Fs/2) fh./(Fs/2)];
[b,a]=butter(3,wp,'bandpass');
outhi = filter(b,a,data);
axis tight
% downsampling
outhi = resample(double(outhi),30,round(Fs));
outhi=outhi(500:end);
plot(outhi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_size = 500;
data_len = 6000;
num_cols = floor((length(outhi) - data_len) / step_size) + 1;
data_t = zeros(data_len, num_cols);
for i = 1:num_cols
    start_index = (i-1) * step_size + 1;
    end_index = start_index + data_len - 1;
    data_t(:, i) = outhi(start_index:end_index);
end
[~,iter]=size(data_t);
toc
%%
tic
ATTRIBUTES=zeros(iter,20);  %iter is the total number of data (all data have been preprocessed)
for i=1 : iter
    signal_dataZ=data_t(:,i);
    Fs=30;
    SIG=signal_dataZ;
    dataZ_nor=SIG/max(abs(SIG));
    %---------A(1-3)---------------
    L=length(SIG);
    N=100;
        for j=1:1:L-N
            iUy(j,1)=mean(abs(dataZ_nor(j:(j+N-1)))); %MA
            iUy2(j,1)=mean(dataZ_nor(j:(j+N-1)).^2);
            SI(j,1)=sqrt((iUy2(j)-iUy(j).^2)/iUy(j).^2);
        end
    uSI = mean(SI);%uSI
    SIR = max(SI)./uSI;%A(1)
    MAR = max(iUy)./SIR; %A(2)
    CC=conv(SI,iUy);
    LL=(max(CC)/mean(CC))./std(CC);
    NN=mean(CC)./std(CC);
    LAMDA=MAR./std(SI);%A(3)

    %---------A(4)---------------Mean
    X_ba=mean(abs(dataZ_nor));%A(4)
    X_nor=SIG./std(SIG);

    %---------A(5)---------------Kurtosis
    KurtoSig=kurtosis(dataZ_nor);%A(5)

    %---------A(6)---------------H
    H_shannon=-sum(dataZ_nor.*log2(dataZ_nor));%Shannon Entropy
    H_shannon=abs(H_shannon);
    %---------A(7)---------------std(AIC)
    [Rind,RAic] = M_aic(dataZ_nor); %M-AIC
    STD_AIC=std(RAic(5:5995));% 

    %---------A(8)--------------- 
    F_W1 = sum(abs(dataZ_nor))/length(SIG);
    F_W2 = sqrt(sum(dataZ_nor.^2)/length(SIG));
    W = F_W2./F_W1;%

    %---------A(9-11)---------------DISTQ2Q1  Spectral centroid，SC
    NyF=Fs/2;% Nyquist Frequency
    n=2*Fs;
    Freq1=linspace(0,1,n/2)*(Fs/2);% Frequency array 
    SpecWdow=20; % Length of th windos used to compute the spectrogram (samples)
    noverlap=0.90*SpecWdow; % overlap of the moving window
    SPEC=filter(ones(1,10)./10,1,(abs(spectrogram(SIG,SpecWdow,noverlap,Freq1,Fs))),[],1); % Compute and smooth Spectro
    Xcentroid=sum(repmat(1:size(SPEC,1),size(SPEC,2),1)'.*SPEC,1)./...
        sum(SPEC,1); % Compute centroid of each DFT in Spec
        for kk=1:round(size(SPEC,2))
            Xquart1(kk)=sum((1:length(SPEC(1:round(Xcentroid(kk)),kk)))'.*SPEC(1:round(Xcentroid(kk)),kk),1)./...
                sum(SPEC(1:round(Xcentroid(kk)),kk)); % Compute Q1 of each DFT in Spec
        end
        for kk2=1:round(size(SPEC,2))
            Xquart3(kk2)=round(Xcentroid(kk2))+sum((1:length(SPEC(round(Xcentroid(kk2)):end,kk2)))'...
                .*SPEC(round(Xcentroid(kk2)):end,kk2),1)./sum(SPEC(round(Xcentroid(kk2)):end,kk2)); % Compute Q3 of each DFT in Spec
        end
    FREQCENTER=Freq1(round(Xcentroid)); % Curve of the centroid position over time 
    FREQQ1=Freq1(round(Xquart1));
    FREQQ3=Freq1(round(Xquart3));
    DISTQ2Q1=mean(abs(FREQCENTER-FREQQ1)); % Distance Q2 curve to Q1 curve
    F_DISTQ3Q2=mean(abs(FREQQ3-FREQCENTER)); % Distance Q3 curve to Q2 curve
    F_DISTQ3Q1=mean(abs(FREQQ3-FREQQ1)); % Distance Q3 curve to Q1 curve

    %---------A(12-15)---------------F_ES(Energy between 0.1 and 1 Hz)
    NyF=Fs/2;% Nyquist Frequency
    DATAF=cell(1,5);
    % DATAF=zeros(length(SIG),5);
    %     [Fa,Fb]=butter(3,[0.1/30 1/30],'bandpass'); % Butterworse Filter with normalized frequency, 
    %     DATAF{1}=filtfilt(Fa,Fb,double(SIG)); % acausal filterting
    %     F_ES=log10(trapz(abs(hilbert(DATAF{1})))); % compute energy in several frequency band

    FilterFrequencyI=[0.1,  1,  5,  10]; % to test
    FilterFrequencyE=[  1,  5,  10, 14]; % to test
        for jj=1:length(FilterFrequencyI) 
            [Fa,Fb]=butter(3,[FilterFrequencyI(jj)/NyF FilterFrequencyE(jj)/NyF],'bandpass'); % Butterworse Filter with normalized frequency, 1=Nyquist
            DATAF{jj}=filtfilt(Fa,Fb,double(SIG)); % acausal filterting
            F_ES(jj)=log10(trapz(abs(hilbert(DATAF{jj})))); % compute energy in several frequency band
        end
    % Ew    

    %---------A(16)---------------Margin factor Detection of shock signals in the signal
    kf=(max(SIG)-min(SIG))./mean(sqrt(abs(SIG)))^2;

    %---------A(17)---------------Spectrum Mean
    FFT_sig= fft(SIG(1:6000));
    A12=mean(abs(FFT_sig(1:6000/2))./max(abs(FFT_sig(1:6000/2))));

    %---------A(18)---------------Entropy of a permutation
    x = SIG;
    a=1;
        for j=1:1:5950
            X= x(j:j+50-1)';
            [pe hist] = pec(X, 4, 5);
            plist(a) = pe;
            a=a+1;
        end
    A13 = max(plist)-min(plist);
    %---------A(19)---------------frequency-wave number 
    fk = abs(fft2(X_nor));
    fk=std(fk)./mean(fk);

    % end
    ATTRIBUTES(i,1) = SIR;
    ATTRIBUTES(i,2) = MAR;
    ATTRIBUTES(i,3) = LAMDA;
    ATTRIBUTES(i,4) = X_ba;
    ATTRIBUTES(i,5) = KurtoSig;
    ATTRIBUTES(i,6) = H_shannon;
    ATTRIBUTES(i,7) = STD_AIC;
    ATTRIBUTES(i,8) = W;
    ATTRIBUTES(i,9) = DISTQ2Q1;
    ATTRIBUTES(i,10) =F_DISTQ3Q2;
    ATTRIBUTES(i,11) = F_DISTQ3Q1;
    ATTRIBUTES(i,12) = F_ES(1);
    ATTRIBUTES(i,13) = F_ES(2);
    ATTRIBUTES(i,14) = F_ES(3);
    ATTRIBUTES(i,15) = F_ES(4);
    ATTRIBUTES(i,16) = kf;
    ATTRIBUTES(i,17) = A12;
    ATTRIBUTES(i,18) = A13;
    ATTRIBUTES(i,19) = fk;
    ATTRIBUTES(i,20) = NN;



    if mod(i,10)  == 0
        i                  % Output convergence progress
    end


end
% save('...\name.mat ','ATTRIBUTES');
% clear ATTRIBUTES
toc
%% predict data plot
% Read the 'ps_input' first.  
load('Example_models.mat'); % Load the model trained in the 'model_training' session
tic
    Xtrain = mapminmax('apply', ATTRIBUTES', ps_input );
    [trainPred, Vote1] = classRF_predict(Xtrain', model);% 'model': the model you've finished training

    %-----------figure----------------

% figure('Position',[500,200,1000,550])
color1 = [120 000 001]./255;% Burgundy
color2 = [193 018 033]./255;% Big Red
color3 = [000 047 073]./255;% Navy
color4 = [102 155 187]./255;% Indigo
color5 = [254 240 213]./255;% Flesh
color6 = [242,242,230]./255;% Flesh
color=[[120 000 001]./255;[193 018 033]./255;[000 047 073]./255;[102 155 187]./255;[254 240 213]./255;[242,242,230]./255];
type='compact';

    ax1=subplot(2,1,1);hold on
    N=1:length(outhi);%filter_dataZ data
    box on;
    plot(N,outhi,'Color',color(3,:),'linewidth',1);
    axis tight
    set(gca,'FontSize',20,'Fontname', 'Times New Roman','ytick',[],'xtick',[]);
    set(gca,'Position',[.1 .55  0.87 .4]);

ax2=subplot(2,1,2);hold on
    for i = 1:iter
        if trainPred(i)== 1 %&& trainPred(i+3) ==1 
        xline(ax2,i,'LineWidth',0.2,'LineStyle','-','Color',color(2,:));
    %     elseif trainPred(i)== 2 %&& trainPred(i+5) ==1 %&& trainPred(i+9) ==1
    %     xline(ax2,i,'LineWidth',0.2,'LineStyle','-','Color',color(4,:))    
        else
        xline(ax2,i,'LineWidth',0.2,'LineStyle','-','Color',color(6,:));
        end 

    end
toc
ax2.TickDir = 'out'; % Set the scale direction to face outward
ax2.XAxisLocation = 'bottom'; % Set the x-axis position at the bottom
axis tight
set(gca,'FontSize',20,'Fontname', 'Times New Roman','ytick',[]);
set(gca,'FontSize',20,'Fontname', 'Helvetica');
set(gca,'Position',[.1 .15  0.87 .4]);

