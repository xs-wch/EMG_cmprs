%%%% for test wavmat
%%%% dwt with matrix multiplication form
clear all
close all
%%%% setup for wavelab toolbox
% eeglab
%  
% addpath ./wavelab/Wavelab850
% WavePath
% 
% addpath ./package/utilities
%%%% load data
EMG = loadcnt('/home/chaohua/Documents/EMG/emg/basic/aslw.cnt','dataformat','int32');
fs = 1000;
sig = EMG.data(1,:)';
%%%% filter - hightpass: pass band 6Hz-0.5*fs
fp1 = 6*2/fs;fst1 = 4*2/fs;
B1 = firls(120,[0 fst1 fp1 1],[0 0 1 1]);
%fvtool(B1);
sig = filter(B1,1,sig);
% [B2,A2] = iirnotch(50*2/fs,3*2/fs);
% fvtool(B2,A2);
% sig = filter(B2,A2,sig);
% figure
% subplot(121)
% plot(sig);
% subplot(122)
% plot(abs(fft(sig)));

testsig = sig(6000:6000+1024*19-1);
% figure
% subplot(1,2,1)
% plot(testsig);
% subplot(1,2,2)
% plot(abs(fft(testsig)))

len = 256;
N = log2(len);
%%%% select wavelet: orthognal wavelet basis used  
W = WavMat(MakeONFilter('Vaidyanathan'),len,N);%(db10)
%figure
%imagesc(W*W');
% rng(78)
% wt = rand(len,1);
% I = wt < 0.8;
% sum(I)
% wt(I) = 0;
% testsig = W'*wt;
%testsig = sin(2*pi*40/fs*[0:len-1]+pi/3).*sin(2*pi*5/fs*[0:len-1]+pi/5);%+0.3*randn(1,len);
%testsig = testsig';
wt = W*testsig(1:len);
figure;bar(wt);

%%%% encode
[CS_data,sensematrix] = EMG_CS_encode(testsig,len,128);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% transmit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% decode
recover_data = EMG_CS_decode_EB(CS_data,sensematrix,W');

figure
subplot(2,1,1)
plot(testsig);
xlabel('time(ms)')
ylabel('amplitude(\muV)')
title('original')
subplot(2,1,2)
plot(recover_data);
xlabel('time(ms)')
ylabel('amplitude(\muV)')
title('recovered:compress ratio 512->256')

