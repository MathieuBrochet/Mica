
clear; 
close all; 
clc;


%% Load a signal
[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = signal.ecg; % Your ecg data
Fs = signal.Fs; % Sampling frequency
N = size(data,2); % Data length
time_axis = (1:N)/Fs;

%% Band pass filter 
%% creation of low pass filter 
Ts = 1 / Fs; 

H1_num = [1 0 0 0 0 0 -2 0 0 0 0 0 1]; %numerator of the low pass filter 

H1_den = [1 -2 1 0 0 0 0 0 0 0 0 0 0]; %denominator of the low pass filter 
h1 = filter(H1_num, H1_den,data); % out signal of the low pass filter 

%plot(h1);

%plot(h1);
%% creation of high pass filter 

H2_num = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; % numerator  of the high pass filter 

H2_den = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % denominator of the high pass filter 

h2 = filter(H2_num,H2_den,h1); % out signal of the high pass filter 


%% Derivative 
Deriv_num = [-1 -2 0 2 1]; % numerator of derivative

h3 = filter(Deriv_num,8*Ts,h2); % derivative filter 



%plot(h3);

%% | |^2

h4 = abs(h3).^2; % taking the absolute value of the h3 

%plot(h4);



%% Moving Window Integration 

h = ones(1,21)/21; % creation of the window with 20 samples per 0.1s (QRS complex)

h5 = conv(h4,h); % apply the window to the signal h4 
delay = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0 1]; 
data = conv(delay,data); % delay = 5+16+2+10



% after the MWI, the time duration of the rising edge is 20 samples,
% corresponding to 0.1s. So it corresponds to the signal QRS.


%% thresholding 

tresh = 1e13;
p = (h5 > tresh);

figure(1);
plot(h5);
hold all;
plot(1e13*p);

% max_h5 = max(h5);
% mean5 = mean (h5 );
% poss_reg =(h5 > mean5*max_h5)';
%  
% left = find(diff([0 poss_reg])==1);
% right = find(diff([poss_reg 0])==-1);
%  
% left=left-(6+16);  % cancle delay because of LP and HP
% right=right-(6+16);% cancle delay because of LP and HP
%  
% for i=1:length(left)
%     [R_value(i) R_loc(i)] = max( x1(left(i):right(i)) );
%     R_loc(i) = R_loc(i)-1+left(i); % add offset
%  
%     [Q_value(i) Q_loc(i)] = min( x1(left(i):R_loc(i)) );
%     Q_loc(i) = Q_loc(i)-1+left(i); % add offset
%  
%     [S_value(i) S_loc(i)] = min( x1(left(i):right(i)) );
%     S_loc(i) = S_loc(i)-1+left(i); % add offset
%  
% end
%  
% % there is no selective wave
% Q_loc=Q_loc(find(Q_loc~=0));
% R_loc=R_loc(find(R_loc~=0));
% S_loc=S_loc(find(S_loc~=0));
%  
% figure
% subplot(2,1,1)
% title('ECG Signal with R points');
% plot (t,x1/max(x1) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o');
% legend('ECG','R','S','Q');
% subplot(2,1,2)
% plot (t,x1/max(x1) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o');
% xlim([1 3])


%% detection of maxima





