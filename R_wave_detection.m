
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

tresh = 3e13;
p = (h5 > tresh);

% figure(1);
% plot(h5);
% hold all;
% plot(1e13*p);
% 


%% detection of maxima

mult = h5.*p;
 
mult1 = data(1:length(p)).*p;


% figure(2);
% plot(1e11*data);
% hold all;
% plot(mult);

[pks, abscisse] = findpeaks(mult1);
pks(pks < 0) = 0;


R_value = abscisse(pks > 0); % value of R peak abcsissa 
% R_value is the abscissa of R peaks and pks are the ordinates of R peaks



%% Q and S wave detection


% figure(2);
% plot(data);

deriv_data = diff(data);

% figure(3);
% plot(deriv_data);

[pks1, abscisse1] = findpeaks(-data); % display the pks1 and their abcissa 

for i=1:length(R_value) 
    j = 1;
    while abscisse1(j) < R_value(i)
        j=j+1;
    end;
    S_value(i) = abscisse1(j);
    Q_value(i) = abscisse1(j-1);
end;


%% P and T wave detection



G1 = [ 1 0 0 0 0 0 -1 ];   % delay = 3

g1 = filter(G1,1,data);   % First filter


G2_num = [ 1 0 0 0 0 0 0 0 -1 ];
G2_den = [ 1 0 -1 ];                 % delay = 3

g2 = filter(G2_num, G2_den, g1);  % Second filter



delay2 = [0 0 0 0 0 0 1];
data = conv(data,delay2);   % correction of the 2 filters delays

% figure(5);
% plot(data);
% hold all; 
% plot(g1);
% 
figure(6);
plot(data);
hold all; 
plot(g2);
grid on; 



for p=1:length(R_value)-1
    k=1;
    j=R_value(p);
    while j < R_value(p)+0.7*(R_value(p+1)-R_value(p))
        if g2(j)*g2(j+1) < 0
            t_vect(k)=j;
            k=k+1;
        end
        j=j+1;
    end;
    %T_value(i)=max(t_vect);
end;



%% 




