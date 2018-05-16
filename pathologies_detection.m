
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

R_value = R_wave_detection( data, Fs);
%% Rythm cardiac

sumR = 0;

for i=1:length(R_value)-1
   R_R(i) = R_value(i+1)-R_value(i);
   sumR = sumR + R_R(i);
end



rythm = sumR / (length(R_value)-1);

