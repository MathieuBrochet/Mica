%% Band pass filter 
%% creation of low pass filter 
Ts = 1 / Fs; 

H1_num = [1 0 0 0 0 0 -2 0 0 0 0 0 1]; %numerator of the low pass filter 

H1_den = [1 -2 1 0 0 0 0 0 0 0 0 0 0]; %denominator of the low pass filter 

h1 = filter(H1_num, H1_den,ecg); % out signal of the low pass filter 

plot(h1);

%plot(h1);
%% creation of high pass filter 

H2_num = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; % numerator  of the high pass filter 

H2_den = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % denominator of the high pass filter 

h2 = filter(H2_num,H2_den,h1); % out signal of the high pass filter 

%% Derivative 
Vect_caus = [ 0 0 0 0 1 ];
Deriv_num = [-1 -2 0 2 1]; % numerator of derivative

Deriv_num = conv(Vect_caus,Deriv_num); % we transform the signal in a causal signal thank to Vect_caus 

h3 = filter(Deriv_num,8*Ts,h2); % derivative filter 



%plot(h3);

%% | |^2

h4 = abs(h3); % taking the absolute value of the h3 

%plot(h4);



%% Moving Window Integration 

h = ones(1,30)/30; % creation of the window 

h5 = conv(h4,h); % apply the window to the signal h4 

h5 = h5 / max(abs(h5)); % normalisation of the signal h5 

%plot(h5);
%% thresholding 







%% detection of maxima 