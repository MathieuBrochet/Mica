
% R_wave_detection.m :
% filter between pass-band filter and detection of maxima
% detection of R_value
% detection of P & T value
% detection of Q & S value
% displaying all signal during the fitrage


function [R_value, Q_value, S_value, P_value, T_value,tresh] = R_wave_detection( data, Fs)
%% Band pass filter
%% creation of low pass filter
Ts = 1 / Fs;  % definition of sample period

H1_num = [1 0 0 0 0 0 -2 0 0 0 0 0 1]; %numerator of the low pass filter

H1_den = [1 -2 1 0 0 0 0 0 0 0 0 0 0]; %denominator of the low pass filter
h1 = filter(H1_num, H1_den,data); % out signal of the low pass filter

%% creation of high pass filter

H2_num = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; % numerator  of the high pass filter

H2_den = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % denominator of the high pass filter

h2 = filter(H2_num,H2_den,h1); % out signal of the high pass filter


%% Derivative
Deriv_num = [-1 -2 0 2 1]; % numerator of derivative

h3 = filter(Deriv_num,8*Ts,h2); % derivative filter

%% | |^2

h4 = abs(h3).^2; % taking the absolute value of the h3

%% Moving Window Integration

h = ones(1,21)/21; % creation of the window with 20 samples per 0.1s (QRS complex)

h5 = conv(h4,h); % apply the window to the signal h4
delay = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0 1]; % delay = 5+16+2+10
data = conv(delay,data); % to cunter the delay
% after the MWI, the time duration of the rising edge is 20 samples,
% corresponding to 0.1s. So it corresponds to the signal QRS.


%% thresholding
%  we choose to threshold the signal with the average of the signal after
%  the moving window integration
h5_moy =0 ;
for i=1:length(h5)
    h5_moy = h5_moy +h5(i);
end


h5_moy = h5_moy / length(h5);
tresh = h5_moy;   % treshold is the average of the MWI

p = (h5 > tresh); % p is the vector which contain all the value of h5 superior to tresh







%% detection of maxima
mult1 = data(1:length(p)).*p;

[pks, abscisse] = findpeaks(mult1); % find the peaks R : pks , and the abcisse of the peaks R
pks(pks < 0) = 0;


R_value = abscisse(pks > 0); % value of R peak abcsissa
% R_value is the abscissa of R peaks and pks are the ordinates of R peaks



%% Q and S wave detection

[pks1, abscisse1] = findpeaks(-data); % display the pks1 of -data and their abcissa, to
% consider Q and S waves as local maximum

for i=1:length(R_value)
    j = 1;
    while abscisse1(j) < R_value(i)
        j=j+1;
    end;
    S_value(i) = abscisse1(j); % max after R_value is S_value
    Q_value(i) = abscisse1(j-1);% max before R_value is Q_value
end;


%% P and T wave detection

% first filter for P & T value
G1 = [ 1 0 0 0 0 0 -1 ];   % first filter gives a delay = 3
g1 = filter(G1,1,data);


% second filter for P & T value
G2_num = [ 1 0 0 0 0 0 0 0 -1 ]; % numerator
G2_den = [ 1 0 -1 ];             % denominator
g2 = filter(G2_num, G2_den, g1); % second filter gives a delay =


% correction of the 2 filters delays
delay2 = [0 0 0 0 0 0 1];    % delay = 7
data2 = conv(data,delay2);   % cunter the delay



% case of T_values
for p=1:length(R_value)-1
    k=1;
    j=R_value(p);
    t_vect=[];
    while j < R_value(p)+0.7*(R_value(p+1)-R_value(p)) % value between R and 0.7*R-R'
        if g2(j)*g2(j+1) < 0 % if two folowing sample have different sign
            t_vect(k)=j;   % t_vect contain all the zeros of g2
            k=k+1;
        end
        j=j+1;
    end;
    max_ord = max(data(t_vect(2:length(t_vect))));
    for m = 1:length(t_vect) % max value in data of the zeros
        if data(t_vect(m))==max_ord
            T_value(p)=t_vect(m); %abcissa of T_value
        end
    end
end



% case of P values
for p=2:length(R_value)
    k=1;
    j=R_value(p);
    p_vect=[];
    while j > R_value(p)-0.3*(R_value(p)-R_value(p-1)) % value between R and 0.7*R-R'
        if g2(j)*g2(j+1) < 0 % if two following sample have different sign
            p_vect(k) = 0;
            p_vect(k)=j;   % t_vect contain all the zeros of g2
            k=k+1;
        end
        j=j-1;
    end
    max_ordp = max(data(p_vect(2:length(p_vect))));
    for m = 1:length(p_vect)                            % max value in data of the zeros
        if data(p_vect(m))==max_ordp
            P_value_int(p)=p_vect(m);
        end
    end
    P_value = P_value_int(2:length(P_value_int)); %abcissa of P_value
end


%% Display the signal

% display h1

%figure(1);
%plot(h1);
%title('Signal after low-pass filter');
%xlabel('Time');
%ylabel('signal h1');

% display h2

%figure(2);
%plot(h2);
%title('Signal after high-pass filter');
%xlabel('Time');
%ylabel('signal h2');

% display h3

%figure(3);
%plot(h3);
%title('Signal after derivative filter');
%xlabel('Time');
%ylabel('signal h3');

% display h4

%figure(4);
%plot(h4);
%title('Signal after square module filter');
%xlabel('Time');
%ylabel('signal h4');

%display h5 with p

% figure(5);
% plot(h5);
% hold all;
% plot(1e13*p);
%title('Superposition of h5 and p');
%xlabel('Time');
%ylabel('Amplitude of h5 and 10^13*p');

%display g1 (first filter of P & T detection method) with the data

% figure(6);
% plot(data);
% hold all;
% plot(g1);
%title('Superposition of g1 and data');
%xlabel('Time');
%ylabel('Amplitude of g1 and data');


%display g2 (second filter of P & T detection method) with the data

% figure(7);
% plot(data);
% hold all;
% plot(g2);
% grid on;
%title('Superposition of g2 and data');
%xlabel('Time');
%ylabel('Amplitude of g2 and data');


end % the function enable to have access to variable that we needs in pathologies_detection.m