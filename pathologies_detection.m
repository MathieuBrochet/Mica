function [bpm,perc_of_p_value_AF,perc_sample_brady,perc_sample_tachy,percent_of_extopic_beat,gamma] = pathologies_detection( data, Fs)



%% Load a signal

% add R,Q,S,P,T_value to file pathologies_detection.m
[R_value, Q_value, S_value, P_value, T_value] = R_wave_detection( data, Fs);


%% Rythm cardiac

delay = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0 1]; 
data = conv(delay,data);
% we have to insert the delay an other time because it load the data from
% the file .mat and not with the correct delay 

sumR = 0;
for i=1:length(R_value)-1
   R_R(i) = R_value(i+1)-R_value(i);
   sumR = sumR + R_R(i);
end
R_average = sumR / (length(R_value)-1); % calcul of the R average of the signal 

%calcul du bpm 
total_ech_nb = length(data);  % total number of sample in data

%on a f_ech =200 ech/s
time_of_ech_total = total_ech_nb/Fs;  % time corresponding to all those sample
nb_bat_total = length(R_value); % number of beat on the signal 
nb_bat_par_sec = length(R_value)/ time_of_ech_total; % number of beat per second
bpm = nb_bat_par_sec *60; % number of beat per minute




%% Ectopic beat 
% delta(n) - delta(n-1) , cunter_e_beat compte le nombre de ectopic beat

Ectopic_beat=zeros();
k=1;
for j=1: length(R_R)-1
    Ect1(j) = abs(R_R(j+1)-R_R(j)); % minus between two consecutives segment  R_R
    if Ect1(j) >= R_average % condition for an ectopic beat anomaly
        Ectopic_beat(k) = Ect1(j); % vector which contain all Ectopic position
        k=k+1;
    end
end
cunter_e_beat= length(Ectopic_beat);  % number of ectopic beat
percent_of_extopic_beat = cunter_e_beat *100 / length(R_R); % percentage of ectopic beat 


%%  Atrial Fibrillation 
 % the autocovariance of a white gaussian noise is a Dirac, cause all of
 % those sample have no corrélation betwwen them.
 %If the autocariance is similar to white gaussian noise, that is to say it
 %is represent by a Dirac, the ecg is atrial fibrillate
N=length(R_R);
auto_cov_fibr_k1=0;
auto_cov_fibr_k2=0;
for k=1:N-1
    R_k=0;
    for n=1:N-k-1
        auto_cov_fibr_k1=(R_R(n+k)-R_average);
        auto_cov_fibr_k2=(R_R(n) - R_average);
        auto_cov_fibr = auto_cov_fibr_k1*auto_cov_fibr_k2;
        R_k = R_k + auto_cov_fibr;  
    end
     gamma(k)=R_k/(N-k-1); % contain all the covariance of the signal 
end
 
% It is the first condition to diagnosticate an Atrial fibrillation
 
%% detection du % de pic P valable

P_value_real = P_value( P_value >0); % consider the P_value not equal to 0
R_value_real = R_value ( R_value >0);% consider the R_value not equal to 0

P_ord = data(P_value_real); % P ordonate
R_ord = data (R_value_real); % R_ordonate
P_available = zeros(); % going to contain the existing P_values


for i=1:length(P_value_real)
   if (abs(R_ord(i)/P_ord(i))) > 5 && ((R_ord(i)/P_ord(i))) < 400 % on a 'normal' ecg, the quotient of 
       % R_ord/P_ordi is 95% of the time between 5 and 400. So if a
       % frequent occurence of P_value disrespect that, we consider they
       % are not available
       P_available(i) = abs(P_ord(i));
       frac(i)=abs(R_ord(i)/P_ord(i));
   else
       P_available(i) = 0 ;
       frac(i)=abs(R_ord(i)/P_ord(i));
   end
end

P_available_ss_zero = P_available(P_available >0); % ordonate of P values available

nb_p_values = length(P_value_real); %number of all the 'p_value' with ones which are false
nb_p_values_real = length(P_available_ss_zero); %number of available P_value
perc_of_p_value_AF = (nb_p_values_real/nb_p_values)*100; % percentage of P_value on the signal 

% that is the second condition to diagnosticate an Atrial fibrillation 

 %% Ventricular fibrillation 
% absence of traditional P Q R S T waves , 

% import les Q S T values et comparé leur valeur à un seuil pour dire qu'il
% y a une abscence de P,Q,R,S,T


% % comparaison to a pure sine, with

% %creation of a pure sine 
% % f= 240 to 600 bpm 

N_value = length(data);
for f=4:10 % bpm 240 to 600
    comparaison_signal = sin(2*pi*f*(0:N_value-1)/Fs); % pure sine
end

% % method de comparaison avec autocov de ecg, et auto cov de sine  

y=autocorr(comparaison_signal); 
  
%  %determienr une periode d'analyse du signal  pour la parole c'est 10-20
%  %ms pour 20-20kHz 


%% Releving interessant part of the signal with a Window 

decalage_sample=200; % sample offset
partie_entiere = floor(length(data)/2);
BPM =zeros(1,partie_entiere);  % vector which will contain all the bpm calculated with all the window 
n=0; %initialisation of n
window = zeros(1,length(data)); % initialisation of the window 
for m=1:decalage_sample:(length(data)/2) % calcul of bpm all 200 sample , cause under 200 sample offset there is no significative difference on the bpm 
    n = n+decalage_sample;
    for k=n:((length(data)/2) + n-1) 
        window(k) = 1; % creation of the window , augmention by one the k'th sample to move the window on the right way
    end
    for k=1:n
        window(k) = 0; % diminution of the k th sample to move the window on the right 
    end 
    data_value_window = zeros(length(data),1);


    for p=1:length(data)
        data_value_window(p)= window(p).*data(p); % on recupere seulement les r values de la fenetre 
    end 

% Find the R peaks 

    new_vect_window_R = zeros();
    for k=1:length(R_value)
        new_vect_window_R(k)=data_value_window(R_value(k)); %contient les pics R sur la fenetre 
    end
  
    R_window = new_vect_window_R(new_vect_window_R > 0);
% Calcul of the bpm with all the window 

    total_ech_nb_window = length(data)/2;
% We have 200 sample per second  
    time_of_ech_total_window = total_ech_nb_window/Fs; 
    nb_bat_total_window = length(R_window);
    nb_bat_par_sec_window = length(R_window)/ time_of_ech_total_window; 
    bpm_window = nb_bat_par_sec_window *60;  
    BPM(m) = bpm_window;
end 

BPM_sans_zero = BPM( BPM >0); % value of the bpm for all the window , with offset of 200 sample

BPM_average = 0; %initialisaiton 

for i = 1:length(BPM_sans_zero)
    BPM_average = BPM_average + BPM_sans_zero(i) ; 
end

BPM_average = BPM_average / length(BPM_sans_zero); % Bpm average of all the bpm calculated 
% on the previous lines

%% Bradycardia 
% detection of bradycardia with the bpms on BPM_sans_zero which are under
% 60
cunter_brady = 0; 
for i=1:length(BPM_sans_zero)
    if bpm < 60
        cunter_brady= cunter_brady +1; % cunter of bpm under 60
    end
end
perc_sample_brady = (cunter_brady / length(BPM_sans_zero))*100;
%% Tachycardia 
% detection of bradycardia with the bpms on BPM_sans_zero which are above
% 100
cunter_tachy = 0 ;
for i=1:length(BPM_sans_zero)
    if BPM_sans_zero(i) > 100 
        cunter_tachy = cunter_tachy +1; %cunter of bpm abode 100
    end
end

perc_sample_tachy = (cunter_tachy / length(BPM_sans_zero))*100;
end