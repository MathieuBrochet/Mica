
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

R_average = sumR / (length(R_value)-1); % probleme d'unit� 

%calcul du bpm 

total_ech_nb = length(data);

% on a f_ech =200 ech/s

time_of_ech_total = total_ech_nb/Fs; 

nb_bat_total = length(R_value);
nb_bat_par_sec = length(R_value)/ time_of_ech_total; 

bpm = nb_bat_par_sec *60; % a peu pr�s correct 


%% Bradycardia 
%under 60 bpm // bonus : faire un systeme qui note l'importance de la
%tachycharide en fonction du bpm qui se r�duit 

% if bpm < R_average 
%     % interface graphique qui le classifie en Bradycardia 
% end


%% Tachycardia 
%above 100 bpm // bonus faire pareil 


% if bpm > 100 
%     % interface graphique qui le classifie en Tachycardia
% end


%% Ectopic beat 
% delta(n) - delta(n-1) , cunter_e_beat compte le nombre de ectopic beat

k=1;
for j=1: length(R_R)-1
    Ect1(j) = abs(R_R(j+1)-R_R(j));
    if Ect1(j) >= R_average
        Ectopic_beat(k) = Ect1(j);
        k=k+1;
        
    end
end
cunter_e_beat= length(Ectopic_beat);
percent_of_extopic_beat = cunter_e_beat *100 / length(R_R);


%% Fibrillation 
 % on aura le 1 echantillon de covariance , il faut maintenant boucler sur �a 
N=length(R_R);
auto_cov_int_k1=0;
auto_cov_int_k2=0;
k=1;
for k=1:N-1
    R_k=0;
    for n=1:N-k-1
        auto_cov_int_k1=(R_R(n+k)-R_average);
        auto_cov_int_k2=(R_R(n) - R_average);
        auto_cov_int = auto_cov_int_k1*auto_cov_int_k2;
        R_k = R_k + auto_cov_int; % auto_covariance interm�diaire  
    end
     gamma(k)=R_k/(N-k-1);
 end
 %autocovariance d'un bruit blanc = dirac , grosse fluctuation � un moment
 % irr�gularit� 
 % ici normale 1 ne poss�de pas une auto covariance qui ressemble � un
 % dirac 


% importer les P_values, si
% length(P_value) ont des toutes petits values ( seuil � d�terminer) alors
% on peut dire qu'il y a abscence de P_value et donc par suite si 



%% Ventricular fibrillation 
%
