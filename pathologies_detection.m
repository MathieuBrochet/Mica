function [bpm,perc_of_p_value_AF,perc_sample_brady,perc_sample_tachy,percent_of_extopic_beat] = pathologies_detection( data, Fs)



%% Load a signal


% [file,path] = uigetfile('*.mat', 'rt');
% signal = load(fullfile(path, file));
% data = signal.ecg; % Your ecg data
% Fs = signal.Fs; % Sampling frequency
% N = size(data,2); % Data length
% time_axis = (1:N)/Fs;
 [R_value, Q_value, S_value, P_value, T_value] = R_wave_detection( data, Fs); % récupération des valeurs PQRST précédentes
% 
% 

%% Rythm cardiac


delay = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0 1]; 
data = conv(delay,data);
% obligé de remettre le delais car 
% on load de nouveau le signal data avec un delais par rapport à nos signaux de filtre 

sumR = 0;
for i=1:length(R_value)-1
   R_R(i) = R_value(i+1)-R_value(i);
   sumR = sumR + R_R(i);
end
R_average = sumR / (length(R_value)-1); % probleme d'unité 

%calcul du bpm 
total_ech_nb = length(data);

%on a f_ech =200 ech/s
time_of_ech_total = total_ech_nb/Fs; 
nb_bat_total = length(R_value);
nb_bat_par_sec = length(R_value)/ time_of_ech_total; 
bpm = nb_bat_par_sec *60; % correct 
% 



%% Ectopic beat 
% delta(n) - delta(n-1) , cunter_e_beat compte le nombre de ectopic beat
Ectopic_beat=zeros();
k=1;
for j=1: length(R_R)-1
    Ect1(j) = abs(R_R(j+1)-R_R(j));
    if Ect1(j) >= R_average
        Ectopic_beat(k) = Ect1(j);
        k=k+1;
        
    end
end
cunter_e_beat= length(Ectopic_beat);  % le nombre de ectopic beat 
percent_of_extopic_beat = cunter_e_beat *100 / length(R_R); % le pourcentage d'ectopic beat 


%%  Atrial Fibrillation 
N=length(R_R);
auto_cov_fibr_k1=0;
auto_cov_fibr_k2=0;
for k=1:N-1
    R_k=0;
    for n=1:N-k-1
        auto_cov_fibr_k1=(R_R(n+k)-R_average);
        auto_cov_fibr_k2=(R_R(n) - R_average);
        auto_cov_fibr = auto_cov_fibr_k1*auto_cov_fibr_k2;
        R_k = R_k + auto_cov_fibr; % auto_covariance intermédiaire  
    end
     gamma(k)=R_k/(N-k-1);
 end
%  autocovariance d'un bruit blanc = dirac , grosse fluctuation à un moment
%  irrégularité 
%  ici normale 1 ne possède pas une auto covariance qui ressemble à un
%  dirac , car les signaux ecarté d'une periode se ressemble --> ils
%  possède un lien important entre eux = autocorélation qui ne fluctue pas
%  trop 
%  par ailleur si un signal à une autocovariance de type dirac, alors le
%  signal a pas de lien entre les signaux qui se suivent et donc, le signal
%  fluctue beaucoup entre chaque batement, pas de lien == fibrillation
%  =malade 
% % importer les P_values, si
% % length(P_value) ont des toutes petits values ( seuil à déterminer) alors
% % on peut dire qu'il y a abscence de P_value et donc par suite si 

 
%% detection du % de pic P valable
P_value_real = P_value( P_value >0);
R_value_real = R_value ( R_value >0);

P_ord = data(P_value_real);
R_ord = data (R_value_real);
P_available = zeros();


for i=1:length(P_value_real)
   if (abs(R_ord(i)/P_ord(i))) > 5 && ((R_ord(i)/P_ord(i))) < 400
       P_available(i) = abs(P_ord(i));
       frac(i)=abs(R_ord(i)/P_ord(i));
   else
       P_available(i) = 0 ;
       frac(i)=abs(R_ord(i)/P_ord(i));
   end
end

P_available_ss_zero = P_available(P_available >0);

nb_p_values = length(P_value_real);
nb_p_values_real = length(P_available_ss_zero);
perc_of_p_value_AF = (nb_p_values_real/nb_p_values)*100;

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


%% releving interessant part of the signal 
pas_window=200;
partie_entiere = floor(length(data)/2);
BPM =zeros(1,partie_entiere); % vecteur qui contiendra les bpm pour la fenetre quand elle se décale 
n=0;
window = zeros(1,length(data));
for m=1:pas_window:(length(data)/2) % je calcul la bpm toute es 58 fenetre , sinon trop long et peu de de variation si on 
    % décale la fenetre d'un echantilllon 
% % window 
% j'ai choisi 58 car diviseur de 149 119 pile poil ( jusqu'ou va m)
    n = n+pas_window;
    for k=n:((length(data)/2) + n-1) % je cree ma fenetre en l'augmentant de 1 vers la droite a chaue itt de m
        window(k) = 1;
    end
    for k=1:n
        window(k) = 0; 
    end % je met le derniere element de ma fenetre à 0 (reduction
    % à gauche ) comme si elle se déplacer à droite 
    data_value_window = zeros(length(data),1);


    for p=1:length(data)
        data_value_window(p)= window(p).*data(p); % on recupere seulement les r values de la fenetre 
    end 

% %% trouver le R pic de data_value_window 

    new_vect_window_R = zeros();

% plot(data);
% hold all; 
% plot (window);
% plot(data_value_window);

    for k=1:length(R_value)
        new_vect_window_R(k)=data_value_window(R_value(k)); %contient les pics R sur la fenetre 
    end
  
    R_window = new_vect_window_R(new_vect_window_R > 0);
% % %calcul du bpm window 

    total_ech_nb_window = length(data)/2;

% % % on a f_ech =200 ech/s
 
    time_of_ech_total_window = total_ech_nb_window/Fs; 


    nb_bat_total_window = length(R_window);
    nb_bat_par_sec_window = length(R_window)/ time_of_ech_total_window; 
    bpm_window = nb_bat_par_sec_window *60; % correct 
    BPM(m) = bpm_window;
end 

BPM_sans_zero = BPM( BPM >0); % valeur des bpm du signal en entrée pour une porte décalé de 58 à chaque fois 

BPM_average = 0; %initialisation du vecteur qui contient les bpm de chaque fenetre 

for i = 1:length(BPM_sans_zero)
    BPM_average = BPM_average + BPM_sans_zero(i) ; 
end

BPM_average = BPM_average / length(BPM_sans_zero); % on calcul la moyenne des bpm du signal fenetré 
% on devrait retrouver la meme moyenne que bpm 

%% Bradycardia 
% under 60 bpm // bonus : faire un systeme qui note l'importance de la
% tachycharide en fonction du bpm qui se réduit 

cunter_brady = 0;
for i=1:length(BPM_sans_zero)
    if bpm < 60
        cunter_brady= cunter_brady +1;
    end
end
perc_sample_brady = (cunter_brady / length(BPM_sans_zero))*100;
%% Tachycardia 
% above 100 bpm // bonus faire pareil 
cunter_tachy = 0 ;
for i=1:length(BPM_sans_zero)
    if BPM_sans_zero(i) > 100 
        cunter_tachy = cunter_tachy +1;
    end
end

perc_sample_tachy = (cunter_tachy / length(BPM_sans_zero))*100;
end