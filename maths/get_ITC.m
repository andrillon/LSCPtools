%%%% Computing ITC
%%%% Author:    Thomas Andrillon (thomas.andrillon@ens.fr)
%%%% Created:   26-11-13

% % From Romain Grandchamp
% % 1) Tu fais ton temps-fréquence avec timefreq
% % Exemple
% % elec_number = 64;
% % 
% % for elec = 1:elec_number
% % 
% % [tf(elec,:,:,:), freqs, times]  = timefreq(squeeze(EEG.data(
% % elec,:,:)), 128,'ntimesout', 360, 'cycles', [3 0.5],'freqs',[2 50],'nfreqs',96);
% % 
% % end;
% % 
% % Fait pour toutes les électrodes (elec), fréquences de 2 à 50 Hz et avec un sampling rate de 128Hz.  'cycles' [3 0.5], ça veut dire que c'est une 'wavelet'  qui commence avec 3 cycles à 2 Hz et qui augmente de façon linéaire jusqu'à 37.5 cycles à 50 Hz et un nombre de fréquences selon l'intervalle de fréquences que tu veux analyser ('nfreqs' 96 était une bonne résolution pour 2 à 50Hz).
% % 
% % 'ntimesout' correspond à tes output times, cad le nombre de timepoints que tu veux avoir dans ton tf. En fait, ici on a choisi 360 timepoints pour une durée de 120 secondes (époques de 120 secondes), donc des steps de ~33 ms. Cette résolution temporelle dépend de ce que tu veux regarder. 33ms peut être un peu trop long...après, ça dépend aussi de la puissance de ton ordi - plus de timesteps tu choisis, plus ca va être long à calculer.
% % Dans ton cas, si tu prends des epochs de 2s échantillonnées à 1000Hz, tu peux prendre 200 timepoints en sortie je pense, ca te fera des steps à 10ms.
% % 
% % 
% % 2) Tu moyennes ton tf => tfmean
% % Le tfmean est juste une variable, cad j'ai d'abord juste moyenné sur les essais:
% % 
% % tfmean = mean(abs(tf),4);
% % 
% % 3) Tu normalise ton tfmean (si tu veux)
% % tfnorm = repmat(mean(tfmean,3),[1 1 size(tfmean,3)]);
% % tfmean = tfmean ./ tfnorm;
% % 
% % 4) Et maintenant tu plottes tes résultats (en moyennant sur les électrodes que tu trouves intéressantes)
% % myelecs = [30 63 62 29 64 31 57 59 60];


function [tfmean_pow, tfmean_itc, freq, temps]=get_ITC(data,freqband,timewin,elecs,sampling)

ntimesout=size(data,2)/timewin;
nelec=0;
for elec = elecs
    nelec=nelec+1;
    [tf_pow(nelec,:,:,:), freq, temps, tf_itc(nelec,:,:,:)]  = timefreq(squeeze(data(elec,:,:)), sampling,'ntimesout', ntimesout, 'cycles', [3 0.5],'freqs',freqband,'nfreqs',96);
end;

% average across trials
tfmean_pow = mean(abs(tf_pow),4);
tfmean_itc = mean(abs(tf_itc),4);

% normalize
tfnorm_pow = repmat(mean(tfmean_pow,3),[1 1 size(tfmean_pow,3)]);
tfmean_pow = tfmean_pow ./ tfnorm_pow;

tfnorm_itc = repmat(mean(tfmean_itc,3),[1 1 size(tfmean_itc,3)]);
tfmean_itc = tfmean_itc ./ tfnorm_itc;

end