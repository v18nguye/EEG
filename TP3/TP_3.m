clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%%
%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%%
%%% 2.1 Comparaison qualitative
% 1. RSB = [0.1, 1, 10]
SNR = 10;
alpha_sissy = 10;
lamb_sissy = 10;
lamb_mne = 100;

%generate noisy data according to given SNR
X=Xs+ 1/sqrt(SNR)*Noise;

% MNE
S_h_mne = MNE(X(:,id),G,lamb_mne);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h_mne);
%title(['Esitmated sources by MNE'],'FontSize',18); axis off;
title(['Esitmated sources by MNE (lambda = ',num2str(lamb_mne),')'],'FontSize',18); axis off;

% SISSY
T = variation_operator(mesh,'face');
S_h_sissy = SISSY(X(:,id),G,T,lamb_sissy,alpha_sissy);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h_sissy);
%title(['Esitmated sources by SISSY'],'FontSize',18); axis off;
title(['Esitmated sources by SISSY (lambda = ',num2str(lamb_sissy),')'],'FontSize',18); axis off;

%%
%%% 2.1 Comparaison qualitative
% 2. Bruit spatialement corrélé, SNR = [0.1 1 10]

%# Génerer un bruit spatialement corrélé
Snoise = zeros(19626,200);
mask =  S == 0;
Snoise(mask) = rand(size(S(mask),1),1);
Noise = G*Snoise;

alpha_sissy = 0.1;
lamb_sissy = 10;
lamb_mne = 3.5;

SNR = 10;

%generate noisy data according to given SNR
X=Xs+ 1/sqrt(SNR)*Noise;

% MNE
S_h_mne = MNE(X(:,id),G,lamb_mne);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h_mne);
title(['Esitmated sources by MNE (lambda = ',num2str(lamb_mne),')'],'FontSize',18); axis off;

% SISSY
T = variation_operator(mesh,'face');
S_h_sissy = SISSY(X(:,id),G,T,lamb_sissy,alpha_sissy);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h_sissy);
title(['Esitmated sources by SISSY (lambda = ',num2str(lamb_sissy),')'],'FontSize',18); axis off;
%%
%%% 2.2 Comparer quantitative
% 3. DLEs
SNR = 1;
T = variation_operator(mesh,'face');
alpha_sissy = 0.1;
lamb_sissy = 10;
lamb_mne = 3.5;
r_grid = mesh.f;

%# définir une seuil pour la source estimé
seuil = 0.1;
idx = sum(abs(S) > 0, 2) ~= 0;

%#
N = 10;
dle_mne = zeros(N,1);
dle_sissy = zeros(N,1);
%snr_matrice = ones(N,1)*SNR;

%# DLEs for 10 ensembles des données
for i = 1:N
    
    %# generer le bruit spatialement corrélé 
    Snoise = zeros(19626,200);
    mask =  S == 0;
    Snoise(mask) = rand(size(S(mask),1),1);
    Noise = G*Snoise;
    
    %# generate noisy data according to given SNR
    X=Xs+ 1/sqrt(SNR)*Noise;
    
    %# MNE
    S_mne = MNE(X,G,lamb_mne);
    %## seuiller
    S_mne(abs(S_mne)<seuil) = 0;
    idx_mne = sum(S_mne > 0, 2) ~= 0;
    dle_mne(i,1) = DLE(idx,idx_mne,r_grid);
    
    %# SISSY
    S_sissy = SISSY(X,G,T,lamb_sissy,alpha_sissy);
    S_sissy(abs(S_sissy)<seuil) = 0;
    idx_sissy = sum(S_sissy > 0, 2) ~= 0;
    dle_sissy(i,1) = DLE(idx,idx_sissy,r_grid);
    
end

figure;
boxplot(dle_mne);
xlabel('SNR')
ylabel('DLE')
title('Boxplot MNE')

figure;
boxplot(dle_sissy);
xlabel('SNR')
ylabel('DLE')
title('Boxplot SISSY')

%%
%%% 2.2 Comparer quantitative
% 4. DLEs sur les différentes valeurs du SNR
SNRs = [0.1,1,3,7,10];

T = variation_operator(mesh,'face');
alpha_sissy = 0.1;
lamb_sissy = 10;
lamb_mne = 3.5;
r_grid = mesh.f;

%# définir une seuil pour la source estimé
seuil = 0.1;
idx = sum(abs(S) > 0, 2) ~= 0;

%#
N = 10;
dle_mne = zeros(5,N);
dle_sissy = zeros(5,N);

for j = 1:5
    
    SNR = SNRs(j);
    for i = 1:10
        %# generer le bruit spatialement corrélé 
        Snoise = zeros(19626,200);
        mask =  S == 0;
        Snoise(mask) = rand(size(S(mask),1),1);
        Noise = G*Snoise;

        %# generate noisy data according to given SNR
        X=Xs+ 1/sqrt(SNR)*Noise;

        %# MNE
        S_mne = MNE(X,G,lamb_mne);
        %## seuiller
        S_mne(abs(S_mne)<seuil) = 0;
        idx_mne = sum(S_mne > 0, 2) ~= 0;
        dle_mne(j,i) = DLE(idx,idx_mne,r_grid);

        %# SISSY
        S_sissy = SISSY(X,G,T,lamb_sissy,alpha_sissy);
        S_sissy(abs(S_sissy)<seuil) = 0;
        idx_sissy = sum(S_sissy > 0, 2) ~= 0;
        dle_sissy(j,i) = DLE(idx,idx_sissy,r_grid);
        
    end
end

%#
dle_mne_avg = mean(dle_mne,2);
dle_sissy_avg = mean(dle_sissy,2);

figure;
plot(SNRs,dle_mne_avg)
xlabel('SNR')
ylabel('DLE_avgs')

title('DLE Average')
hold on
plot(SNRs,dle_sissy_avg)

legend('MNE','SISSY')
hold off