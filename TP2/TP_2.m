clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('original source configuration: two source regions','FontSize',18); axis off;
%%
%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%signal to noise ratio
SNR=1;

%generate noisy data according to given SNR
X=Xs+ 1/sqrt(SNR)*Noise;

%visualize data (for a reduced number of sensors whose indices are 
%specified by idx_electrodes)
plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
title('noisy EEG data','FontSize',18);

%% 3.2. Studying the MNE algorithme

%%
%%% 3.2.1 Test the MNE case where lambda = 1
lambda = 1;
%S_h = MNE(X,G,lambda);
S_h = MNE(X(:,id),G,lambda);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h);
title(['Esitmated sources by MNE (lambda = ',num2str(lambda),')'],'FontSize',18); axis off;

% Remarks:
%   The algorithm detected many sources, less accurate method, need to tune the
%   regularization parameter.

%%
%%% 3.2.2 Vary the lambda from 0.1 1000
lambda = [0.1, 10, 100, 500, 1000];

for i = 1:5

    S_h = MNE(X(:,id),G,lambda(i));

    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h);
    title(['Esitmated sources by MNE (lambda = ',num2str(lambda(i)),')'],'FontSize',18); axis off;
end

% Remarks:
%   During the simulations, we observed that when the value of the
%   regularization parameter increases, then the algorithm tends to
%   penalize more, and give more sparse results. That's why we observed
%   smaller number of estimated sources when increasing the lambda.

%%
%%% 3.2.3 Vary the SNR and the lambda
SNR = [0.1, 7, 10];
lambda = [0.1, 500, 1000];

for i = 1:3
    for j = 1:3
        X=Xs+1/sqrt(SNR(j))*Noise;      
        S_h = MNE(X(:,id),G,lambda(i));
        figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h);
        title(['Esitmated sources by MNE SNR = ',num2str(SNR(j)),'(lambda = ',num2str(lambda(i)),')'],'FontSize',18); axis off;
    end
end

% Remarks:
%   The we have a higner SNR, we need a smaller value of the lambda to
%   estimate more accurate number of sources.

%%
%%% 3.2.4 Criteria fixing the value of the lambda.
SNR = 1;
X=Xs+1/sqrt(SNR)*Noise;
N = 30;
lambda = linspace(0.001,50,N);
%%
%%%% a. L-curve criterion
error_MNE = zeros(1,N);
S_norm = zeros(1,N);

for i = 1:N
    
    S_h = MNE(X(:,id),G,lambda(i));
    error_MNE(i) = norm(X(:,id)-G*S_h);
    S_norm(i) = norm(S_h); 
end

figure;
plot(error_MNE,S_norm);
xlabel('||x-As||')
ylabel('||s||')
hold on
scatter(error_MNE,S_norm);
hold on
text(error_MNE,S_norm, string(lambda),'HorizontalAlignment', 'right',...
                'VerticalAlignment', 'bottom');
title("L-curve criterion");

%%
%%%% b. Discrepancy principle
noise_power = norm(1/sqrt(SNR)*Noise(:,id),'fro');
N = 20;
res = zeros(1,N);
lambda = linspace(0.001,1000,N);

for i = 1:N
    S_h = MNE(X(:,id),G,lambda(i));
    res(i) = 2*log(norm(X(:,id)-G*S_h,'fro')) - 2*log(noise_power);
end

figure;
plot(lambda,res);
xlabel('lambda')
%dim = [.2 .5 .3 .3];
ylabel('2*(log(||X-G*S_h||) - log(||n||))')
hold on
scatter(lambda,res);
hold on
text(lambda,res, string(lambda),'HorizontalAlignment', 'right',...
                'VerticalAlignment', 'bottom');
title("L-curve criterion");
%[~,index] = min(res);
%annotation('textbox',dim,'String','min-lambda = '+string(lambda(1,index)),'FitBoxToText','on');
%%
%%%% c. Generalized cross-validation
error_CGV = zeros(1,N);
N = 30;
res = zeros(1,N);
lambda = linspace(0.001,50,N);

[n,~] = size(G);
for i = 1:N
    S_h = MNE(X(:,id),G,lambda(i));
    error_CGV(i) = (norm(X(:,id)-G*S_h, 'fro')^2)/(trace(eye(n)-G*G.'*inv((G*G.' + lambda(i)*eye(n))))^2);
end

figure;
plot(lambda, error_CGV);
hold on
scatter(lambda,error_CGV);
hold on
text(lambda, error_CGV, string(lambda));
xlabel('\lambda');
ylabel('CGV');
title("CGV-curve criterion");

%Remarks:
% After implenting several tests, we realize the approriate range of values of the
% regularization parameter is from 3 to 7. We can choose one of the
% value in this range for the lambda. The CGV criterion is suitable for evaluating precisely
% the performance of the algorithm,while using the discrepancy principle is
% hard to evaluate performances in the context.

%% 3.4 Studying the SISSY algorithme

%%
%%% 3.4.1 Teste the algorithm with the first parameter configuration.

SNR = 10;
X=Xs+1/sqrt(SNR)*Noise;
lambda = 1;
alpha = 0.1;
T = variation_operator(mesh,'face');
S_h = SISSY(X(:,id),G,T,lambda,alpha);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h);
title('Esitmated sources by SISSY ','FontSize',18); axis off;

% Remarks:
%   The algorithm works not well with the chosen parameters.

%%
%%% 3.4.2 Tester the algorithm by varying lambda from 0.1 - 1000
lambda = [0.1, 10, 100, 500, 1000];

for i = 1:5
    S_h = SISSY(X(:,id),G,T,lambda(i),alpha);

    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h);
    title(['Esitmated sources by SISSY (lambda = ',num2str(lambda(i)),')'],'FontSize',18); axis off;
end

% Remarks:
%   The algorithm tends to determine well the initial sources when
%   increasing the value of lambda.

%%
%%% 3.4.3 Tester the algorithm with a chosen lambda, and different values
%%% of the alpha.
alpha = [0.1, 0.3, 0.5, 0.7, 1];
lambda = 1;

for i = 1:5
    S_h = SISSY(X(:,id),G,T,lambda,alpha(i));

    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S_h);
    title(['Esitmated sources by SISSY (alpha = ',num2str(alpha(i)),')'],'FontSize',18); axis off;
end

% Remarks:
%   With the chosen lamda, the algorithm give more sparely
%   solution as the alpha increases. It means the lambda impact to sparsity
%   of the solution

%%
%%% 3.4.4 By observation, we find that the two parameters lambda and alpha
%%% impacts strongly to sparsity of solution of the algorithm. The solution
%%% is more spares when we increase these values of the lambda and the
%%% alpha.


%%
%%% 3.4.5 Tester the algorithm with a chosen lambda, and different values
%%% of the alpha.
