function [theta]=parameter_RAND_IPM(SST, wind, vwind, NOW, param_post, tlimit)

% theta
%   1     2     3      4      5      6    
% [fec, phiA, phi1, Pbreed, recr5, recr6]

% Parameters distributions
beta0_phiA=param_post(:,1);
betaSST=param_post(:,2);
betawind=param_post(:,3);
sigma2_phiA=param_post(:,4); % used to draw yearly random effect
beta0_breed=param_post(:,5);
betavwind=param_post(:,6);
sigma2_breed=param_post(:,7); 
beta0_fec=param_post(:,8);
betaNOW=param_post(:,9);
sigma2_fec=param_post(:,10);
mu_phi1=param_post(:,11);


% Parameter estimation
nsim_par=length(beta0_phiA);

theta=zeros(6,tlimit);

theta(5,:)= 0.22; %recr5
theta(6,:)= 0.32; %recr6


for t=1:tlimit

    index=randi(nsim_par);

    %Draw yearly random effect from normal distribution
    eps_phiA = normrnd(0, sqrt(sigma2_phiA(index))) ;
    eps_breed = normrnd(0, sqrt(sigma2_breed(index))) ;
    eps_fec = normrnd(0, sqrt(sigma2_fec(index))) ;


    theta(3,t)= mu_phi1(index); %juvenile survival

    % stochasticity due to interannual variability
    theta(2,t)= invlogit(beta0_phiA(index) + (betaSST(index) * SST(1,t)) + (betawind(index) * wind(1,t)) ...
        + eps_phiA); %adult survival

    theta(4,t)=invlogit(beta0_breed(index) + (betavwind(index) * vwind(1,t)) ...
        + eps_breed); % breeding proba

    theta(1,t)=invlogit(beta0_fec(index) + (betaNOW(index) * NOW(1,t)) ...
        + eps_fec); % fecondity

end %t