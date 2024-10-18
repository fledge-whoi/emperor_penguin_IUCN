% Code to calculate stable population before calculating generation length

% Global population projection with all scenario up to 2100, without
% movement and density dependance

% LE2 CESM CMIP6 SSP370

% Out : Stable Matrix A of the population projected between 1900 and 1950
% to calculate table population

% Calcul des conditions initiales pour chaque colonie (SAD precedente * comptages bilgecan 2009) 

%% Initiate
clear; close all; clc;

% DIRECTORY to Codes_EP

ordi = "/Users/aliceeparvier/Desktop";

addpath(sprintf('%s/Codes_EP/Codes_CMR/', ordi))

%load seaiceLE.mat %LE2 CESM CMIP6 SSP370
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/LE_CESM2_seaice_std.mat', ordi));

nens=length(ENS);

% Get covariance matrix of demo parameter
% uncertainties on the estimated parameters

nsim=100; %50 %parameter uncertainty & env stochasticity
[bR,bS]=beta; % estimate CMR, return and survival
xA=0.058; % tag loss
[betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);

%% Calculate mean SIC for the 50 ens for each season
% ncol=66;
% SIC=zeros(ncol,4);
% 
% % Choose period 
% tstart=1; % 0 =1900
% tend=51; %50 =1950
% nt=tend-tstart+1;
% 
% for c=1:ncol
%     for si=1:4 %season
%         SICmoy=zeros(nens*nt,1);
%         x=1;
%         for ens=1:nens
%             for t=1:50
%                 SICmoy(x,1)=ENS(ens).SICa(t, si, c);
%                 x=x+1;
%             end
%         end
%         SIC(c,si)=median(SICmoy(:,1));
%     end
% end

%% Take observations of SIC 1979-2018 for each colony

ncol = 66;
nt = 40;
SIC_tot = zeros(ncol, 4, nt);
for col=1:ncol
    SIC_col = readmatrix(sprintf("%s/Codes_EP/Codes_CMR/SIC_obs_allcol/SIC_col%d.csv", ordi, col));
    for season=1:4
        SIC_tot(col, season, :)= SIC_col(:,season);
        SIC_mean = mean(SIC_tot(col,season,1:30));
        SIC_tot(col,season,:) = (SIC_tot(col,season,:) - SIC_mean) ./ SIC_mean;
    end
end

% Calculate anomalies
% for s=1:4
%     SIC_mean = mean(mean(SIC_tot(:,s,:), 3));
%     SIC_tot(:,s,:) = (SIC_tot(:,s,:) - SIC_mean) ./ SIC_mean;
% end

SIC=mean(SIC_tot(:,:,:), 3); %mean on the 40 years


%% LOOP TO CALCULATE STABLE A, SAD

nsim=100

tconv=300;
ltot=zeros(tconv, nsim, ncol);
Atot(nsim, ncol) = {zeros(5)}; 
Vtot(nsim, ncol) = {zeros(5,1)}; 
Wtot(nsim, ncol) = {zeros(1,5)}; 
uftot=zeros(nsim, ncol);
umtot=zeros(nsim, ncol);
thetatot=zeros(9,nsim,ncol);

x=1;

for c=1:ncol
    display(c)

    %theta simulations
    thetasim=0;
    kmaxi=1000; 
    for k=1:kmaxi
        [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
        thetasim=thetasim+parameter_RAND3(xA,SIC(c,:),bR,bS,SICf,SICm,W,1);
    end
    theta=thetasim/kmaxi;

    for s=1:nsim %initial conditions variability
        wvec=simpsamp(5)';
        r=zeros(tconv,1);

        for t=1:tconv
            [Fav Mav]=breederav2(theta,1);

            %calculate demographic parameters'
            uf=min(Fav(:)'*wvec, Mav(:)'*wvec)/(Fav(:)'*wvec);
            um=min(Fav(:)'*wvec, Mav(:)'*wvec)/(Mav(:)'*wvec);

            % project population
            A=popmat(theta,uf,um); % define population matrix
            ltot(t, s, c)=max(eig(A));
            wvec=A*wvec;
            r(t,1)=sum(wvec);  % calculate growth rate
            wvec=wvec/r(t,1);
            if t==300
                uftot(s,c)=uf;
                umtot(s,c)=um;
                thetatot(:,s,c)=theta;
            end
        end % time
        Atot{s,c}=A; %save A stable
        [V,D,W2]=eig(A);
        ivp=find(eig(A)==max(eig(A))); %indice de la vp dominante
        Vtot{s,c}=V(:,ivp); %vecteur propre a droite
        Wtot{s,c}=W2(:,ivp);

        x=x+1;
    end %sim
end %col

% New calculations with SIC observed for each colony
% save("Astable_SICobs.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot") 
% Here simulations and colonies are separated

% Visualize lambda convergence
 figure
 for col=1:66
    plot(mean(ltot(1:150,:,col),2))
    hold on
 end

