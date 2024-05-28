% Global population projection with all scenario up to 2100, without
% movement and density dependance

% LE2 CESM CMIP6 SSP370

% Version ou on garde les 50 ens separes et on fait mediane des simu pour
% chaque ens (pour carrying capacity)

% En sortie de ce code : Matrice A stable projetee entre 1900 et 1950 qui
% sert a calculer le temps de generation

% Calcul des conditions initiales pour chaque colonie (SAD precedente * comptages bilgecan 2009) 

%% Initiate
clear; close all; clc;
addpath('/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/')

%load seaiceLE.mat %LE2 CESM CMIP6 SSP370
load('/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/mat_data/LE_CESM2_seaice_std.mat');

nens=length(ENS);

% Get covariance matrix of demo parameter
% uncertainties on the estimated parameters

nsim=100; %50 %parameter uncertainty & env stochasticity
[bR,bS]=beta; % estimate CMR, return and survival
xA=0.058; % tag loss
[betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);

%% Calculate mean SIC for the 50 ens for each season
ncol=66;
SIC=zeros(ncol,4);

% Choose period 
tstart=1; % 0 =1900
tend=51; %50 =1950
nt=tend-tstart+1;

for c=1:ncol
    for si=1:4 %season
        SICmoy=zeros(nens*nt,1);
        x=1;
        for ens=1:nens
            for t=1:50
                SICmoy(x,1)=ENS(ens).SICa(t, si, c);
                x=x+1;
            end
        end
        SIC(c,si)=median(SICmoy(:,1));
    end
end

%% Take observation of SIC 1979-2018

ncol = 66;
nt = 40;
SIC_tot = zeros(ncol, 4, nt);
for col=1:66
    SIC_col = readmatrix(sprintf("/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/SIC_obs_allcol/SIC_col%d.csv", col));
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


%% BOUCLE POUR CALCULER A STABLE, SAD

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

    % Modify to obtain lambda closer to 1
%     theta(8)=0.87;
%     theta(9)=0.87;

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

% save("Astable1.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot") 
% 1 1900-1950
% 2 1950-2000
% 3 2000-2050 present
% 4 2050-2100

% New calculations with SIC std
% save("Astable_SICstd.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot") 
% Here simulations and colonies are separated

% New calculations with SIC observed for each colony
% save("Astable_SICobs.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot") 
% Here simulations and colonies are separated

% Visualize lambda convergence
 figure
 for col=1:66
    plot(mean(ltot(1:150,:,col),2))
    hold on
 end



%% Calcul des conditions intitiales
%sad=[0.2806    0.1642    0.2713    0.1030    0.1809]; %Stef

% Get covariance matrix of demo parameter
% uncertainties on the estimated parameters

nsim=100; %50 %parameter uncertainty & env stochasticity
[bR,bS]=beta; % estimate CMR, return and survival
xA=0.058; % tag loss
[betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);

load("Astable_SICstd.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot");

% Calcul de la SAD
sad=zeros(1,5);

% Convert Vtot into vector sum=1 and >0

for i=1:66 %col
    for s=1:100 %sim
        if Vtot{s,i}<0
            Vtot{s,i}=-Vtot{s,i};
        end
        Vtot{s,i} = Vtot{s,i} ./ sum(Vtot{s,i});
    end
end

for j=1:5 %class
    sadinter=0;
    for i=1:66 %col
        for s=1:100 %sim
            sadinter = sadinter + Vtot{s,i}(j);
        end
    end
    sad(1,j)=sadinter/(66*100);
end

dataOBSSAT_bilgecan_col=readmatrix("empe_sitesNewNB.csv");
nc=length(dataOBSSAT_bilgecan_col(:,2));
dataOBSSAT_bilgecan_col(65:66,2)=1000; %conservatif

N0col = zeros(5,nc);
for c=1:nc
    N0col(:,c) = sad.*dataOBSSAT_bilgecan_col(c,2);
end

BE=N0col(5,:); %only breeding pairs

% sad [0.2837    0.1578    0.2756    0.0982    0.1848]


%% Main loop for separate colonies

load proportion_ext_event.mat
proportion = proportion_ext_event(6); %LE_CESM2


nt= 51

ncol=66;
tic
SCEN2=struct; % scenario d'extreme event (# 5)
nsimu=100;
y2009=length(1900:2009);
y2018=length(1900:2018);

% perturbation - extreme events
percentSI=13/355; % =3.66% of extreme events counted on historical period
PERTclimat=binornd(1,percentSI,nsimu,nt);

[bR,bS]=beta; % estimate CMR, return and survival
xA=0.058; % tag loss
[betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);


for scenEXT=[1 2 4] % scenario d'extreme event
    Ncol=zeros(5,nt,ncol,nens,nsimu);
    Rcol=zeros(nsimu, nens, nt, ncol); 

    display(scenEXT)

    for c=1:ncol
        display(c)

        for ens=1:nens; % climate ensemble

            % SEA ICE FORECAST ENS(climatic ens 50) matrix SICa
            % time*season*colony

            % Extreme sea ice loss threshold
            threshold_SIC2=prctile(ENS(ens).SICa(y2009:y2018,2,c),percentSI); %laying season
            threshold_SIC4=prctile(ENS(ens).SICa(y2009:y2018,4,c),percentSI); %rearing season

            EXT=find(ENS(ens).SICa(:,2,c)<threshold_SIC2 | ENS(ens).SICa(:,4,c)<threshold_SIC4);
            iEXT=EXT(binornd(1,proportion,length(EXT),1)==1); % proportion of extreme events leading to total breeding failure

            % simulations

            for s=1:nsimu % sampling parameter uncertainty

                % parameter uncertainty
                [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
                
                %for extreme=1:nsimE % sampling environemental stochascity: year to year variability in theta not explained by sea ice.  

                    % POP INITIAL CONDITIONS
                    Ncol(:,1,c,ens,s)=N0col(:,c); 
                    wvec=sad';
                    r=zeros(nt,1);
                    r(1,1)=1;

                    % vital rates sampling stochasticity
                    [THETA]=parameter_RAND3(xA,ENS(ens).SICa(1:nt,:,c),bR,bS,SICf,SICm,W,nt);

                    % Extreme demographic scenario
                    % scenEXT=1 pas dextreme event
                    if scenEXT==2  % historical extreme event -> BS
                        THETA(2,PERTclimat(s,:)==1)=0; %si extreme event, BS=0
                    elseif scenEXT==3 % historical  extreme event -> BS & Sa
                        THETA(2,PERTclimat(s,:)==1)=0; %BS=0
                        THETA(8,PERTclimat(s,:)==1)=0.9*THETA(8,PERTclimat(x,:)==1); 
                        THETA(9,PERTclimat(s,:)==1)=0.9*THETA(9,PERTclimat(x,:)==1); 
                    elseif scenEXT==4 % sea ice loss extreme event-> BS
                        THETA(2,iEXT)=0;
                    elseif scenEXT==5 % sea ice loss extreme event-> BS & Sa
                        THETA(2,iEXT)=0;
                        THETA(8,iEXT)=THETA(8,iEXT).*0.9;
                        THETA(9,iEXT)=THETA(9,iEXT).*0.9;
                    end

                    % POP projection
                    THETA=THETA(:,1:nt);
                    [Fav Mav]=breederav2(THETA,nt); % vector to calculate males and females available to mate

                    for t=2:nt % 1900 to 1950
                        %calculate demographic parameters'
                        uf=min(Fav(:,t)'*Ncol(:,t-1,c,ens,s), Mav(:,t)'*Ncol(:,t-1,c,ens,s))/(Fav(:,t)'*Ncol(:,t-1,c,ens, s));
                        um=min(Fav(:,t)'*Ncol(:,t-1,c,ens, s), Mav(:,t)'*Ncol(:,t-1,c,ens,s))/(Mav(:,t)'*Ncol(:,t-1,c,ens, s));
                        % project population
                        wvec=popmat(THETA(:,t),uf,um)*wvec;
                        Ncol(:,t,c,ens,s)=popmat(THETA(:,t),uf,um)*Ncol(:,t-1,c,ens,s);
                        r(t)=sum(wvec);  % calculate growth rate
                        wvec=wvec/r(t,1);
                    end % time

                    Rcol(s,ens,:,c)=r;

                %end % Envt sto sim
            end % demographic sim
        end  % CLIMATE ENS
    end %col
    SCEN2(scenEXT).Ncol=Ncol;
    SCEN2(scenEXT).GRcol=Rcol;
end %scenEXT

toc

%save('NRBE_temp2.mat', 'BE', 'Ncol', 'Rcol')

%% Calculate mean on scenEXT
Ncol_mean = zeros(5,nt,ncol,nens,nsimu);
Rcol_mean = zeros(nsimu, nens, nt, ncol); 

for scenEXT=[1 2 4]
    Ncol_mean = Ncol_mean + SCEN2(scenEXT).Ncol;
    Rcol_mean = Rcol_mean + SCEN2(scenEXT).GRcol;
end
Ncol_mean = Ncol_mean ./3; %mean on scenEXT
Rcol_mean = Rcol_mean ./3;

%save ('NRBE_new.mat', 'BE', 'Ncol_mean', 'Rcol_mean')


%% Plot breeding pairs
figure
for c=1:66
    clf
    plot(squeeze(Ncol(5,:,c,1:100)))
    pause
end         
