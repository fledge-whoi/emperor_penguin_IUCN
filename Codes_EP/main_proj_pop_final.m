% Intrinsic growth rate projection with all scenario up to 2100
% SCEN(5).GR=temporaire;
%% Initiate
clear; close all; clc;
addpath('/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/')

%load seaiceLE.mat
load('mat_data/LE_CESM2_seaice_std.mat')

%load the proportion of extreme events for scenEXT4 for each climate model
load proportion_ext_event.mat

nsim=10; %50 %parameter uncertainty & env stochasticity
nsimE=10; %50 % environmental varibility not related to SIC

load NewNameColony.mat

% colony
yearstart=1900;
timePOP=yearstart:2100;
nt=length(timePOP);
nc=66; % to be uptdated

% climate
ntime=1900:2100;
nens=50; % climate ensemble

nsimulation=nsim*nsimE; % number of simulation total sampled within nsim*nsimE*nens 

% perturbation - extreme events
percentSI=13/355; % =3.66% of extreme events counted on historical period
PERTclimat=binornd(1,percentSI,nsimulation,nt);
y2009=find(ntime==2009);y2018=find(ntime==2018);


% Get covariance matrix of demo parameter
% uncertainties on the estimated parameters
[bR,bS]=beta; % estimate CMR, return and survival
xA=0.058; % tag loss
[betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);


%% main loop

tic %start clo                                                                                                                                                                                                                                                                                                                            ck
SCEN=struct; % scenario d'extreme event (# 5)

for scenEXT=[1 2 4] % scenario d'extreme event
    Rtot=zeros(nens, nsimulation,nt,nc);

    display(scenEXT)
    
    for c=1:nc % colony
        display(c)

        for ens=1:nens % climate ensemble

            % SEA ICE FORECAST ENS(climatic ens 50) matrix SICa
            % time*season*colony

            % Extreme sea ice loss threshold
            threshold_SIC2=prctile(ENS(ens).SICa(y2009:y2018,2,c),percentSI);
            threshold_SIC4=prctile(ENS(ens).SICa(y2009:y2018,4,c),percentSI);

            EXT=find(ENS(ens).SICa(:,2,c)<threshold_SIC2 | ENS(ens).SICa(:,4,c)<threshold_SIC4);
            proportion = proportion_ext_event(6);
            iEXT=EXT(binornd(1,proportion,length(EXT),1)==1); % likely as not to observe an extreme events

            % simulations
            x=1;
            
            for s=1:nsim*nsimE % sampling parameter uncertainty

                % parameter uncertainty
                [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
           %breeding success, proportion of 1st return, Return, survival
                
                    % POP INITIAL CONDITIONS
                    r=zeros(nt,1);

                    % SAD
                    wvec=simpsamp(5)';

                    % vital rates sampling stochasticity
                    [THETA]=parameter_RAND3(xA,ENS(ens).SICa(:,:,c),bR,bS,SICf,SICm,W,nt);
                    %[THETA]=parameter_DET3(xA,ENS(ens).SICa(:,:,c),bS,SICf,SICm,W,nt);

                    % Extreme demographic scenario
                    % scenEXT=1 pas dextreme event
                    if scenEXT==2  % historical extreme event -> BS
                        THETA(2,PERTclimat(x,:)==1)=0; %si extreme event, BS=0
                    elseif scenEXT==3 % historical  extreme event -> BS & Sa
                        THETA(2,PERTclimat(x,:)==1)=0; %BS=0
                        THETA(8,PERTclimat(x,:)==1)=0.9*THETA(8,PERTclimat(x,:)==1); 
                        THETA(9,PERTclimat(x,:)==1)=0.9*THETA(9,PERTclimat(x,:)==1); 
                    elseif scenEXT==4 % sea ice loss extreme event-> BS
                        THETA(2,iEXT)=0;
                    elseif scenEXT==5 % sea ice loss extreme event-> BS & Sa
                        THETA(2,iEXT)=0;
                        THETA(8,iEXT)=THETA(8,iEXT).*0.9;
                        THETA(9,iEXT)=THETA(9,iEXT).*0.9;
                    end

                    % POP projection

                    [Fav Mav]=breederav2(THETA,nt); % vector to calculate males and females available to mate

                    for t=1:nt % from year start (1920 or 2006) to 2100
                        %calculate demographic parameters'
                        uf=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Fav(:,t)'*wvec);
                        um=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Mav(:,t)'*wvec);
                        % project population
                        wvec=popmat(THETA(:,t),uf,um)*wvec; % define population matrix
                        r(t,1)=sum(wvec);  % calculate growth rate
                        wvec=wvec/r(t,1);
                    end % time

                    Rtot(ens, x,:,c)=r;
                    x=x+1;

                %end % Envt sto sim
            end % demographic sim
        end  % CLIMATE ENS
    end % colony
    SCEN(scenEXT).GR=Rtot;

%     % median of each colony
%     Rcol_med=zeros(3, nt, nc);
%     for c=1:nc
%         Rcol_med(:,:,c) = quantile(Rtot(:,:,c),[0.05,0.5,0.95]);
%     end
%     Rcol_med = permute(Rcol_med, [3,2,1]);
%     %fileName = sprintf('R_temp_scen_%d', scenEXT)
%     %save(fileName, 'Rtot', 'Rcol_med')

    %save('R_scen3','Rtot', 'Rcol_med')
end  % EXT scenario

%save('GR_CESM2_SIC_std', 'SCEN')

%% Fusion 5 scenario
Rtot=SCEN(1).GR;
for scenEXT=2:5
    Rtot=cat(1, Rtot, SCEN(scenEXT).GR);
end


%% FIGURES
load R10000-001.mat
%% figure of pop growth rate for one colony

Rtot=R1;
nsimulation=length(R1(:,1,1))
col=15
figure
plot(Rtot(randi(nsimulation,10,1),:,col)')
hold on
plot(median(Rtot(:,:,col)),'r', LineWidth=2)

%%  after diffusion
%load pop_R10000_diff.mat
col=15
nsimulation=length(NCOL(:,1,1))
figure
plot(NCOL(randi(nsimulation,10,1), :, col)')

%%

scenEXT=5;
col= 54 %LEDD=54 %CROZ=50 %KLOA=21; %39=PG

figure
plot(SCEN(scenEXT).GR(randi(nsimulation,10,1),:,col)')
hold on
plot(median(SCEN(scenEXT).GR(:,:,col)),'r', LineWidth=2)

%visualiser declin
nt=length(median(SCEN(scenEXT).GR(:,:,col)))
    delai=20
    i=20
    while median(median(SCEN(scenEXT).GR(:,i:i+delai,col))) > 1
        i=i+1
    end
    limite_t = i+delai

hold on
yline(1, '--blue' )
xline(limite_t, 'r', {limite_t+1900}, LineWidth=2)

%% FIGURE FOR EACH SCENARIO

figure
for scenEXT=1:6
    subplot(2,3,scenEXT)
    col= 54 %LEDD=54 %CROZ=50 %KLOA=21; %39=PG
    plot(SCEN(scenEXT).GR(randi(nsimulation,10,1),:,col)')
    hold on
    plot(median(SCEN(scenEXT).GR(:,:,col)),'r',LineWidth=2)
    title('Scenario', scenEXT)

    %visualiser declin
    nt=length(median(SCEN(scenEXT).GR(:,:,col)))
    delai=30
    i=20
    while median(median(SCEN(scenEXT).GR(:,i:i+delai,col))) > 1
        i=i+1
    end
    limite_t = i + (delai/2)
    hold on
    yline(1, '--black' )
    xline(limite_t, 'r', {limite_t+1900}, LineWidth=2)
end


%% FIGURE FOR FOUR COLONY
load('/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/N_results/N_CanESM2_5scen.mat')
figure
scenEXT=1
i=1
for col=[54, 50, 21, 39]
    subplot(2,2,i)
    plot(SCEN_N(scenEXT).NCOL(randi(nsimulation,10,1),:,col)')
    hold on
    plot(median(SCEN_N(scenEXT).NCOL(:,:,col)),'r',LineWidth=2)
    title('Colony', col)
    i=i+1

%     %visualiser declin
%     nt=length(median(SCEN_N(scenEXT).GR(:,:,col)))
%     delai=20
%     j=20
%     while median(median(SCEN_N(scenEXT).GR(:,j:j+delai,col))) > 1
%         j=j+1
%     end
%     limite_t = j+delai
    
%     hold on
%     yline(1, '--black' )
%     xline(limite_t, 'r', {limite_t+1900}, LineWidth=2)
end

%% with IC and median

load('R10000-001.mat','Rmedian_comb','R1')
col= 21; %LEDD=54 %CROZ=50 %KLOA=21; %39=PG
couleur=cbrewer('qual', 'Paired', 12);
fontsi=16;

figure
hold on

ttt = [timePOP, fliplr(timePOP)];

inBetween = [Rmedian_comb(col,:,1), fliplr(Rmedian_comb(col,:,3))];
f = fill(ttt, inBetween, couleur(1,:));
set(f,'EdgeColor','none','FaceAlpha', 0.3)
plot(timePOP,Rmedian_comb(col,:,2),'color',couleur(4,:),'linewidth',5) %median
plot(timePOP,Rmedian_comb(col,:,1),'color',couleur(3,:),'linewidth',2)
plot(timePOP,Rmedian_comb(col,:,3),'color',couleur(3,:),'linewidth',2)


