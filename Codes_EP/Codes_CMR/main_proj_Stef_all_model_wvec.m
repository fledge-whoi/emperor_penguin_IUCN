% Intrinsic growth rate projection with all scenario up to 2100

clear; close all; clc;

% DIRECTORY to Codes_EP
ordi = '/Users/aliceeparvier/Desktop';


%% Initiate

%load informations on climate models (years, number of ens)
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi));

%load the proportion of extreme events for scenEXT4 for each climate model
load proportion_ext_event.mat

model= ["LE_CanESM2","LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"] 

for i=[1:3 6 7] %Climate model
    mod=model(i); %name
    display(mod)
    
    % Standardized SIC data (Bilgecan method), 500 and 850km2
    file_name=sprintf(sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis.mat', ordi, mod));

    if i==6
        file_name=sprintf(sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod));
    end

    load(file_name)

    load NewNameColony.mat

    % colony
    yearstart=info_models(i,3);
    yearend=info_models(i,4);
    if i==7
        yearstart=2006;
    end

    timePOP=yearstart:yearend;
    nt=length(timePOP);
    nc=66;

    % climate
    ntime=yearstart:yearend;
    nens=info_models(i,5); % climate ensemble, nens minimum is 5
    nsimulation=100;

    % perturbation - extreme events
    percentSI=13/355; % =3.66% of extreme events counted on historical period
    PERTclimat=binornd(1,percentSI,nsimulation,nt);
    y2009=find(ntime==2009);y2018=find(ntime==2018);


    % Get covariance matrix of demo parameter
    % uncertainties on the estimated parameters
    [bR,bS]=beta; % estimate CMR, return and survival
    xA=0.058; % tag loss
    [betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);


    % MAIN LOOP

    tic %start clock                                                                                                                                                                                                                                                                                                                            ck
    SCEN=struct; 

    for scenEXT=[1 2 4]  %extreme event scenario (# 5)
        display(scenEXT)

        for col=1:nc % colony
            display(col)

            parfor ens=1:nens % climate ensemble
                
                % SEA ICE FORECAST ENS (climatic ens 50) matrix SICa
                % time*season*colony

                % Extreme sea ice loss threshold
                threshold_SIC2=prctile(ENS(ens).SICa(y2009:y2018,2,col),percentSI);
                threshold_SIC4=prctile(ENS(ens).SICa(y2009:y2018,4,col),percentSI);

                EXT=find(ENS(ens).SICa(:,2,col)<threshold_SIC2 | ENS(ens).SICa(:,4,col)<threshold_SIC4);
                proportion = proportion_ext_event(i); %Calculated to have same nb of colonies having total breeding failure 2018-2022 than observed
                iEXT=EXT(binornd(1,proportion,length(EXT),1)==1); % proportion of extreme events leading to total breeding failure

                % SIMULATIONS

                Rtot=zeros(nt, nsimulation);
                WVEC=zeros(5, nt, nsimulation); %to save SAD
                x=1;

                for s=1:nsimulation % sampling parameter uncertainty

                    % parameter uncertainty
                    [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
                    %breeding success, proportion of 1st return, Return, survival

                    % POP INITIAL CONDITIONS
                    r=zeros(nt,1);
                    w=zeros(5, nt);

                    % SAD
                    wvec=simpsamp(5)';
                   

                    % vital rates sampling stochasticity
                    [THETA]=parameter_RAND3(xA,ENS(ens).SICa(:,:,col),bR,bS,SICf,SICm,W,nt);

                    % Extreme demographic scenario
                    % scenEXT=1 no extreme event
                    if scenEXT==2  % historical extreme event -> BS
                        THETA(2,PERTclimat(x,:)==1)=0;
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

                    % Pop projection

                    [Fav Mav]=breederav2(THETA,nt); % calculate males and females available to mate

                    for t=1:nt % from 1900 to 2100
                        %calculate demographic parameters
                        uf=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Fav(:,t)'*wvec);
                        um=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Mav(:,t)'*wvec);

                        % project population
                        wvec=popmat(THETA(:,t),uf,um)*wvec; % define population matrix
                        r(t,1)=sum(wvec);  % calculate growth rate
                        wvec=wvec/r(t,1);
                        w(:,t,1)=wvec;
                    end % time

                    Rtot(:, s)=r; %R for each simulation for each t for one ens
                    WVEC(:,:,s)=w;
                    x=x+1;

                end % demographic sim

                MAT=WVEC(2,:,:)+WVEC(4,:,:)+WVEC(5,:,:); %proportion of matures individuals
                MAT=permute(MAT, [2 3 1]);
                % Save results for each mod, scenEXT, col, sim
                filename=sprintf('%s/Codes_EP/GR_results_tot/GR_results_Stef/%s/GR_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);
                filename2=sprintf('%s/Codes_EP/MAT_results_tot/MAT_results_Stef/%s/MAT_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);
                % CESM1 starting in 1900
                %filename=sprintf('%s/Codes_EP/GR_results_tot/GR_results_Stef/%s/GR_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);

                writematrix(Rtot, filename);
                writematrix(MAT, filename2);

            end  % climate ensemble
        end % colony
    end  % extreme event scenario
end %model





%% MEAN OF ENSEMBLES

%load informations on climate models (years, number of ens)
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi));

%load the proportion of extreme events for scenEXT4 for each climate model
load proportion_ext_event.mat

model= ["LE_CanESM2","LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"] 

for i=[1.3 6 7] %Climate model
    mod=model(i); %name
    display(mod)
    
    % Standardized SIC data (Bilgecan method), 500 and 850km2
    file_name=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis.mat', ordi, mod);

    if i==6 %CESM2
        file_name=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod);
    end
    load(file_name)

    load NewNameColony.mat

    % colony
    yearstart=info_models(i,3);
    yearend=info_models(i,4);
    timePOP=yearstart:yearend;
    nt=length(timePOP);
    nc=66;

    % climate
    ntime=yearstart:yearend;
    nens=info_models(i,5); % climate ensemble, nens minimum is 5
    nsim=100;
    nsimulation=100;

    % perturbation - extreme events
    percentSI=13/355; % =3.66% of extreme events counted on historical period
    PERTclimat=binornd(1,percentSI,nsimulation,nt);
    y2009=find(ntime==2009);y2018=find(ntime==2018);


    % Get covariance matrix of demo parameter
    % uncertainties on the estimated parameters
    [bR,bS]=beta; % estimate CMR, return and survival
    xA=0.058; % tag loss
    [betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);


    % Calculate mean of ensemble
    
    SICa_mean = zeros(nt, 4, nc);
    for ens = 1:nens
        SICa_mean = SICa_mean + ENS(ens).SICa(:,:,:);
    end
    SICa_mean = SICa_mean ./ nens;


    % MAIN LOOP

    tic %start clock                                                                                                                                                                                                                                                                                                                            ck
    SCEN=struct; 

    for scenEXT=[1 2 4]  %extreme event scenario (# 5)
        display(scenEXT)

        for col=1:nc % colony
            display(col)

                
                % SEA ICE FORECAST ENS (climatic ens 50) matrix SICa
                % time*season*colony

                % Extreme sea ice loss threshold
                threshold_SIC2=prctile(SICa_mean(y2009:y2018,2,col),percentSI);
                threshold_SIC4=prctile(SICa_mean(y2009:y2018,4,col),percentSI);

                EXT=find(SICa_mean(:,2,col)<threshold_SIC2 | SICa_mean(:,4,col)<threshold_SIC4);
                proportion = proportion_ext_event(i); %Calculated to have same nb of colonies having total breeding failure 2018-2022 than observed
                iEXT=EXT(binornd(1,proportion,length(EXT),1)==1); % proportion of extreme events leading to total breeding failure

                % SIMULATIONS

                Rtot=zeros(nt, nsimulation);
                WVEC=zeros(5, nt, nsimulation); %to save SAD
                x=1;

                for s=1:nsimulation % sampling parameter uncertainty

                    % parameter uncertainty
                    [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
                    %breeding success, proportion of 1st return, Return, survival

                    % POP INITIAL CONDITIONS
                    r=zeros(nt,1);
                    w=zeros(5, nt);

                    % SAD
                    wvec=simpsamp(5)';

                    % vital rates sampling stochasticity
                    [THETA]=parameter_RAND3(xA,SICa_mean(:,:,col),bR,bS,SICf,SICm,W,nt);

                    % Extreme demographic scenario
                    % scenEXT=1 no extreme event
                    if scenEXT==2  % historical extreme event -> BS
                        THETA(2,PERTclimat(x,:)==1)=0;
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

                    % Pop projection

                    [Fav Mav]=breederav2(THETA,nt); % calculate males and females available to mate

                    for t=1:nt % from 1900 to 2100
                        %calculate demographic parameters
                        uf=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Fav(:,t)'*wvec);
                        um=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Mav(:,t)'*wvec);

                        % project population
                        wvec=popmat(THETA(:,t),uf,um)*wvec; % define population matrix
                        r(t,1)=sum(wvec);  % calculate growth rate
                        wvec=wvec/r(t,1);
                        w(:,t,1)=wvec;
                    end % time

                    Rtot(:, s)=r; %R for each simulation for each t for one ens
                    WVEC(:,:,s)=w;
                    x=x+1;

                end % demographic sim
                
                MAT=WVEC(2,:,:)+WVEC(4,:,:)+WVEC(5,:,:); %proportion of matures individuals
                MAT=permute(MAT, [2 3 1]);
                % Save results for each mod, scenEXT, col, sim
                filename=sprintf('%s/Codes_EP/GR_results_tot/mean_ensemble/GR_results_Stef/%s/GRmean_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);

                filename2=sprintf('%s/Codes_EP/MAT_results_tot/mean_ensemble/MAT_results_Stef/%s/MAT_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);

                writematrix(Rtot, filename);
                writematrix(MAT, filename2);

            
        end % colony
    end  % extreme event scenario
end %model




%% START IN 2009


%load informations on climate models (years, number of ens)
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi));

%load the proportion of extreme events for scenEXT4 for each climate model
load proportion_ext_event.mat

model= ["LE_CanESM2","LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"] 
                                                                                                                                                                                                                                                                                                                          ck

for i=[1:3 6 7] %Climate model
    mod=model(i); %name
    display(mod)
    
    % Standardized SIC data (Bilgecan method), 500 and 850km2
    if i==6
            file_name=sprintf(sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod));
    else
        file_name=sprintf(sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis.mat', ordi, mod));
    end

    % if i==2 % Data for ToE, start in 1900
    %         file_name=sprintf(sprintf('/Users/%s/Desktop/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis1900.mat', ordi, mod));
    % end

    load(file_name)

    load NewNameColony.mat

    % colony
    yearstart=info_models(i,3);
    yearstart_bis = 2009; %Here we want to start in 2009 for all models
    y_bis = length(yearstart:yearstart_bis); % find year 2009 in SIC projections

    yearend=info_models(i,4);
    % if i==7
    %     yearstart=2006;
    % end

    timePOP=yearstart_bis:yearend;
    nt=length(timePOP);
    nc=66;

    % climate
    ntime=yearstart_bis:yearend;
    nens=info_models(i,5); % climate ensemble, nens minimum is 5
    nsim=10;
    nsimE=10;
    nsimulation=nsim*nsimE;

    % perturbation - extreme events
    percentSI=13/355; % =3.66% of extreme events counted on historical period
    PERTclimat=binornd(1,percentSI,nsimulation,nt);
    y2009=find(ntime==2009);y2018=find(ntime==2018);


    % Get covariance matrix of demo parameter
    % uncertainties on the estimated parameters
    [bR,bS]=beta; % estimate CMR, return and survival
    xA=0.058; % tag loss
    [betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);


    % MAIN LOOP

    SCEN=struct; 

    for scenEXT=[1 2 4]  %extreme event scenario (# 5)
        display(scenEXT)

        for col=1:nc % colony
            display(col)

            parfor ens=1:nens % climate ensemble
                
                % SEA ICE FORECAST ENS (climatic ens 50) matrix SICa
                % time*season*colony

                % Extreme sea ice loss threshold
                threshold_SIC2=prctile(ENS(ens).SICa(y2009:y2018,2,col),percentSI);
                threshold_SIC4=prctile(ENS(ens).SICa(y2009:y2018,4,col),percentSI);

                EXT=find(ENS(ens).SICa(y_bis:y_bis+nt-1,2,col)<threshold_SIC2 | ENS(ens).SICa(y_bis:y_bis+nt-1,4,col)<threshold_SIC4);
                proportion = proportion_ext_event(i); %Calculated to have same nb of colonies having total breeding failure 2018-2022 than observed
                iEXT=EXT(binornd(1,proportion,length(EXT),1)==1); % proportion of extreme events leading to total breeding failure

                % SIMULATIONS

                Rtot=zeros(nt, nsimulation);
                WVEC=zeros(5, nt, nsimulation); %to save SAD
                x=1;

                for s=1:nsimulation % sampling parameter uncertainty

                    % parameter uncertainty
                    [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
                    %breeding success, proportion of 1st return, Return, survival

                    % POP INITIAL CONDITIONS
                    r=zeros(nt,1);
                    w=zeros(5, nt);

                    % SAD
                    wvec=simpsamp(5)';

                    % vital rates sampling stochasticity
                    [THETA]=parameter_RAND3(xA,ENS(ens).SICa(y_bis:y_bis+nt-1,:,col),bR,bS,SICf,SICm,W,nt);

                    % Extreme demographic scenario
                    % scenEXT=1 no extreme event
                    if scenEXT==2  % historical extreme event -> BS
                        THETA(2,PERTclimat(x,:)==1)=0;
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

                    % Pop projection

                    [Fav Mav]=breederav2(THETA,nt); % calculate males and females available to mate

                    for t=1:nt % from 1900 to 2100
                        %calculate demographic parameters
                        uf=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Fav(:,t)'*wvec);
                        um=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Mav(:,t)'*wvec);

                        % project population
                        wvec=popmat(THETA(:,t),uf,um)*wvec; % define population matrix
                        r(t,1)=sum(wvec);  % calculate growth rate
                        wvec=wvec/r(t,1);
                        w(:,t,1)=wvec;
                    end % time

                    Rtot(:, s)=r; %R for each simulation for each t for one ens
                    WVEC(:,:,s)=w;
                    x=x+1;

                end % demographic sim

                MAT=WVEC(2,:,:)+WVEC(4,:,:)+WVEC(5,:,:); %proportion of matures individuals
                MAT=permute(MAT, [2 3 1]);
                % Save results for each mod, scenEXT, col, sim
                filename=sprintf('%s/Codes_EP/GR_results_tot/2009/GR_results_Stef/%s/GR_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);
                filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/MAT_results_Stef/%s/MAT_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);

                writematrix(Rtot, filename);
                writematrix(MAT, filename2);

            end  % climate ensemble
        end % colony
    end  % extreme event scenario
end %model


%% MEAN OF ENSEMBLES 2009

%load informations on climate models (years, number of ens)
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi));

%load the proportion of extreme events for scenEXT4 for each climate model
load proportion_ext_event.mat

model= ["LE_CanESM2","LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"] 

for i=[1:3 6 7] %Climate model
    mod=model(i); %name
    display(mod)
    
    % Standardized SIC data (Bilgecan method), 500 and 850km2
    file_name=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis.mat', ordi, mod);

    if i==6 %CESM2
        file_name=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod);
    end
    load(file_name)

    load NewNameColony.mat

    % colony
    yearstart_bis=info_models(i,3);
    yearstart=2009;
    y_bis = length(yearstart_bis:yearstart);

    yearend=info_models(i,4);
    timePOP=yearstart:yearend;
    nt=length(timePOP);
    nc=66;

    % climate
    ntime=yearstart:yearend;
    nens=info_models(i,5); % climate ensemble, nens minimum is 5
    nsim=100;
    nsimulation=100;

    % perturbation - extreme events
    percentSI=13/355; % =3.66% of extreme events counted on historical period
    PERTclimat=binornd(1,percentSI,nsimulation,nt);
    y2009=find(ntime==2009);y2018=find(ntime==2018);


    % Get covariance matrix of demo parameter
    % uncertainties on the estimated parameters
    [bR,bS]=beta; % estimate CMR, return and survival
    xA=0.058; % tag loss
    [betaBS covBS betaPB covPB betaR covR betaS covS W SICf SICm array]=covbeta(nsim);


    % Calculate mean of ensemble
    % Select years 2009-2100
    
    SICa_mean = zeros(nt, 4, nc);
    for ens = 1:nens
        SICa_mean = SICa_mean + ENS(ens).SICa(y_bis:y_bis+nt-1,:,:);
    end
    SICa_mean = SICa_mean ./ nens;


    % MAIN LOOP

    tic %start clock                                                                                                                                                                                                                                                                                                                            ck
    SCEN=struct; 

    for scenEXT=[1 2 4]  %extreme event scenario (# 5)
        display(scenEXT)

        for col=1:nc % colony
            display(col)

                
                % SEA ICE FORECAST ENS (climatic ens 50) matrix SICa
                % time*season*colony

                % Extreme sea ice loss threshold
                threshold_SIC2=prctile(SICa_mean(y2009:y2018,2,col),percentSI);
                threshold_SIC4=prctile(SICa_mean(y2009:y2018,4,col),percentSI);

                EXT=find(SICa_mean(:,2,col)<threshold_SIC2 | SICa_mean(:,4,col)<threshold_SIC4);
                proportion = proportion_ext_event(i); %Calculated to have same nb of colonies having total breeding failure 2018-2022 than observed
                iEXT=EXT(binornd(1,proportion,length(EXT),1)==1); % proportion of extreme events leading to total breeding failure

                % SIMULATIONS

                Rtot=zeros(nt, nsimulation);
                WVEC=zeros(5, nt, nsimulation); %to save SAD
                x=1;

                for s=1:nsimulation % sampling parameter uncertainty

                    % parameter uncertainty
                    [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS);
                    %breeding success, proportion of 1st return, Return, survival

                    % POP INITIAL CONDITIONS
                    r=zeros(nt,1);
                    w=zeros(5,nt);

                    % SAD
                    wvec=simpsamp(5)';

                    % vital rates sampling stochasticity
                    [THETA]=parameter_RAND3(xA,SICa_mean(:,:,col),bR,bS,SICf,SICm,W,nt);

                    % Extreme demographic scenario
                    % scenEXT=1 no extreme event
                    if scenEXT==2  % historical extreme event -> BS
                        THETA(2,PERTclimat(x,:)==1)=0;
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

                    % Pop projection

                    [Fav Mav]=breederav2(THETA,nt); % calculate males and females available to mate

                    for t=1:nt % from 1900 to 2100
                        %calculate demographic parameters
                        uf=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Fav(:,t)'*wvec);
                        um=min(Fav(:,t)'*wvec, Mav(:,t)'*wvec)/(Mav(:,t)'*wvec);

                        % project population
                        wvec=popmat(THETA(:,t),uf,um)*wvec; % define population matrix
                        r(t,1)=sum(wvec);  % calculate growth rate
                        wvec=wvec/r(t,1);
                        w(:,t,1)=wvec;
                    end % time

                    Rtot(:, s)=r; %R for each simulation for each t for one ens
                    WVEC(:,:,s)=w;

                    x=x+1;

                end % demographic sim

                MAT=WVEC(2,:,:)+WVEC(4,:,:)+WVEC(5,:,:); %proportion of matures individuals
                MAT=permute(MAT, [2 3 1]);
                % Save results for each mod, scenEXT, col, sim
                filename=sprintf('%s/Codes_EP/GR_results_tot/2009/mean_ensemble/GR_results_Stef/%s/GRmean_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);

                filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/mean_ensemble/MAT_results_Stef/%s/MAT_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);

                writematrix(Rtot, filename);
                writematrix(MAT, filename2);

        end % colony
    end  % extreme event scenario
end %model

toc
