%% IPM Model, designed by Francesco Ventura
% Final version, save one txt file per ens, per col

%% Main loop - GR

% DIRECTORY to Codes_EP
ordi='/Users/aliceeparvier/Desktop';

% yearstart needs to be one year after the start of your environmental data
% because of the fact that wind vwind year t affect year t+1
yearstart=1921;
yearstop=2100; %2020 %2100
timePOP=yearstart:yearstop;
nt=length(timePOP);

% Parameter posteriors from Francesco
% col 1,2,3,4 beta0 phi, betasst, betawind, sigma2
% col 5,6,7 beta0breed betavwind sigma2breed
% col 8,9,10 beta0fec betaNOW sigma2fec
% col 11 mu phi1
param_post=readmatrix("param_postIPM2.csv");
param_post(1,:)=[];
param_post(:,1)=[];

% Simulation parameters
nsimulation=100;
nens=50;
ncol = 66;

% Nearest open water - optimistic
NOW = zeros(50, nt+1); %2020-2100 to project on 2021-2100
% for t=1:length(NOW(1,:))
%     NOW(:,:)=-0.8;
% end


for col = 57:66
    display(col)

    % import climate data
    file_name_temp = sprintf('IPM_R/All_colonies/env_data_all_col_csv/SST_forecast_1920-2100_col%d.csv', col);
    SST = readmatrix(file_name_temp);
    SST(1,:)=[]; SST(:,1)=[];

    file_name_wind = sprintf('IPM_R/All_colonies/env_data_all_col_csv/wind_forecast_1920-2100_col%d.csv', col);
    wind = readmatrix(file_name_wind);
    wind(1,:)=[]; wind(:,1)=[];

    file_name_vwind = sprintf('IPM_R/All_colonies/env_data_all_col_csv/vwind_forecast_1920-2100_col%d.csv', col);
    vwind = readmatrix(file_name_vwind);
    vwind(1,:)=[]; vwind(:,1)=[];

    % Data start in 2020, projections in 2021
    % Data start in 2014
    SST(:,1:length(1979:yearstart)-2)=[];
    vwind(:,1:length(1979:yearstart)-2)=[];
    wind(:,1:length(1979:yearstart)-2)=[];

    % Main loop - GR

    for ens=1:50

        Rtot=zeros(nt, nsimulation);
        WVEC=zeros(7, nt, nsimulation); %to save SAD

        for s=1:nsimulation

            % POP INITIAL CONDITIONS
            r=zeros(nt,1);
            w=zeros(7, nt);

            %    1        2    3   4    5   6  7
            % Breeders  Non-B  J1  J2  J3  J4  J5
            wvec=simpsamp(7)';

            % theta
            %   1     2     3      4      5      6
            % [fec, phiA, phiJ, Pbreed, recr5, recr6]
            [THETA]=parameter_RAND_IPM(SST(ens,2:end), wind(ens,1:end-1), vwind(ens,1:end-1), NOW(ens,2:end), param_post, nt);

            % POP PROJECTION

            for t=1:nt % from year start 2021 to 2099
                wvec=popmatIPM(THETA(:,t))*wvec; % define population matrix
                r(t,1)=sum(wvec);  % calculate growth rate
                wvec=wvec/r(t,1);
                w(:,t,1)=wvec;
            end % time

            Rtot(:, s)=r;
            WVEC(:,:,s)=w;

        end %sim

        MAT=WVEC(1,:,:)+WVEC(2,:,:); %proportion of matures individuals
        MAT=permute(MAT, [2 3 1]);

        % Save results for each ens, col, sim
        filename=sprintf('%s/Codes_EP/GR_results_tot/GR_results_IPM/CESM2_1920/GR_IPM_CESM2_ens%d_col%d.txt', ordi, ens, col);
        filename2=sprintf('%s/Codes_EP/MAT_results_tot/MAT_results_IPM/CESM2_1920/MAT_IPM_CESM2_ens%d_col%d.txt', ordi, ens, col);

        % Folders CESM2_1980, CESM2_1920, CESM2 starts in 2021
        writematrix(Rtot, filename); %(nt, nsim)
        writematrix(MAT, filename2);

    end %ens

end %col





%% Main loop from 2009

% colony
yearstart=2009; 
yearstop=2100;
timePOP=yearstart:yearstop;
nt=length(timePOP);

% Parameter posteriors from Francesco
% col 1,2,3,4 beta0 phi, betasst, betawind, sigma2
% col 5,6,7 beta0breed betavwind sigma2breed
% col 8,9,10 beta0fec betaNOW sigma2fec
% col 11 mu phi1
param_post=readmatrix("param_postIPM2.csv");
param_post(1,:)=[];
param_post(:,1)=[];

% Simulation parameters
nsimulation=100;
nens=50;
ncol = 66;

% optimistic
NOW = zeros(50, nt+1); %2020-2100 to project on 2021-2100
% for t=1:length(NOW(1,:))
%     NOW(:,:)=-0.8;
% end


for col = 1:66

    display(col)
    file_name_temp = sprintf('IPM_R/All_colonies/env_data_all_col_csv/SST_forecast_col%d.csv', col);
    SST = readmatrix(file_name_temp);
    SST(1,:)=[]; SST(:,1)=[];

    file_name_wind = sprintf('IPM_R/All_colonies/env_data_all_col_csv/wind_forecast_col%d.csv', col);
    wind = readmatrix(file_name_wind);
    wind(1,:)=[]; wind(:,1)=[];

    file_name_vwind = sprintf('IPM_R/All_colonies/env_data_all_col_csv/vwind_forecast_col%d.csv', col);
    vwind = readmatrix(file_name_vwind);
    vwind(1,:)=[]; vwind(:,1)=[];

    % Data start in 2020, projections in 2021
    % Data start in 2014
    SST(:,1:length(1979:yearstart)-2)=[];
    vwind(:,1:length(1979:yearstart)-2)=[];
    wind(:,1:length(1979:yearstart)-2)=[];

    % Main loop - GR

    for ens=1:50

        Rtot=zeros(nt, nsimulation);
        WVEC=zeros(7, nt, nsimulation); %to save SAD

        for s=1:nsimulation

            % POP INITIAL CONDITIONS
            r=zeros(nt,1);
            w=zeros(7, nt);

            %    1        2    3   4    5   6  7
            % Breeders  Non-B  J1  J2  J3  J4  J5
            wvec=simpsamp(7)';

            % theta
            %   1     2     3      4      5      6
            % [fec, phiA, phiJ, Pbreed, recr5, recr6]
            [THETA]=parameter_RAND_IPM(SST(ens,2:end), wind(ens,1:end-1), vwind(ens,1:end-1), NOW(ens,2:end), param_post, nt);

            % POP PROJECTION

            for t=1:nt % from year start 2021 to 2099
                wvec=popmatIPM(THETA(:,t))*wvec; % define population matrix
                r(t,1)=sum(wvec);  % calculate growth rate
                wvec=wvec/r(t,1);
                w(:,t,1)=wvec;
            end % time

            Rtot(:, s)=r;
            WVEC(:,:,s)=w;

        end %sim

        MAT=WVEC(1,:,:)+WVEC(2,:,:); %proportion of matures individuals
        MAT=permute(MAT, [2 3 1]);

        % Save results for each ens, col, sim
        filename=sprintf('%s/Codes_EP/GR_results_tot/2009/GR_results_IPM/CESM2/GR_IPM_CESM2_ens%d_col%d.txt', ordi, ens, col);
        filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/MAT_results_IPM/CESM2/MAT_IPM_CESM2_ens%d_col%d.txt', ordi, ens, col);
        
        writematrix(Rtot, filename); %(nt, nsim)
        writematrix(MAT, filename2);
        
    end %ens

end %col


%% MEAN OF ENSEMBLES

% colony
yearstart=2009; %1980 %2009
yearstop=2100; %2020 %2100
timePOP=yearstart:yearstop;
nt=length(timePOP);

% Parameter posteriors from Francesco
% col 1,2,3,4 beta0 phi, betasst, betawind, sigma2
% col 5,6,7 beta0breed betavwind sigma2breed
% col 8,9,10 beta0fec betaNOW sigma2fec
% col 11 mu phi1
param_post=readmatrix("param_postIPM2.csv");
param_post(1,:)=[];
param_post(:,1)=[];

% Simulation parameters
nsimulation=100;
nens=50;
ncol = 66;

% optimistic
NOW = zeros(50, nt+1); %2020-2100 to project on 2021-2100
% for t=1:length(NOW(1,:))
%     NOW(:,:)=-0.8;
% end


for col = 1:66

    display(col)
    file_name_temp = sprintf('IPM_R/All_colonies/env_data_all_col_csv/SST_forecast_col%d.csv', col);
    SST = readmatrix(file_name_temp);
    SST(1,:)=[]; SST(:,1)=[];

    file_name_wind = sprintf('IPM_R/All_colonies/env_data_all_col_csv/wind_forecast_col%d.csv', col);
    wind = readmatrix(file_name_wind);
    wind(1,:)=[]; wind(:,1)=[];

    file_name_vwind = sprintf('IPM_R/All_colonies/env_data_all_col_csv/vwind_forecast_col%d.csv', col);
    vwind = readmatrix(file_name_vwind);
    vwind(1,:)=[]; vwind(:,1)=[];

    % Data start in 2020, projections in 2021
    % Data start in 2014
    SST(:,1:length(1979:yearstart)-2)=[];
    vwind(:,1:length(1979:yearstart)-2)=[];
    wind(:,1:length(1979:yearstart)-2)=[];

    % Main loop - GR


    % Calculate mean of ensemble
    
    SST_mean = zeros(1,nt+1);
    for ens = 1:nens
        SST_mean = SST_mean + SST(ens,:);
    end
    SST_mean = SST_mean ./ nens;


    vwind_mean = zeros(1, nt+1);
    for ens = 1:nens
        vwind_mean = vwind_mean + vwind(ens,:);
    end
    vwind_mean = vwind_mean ./ nens;


    wind_mean = zeros(1, nt+1);
    for ens = 1:nens
        wind_mean = wind_mean + wind(ens,:);
    end
    wind_mean = wind_mean ./ nens;


        Rtot=zeros(nt, nsimulation);
        WVEC=zeros(7, nt, nsimulation); %to save SAD

        for s=1:nsimulation

            % POP INITIAL CONDITIONS
            r=zeros(nt,1);
            w=zeros(7, nt);
            
            %    1        2    3   4    5   6  7
            % Breeders  Non-B  J1  J2  J3  J4  J5
            wvec=simpsamp(7)';

            % theta
            %   1     2     3      4      5      6
            % [fec, phiA, phiJ, Pbreed, recr5, recr6]
            [THETA]=parameter_RAND_IPM(SST(1,2:end), wind(1,1:end-1), vwind(1,1:end-1), NOW(1,2:end), param_post, nt);

            % POP PROJECTION

            for t=1:nt % from year start 2021 to 2099
                wvec=popmatIPM(THETA(:,t))*wvec; % define population matrix
                r(t,1)=sum(wvec);  % calculate growth rate
                wvec=wvec/r(t,1);
                w(:,t,1)=wvec;
            end % time

            Rtot(:, s)=r;
            WVEC(:,:,s)=w;

        end %sim

        MAT=WVEC(1,:,:)+WVEC(2,:,:); %proportion of matures individuals
        MAT=permute(MAT, [2 3 1]);

        % Save results for each ens, col, sim
        filename=sprintf('%s/Codes_EP/GR_results_tot/2009/mean_ensemble/GR_results_IPM/CESM2/GRmean_IPM_CESM2_col%d.txt', ordi, col);
        filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/mean_ensemble/MAT_results_IPM/CESM2/MAT_IPM_CESM2_col%d.txt', ordi, col);
        
        writematrix(Rtot, filename); %(nt, nsim)
        writematrix(MAT, filename2);


end %col


