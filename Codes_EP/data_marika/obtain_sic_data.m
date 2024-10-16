% In this file, code to transform ncdf files into matlab files, to do the
% transformation to have the same variance and mean, to add 20 years to the
% CESM1 model for the time of emergence project.

% DIRECTORY to Codes_EP
ordi = "/Users/aliceeparvier/Desktop"

%% Transform Marika's data into matlab files

% No transformation of the data here (variance and mean)
% Non-breeding season = 850 km2 area
% Other seasons = 500 km2 area

addpath(sprintf('%s/Codes_EP/data_marika', ordi));

model= ["CanESM2", "CESM1-CAM5", "CSIRO-Mk3-6-0", "GFDL-CM3", "MPI-ESM"];
n_mod=length(model); % #models

info_models=zeros(n_mod,6); %table with col: mod name, n_year, year_start, year_end, n_memb, col indexes 1->all good

i=1;

for mod=[model]

    % Import marika's data
    fname1 = sprintf('%s_colony_sic_%s.nc', mod, 'poly500'); %sic data 500km2
    fname2 = sprintf('%s_colony_sic_%s.nc', mod, 'poly850'); %850km2
    fname3 = sprintf('tas_regavg_NH_%s.nc', mod); %temp data

    % sic : month, year, ens, col
    [year1, num1, sic1]=read_climat_nc(fname1); %read from Marika's nc files
    [year2, num2, sic2]=read_climat_nc(fname2);

    % temp : month, year, ens
    [year, temp]=read_temp_nc(fname3);
    

    % Create mat data
    nt=length(year1);
    nens=length(temp(1,1,:));
    y1979=nt-121;
    y2007=y1979+28;

    ENS=struct;
    for ens=1:nens
        ENS(ens).SIC=zeros(nt, 4, 66); %(nt, season, col)
        ENS(ens).SICa=zeros(nt, 4, 66);

        %temp
        ENS(ens).temp=zeros(nt, 4);
        ENS(ens).tempa=zeros(nt, 4);
    end %ens

    seasons={1:3, 4:5, 6:7, 8:12}; %from months to 4 seasons
    for s = 1:4 %season

        if s==1 %Non breeding season -> 750km

            for ens=1:nens
                sic_season2=mean(sic2(seasons{s}, :, ens, :)); %mean on months to have one value for the season 850km2
                ENS(ens).SIC(:, s, :)=sic_season2;
                ENS(ens).SICa(:,s, :)=(sic_season2-mean(sic_season2(1,y1979:y2007,1, :)))./mean(sic_season2(1,y1979:y2007,1, :)); %anomalies, SICmean=1979-2007

                %temperature
                temp_season=mean(temp(seasons{s}, :, ens));
                ENS(ens).temp(:, s)=temp_season;
                ENS(ens).tempa(:,s)=(temp_season-mean(temp_season(1,y1979:y2007 ,1)))./mean(temp_season(1,y1979:y2007 ,1));
            end %ens


        else %other 3 seasons
            for ens=1:nens
                sic_season=mean(sic1(seasons{s}, :, ens, :)); %mean on months to have one value for the season
                ENS(ens).SIC(:, s, :)=sic_season;
                ENS(ens).SICa(:,s, :)=(sic_season-mean(sic_season(1,y1979:y2007,1, :)))./mean(sic_season(1,y1979:y2007,1, :)); %anomalies, SICmean=1979-2007

                %temperature
                temp_season=mean(temp(seasons{s}, :, ens));
                ENS(ens).temp(:, s)=temp_season;
                ENS(ens).tempa(:,s)=(temp_season-mean(temp_season(1,y1979:y2007 ,1)))./mean(temp_season(1,y1979:y2007 ,1));
            end %ens

        end

        end %season

        % Verify colonies indexes between Marika's data and Stef indexes
        colind = 0;
        v= zeros(66,1);
        v(1:66)=(1:66);
        if num1==v %if same indexes
            colind=1;
        end

        % Info model
        info_models(i,1)=i; % #model
        info_models(i,2)=nt; % ntime
        info_models(i,3)=year(1); %year start
        info_models(i,4)=year(end); %year end
        info_models(i,5)=nens; %n ens members
        info_models(i,6)=colind;%1 if good indexes

    % Save the mat data
    mfile = sprintf('%s/Codes_EP/Codes_Stef/mat_data/not_std/LE_%s_seaice_bis.mat', ordi, mod); %create the mat file name, data not standardized
    save(mfile, 'ENS');
    i=i+1;
    
end %model

 %save("/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/mat_data/infos_models.mat", 'info_models')



%% Create the mat file for PARIS AGREEMENT 66 COLONIES
% Read Marika's data 2006-2100

% To calculate the anomalies, we take the median of the 50 ens mean
% 1979-2007 from the CESM2-LE (the Paris agreement runs were initialized
% from 2006 CESM2-LE).


load(sprintf('%s/Codes_EP/Codes_Stef/mat_data/not_std/LE_CESM1-CAM5_seaice_bis.mat', ordi))
ENS_hist=ENS;
nens_hist=40;
nt_hist=181;
ncol=66;
y1979=nt_hist-121;
y2007=y1979+29;

mod= ["LE_CESM1-CAM5-PARIS2"];


% Import marika's data
fname1 = sprintf('%s_colony_sic_%s.nc', mod, 'poly500'); %sic data 500km2
fname2 = sprintf('%s_colony_sic_%s.nc', mod, 'poly850'); %850km2

% sic : month, year, ens, col
[year1, num1, sic1]=read_climat_nc(fname1); %read from Marika's nc files
[year2, num2, sic2]=read_climat_nc(fname2);


% Create mat data
nt=length(year1);
nens=length(sic1(1,1,:,1));


ENS=struct;
for ens=1:nens
    ENS(ens).SIC=zeros(nt, 4, 66); %(nt, season, col)
    ENS(ens).SICa=zeros(nt, 4, 66);

end %ens

seasons={1:3, 4:5, 6:7, 8:12}; %from months to 4 seasons

for s = 1:4 %season
    SIC_hist = zeros(ncol, nens_hist);
    for ens=1:nens_hist
        SIC_hist(:, ens)= mean(ENS_hist(ens).SIC(y1979:y2007, s, :));
    end
    historical_mean = median(SIC_hist(:,:),2);

    if s==1 %Non breeding season -> 750km
        
        for ens=1:nens
            sic_season2=mean(sic2(seasons{s}, :, ens, :)); %mean on months to have one value for the season 850km2
            ENS(ens).SIC(:, s, :)=sic_season2;
            for col=1:66
                ENS(ens).SICa(:,s, col)=(ENS(ens).SIC(:,s,col)-historical_mean(col))./historical_mean(col); %anomalies, SICmean=1979-2007
            end
        end %ens


    else %other 3 seasons
        for ens=1:nens
            sic_season=mean(sic1(seasons{s}, :, ens, :)); %mean on months to have one value for the season
            ENS(ens).SIC(:, s, :)=sic_season;
            for col=1:66
                ENS(ens).SICa(:,s, col)=(ENS(ens).SIC(:,s,col)-historical_mean(col))./historical_mean(col); %anomalies, SICmean=1979-2007
            end
        end %ens

    end

end %season


% Save the mat data
mfile = sprintf('/Users/aliceeparvier/Desktop/Codes_EP/Codes_Stef/mat_data/not_std/%s_seaice_bis.mat', mod); %create the mat file name
%save(mfile, 'ENS');


% -------------------------------------------------------------------------
%% VARIANCE CORRECTION - STANDARDIZATION OF THE DATA

% modify the data with Bilgecan method to have same variance and mean than observations
% (sigma_obs_SIC / sigma_model_SIC) * (SIC[i,j] - mu_model_SIC) + mu_obs_SIC
% Out : mat_data/LE_CESM2_seaice_std.mat


%% LE CESM2

% Observations
mu_obs = zeros(66,4);
sigma_obs = zeros(66,4);
y1979=1;
y2007=29;
for col=1:66
    file=sprintf('%s/Codes_EP/Codes_Stef/SIC_obs_allcol/SIC_col%d.csv', ordi, col);
    SIC = readmatrix(file);
    mu_obs(col,:) = mean(SIC(y1979:y2007, :)); %mean on 1980-2007 for each season each colony
    sigma_obs(col,:) = std(SIC(y1979:y2007, :));
end


% CESM2 LE
load(sprintf('%s/Codes_EP/Codes_Stef/seaiceLE.mat', ordi);
ENS_old = ENS;

% Create mat data
nt=201;
nens=50;
y1979=nt-121;
y2007=y1979+29;
ENS=struct;

for s = 1:4 %season
    for ens=1:nens
        for col=1:66
            mu_model = mean(ENS_old(ens).SIC(y1979:y2007,s,col));
            sigma_model = std(ENS_old(ens).SIC(y1979:y2007,s,col));

            % TRANSFORMATION
            ENS(ens).SIC(:,s,col) = ((sigma_obs(col,s) / sigma_model) * (ENS_old(ens).SIC(:,s,col) - mu_model)) + mu_obs(col,s);

            % Make sure the data is between 0 and 1
            ENS(ens).SIC(ENS(ens).SIC(:,s,col) >1, s, col)=1;
            ENS(ens).SIC(ENS(ens).SIC(:,s,col) <0, s, col)=0;
        end %col

        % ANOMALIES
        ENS(ens).SICa(:, s, :)=(ENS(ens).SIC(:,s,:)-mean(ENS(ens).SIC(y1979:y2007,s,:)))./mean(ENS(ens).SIC(y1979:y2007,s,:)); %anomalies, SICmean=1979-2007

    end %ens
end %season

mfile = sprintf('%s/Codes_EP/Codes_Stef/mat_data/LE_CESM2_seaice_std.mat', ordi); %create the mat file name
%save(mfile, 'ENS');



%% The other 3 models CanESM, CESM1, CSIRO

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];

mu_obs = zeros(66,4);
sigma_obs = zeros(66,4);
y1979=1;
y2007=29;

for col=1:66
    file=sprintf('%s/Codes_EP/Codes_Stef/SIC_obs_allcol/SIC_col%d.csv', ordi, col);
    SIC = readmatrix(file);
    mu_obs(col,:) = mean(SIC(y1979:y2007, :)); %mean on 1980-2007 for each season each colony
    sigma_obs(col,:) = std(SIC(y1979:y2007, :));
end


for i=1 %models

    mod=model(i)
    filename=sprintf('%s/Codes_EP/Codes_Stef/mat_data/not_std/%s_seaice_bis.mat', ordi, mod);
    load(filename)

    ENS_old = ENS;

    % Create mat data
    nt=info_models(i, 2);
    nens=info_models(i, 5);
    y1979=nt-121;
    y2007=y1979+29;

    ENS=struct;

    for s = 1:4 %season
        for ens=1:nens
            for col=1:66
                mu_model = mean(ENS_old(ens).SIC(y1979:y2007,s,col));
                sigma_model = std(ENS_old(ens).SIC(y1979:y2007,s,col));

                % TRANSFORMATION
                ENS(ens).SIC(:,s,col) = ((sigma_obs(col,s) / sigma_model) * (ENS_old(ens).SIC(:,s,col) - mu_model)) + mu_obs(col,s);
                ENS(ens).SIC(ENS(ens).SIC(:,s,col) >1, s, col)=1;
                ENS(ens).SIC(ENS(ens).SIC(:,s,col) <0, s, col)=0;
            end %col

            % ANOMALIES
            ENS(ens).SICa(:, s, :)=(ENS(ens).SIC(:,s,:)-mean(ENS(ens).SIC(y1979:y2007,s,:)))./mean(ENS(ens).SIC(y1979:y2007,s,:)); %anomalies, SICmean=1979-2007

        end %ens
    end %season

    % bis = mix between 500 and 850 km2
    mfile = sprintf('%s/Codes_EP/Codes_Stef/mat_data/%s_seaice_std_bis.mat', ordi, mod); %create the mat file name
    save(mfile, 'ENS');

end %model


%% Finally, for the Paris Agreement model
 mod="LE_CESM1-CAM5-PARIS2";

 mu_obs = zeros(66,4);
 sigma_obs = zeros(66,4);
 y1979=1;
 y2007=29;

 for col=1:66
     file=sprintf('%s/Codes_EP/Codes_Stef/SIC_obs_allcol/SIC_col%d.csv', ordi, col);
     SIC = readmatrix(file);
     mu_obs(col,:) = mean(SIC(y1979:y2007, :)); %mean on 1980-2007 for each season each colony
     sigma_obs(col,:) = std(SIC(y1979:y2007, :));
 end


 load('%s/Codes_EP/Codes_Stef/mat_data/not_std/LE_CESM1-CAM5_seaice_bis.mat', ordi)
 ENS_hist=ENS;

 load('%s/Codes_EP/Codes_Stef/mat_data/LE_CESM1-CAM5_seaice_std_bis.mat', ordi)
 ENS_hist_std=ENS;

 filename=sprintf('%s/Codes_EP/Codes_Stef/mat_data/not_std/%s_seaice_bis.mat', ordi, mod);
 load(filename)

 ENS_old = ENS;

 % Create mat data
 nt=95;
 nens=11;

 nt_hist=181;
 nens_hist=40;
 y1979 = nt_hist - 121;
 y2007 = y1979 + 29;

 ENS=struct;

 for s = 1:4 %season

     % Calculate mu_model and sigma_model with CESM1-CAM5 (median of ens)
     mu_model=zeros(4, ncol);
     for col=1:ncol
         mu=zeros(nens_hist,1);
         sigma=zeros(nens_hist,1);
         for ens=1:nens_hist
             mu(ens) = mean(ENS_hist(ens).SIC(y1979:y2007,s,col));
             sigma(ens) = std(ENS_hist(ens).SIC(y1979:y2007,s,col));
         end
         mu_model(s,col) = median(mu);
         sigma_model(s,col) = median(sigma);
     end

     % Calculate historical mean std to calculate the anomalies
    SIC_hist_std = zeros(ncol, nens_hist);
    for ens=1:nens_hist
        SIC_hist_std(:, ens)= mean(ENS_hist_std(ens).SIC(y1979:y2007, s, :));
    end
    historical_mean = median(SIC_hist_std(:,:),2);

     for ens=1:nens
         for col=1:66

             % TRANSFORMATION
             ENS(ens).SIC(:,s,col) = ((sigma_obs(col,s) / sigma_model(s,col)) * (ENS_old(ens).SIC(:,s,col) - mu_model(s,col))) + mu_obs(col,s);
             ENS(ens).SIC(ENS(ens).SIC(:,s,col) >1, s, col)=1;
             ENS(ens).SIC(ENS(ens).SIC(:,s,col) <0, s, col)=0;

         end %col

         % ANOMALIES
         for col=1:ncol
            ENS(ens).SICa(:, s, col)=(ENS(ens).SIC(:,s,col)-historical_mean(col))./historical_mean(col); %anomalies, SICmean=1979-2007
         end
     end %ens
 end %season

 % bis = mix between 500 and 850 km2
 mfile = sprintf('%s/Codes_EP/Codes_Stef/mat_data/%s_seaice_std_bis.mat', ordi, mod); %create the mat file name
 save(mfile, 'ENS');





% ------------------------------------------------------------
% Time of emergence project, we needed the CESM1 model to start in 1900, so
% copied the data from 1920 to 1940 to 1900 to 1920

%% Add 20 years to the CESM1 model so that it starts in 1900 (data from 1900 to 1919 will be that of 1920-1940)

load(sprintf("%s/Codes_EP/Codes_Stef/mat_data/LE_CESM1-CAM5_seaice_std_bis.mat", ordi));

yearstart=1920;
yearend=2100;
nt=length(yearstart:yearend);
nens=40;

ENS_1900=struct; % New data

for ens=1:nens
    for s=1:4 %season
            ENS_1900(ens).SIC(1:20, s, :) = ENS(ens).SIC(1:20, s, :);
            ENS_1900(ens).SICa(1:20, s, :) = ENS(ens).SICa(1:20, s, :);
    end

    % Add 1920-2100 data
    ENS_1900(ens).SIC(21:201, :, :) = ENS(ens).SIC(:, :, :);
    ENS_1900(ens).SICa(21:201, :, :) = ENS(ens).SICa(:, :, :);
end

ENS = ENS_1900;

%save(sprintf("%s/Codes_EP/Codes_Stef/mat_data/LE_CESM1-CAM5_seaice_std_bis1900.mat", ordi),'ENS',)