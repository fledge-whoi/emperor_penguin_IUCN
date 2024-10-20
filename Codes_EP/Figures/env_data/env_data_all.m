%% Figures to visualize the environmental data

% DIRECTORY to Codes_EP
ordi = "/Volumes/My Passport";

load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/SICaOBS.mat', ordi)) %1980-2011
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi))

%% Sea ice concentratino - Yearly estimates
% Figure mean of ens, mean of colonies for all years and 4 separate seasons 
% Change SIC or SICa depending on what you want 

couleur = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
n_mod=length(model);

fig = figure;
set(fig, 'Position', [0, 0, 1150, 1600]);

sea={'Non breeding season', 'Laying season', 'Incubation', 'Rearing season'};
tiledlayout(2,2);
leg = ['a' 'b' 'c' 'd'];

for s=1:4 % season
    nexttile

    for i=[1:3 6:7]
        mod=model(i);

        if i<4 | i>6
            filename=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis.mat', ordi, mod);
        else
            filename=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod);
        end

        load(filename)
        yearstart = info_models(i,3); yearstop = info_models(i,4);
        if i==7
            yearstart = 2006;
        end
        timePOP   = yearstart:yearstop;

        ice_mean=zeros(info_models(i,5),length(timePOP)); %(ens, nt)

        for ens=1:length(ENS)
            SIC=ENS(ens).SICa; % standardized data

            %ice_mean(ens,:)=nanmean(SIC(:,s,:),3); %mean on colonies for each ens for each year
            ice_mean(ens,:)=SIC(:,s,39); %Pointe Géologie

        end

        %ylim([-1,1.5])
        plot(timePOP, nanmean(ice_mean(:,:)), 'linewidth', 1.5) %mean on ens
        hold on

    end %model

    %plot(1980:2011, SICa(:,s,39),'k','linewidth', 1.5) %observed mean on colonies

    title(sea{s}, 'fontsize', 24, 'fontweight', 'bold');
    legend({'CanESM2', 'CESM1', 'CSIRO', 'CESM2', 'Paris 2°C'}, 'Location','southwest','fontsize', 16)
    ylabel('SICa', 'fontsize', 24, 'fontweight', 'bold')
    xlabel('Years', 'fontsize', 24, 'fontweight', 'bold')
    set(gca, 'fontsize', 24, 'fontweight', 'bold')
    text(-0.2, 1.1, leg(s), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

end %season




%% Figure raw data IPM SAT

fig = figure;
set(fig, 'Position', [0, 0, 1150, 800]);

leg = ['a' 'b' 'c' 'd'];

% SAT
% not std just transformed
nexttile
timePOP=1909:2099;
col=39;

load(sprintf('%s/Codes_EP/Codes_SAT/SIC_proj_trans_total.mat', ordi)); %from obtain_SIC_proj_data

plot(timePOP, mean(SIC_proj_total(:,col,:), 3),'linewidth',1.5, 'color', [0, 0.4470, 0.7410])
xlabel('Years', 'FontSize', 24)
ylabel('SIC', 'FontSize', 24)
title("SAT", 'FontSize', 24, 'FontWeight', 'bold');
set(gca, 'FontSize', 24, 'FontWeight', 'bold');
text(-0.2, 1.1, leg(1), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');


%IPM

timePOP=1920:2100;

% Environmental data SIMULATED from Francesco
% 50 ensembles for each year 1979-2100
col = 39

file_name_temp = sprintf('%s/Codes_EP/Codes_IPM/IPM_R/All_colonies/SST_forecast_1920-2100_col_%d.csv', ordi, col)
SSTsim = readmatrix(file_name_temp);
SSTsim(1,:)=[]; SSTsim(:,1)=[]; 

file_name_wind = sprintf('%s/Codes_EP/Codes_IPM/IPM_R/All_colonies/wind_forecast_1920-2100_col_%d.csv', ordi, col)
windsim = readmatrix(file_name_wind);
windsim(1,:)=[]; windsim(:,1)=[];

file_name_vwind = sprintf('%s/Codes_EP/Codes_IPM/IPM_R/All_colonies/vwind_forecast_1920-2100_col_%d.csv', ordi, col)
vwindsim = readmatrix(file_name_vwind);
vwindsim(1,:)=[]; vwindsim(:,1)=[];


envdata={SSTsim, vwindsim, windsim};
names={'SST', 'vWind', 'Wind'}
nexttile
i=1

env=envdata{i};
y2=plot(timePOP, mean(env(:,:)),'linewidth',1.5, 'color', [0, 0.4470, 0.7410]);
hold on
xlabel('Years', 'FontSize', 24)
ylabel(names{i}, 'FontSize', 24)
title("IPM", 'FontSize', 24, 'FontWeight', 'bold');
set(gca, 'FontSize', 24, 'FontWeight', 'bold');
text(-0.2, 1.1, leg(2), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

% Wind
hold on

for i=2:3
    nexttile
    env=envdata{i};
    
    if i==2
        y2=plot(timePOP, mean(env(:,:)), 'linewidth',1.5, 'color', [0, 0.4470, 0.7410]);
    else
        y2=plot(timePOP, mean(env(:,:)), 'linewidth',1.5);
    end
    
    xlabel('Years', 'FontSize', 24)
    ylabel(names{i}, 'FontSize', 24)
    title("IPM", 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 24, 'FontWeight', 'bold');
    hold on
    %legend({'Wind', 'vWind'}, 'Location','southwest','fontsize', 16)
    text(-0.2, 1.1, leg(i+1), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');
end



%% Figure Anomalies
couleur = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];


% Panel a : CMR 4 SEASONS SICa

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
n_mod=length(model);

fig = figure;
set(fig, 'Position', [0, 0, 1150, 800]);
sea={'Non breeding season', 'Laying season', 'Incubation', 'Rearing season'};
tiledlayout(2,2);
leg = ['a' 'b' 'c' 'd'];

nexttile
for s=1:4 % season
    
    for i=6
        mod=model(i);

        if i<4 | i>6
            filename=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std_bis.mat', ordi, mod);
        else
            filename=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod);
        end

        load(filename)
        yearstart = info_models(i,3); yearstop = info_models(i,4);
        if i==7
            yearstart = 2006;
        end
        timePOP   = yearstart:yearstop;

        ice_mean=zeros(info_models(i,5),length(timePOP)); %(ens, nt)

        for ens=1:length(ENS)
            SIC=ENS(ens).SICa; % standardized data

            %ice_mean(ens,:)=nanmean(SIC(:,s,:),3); %mean on colonies for each ens for each year
            ice_mean(ens,:)=SIC(:,s,39); %Pointe Géologie

        end

        %ylim([-1,1.5])
        plot(timePOP, nanmean(ice_mean(:,:)), 'linewidth', 1.5) %mean on ens

        %yline(0, 'k--')
        hold on

    end %model

    %plot(1980:2011, SICa(:,s,39),'k','linewidth', 1.5) %observed mean on colonies

    title('CMR', 'fontsize', 24, 'fontweight', 'bold');
    legend({'Non-breeding', 'Laying', 'Incubation', 'Rearing'}, 'Location','southwest','fontsize', 16)
    ylabel('SICa', 'fontsize', 24, 'fontweight', 'bold')
    xlabel('Years', 'fontsize', 24, 'fontweight', 'bold')
    set(gca, 'fontsize', 24, 'fontweight', 'bold')
    hold on
    text(-0.2, 1.1, leg(1), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');
end %season

%nexttile

% % SAT
% % not std just transformed
% nexttile
% timePOP=1909:2099;
% col=39;
% 
% load(sprintf('%s/Codes_EP/Codes_SAT/SIC_proj_trans_total.mat', ordi)); %from obtain_SIC_proj_data
% 
% plot(timePOP, mean(SIC_proj_total(:,col,:), 3),'linewidth',1.5, 'color', [0, 0.4470, 0.7410])
% xlabel('Years', 'FontSize', 24)
% ylabel('SICa', 'FontSize', 24)
% title("SAT", 'FontSize', 24, 'FontWeight', 'bold');
% set(gca, 'FontSize', 24, 'FontWeight', 'bold');


%IPM

timePOP=1920:2100;

% Environmental data SIMULATED from Francesco
% 50 ensembles for each year 1979-2100
col = 39

file_name_temp = sprintf('%s/Codes_EP/Codes_IPM/IPM_R/All_colonies/env_data_all_col_csv/SST_forecast_1920-2100_col%d.csv', ordi, col)
SSTsim = readmatrix(file_name_temp);
SSTsim(1,:)=[]; SSTsim(:,1)=[]; 

file_name_wind = sprintf('%s/Codes_EP/Codes_IPM/IPM_R/All_colonies/env_data_all_col_csv/wind_forecast_1920-2100_col%d.csv', ordi, col)
windsim = readmatrix(file_name_wind);
windsim(1,:)=[]; windsim(:,1)=[];


file_name_vwind = sprintf('%s/Codes_EP/Codes_IPM/IPM_R/All_colonies/env_data_all_col_csv/vwind_forecast_1920-2100_col%d.csv', ordi, col)
vwindsim = readmatrix(file_name_vwind);
vwindsim(1,:)=[]; vwindsim(:,1)=[];


envdata={SSTsim, vwindsim, windsim};
names={'SST anomalies', 'vWind anomalies', 'Wind anomalies'}
nexttile
i=1

env=envdata{i};
y2=plot(timePOP, mean(env(:,:)),'linewidth',1.5, 'color', [0, 0.4470, 0.7410]);
hold on
xlabel('Years', 'FontSize', 24)
ylabel(names{i}, 'FontSize', 24)
title("IPM", 'FontSize', 24, 'FontWeight', 'bold');
set(gca, 'FontSize', 24, 'FontWeight', 'bold');
text(-0.2, 1.1, leg(2), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

% Wind
hold on


for i=2:3
    nexttile
    env=envdata{i};
    
    if i==2
        y2=plot(timePOP, mean(env(:,:)), 'linewidth',1.5, 'color', [0, 0.4470, 0.7410]);
    else
        y2=plot(timePOP, mean(env(:,:)), 'linewidth',1.5);
    end
    
    xlabel('Years', 'FontSize', 24)
    ylabel(names{i}, 'FontSize', 24)
    title("IPM", 'FontSize', 24, 'FontWeight', 'bold');
    set(gca, 'FontSize', 24, 'FontWeight', 'bold');
    hold on
    %legend({'Wind', 'vWind'}, 'Location','southwest','fontsize', 16)
    text(-0.2, 1.1, leg(i+1), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');
end

