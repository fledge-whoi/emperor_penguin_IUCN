%% Figure to plot the different ecological models

% Figures for : 
% ecological models + eco-ensemble from 1900, 1950 or 2009
% internal variability of ecological models


% We plot the mean of the demographic simulations after doing the mean of
% the ensemble for each simulation

% DIRECTORY to Codes_EP
ordi="/Users/aliceeparvier/Desktop"

% Load info models
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi))

% Load observed data
dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));

% Climate model, scenario SSP370
mod = "LE_CESM2";

% Ecological models
eco_mod = ["CMR", "IPM", "Sat"];

% Model parameters
years_start = [1900 1921 1909];
years_end = [2100 2100 2100];
nens=50;
nc=66;

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

% Extreme event scenario for CMR model = 1, 2 and 4 gathered 
scenEXT = [1, 2, 4];


%% 3 models 1950-2100 and panel D: relative change - FIGURE 1

% Climate model, scenario SSP370
mod = "LE_CESM2";

% Ecological models
eco_mod = ["CMR", "IPM", "Sat"];

% Model parameters
years_start = [1900 1921 1909];
years_end = [2100 2100 2100];
nens=50;
nc=66;

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

% Extreme event scenario for CMR model = 1, 2 and 4 gathered 
scenEXT = [1, 2, 4];


leg = ['a' 'b' 'c'];

% Create a figure with four subplots
fig = figure;
set(fig, 'Position', [80, 0, 1300, 1100]);


y_tot=zeros(4,1); %for the legend
for m = [1:3] % Ecological model
    subplot(2,2,m)

    e_mod=eco_mod(m);

    % Model parameters
    yearstart = years_start(m);
    yearstop = years_end(m);
    timePOP   = yearstart:yearstop;
    
    nsim = 100;
    
    y1950=length(yearstart:1950);
    timePOP   = 1950:yearstop;
    nt        = length(timePOP);

    if m==1 %CMR
        Ntot = zeros(nt, 3*nsim*nens);
    else
        Ntot = zeros(nt, nsim*nens);
    end


    x=1;
    for ens = 1:nens

        if m == 1 %CMR : We select nsim simulations among 3*nsim (3 scenEXT)

            scen=1;
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
            N = readmatrix(file_name); % N = (nt, nsim)
            N=N(y1950:end, :);
            N_old=N;

            for scen=[2 4] % Extreme event scenario

                file_name = sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
                N = readmatrix(file_name); % N = (nt, nsim)
                N=N(y1950:end, :);
                N_old = cat(2, N_old, N);
               
            end %scenEXT

            %N = N_old(:, randi(3*nsim, nsim,1));
            % We have N=(nt, nsim) with nsim random mix of 3 scenEXT
            N=N_old;

        elseif m==2
            % for IPM, BE=1.3 BE and Kmax=4.5
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
            N=N(y1950:end, :);
        else %Sat
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/LE_CESM2/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, e_mod, mod, ens));
            N=N(y1950:end, :);
        end

        if m==1
            Ntot(:,x:x+3*nsim-1) = N;
            x=x+3*nsim;
        else %IPM
            Ntot(:, x:x+nsim-1) = N;
            x=x+nsim;
        end

    end %ens

    Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % PLOT

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot 3 simulations
    plot(timePOP, Ntot(randi(nens*nsim, 3,1),:),'linewidth',1)

    hold on

    % Plot median and 95% quantiles
    y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur(m,:),'linewidth',5);
    y_tot(m)=y;
    hold on
    plot(timePOP,Ntot_med(1,:),'-','color', couleur(m,:),'linewidth',2)
    hold on
    plot(timePOP, Ntot_med(3,:),'-','color', couleur(m,:),'linewidth',2)

    %   y_tot(3) = plot(2009:2018, dataOBSSAT_bilgecan(:,2), 'k','linewidth',5);
    %     plot(2009:2018, dataOBSSAT_bilgecan(:,4), 'k:','linewidth',2)
    %     plot(2009:2018, dataOBSSAT_bilgecan(:,3), 'k:','linewidth',2)

    hold on

    fontsi = 20;

    % y_tot(4) = plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylabel('Nb of mature individuals', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 500000]);
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ax=gca;
    ax.YColor = [0 0 0];
    text(-0.2, 1.1, leg(m), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');


end % ecological model


%% Panel D
% Figure all ecological models

subplot(2,2,4)
yearstart=1950;
yearstop=2100;
timePOP=1950:2100;
nt=length(timePOP);
y1950=length(1921:1950);
% Import eco-ensemble results

%ens_results=readmatrix(sprintf('%s/Figures/Fig_ecomod_mature/ensamble_forecasts_slopeweights_v2.csv', ordi));
ens_med = permute(ens_results(y1950:end,3), [2 1]);
ens_inf = permute(ens_results(y1950:end,4), [2 1]);
ens_sup = permute(ens_results(y1950:end,5), [2 1]);

% Plot mean, 95% interval, internal variability
couleur = [0.3 0.3 0.6];
ttt = [timePOP, fliplr(timePOP)];

inBetween = [ens_inf, fliplr(ens_sup)];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.1)
hold on

% Plot median and 95% quantiles
y=plot(timePOP, ens_med,'-', 'color', couleur,'linewidth',7);
y_tot(m)=y;
hold on
plot(timePOP, ens_inf,'-','color', couleur,'linewidth',2)
hold on
plot(timePOP, ens_sup,'-','color', couleur,'linewidth',2)

obs=(dataOBSSAT_bilgecan(:,2)-mean(dataOBSSAT_bilgecan(:,2)))./mean(dataOBSSAT_bilgecan(:,2)) *100;
yobs = plot(2009:2018, obs,'-k','linewidth',7);
hold on

fontsi = 20;

xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
ylabel('% change', 'FontSize', fontsi, 'FontWeight', 'bold');
%ylim([0 300000]);
title('Eco-ensemble');

set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
text(-0.2, 1.1, 'd', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');



%% Three models separated + panel d: percent change of three models gathered, 1900-2100

leg = ['a' 'b' 'c'];

% Create a figure with four subplots
fig = figure;
set(fig, 'Position', [80, 0, 1300, 1100]);


y_tot=zeros(4,1); %for the legend
for m = [1:3] % Ecological model
    subplot(2,2,m)

    e_mod=eco_mod(m);

    % Model parameters
    yearstart = years_start(m);
    yearstop = years_end(m);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nsim = 100;

    if m==1 %CMR
        Ntot = zeros(nt, 3*nsim*nens);
    else
        Ntot = zeros(nt, nsim*nens);
    end


    x=1;
    for ens = 1:nens

        if m == 1 %CMR : We select nsim simulations among 3*nsim (3 scenEXT)

            scen=1;
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
            N = readmatrix(file_name); % N = (nt, nsim)
            N_old=N;

            for scen=[2 4] % Extreme event scenario

                file_name = sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
                N = readmatrix(file_name); % N = (nt, nsim)
                N_old = cat(2, N_old, N);
               
            end %scenEXT

            %N = N_old(:, randi(3*nsim, nsim,1));
            % We have N=(nt, nsim) with nsim random mix of 3 scenEXT
            N=N_old;

        elseif m==2
            % for IPM, BE=1.3 BE and Kmax=4.5
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
        else %Sat
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/LE_CESM2/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, e_mod, mod, ens));

        end

        if m==1
            Ntot(:,x:x+3*nsim-1) = N;
            x=x+3*nsim;
        else %IPM
            Ntot(:, x:x+nsim-1) = N;
            x=x+nsim;
        end

    end %ens

    Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % PLOT

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot the ensembles (one ens = mean of 100 simulations)
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median and 95% quantiles
    y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur(m,:),'linewidth',5);
    y_tot(m)=y;
    hold on
    plot(timePOP,Ntot_med(1,:),'-','color', couleur(m,:),'linewidth',2)
    hold on
    plot(timePOP, Ntot_med(3,:),'-','color', couleur(m,:),'linewidth',2)

    %   y_tot(3) = plot(2009:2018, dataOBSSAT_bilgecan(:,2), 'k','linewidth',5);
    %     plot(2009:2018, dataOBSSAT_bilgecan(:,4), 'k:','linewidth',2)
    %     plot(2009:2018, dataOBSSAT_bilgecan(:,3), 'k:','linewidth',2)

    hold on

    fontsi = 20;

    % y_tot(4) = plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylabel('Nb of mature individuals', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 500000]);
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ax=gca;
    ax.YColor = [0 0 0];
    text(-0.2, 1.1, leg(m), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

end % ecological model


% Panel D
% Figure all ecological models

subplot(2,2,4)
yearstart=1921;
yearstop=2100;
timePOP=1921:2100;
nt=length(timePOP);
% Import eco-ensemble results

%ens_results=readmatrix(sprintf('%s/Figures/Fig_ecomod_mature/ensamble_forecasts_slopeweights_v2.csv', ordi));
ens_med = permute(ens_results(:,3), [2 1]);
ens_inf = permute(ens_results(:,4), [2 1]);
ens_sup = permute(ens_results(:,5), [2 1]);

% Plot mean, 95% interval, internal variability
couleur = [0.3 0.3 0.6];
ttt = [timePOP, fliplr(timePOP)];

inBetween = [ens_inf, fliplr(ens_sup)];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.1)
hold on

% Plot median and 95% quantiles
y=plot(timePOP, ens_med,'-', 'color', couleur,'linewidth',7);
y_tot(m)=y;
hold on
plot(timePOP, ens_inf,'-','color', couleur,'linewidth',2)
hold on
plot(timePOP, ens_sup,'-','color', couleur,'linewidth',2)

obs=(dataOBSSAT_bilgecan(:,2)-mean(dataOBSSAT_bilgecan(:,2)))./mean(dataOBSSAT_bilgecan(:,2)) *100;
yobs = plot(2009:2018, obs,'-k','linewidth',7);
hold on

fontsi = 20;

xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
ylabel('% change', 'FontSize', fontsi, 'FontWeight', 'bold');
%ylim([0 300000]);
title('Eco-ensemble');

set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
text(-0.2, 1.1, 'd', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

%% Fig supp Hindcast ecomod relative change (1950)

% Climate model, scenario SSP370
mod = "LE_CESM2";

% Ecological models
eco_mod = ["CMR", "IPM", "Sat"];

% Model parameters
years_start = [1900 1921 1909];
years_end = [2100 2100 2100];
nens=50;
nc=66;

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

% Extreme event scenario for CMR model = 1, 2 and 4 gathered 
scenEXT = [1, 2, 4];


leg = ['a' 'b' 'c'];

% Create a figure with four subplots
fig = figure;
set(fig, 'Position', [80, 0, 1300, 1100]);


y_tot=zeros(4,1); %for the legend
for m = [1:3] % Ecological model
    subplot(2,2,m)

    e_mod=eco_mod(m);

    % Model parameters
    yearstart = years_start(m);
    yearstop = years_end(m);
    timePOP   = yearstart:yearstop;

    y2009 = length(1950:2009);
    
    nsim = 100;
    
    y1950=length(yearstart:1950);
    timePOP   = 1950:yearstop;
    nt        = length(timePOP);

    if m==1 %CMR
        Ntot = zeros(nt, 3*nsim*nens);
        Xtot = zeros(nt, 3*nens*nsim);
    else
        Ntot = zeros(nt, nsim*nens);
        Xtot = zeros(nt, nens*nsim);
    end


    x=1;
    for ens = 1:nens

        if m == 1 %CMR : We select nsim simulations among 3*nsim (3 scenEXT)

            scen=1;
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
            N = readmatrix(file_name); % N = (nt, nsim)
            N=N(y1950:end, :);
            N_old=N;

            for scen=[2 4] % Extreme event scenario

                file_name = sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
                N = readmatrix(file_name); % N = (nt, nsim)
                N=N(y1950:end, :);
                N_old = cat(2, N_old, N);
               
            end %scenEXT

            %N = N_old(:, randi(3*nsim, nsim,1));
            % We have N=(nt, nsim) with nsim random mix of 3 scenEXT
            N=N_old;

        elseif m==2
            % for IPM, BE=1.3 BE and Kmax=4.5
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
            N=N(y1950:end, :);
        else %Sat
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/LE_CESM2/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, e_mod, mod, ens));
            N=N(y1950:end, :);
        end

        if m==1
            Ntot(:,x:x+3*nsim-1) = N;
            x=x+3*nsim;
        else %IPM
            Ntot(:, x:x+nsim-1) = N;
            x=x+nsim;
        end

    end %ens


    Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    Xtot = (Ntot - mean(Ntot(:, y2009:y2009+10),2))./mean(Ntot(:, y2009:y2009+10),2);
    Xtot_med = quantile(Xtot, [0.025 0.5 0.975]);

    % PLOT

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Xtot_med(1,:), fliplr(Xtot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot the ensembles (one ens = mean of 100 simulations)
    plot(timePOP, Xtot(randi(length(Xtot(:,1)),10,1),:),'linewidth',1)

    hold on

    % Plot median and 95% quantiles
    y=plot(timePOP, Xtot_med(2,:),'-', 'color', couleur(m,:),'linewidth',5);
    y_tot(m)=y;
    hold on
    plot(timePOP,Xtot_med(1,:),'-','color', couleur(m,:),'linewidth',2)
    hold on
    plot(timePOP, Xtot_med(3,:),'-','color', couleur(m,:),'linewidth',2)

    %   y_tot(3) = plot(2009:2018, dataOBSSAT_bilgecan(:,2), 'k','linewidth',5);
    %     plot(2009:2018, dataOBSSAT_bilgecan(:,4), 'k:','linewidth',2)
    %     plot(2009:2018, dataOBSSAT_bilgecan(:,3), 'k:','linewidth',2)

    hold on

    fontsi = 20;

    obs=(dataOBSSAT_bilgecan(:,2)-mean(dataOBSSAT_bilgecan(:,2)))./mean(dataOBSSAT_bilgecan(:,2));

    y_tot(4) = plot(2009:2018, obs,'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylabel('% change', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([-1 1]);
    xlim([1950 2100]);

    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ax=gca;
    ax.YColor = [0 0 0];
    text(-0.2, 1.1, leg(m), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

   
end % ecological model



% Panel D
subplot(2,2,4)
yearstart=1950;
yearstop=2100;
timePOP=1950:2100;
nt=length(timePOP);
y1950=length(1921:1950);
% Import eco-ensemble results

%ens_results=readmatrix(sprintf('%s/Figures/Fig_ecomod_mature/ensamble_forecasts_slopeweights_v2.csv', ordi));
ens_med = permute(ens_results(y1950:end,3), [2 1]);
ens_inf = permute(ens_results(y1950:end,4), [2 1]);
ens_sup = permute(ens_results(y1950:end,5), [2 1]);

% Plot mean, 95% interval, internal variability
couleur = [0.3 0.3 0.6];
ttt = [timePOP, fliplr(timePOP)];

inBetween = [ens_inf, fliplr(ens_sup)];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.1)
hold on

% Plot median and 95% quantiles
y=plot(timePOP, ens_med,'-', 'color', couleur,'linewidth',7);
y_tot(m)=y;
hold on
plot(timePOP, ens_inf,'-','color', couleur,'linewidth',2)
hold on
plot(timePOP, ens_sup,'-','color', couleur,'linewidth',2)

obs=(dataOBSSAT_bilgecan(:,2)-mean(dataOBSSAT_bilgecan(:,2)))./mean(dataOBSSAT_bilgecan(:,2)) *100;
yobs = plot(2009:2018, obs,'-k','linewidth',7);
hold on

fontsi = 20;

xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
ylabel('% change', 'FontSize', fontsi, 'FontWeight', 'bold');
%ylim([0 300000]);
title('Eco-ensemble');

set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
text(-0.2, 1.1, 'd', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');





%% Figure from 2009

years_start = [2009 2009 2009];
years_end = [2100 2100 2099];

fig = figure;
set(fig, 'Position', [100, 100, 1600, 1100]);

y_tot=zeros(4,1); %for the legend


for m = [1:3] % Ecological model
    subplot(2,2,m)
    e_mod=eco_mod(m);

    % Model parameters
    yearstart = years_start(m);
    yearstop = years_end(m);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);

    if m==1 %CMR has 3 scenario
        N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, 1, 1));
        nsim = length(N(1,:));
    elseif m==2 % IPM

        N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, 1));
        nsim = length(N(1,:));

    elseif m==3 %Sat, no metapop model
        N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_LE_CESM2_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, 1));
        nsim = length(N(1,:));

    end

    if m==1 %CMR
        Ntot = zeros(nt, 3*nsim*nens);
    else
        Ntot = zeros(nt, nsim*nens);
    end


    x=1;
    for ens = 1:nens

        if m == 1 %CMR : We select nsim simulations among 3*nsim (3 scenEXT)

            scen=1;
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
            N = readmatrix(file_name); % N = (nt, nsim)
            N_old=N;

            for scen=[2 4] % Extreme event scenario

                file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens, scen);
                N = readmatrix(file_name); % N = (nt, nsim)
                N_old = cat(2, N_old, N);
               
            end %scenEXT

            %N = N_old(:, randi(3*nsim, nsim,1));
            % We have N=(nt, nsim) with nsim random mix of 3 scenEXT
            N=N_old;

        elseif m==2
            % for IPM, BE=1.3 BE and Kmax=4.5
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
        else %Sat
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
             
        end         

        if m==1
            Ntot(:,x:x+3*nsim-1) = N; 
            x=x+3*nsim;
        else %IPM
            Ntot(:, x:x+nsim-1) = N; 
            x=x+nsim;
        end

    end %ens

    Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % PLOT

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot the ensembles (one ens = mean of 100 simulations)
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median and 95% quantiles
    y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur(m,:),'linewidth',5);
    y_tot(m)=y;
    hold on
    plot(timePOP,Ntot_med(1,:),'-','color', couleur(m,:),'linewidth',2)
    hold on
    plot(timePOP, Ntot_med(3,:),'-','color', couleur(m,:),'linewidth',2)

    hold on

    fontsi = 24;

    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    xlim([2000 2100])
    ylabel('Nb of mature individuals', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 400000]);
    if m==1
        title('SPCMR');
        %print('Fig_CMR_mature', '-dpng', '-r300')
    elseif m==2
        title('SPIPM');
        %print('Fig_IPM_mature', '-dpng', '-r300')
    else
        title('SCDSAT');
        %print('Fig_SAT_mature', '-dpng', '-r300')
    end
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ax=gca;
    ax.YColor = [0 0 0];

end % eco model




% Figure all ecological models
subplot(2,2,4)
nens=50;
nsim=100;
yearstart=2009;
yearstop=2099;
timePOP=2009:2099;
nt=length(timePOP);

% SIMULATIONS

Ntot_3gathered = zeros(3*nens*nsim, nt);
x2=1;

    for m=[1 2 3]

        e_mod=eco_mod(m)
        if m==1 %CMR
            time=1:91; %subset years
            scenEXT=[1 2 4];
        elseif m==2 %IPM
            time=1:91;
            %nt=180;
            scenEXT=1;
        else %SAT
            time=1:91;
            %nt=191;
            scenEXT=1;
        end

        Ncol = zeros(nt, length(scenEXT)*nens*nsim);

        x=1;
        
        for scen=scenEXT % Extreme event scenario
            
            for ens = 1:nens

                if m==1
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_CMR/LE_CESM2/Ntot_%s_LE_CESM2_ens%d_scen%d.txt", ordi, e_mod, e_mod, ens, scen);
                elseif m==2
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_IPM/LE_CESM2/Ntot_%s_LE_CESM2_ens%d.txt", ordi, e_mod, e_mod, ens);
                else %Sat
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_Sat/LE_CESM2/Ntot_%s_LE_CESM2_ens%d.txt", ordi, e_mod, e_mod, ens);
                end
                N = readmatrix(file_name); %(nt, nsim)

                Ncol(:,x:x+nsim-1) = Ncol(:,x:x+nsim-1) + N(time,:);
                x=x+nsim;
            end %ens

        end %scenEXT

        if m==1
            Ncol = Ncol(:, randi(3*nens*nsim, nens*nsim, 1));
        end

        Ncol = permute(Ncol, [2 1]); %(sim*ens, nt)

        Ntot_3gathered(x2:(x2+nens*nsim-1), :) = Ncol;
        x2=x2+(nens*nsim);

    end %eco_mod

Ntot_med=quantile(Ntot_3gathered, [0.025 0.5 0.975]);


% PLOT


couleur = [0.3 0.3 0.6];

% Plot mean, 95% interval, internal variability

ttt = [timePOP, fliplr(timePOP)];

inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.1)
hold on

% Plot the ensembles (one ens = mean of 100 simulations)
%plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

hold on

% Plot median and 95% quantiles
y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur,'linewidth',7);
y_tot(m)=y;
hold on
plot(timePOP,Ntot_med(1,:),'-','color', couleur,'linewidth',2)
hold on
plot(timePOP, Ntot_med(3,:),'-','color', couleur,'linewidth',2)

yobs = plot(2009:2018, dataOBSSAT_bilgecan(:,2), 'k','linewidth',7);
hold on

fontsi = 20;

xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
ylabel('Nb of mature individuals', 'FontSize', fontsi, 'FontWeight', 'bold');
ylim([0 400000]);
title('Three models combined');

set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
text(-0.2, 1.1, 'd', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');
%print('Fig_ecomod_gathered', '-dpng', '-r300')

legend(yobs, 'Obs')


%% Figure internal variability of ecological models

leg = ['d' 'a' 'b' 'c'];

% Create a figure with four subplots
fig = figure;
set(fig, 'Position', [80, 0, 1300, 500]);

tiledlayout(1, 2);

%load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi))
%dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));

% Climate model, scenario RCP8.5
mod = "LE_CESM2";

% Ecological models
eco_mod = ["CMR", "IPM", "Sat"];


choix = 2;
if choix==1
    years_start = [1900 1980 1909];
else
    years_start = [2009 2009 2009];
end

% Model parameters
years_end = [2100 2100 2100];
nens=50;
nc=66;
scenEXT = [1, 2, 4];

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];


y_tot=zeros(4,1); %for the legend

for m = [2 3] % Ecological model
    nexttile

    e_mod=eco_mod(m);

    % Model parameters
    yearstart = years_start(m);
    yearstop = years_end(m);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nsim = 100;

    if m==1 %CMR
        Ntot = zeros(nt, 3*nsim*nens);
    else
        Ntot = zeros(nt, nsim*nens);
    end


    x=1;
    for ens = 1:nens

        if m==2
            % for IPM, BE=1.3 BE and Kmax=4.5
            if choix==1
                N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
            else
                N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
            end
        else %Sat
            if choix ==1
                N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/no_disp/Ntot_%s_LE_CESM2_ens%d.txt", ordi, e_mod, e_mod, e_mod, ens));
            else
                N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/LE_CESM2/Ntot_%s_LE_CESM2_ens%d.txt", ordi, e_mod, e_mod, e_mod, ens));
            end
        end

        if m==1
            Ntot(:,x:x+3*nsim-1) = N;
            x=x+3*nsim;
        else %IPM
            Ntot(:, x:x+nsim-1) = N;
            x=x+nsim;
        end

    end %ens

    Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % ADD INTERNAL VARIABILITY
   
    if m==2
        Ntot2 = zeros(nt, nsim);
        if choix==1 %1909
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/mean_ensemble/N_results_%s/Ntot_%s/%s/Ntot_%s_%s.txt", ordi, e_mod, e_mod, mod, e_mod, mod);
        else %2009
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_%s/Ntot_%s/%s/Ntot_%s_%s.txt", ordi, e_mod, e_mod, mod, e_mod, mod);
        end
        N2 = readmatrix(file_name); %(nt, nsim)
        Ntot2(:,:) = N2(:,:);
    else m==3
        Ntot2 = zeros(nt, nsim);
        if choix==1 %1909
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/mean_ensemble/N_results_%s/Ntot_%s/%s/Ntot_%s_%s.txt", ordi, e_mod, e_mod, mod, e_mod, mod);
        else %2009
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_%s/Ntot_%s/%s/Ntot_%s_%s.txt", ordi, e_mod, e_mod, mod, e_mod, mod);
        end
        N2 = readmatrix(file_name); %(nt, nsim)
        Ntot2(:,:) = N2(:,:);
    end

    Ntot2 = permute(Ntot2, [2 1]); %(sim, nt)
    Ntot_med2 = quantile(Ntot2, [0.025 0.5 0.975]);


    % PLOT

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot median and 95% quantiles
    plot(timePOP, Ntot_med(2,:),'-', 'color', couleur(m,:),'linewidth',6);
    y_tot(m)=y;
    hold on
    plot(timePOP,Ntot_med(1,:),'-','color', couleur(m,:),'linewidth',2)
    hold on
    plot(timePOP, Ntot_med(3,:),'-','color', couleur(m,:),'linewidth',2)

    hold on

    % Plot mean of ensembles and quantiles
    y=plot(timePOP,Ntot_med2(2,:),'color', [0.6 0.6 0.6],'linewidth',4);
    hold on
    plot(timePOP,Ntot_med2(1,:),'color', [0.6 0.6 0.6],'linewidth',2);
    hold on
    plot(timePOP,Ntot_med2(3,:),'color', [0.6 0.6 0.6],'linewidth',2);


    fontsi = 24;

    %y_tot(4) = plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    text(-0.15, 1.08, leg(m), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

    ylabel('Nb of mature individuals', 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 350000]);
    xlim([2009 2100])
    titre = {'CMR', 'SPIPM', 'SCDSAT'};
    title(titre{m});

    legend(y, "mean of ensembles");


end % climate model

