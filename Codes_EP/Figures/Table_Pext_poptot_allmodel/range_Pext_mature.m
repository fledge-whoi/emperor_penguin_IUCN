%% CALCULATE ALL DECLINE PROBABILITIES FOR ALL THE MODELS

% Plot figures with proba ext, at the end of the code we create a table
% with all Pext calculated (Ecological models, climate models, climate
% scenarios, extreme event scenarios, internal variability)

ordi="/Users/aliceeparvier/Desktop";

% Load info models
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi))

% Load observed data
dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));

time_present= 2024; %2024
GL = 16.4; % Generation length = average age of parents in the current cohort
year_eval = round(time_present + 3*GL);


%% ECOLOGICAL MODELS - 3 separate + gathered
% Figure and Pext_eco_mod table

fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);
tiledlayout(2,2);

% Climate model, scenario SSP370
mod = "LE_CESM2";

% Ecological models
eco_mod = ["CMR", "IPM", "Sat"];

% Model parameters
years_start = [2009 2009 2009];
years_end = [2100 2100 2100];
nens=50;
nc=66;

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.7 0.7 0.7];

% Extreme event scenario for CMR model = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

y_tot=zeros(4,1); %for the legend

Pext_eco_mod = zeros(4,4); %3 models + gathered

for m = [1:3] % Ecological model
    nexttile

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

    % Parameters for the extinction probability
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %mean 2009-2018
    y_eval = year_eval - yearstart;

    %pop_present = mean(dataOBSSAT_bilgecan(:,2)); %predicted current population

    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for i=1:length(Ntot(:,1)); %ensemble
            if Ntot(i, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_eco_mod(m,:) = proba_tot;

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

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)



    fontsi = 24;
    y_tot(4) = plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    xlim([2000 2100])
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 350000])
    if m==3 %SAT
        title('SAT');
        ylabel('Adults in attendance', 'FontSize', fontsi, 'FontWeight', 'bold');
    elseif m==1
        title('CMR')
        ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
    else
        title('IPM')
        ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
    end

end % climate model

% Ecological models gathered

nexttile

% Model parameters
nens=50;
ncol=66;
yearstart = 2009;
yearstop = 2100;
timePOP   = yearstart:yearstop;
nt        = length(timePOP);
nsim=100;
Ntot = zeros(nt, (3*nsim*nens));
% We want as much CMR simulations as IPM and SAT so we select nsim*nens
% simulations (among the three scenEXT)

x=1;
for m = [1:3] % Ecological model

    e_mod=eco_mod(m);

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

            N=N_old(1:nt, randi(3*nsim, nsim, 1)); %we select 1/3 of simulations

        elseif m==2
            % for IPM, BE=1.3 BE and Kmax=4.5
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
        else %Sat
            N=readmatrix(sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, e_mod, e_mod, mod, e_mod, mod, ens));
        end

        if m==1
            Ntot(:,x:x+nsim-1) = N;
            x=x+(nsim);
        else %IPM Sat
            Ntot(:,x:x+nsim-1) = N;
            x=x+nsim;
        end

    end %ens

end

Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

% CALCULATE PEXT
y_eval = year_eval - yearstart;
pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)
threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
threshold_endan = 0.5 * pop_present; % Endangered
threshold_critend = 0.2 * pop_present; % Critically endangered
thr_ext = 0.05 * pop_present;

% Probability to go under the threshold
proba_tot=[0 0 0 0];
j=1;
for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
    proba = 0;
    under_thr = 0;
    for i=1:length(Ntot(:,1)); %total simulations
        if Ntot(i, y_eval) < thr;
            under_thr = under_thr + 1;
            proba = under_thr / length(Ntot(:,1));
        end
    end
    proba_tot(j)=proba;
    j=j+1;
end

Pext_eco_mod(4,:) = proba_tot;

% PLOT

% Plot mean, 95% interval, internal variability
couleur = [0, 0.4470, 0.7410];
ttt = [timePOP, fliplr(timePOP)];

inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.05)
hold on

% Plot 5 simulations
plot(timePOP, Ntot(randi(15000, 5,1 ),:),'color', [0.6 0.6 0.6],'linewidth',0.6)

hold on
% Plot median and 95% quantiles
y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur,'linewidth',5);
y_tot(m)=y;
hold on
plot(timePOP,Ntot_med(1,:),'-','color', couleur,'linewidth',2)
hold on
plot(timePOP, Ntot_med(3,:),'-','color', couleur,'linewidth',2)

hold on

fontsi = 24;
plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);


% Plot Pext and thr
coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl1.Color=coul_P{1};
yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl2.Color=coul_P{2};
yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl3.Color=coul_P{3};
yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl4.Color=coul_P{4};
yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
xlim([2000 2100])
ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
ylim([0 350000]);
title('All ecological models', 'FontSize', fontsi, 'FontWeight', 'bold');
set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');


%save as high resolution
%print('Fig_range_ecomod_mature', '-dpng', '-r300')



%% CLIMATE SCENARIOS - CESM-LE

fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);
tiledlayout(2,2)

% Climate scenario
% RCP8.5, SSP370, Paris Agreement 2degrees
model = ["LE_CESM1-CAM5", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
model_name = ["RCP8.5" "SSP370" "Paris Agreement 2 degrees"];
indice_mod = [2 6 7];


% Ecological model = Stef's model, CMR
eco_mod = "CMR";

% one color per model + one color for ensembles (grey)
couleur = [0.8500 0.3250 0.0980;
    0.4940 0.1840 0.5560;
    0 0.4470 0.7410;
    0.7 0.7 0.7];

% Extreme event scenario = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

% TOTAL POPULATION - 3 CLIMATE SCENARIO + SCENEXT

fontsi=24;

y_tot = zeros(3,1); %for legend

Pext_clim_scen = zeros(4,4);

for m=1:length(model) %Climate scenario
    nexttile
   
    mod=model(m);

    % Model parameters
    i_mod=indice_mod(m);
    yearstart = 2009; yearstop = info_models(i_mod,4);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nc = 66;
    nens = info_models(i_mod,5);
    nsim = 100;

    Ntot = zeros(nt, nens*nsim*3); % To plot all scenEXT, all sim all ens

    x=1;
    for scen=scenEXT % Extreme event scenario
        for ens = 1:nens

            % load Ntot (nt, nsim) for each scenEXT, ens
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            N = readmatrix(file_name); %(nt, nsim)

            Ntot(:,x:x+nsim-1) = N(:,:);
            x=x+nsim;
        end %ens

    end %scenEXT

    Ntot = permute(Ntot, [2 1]); %(ens, nt)

    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % CALCULATE PEXT
    y_eval = year_eval - yearstart;
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)
    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for i=1:length(Ntot(:,1)); %total simulations
            if Ntot(i, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_clim_scen(m, :) = proba_tot;

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)

    % Plot the 50 ens
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median
    y= plot(timePOP,Ntot_med(2,:),'color',couleur(m,:),'linewidth',5);
    y_tot(m)=y;

    % Plot observations
    plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

    hold on

    xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'bold');
    ylabel('Abundance', 'FontSize',fontsi, 'Fontweight', 'bold');
    ylim([0 300000]);
    set(gca, 'FontSize', fontsi, 'Fontweight', 'bold');
    if m==1
        title('RCP8.5', 'Fontsize', fontsi, 'Fontweight', 'bold');
    elseif m==2
        title('SSP370', 'Fontsize', fontsi, 'Fontweight', 'bold');
    elseif m==3
        title('Paris 2 degrees', 'Fontsize', fontsi, 'Fontweight', 'bold');
    end

end % climate scenario

% THREE SCENARIOS GATHERED

nexttile

model = ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
eco_mod='CMR';
yearstart=2009;
yearstop=2099;
timePOP=2009:2099;
nt=length(timePOP);

% nens_tot= 50+40+30+20+100+50+11;
nens_tot= 40+50+11; %CESM1, CEMS2, PARIS

N_all = zeros(nens_tot*nsim*3, nt);

x2=1;

for m=[2 6 7]
    mod=model(m);
    nens=info_models(m,5);

    Ntot = zeros(nt, 3*nens*nsim);
    x=1;
    for scen=scenEXT % Extreme event scenario

        for ens = 1:nens

            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            N = readmatrix(file_name); %(nt, nsim)

            if m==5
                Ntot(:,x:x+nsim-1) = N(:,:);
            else
                Ntot(:,x:x+nsim-1) = N(1:end-1,:);
            end

            x=x+nsim;
        end %ens

    end %scenEXT

    Ntot = permute(Ntot, [2 1]); %(sim*ens, nt)

    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    N_all(x2:x2+(3*nens*nsim)-1, :) = Ntot;

    x2=x2+(3*nens*nsim);

end %model

Ntot_med = quantile(N_all, [0.025 0.5 0.975]);
Ntot = N_all;

% CALCULATE PEXT
    y_eval = year_eval - yearstart;
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)
    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for i=1:length(Ntot(:,1)); %total simulations
            if Ntot(i, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_clim_scen(4,:) = proba_tot;

    % Plot mean, 95% interval, internal variability

    couleur = [.46 .7 .18];

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur);
    set(f,'EdgeColor','none','FaceAlpha', 0.1)

    % Plot the 50 ens
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median
    y= plot(timePOP,Ntot_med(2,:),'color',couleur,'linewidth',5);

    % Plot observations
    plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

    hold on

    xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'bold');
    ylabel('Abundance', 'FontSize',fontsi, 'Fontweight', 'bold');
    ylim([0 300000]);
    title('3 Scenarios gathered', 'Fontsize', fontsi, 'Fontweight', 'bold');
    set(gca, 'FontSize', fontsi, 'Fontweight', 'bold');

%print('Fig_range_clim_scen_mature', '-dpng', 'r300')


%% CLIMATE MODELS - RCP8.5

fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);
tiledlayout(2,2)

% Climate scenario
% RCP8.5, SSP370, Paris Agreement 2degrees
model = ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0"];
model_name = ["CanESM2 RCP8.5" "CESM1-CAM5 RCP8.5" "CSIRO RCP8.5"];
indice_mod = [1:3];


% Ecological model = Stef's model, CMR
eco_mod = "CMR";

couleur = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880;
    0.8 0.8 0.8];

% Extreme event scenario = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

% TOTAL POPULATION - 3 CLIMATE SCENARIO + SCENEXT

fontsi=24;

y_tot = zeros(3,1); %for legend

Pext_clim_mod = zeros(4,4);

for m=1:length(model) %Climate scenario
    nexttile
   
    mod=model(m);

    % Model parameters
    i_mod=indice_mod(m);
    yearstart = 2009; yearstop = info_models(i_mod,4);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nc = 66;
    nens = info_models(i_mod,5);
    nsim = 100;

    Ntot = zeros(nt, nens*nsim*3); % To plot all scenEXT, all sim all ens

    x=1;
    for scen=scenEXT % Extreme event scenario
        for ens = 1:nens

            % load Ntot (nt, nsim) for each scenEXT, ens
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            N = readmatrix(file_name); %(nt, nsim)

            Ntot(:,x:x+nsim-1) = N(:,:);
            x=x+nsim;
        end %ens

    end %scenEXT

    Ntot = permute(Ntot, [2 1]); %(ens, nt)

    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % CALCULATE PEXT
    y_eval = year_eval - yearstart;
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)
    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for i=1:length(Ntot(:,1)); %total simulations
            if Ntot(i, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_clim_mod(m, :) = proba_tot;

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)

    % Plot the 50 ens
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median
    y= plot(timePOP,Ntot_med(2,:),'color',couleur(m,:),'linewidth',5);
    y_tot(m)=y;

    % Plot observations
    plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

    hold on

    xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'bold');
    ylabel('Abundance', 'FontSize',fontsi, 'Fontweight', 'bold');
    ylim([0 300000]);
    set(gca, 'FontSize', fontsi, 'Fontweight', 'bold');

    title(model_name(m), 'Fontsize', fontsi, 'Fontweight', 'bold');

end % climate scenario

% THREE MODELS GATHERED

nexttile

model = ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
eco_mod='CMR';
yearstart=2009;
yearstop=2099;
timePOP=2009:2099;
nt=length(timePOP);

% nens_tot= 50+40+30+20+100+50+11;
nens_tot= 50+40+30; 

N_all = zeros(nens_tot*nsim*3, nt);
x2=1;

for m=[1 2 3]
    mod=model(m);
    nens=info_models(m,5);

    Ntot = zeros(nt, 3*nens*nsim);
    x=1;
    for scen=scenEXT % Extreme event scenario

        for ens = 1:nens

            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            N = readmatrix(file_name); %(nt, nsim)

            if m==5
                Ntot(:,x:x+nsim-1) = N(:,:);
            else
                Ntot(:,x:x+nsim-1) = N(1:end-1,:);
            end

            x=x+nsim;
        end %ens

    end %scenEXT

    Ntot = permute(Ntot, [2 1]); %(sim*ens, nt)

    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    N_all(x2:x2+(3*nens*nsim)-1, :) = Ntot;

    x2=x2+(3*nens*nsim);

end %model

Ntot_med = quantile(N_all, [0.025 0.5 0.975]);
Ntot = N_all;

% CALCULATE PEXT
    y_eval = year_eval - yearstart;
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)
    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for i=1:length(Ntot(:,1)); %total simulations
            if Ntot(i, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_clim_mod(4,:) = proba_tot;
    
    % Plot mean, 95% interval, internal variability

    couleur = [.46 .7 .18];

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur);
    set(f,'EdgeColor','none','FaceAlpha', 0.1)

    % Plot the 50 ens
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median
    y= plot(timePOP,Ntot_med(2,:),'color',couleur,'linewidth',5);

    % Plot observations
    plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

    hold on

    xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'bold');
    ylabel('Abundance', 'FontSize',fontsi, 'Fontweight', 'bold');
    ylim([0 300000]);
    title('3 climate models gathered', 'Fontsize', fontsi, 'Fontweight', 'bold');
    set(gca, 'FontSize', fontsi, 'Fontweight', 'bold');


%print('Fig_range_clim_mod_mature', '-dpng', '-r300')



%% EXTREME EVENT SCENARIO

fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);
tiledlayout(2,2)

% Climate model, scenario SSP370
mod = "LE_CESM2";
eco_mod="CMR";

% Model parameters
i_mod=6;
yearstart = 2009; yearstop = info_models(i_mod,4);
timePOP   = yearstart:yearstop;
nt        = length(timePOP);
nc = 66;
nens = info_models(i_mod,5);
nsim = 100;

i_scen =1;
i=1;

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.7 0.7 0.7];

% Extreme event scenario for CMR model = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

Pext_scenEXT = zeros(4,4); %3 models + gathered


for scen=scenEXT
    nexttile

    Ntot = zeros(nt, nens*nsim);

    x=1;
    for ens = 1:nens

        file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
        N = readmatrix(file_name); %(nt, nsim)

        Ntot(:,x:x+nsim-1) = Ntot(:,x:x+nsim-1) + N(:,:);
        x=x+nsim;
    end %ens


    Ntot = permute(Ntot, [2 1]); %(ens, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    Ntot_allscen(i_scen:i_scen+(nsim*nens)-1,:) = Ntot;

    % Parameters for the extinction probability
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %mean 2009-2018
    y_eval = year_eval - yearstart;

    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for k=1:length(Ntot(:,1)); %ensemble
            if Ntot(k, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_scenEXT(i,:) = proba_tot;

    % PLOT
    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(i,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    hold on

    % Plot median and 95% quantiles
    y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur(i,:),'linewidth',5);
    y_tot(m)=y;
    hold on
    plot(timePOP,Ntot_med(1,:),'-','color', couleur(i,:),'linewidth',2)
    hold on
    plot(timePOP, Ntot_med(3,:),'-','color', couleur(i,:),'linewidth',2)
    hold on

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

    fontsi = 24;
    y_tot(4) = plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    xlim([2000 2100])
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 300000])
    if scen==1 %SAT
        title('scen 1');
        ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
    elseif scen==2
        title('scen 2')
        ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
    else
        title('scen 3')
        ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
    end

    i_scen = i_scen + (nsim*nens);
    i=i+1;

end % scenEXT

nexttile

Ntot=Ntot_allscen;       
Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

% CALCULATE PEXT
y_eval = year_eval - yearstart;
pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)
threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
threshold_endan = 0.5 * pop_present; % Endangered
threshold_critend = 0.2 * pop_present; % Critically endangered
thr_ext = 0.05 * pop_present;

% Probability to go under the threshold
proba_tot=[0 0 0 0];
j=1;
for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
    proba = 0;
    under_thr = 0;
    for i=1:length(Ntot(:,1)); %total simulations
        if Ntot(i, y_eval) < thr;
            under_thr = under_thr + 1;
            proba = under_thr / length(Ntot(:,1));
        end
    end
    proba_tot(j)=proba;
    j=j+1;
end

Pext_scenEXT(4,:) = proba_tot;

% PLOT

% Plot mean, 95% interval, internal variability
couleur = [0, 0.4470, 0.7410];
ttt = [timePOP, fliplr(timePOP)];

inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.05)
hold on

% Plot 5 simulations
plot(timePOP, Ntot(randi(15000, 5,1 ),:),'color', [0.6 0.6 0.6],'linewidth',0.6)

hold on
% Plot median and 95% quantiles
y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur,'linewidth',5);
y_tot(m)=y;
hold on
plot(timePOP,Ntot_med(1,:),'-','color', couleur,'linewidth',2)
hold on
plot(timePOP, Ntot_med(3,:),'-','color', couleur,'linewidth',2)

hold on

fontsi = 24;
plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);


% Plot Pext and thr
coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl1.Color=coul_P{1};
yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl2.Color=coul_P{2};
yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl3.Color=coul_P{3};
yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
yl4.Color=coul_P{4};
yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)

xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
xlim([2000 2100])
ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');
ylim([0 300000]);
title('All ecological models', 'FontSize', fontsi, 'FontWeight', 'bold');
set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');


%% Internal variability

fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);
tiledlayout(1,2);

% Climate model, scenario SSP370
mod = "LE_CESM2";

% Ecological models
e_mod = "CMR";

yearstart=2009;
yearstop=2100;
nc=66;
timePOP   = yearstart:yearstop;
nt        = length(timePOP);

% one color per model + one color for ensembles (grey)
couleur = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.7 0.7 0.7];

% Extreme event scenario for CMR model = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

Pext_intern = zeros(2,4);

for int=1:2
    nexttile

    if int==1
        nens=50;
    else
        nens=1;
    end

    Ntot = zeros(nt, 3*nsim*nens);

    if int==1 % internal var
        x=1;
        for scen=scenEXT % Extreme event scenario
            for ens = 1:nens

                file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
                N = readmatrix(file_name); %(nt, nsim)

                Ntot(:,x:x+nsim-1) = N(:,:);
                x=x+nsim;

            end %ens

        end %scenEXT

    else %no internal var
        x2=1;
        for scen=scenEXT % Extreme event scenario

            % MEAN OF ENSEMBLES

            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, scen);
            N2 = readmatrix(file_name); %(nt, nsim)

            Ntot(:,x2:x2+nsim-1) = N2(:,:);
            x2=x2+nsim;

        end %scenEXT
    end


    Ntot = permute(Ntot, [2 1]); %(nens*nsim, nt)
    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % Parameters for the extinction probability
    pop_present = mean(dataOBSSAT_bilgecan(:,2)); %mean 2009-2018
    y_eval = year_eval - yearstart;

    threshold_vul = 0.7 *  pop_present; % Vulnerable criteria
    threshold_endan = 0.5 * pop_present; % Endangered
    threshold_critend = 0.2 * pop_present; % Critically endangered
    thr_ext = 0.05 * pop_present;

    % Probability to go under the threshold
    proba_tot=[0 0 0 0];
    j=1;
    for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
        proba = 0;
        under_thr = 0;
        for i=1:length(Ntot(:,1)); %ensemble
            if Ntot(i, y_eval) < thr;
                under_thr = under_thr + 1;
                proba = under_thr / length(Ntot(:,1));
            end
        end
        proba_tot(j)=proba;
        j=j+1;
    end

    Pext_eco_mod(m,:) = proba_tot;

    % PLOT

    % Plot mean, 95% interval, internal variability

    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(int,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot the ensembles (one ens = mean of 100 simulations)
    %plot(timePOP, Ntot(:,:),'color', couleur(6,:),'linewidth',1)

    hold on

    % Plot median and 95% quantiles
    y=plot(timePOP, Ntot_med(2,:),'-', 'color', couleur(int,:),'linewidth',5);
    y_tot(m)=y;
    hold on
    plot(timePOP,Ntot_med(1,:),'-','color', couleur(int,:),'linewidth',2)
    hold on
    plot(timePOP, Ntot_med(3,:),'-','color', couleur(int,:),'linewidth',2)
    hold on

    % Plot Pext and thr
    coul_P={'#FFBF00', '#FE8000', '#FF0000', '#620000'};
    yl1=yline(threshold_vul,'-', proba_tot(1)*100, 'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl1.Color=coul_P{1};
    yl2=yline(threshold_endan, '-', proba_tot(2)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl2.Color=coul_P{2};
    yl3=yline(threshold_critend,'-', proba_tot(3)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl3.Color=coul_P{3};
    yl4=yline(thr_ext,'-', proba_tot(4)*100,'linewidth',1.7,'FontSize', 17, 'LabelHorizontalAlignment', 'Left');
    yl4.Color=coul_P{4};
    yline(pop_present, 'k--', 'N0', 'linewidth',1,'FontSize', 20)
    xline(year_eval,'k', year_eval,'linewidth',1, 'LabelOrientation', 'horizontal','FontSize', 20)


    fontsi = 24;
    y_tot(4) = plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);
    xlabel('Years', 'FontSize', fontsi, 'FontWeight', 'bold');
    xlim([2009 2100])
    set(gca, 'FontSize', fontsi, 'FontWeight', 'bold');
    ylim([0 300000])
    if int==1
        title('Internal var');
    elseif int==2
        title('No internal var')
    end
    ylabel('Abundance', 'FontSize', fontsi, 'FontWeight', 'bold');

Pext_intern(int,:) = proba_tot;


end % internal var or no




%% CREATE TABLE

Pext_tot=zeros(18,4);
Pext_tot(1:4,:) = Pext_eco_mod;
Pext_tot(5:8,:) = Pext_clim_scen;
Pext_tot(9:12,:) = Pext_clim_mod;
Pext_tot(13:16, :) = Pext_scenEXT;
Pext_tot(17:18,:) = Pext_intern;
Pext_tot = Pext_tot*100;
Pext_tot = round(Pext_tot);

writematrix(Pext_tot, 'All_Pext.xlsx');

table = array2table(Pext_tot);
table.Properties.VariableNames = {'Vulnerable', 'Endangered', 'Critically', 'Extinct'};
table.Properties.RowNames = {'CMR', 'IPM', 'SAT','Ecological models gathered', 'RCP8.5', 'SSP370', 'Paris2', 'Climate scenario gathered', 'CanESM', 'CESM1', 'CSIRO', 'Climate models gathered', 'scen1', 'scen2', 'scen3', 'gathered', 'internal var', 'no int var'};
T=table;

writetable(table,'All_Pext.csv')

%create a new figure.
fig = uifigure;
set(fig, 'Position', [100, 100, 700, 500]);
%create a uitable
uit = uitable(fig, "Data",T, 'Position', [20 20 660 480]);

%exportapp(fig, 'table.jpg')

% Capture the table window as an image
frame = getframe(fig); % Get the current figure window
imageData = frame2im(frame); % Convert the frame to an image
% Save the image as a PNG file
imwrite(imageData, 'Pext_table.png');








