%% Figure climate models and internal variability
% For CMR model and RCP8.5 climate scenario

% Ecological model = Stef's model, CMR
eco_mod = "CMR";

% Load info models
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi))
% Climate scenario = RCP8.5
model = ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
model_name = ["CanESM", "CESM", "CSIRO", "GFDL-CM3", "MPI-ESM","CESM2", "CESM1-CAM5-PARIS2"];
nmod = length(model);

% one color per model + one color for ensembles (grey)
couleur = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0 0 0; 0 0 0; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880;
    0.8 0.8 0.8];

% Extreme event scenario = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

dataOBSSAT_bilgecan=readmatrix(sprintf("%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv", ordi));


%% FIGURE WITH THE THREE MODELS - from 2009

leg = ['b' 'c' 'd'];

fig = figure;
set(fig, 'Position', [100, 100, 1600, 1100]);

y_tot = zeros(3,1);

subplot(2,2,1)

for m=1:3

    mod=model(m);
    y_tot2=zeros(2,1);

    % Model parameters
    yearstart = 2009; yearstop = info_models(m,4);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nc = 66;
    nens = info_models(m,5);
    nsim = 100;

    Ntot = zeros(nt, 3*nens*nsim);

    x=1;
    for scen=scenEXT % Extreme event scenario
        for ens = 1:nens

            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            N = readmatrix(file_name); %(nt, nsim)

            Ntot(:,x:x+nsim-1) = N(:,:); 
            x=x+nsim;

        end %ens

    end %scenEXT

    Ntot = permute(Ntot, [2 1]); %(ens, nt)

    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % PLOT

    tt = [timePOP, fliplr(timePOP)];
    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on
    y=plot(timePOP, Ntot_med(2,:),'color',couleur(m,:),'linewidth',5);
    y_tot(m)=y;
    hold on

end % climate model


% PARAMETERS FOR FIRST TILE

hold on
xlabel('Years', 'FontSize', 24, 'Fontweight', 'bold');
ylabel('Nb of mature individuals', 'FontSize', 24, 'Fontweight', 'bold');
ylim([0 300000]);
title('Climate models');
set(gca, 'FontSize', 24, 'Fontweight', 'bold');
xlim([yearstart yearstop]);
legend(y_tot, model_name(1:3), 'Location', 'Southwest', 'Fontsize', 16);
text(-0.2, 1.1, 'a', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');



% FIGURE CLIMATE MODELS AND INTERNAL VARIABILITY

eco_mod = 'CMR';

% Total Population - Each model on one figure

model = ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0"];
model_name = ["CanESM", "CESM", "CSIRO"];
label = {"Clim_mod_CanESM2", "Clim_mod_CESM1", "Clim_mod_CSIRO"};

% one color per model + one color for ensembles (grey)
couleur = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880;
    0.8 0.8 0.8];


y_tot = zeros(3,1);

for m=1:3
    subplot(2,2,m+1)

    mod=model(m);
    y_tot2=zeros(2,1);

    % Model parameters
    yearstart = 2009; yearstop = info_models(m,4);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nc = 66;
    nens = info_models(m,5);
    nsim = 100;

    Ntot = zeros(nt, 3*nens*nsim);
    Ntot2 = zeros(nt, 3*nsim);

    x=1;
    x2=1;
    for scen=scenEXT % Extreme event scenario
        for ens = 1:nens

            %file_name = sprintf("/Users/aliceeparvier/Desktop/Codes_EP/N_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens, scen);
            N = readmatrix(file_name); %(nt, nsim)

            Ntot(:,x:x+nsim-1) = N(:,:); 
            x=x+nsim;

        end %ens

            % MEAN OF ENSEMBLES

            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_scen%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, scen);
            N2 = readmatrix(file_name); %(nt, nsim)

            Ntot2(:,x2:x2+nsim-1) = N2(:,:);
            x2=x2+nsim;

    end %scenEXT

    Ntot = permute(Ntot, [2 1]); %(ens, nt)

    Ntot_med = quantile(Ntot, [0.025 0.5 0.975]);

    % MEAN OF ENSEMBLES

    Ntot2 = permute(Ntot2, [2 1]); %(sim, nt)

    Ntot_med2 = quantile(Ntot2, [0.025 0.5 0.975]);

    % PLOT

    ttt = [timePOP, fliplr(timePOP)];
    inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot 5 simulations
    plot(timePOP, Ntot(randi(3*nens*nsim, 5, 1),:),'color', couleur(6,:),'linewidth',1)
    hold on

    % Plot median and quantiles 
    plot(timePOP,Ntot_med(2,:),'color',couleur(m,:),'linewidth',5)
    hold on
    y = plot(timePOP,Ntot_med(1,:),'color',couleur(m,:),'linewidth',2);
    %y_tot2(1)=y;
    hold on
    plot(timePOP,Ntot_med(3,:),'color',couleur(m,:),'linewidth',2)
    hold on

    % Plot mean of ensembles and quantiles
    a=plot(timePOP,Ntot_med2(2,:), '--','color',couleur(m,:),'linewidth',5, 'MarkerSize', 12);
    a.MarkerIndices = 1:10:length(Ntot_med2(2,:));
    hold on
    a=plot(timePOP,Ntot_med2(1,:),'--','color',couleur(m,:),'linewidth',2, 'MarkerSize', 12);
    a.MarkerIndices = 1:10:length(Ntot_med2(2,:));
    y_tot2(2)=a;
    hold on
    a=plot(timePOP,Ntot_med2(3,:),'--','color',couleur(m,:),'linewidth',2, 'MarkerSize', 12);
    a.MarkerIndices = 1:10:length(Ntot_med2(2,:));

    hold on

    % FIGURE PARAMETERS
    xlabel('Years', 'FontSize', 24, 'Fontweight', 'bold');
    ylabel('Nb of mature individuals', 'FontSize', 24, 'Fontweight', 'bold');
    ylim([0 300000]);
    xlim([yearstart 2100])
    title(model_name(m));
    set(gca, 'FontSize', 24, 'Fontweight', 'bold');
    legend(y_tot2(2), 'Mean of ensembles','Location', 'Southwest', 'Fontsize', 16)   


    text(-0.2, 1.1, leg(m), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');

end % climate model


