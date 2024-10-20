%% FIGURE CLIMATE SCENARIO
% Plot 3 climate scenarios and 3 extreme event scenarios
% Figure only for the CMR model

% DIRECTORY to Codes_EP
ordi = "";

% Climate scenario 
% RCP8.5, SSP370, Paris Agreement 2degrees
model = ["LE_CESM1-CAM5", "LE_CESM2", "LE_CESM1-CAM5-PARIS2"];
model_name = ["RCP8.5" "SSP370" "Paris 2°C"];
indice_mod = [2 6 7];

% Load info models
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi));

% Ecological model = Stef's model, CMR
eco_mod = "CMR";

% one color per model + one color for ensembles (grey)
couleur = [0.8500 0.3250 0.0980; 
    0.4940 0.1840 0.5560;
    0 0.4470 0.7410;
    0.7 0.7 0.7];

symbols= {'-*' '-s' '-o'};

% Extreme event scenario = 1, 2 and 4 gathered
scenEXT = [1, 2, 4];

% TOTAL POPULATION - 3 CLIMATE SCENARIO + SCENEXT

fontsi=24;

fig = figure;
set(fig, 'Position', [100, 100, 1600, 1100]);


subplot(2,2,1)
y_tot = zeros(3,1); %for legend 

for m=1:length(model) %Climate scenario
    
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

    hold on

end % climate scenario

legend(y_tot, model_name, 'Location', 'Southwest', 'FontSize', 14);
xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'bold');
ylabel('Nb of mature individuals', 'FontSize',fontsi, 'Fontweight', 'bold');
ylim([0 300000]);
xlim([2009 2100]);
title('Climate scenarios', 'Fontsize', fontsi, 'Fontweight', 'bold');
set(gca, 'FontSize', fontsi, 'Fontweight', 'bold');
text(-0.2, 1.1, 'a', 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');



% EXTREME EVENT SCENARIO
leg = ['b' 'c' 'd'];
label = {'RCP', 'SSP', 'Paris'}; %label file

for m=1:length(model) %Climate model

    subplot(2,2,m+1)
   
    mod=model(m);

    % Model parameters
    i_mod=indice_mod(m);
    yearstart = 2009; yearstop = info_models(i_mod,4);
    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    nc = 66;
    nens = info_models(i_mod,5);
    nsim = 100;

    i_scen =1;
    i=1;
    y_tot=zeros(3,1); %for legend figures

    Ntot_allscen = zeros(3*nens*nsim, nt);

    for scen=scenEXT

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

        hold on

        % Plot median of simulations for one scenEXT
        y = plot(timePOP, Ntot_med(2,:), symbols{i}, 'MarkerSize', 10, 'color', couleur(m,:),'linewidth',2);
        y.MarkerIndices = 1:10:length(Ntot_med(2,:));
        y_tot(i) = y;
        hold on    

        i_scen = i_scen + (nsim*nens);
        i=i+1;

    end %scenEXT

    % Plot 95% confidence interval for all scen, all sim, all ens

    Ntot_allscen = mean(Ntot_allscen(:,:,:), 3);
    Ntot_allscen_med = quantile(Ntot_allscen, [0.025 0.5 0.975]);
    ttt = [timePOP, fliplr(timePOP)];

    inBetween = [Ntot_allscen_med(1,:), fliplr(Ntot_allscen_med(3,:))];
    f = fill(ttt, inBetween, couleur(m,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)

    legend(y_tot, 'scen 1', 'scen 2', 'scen 3', 'Location', 'Southwest', 'Fontsize', 14);

    xlabel('Years', 'Fontsize', fontsi, 'Fontweight', 'bold');
    ylabel('Nb of mature individuals', 'Fontsize', fontsi, 'Fontweight', 'bold');
    ylim([0 300000]);
    title(model_name(m), 'Fontsize',fontsi, 'Fontweight', 'bold');
    set(gca, 'FontSize', fontsi, 'Fontweight', 'bold');
    text(-0.2, 1.1, leg(m), 'Units', 'normalized', 'FontSize', 24, 'FontWeight', 'bold');
    xlim([2009 2100]);

     box on

end % climate scenario
