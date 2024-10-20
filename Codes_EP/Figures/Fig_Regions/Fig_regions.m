%% Figure regions
% Plot Regions (2 types of regions)

%choose path to Codes_EP
ordi = 'D:/Documents/alice'; 

%% Calculate colony sizes
ncol=66;
nens=50;
nsim=100;
yearstart=2009;
yearstop=2100;
timePOP=yearstart:yearstop;
nt=length(timePOP);
scenEXT=[1 2 4];
mod="LE_CESM2";

eco_mod="CMR";

Ntot_allcol = zeros(3*nens*nsim, nt, ncol);

for col=1:ncol

    display(col)

    Ncol = zeros(nt, 3*nens*nsim);

    x=1;
    for scen=scenEXT % Extreme event scenario
        
        for ens = 1:nens

            %load file with new K for each col
            file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_scen%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, scen,col);
            N = readmatrix(file_name); %(nt, nsim)

            Ncol(:,x:x+nsim-1) = Ncol(:,x:x+nsim-1) + N(:,:);
            x=x+nsim;
        end %ens

    end %scenEXT

    Ncol = permute(Ncol, [2 1]); %(sim*ens, nt)

    Ntot_allcol(:,:,col) = Ncol;

end %col

% save('Ntot_allcol_july.mat', 'Ntot_allcol') 


% -------------------------------------------------------------------------
%% REGIONS CCAMLR

load Ntot_allcol_july.mat

% Calculate Region sizes

% Atlantic ocean
%1 48.1
%2 48.5
%3 48.6
% Indian Ocean
%4 58.4.2
%5 58.4.1
% Pacific Ocean
%6 88.1
%7 88.2
%8 88.3

REGIONS=[1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 1 2];
numRegions=max(REGIONS);

regionalSum = zeros(3*nens*nsim, nt, numRegions);

x=1;
for scenEXT = [1 2 4]
    for ens=1:nens
        for sim=1:nsim

            for re=1:numRegions
                regionalSum(x, :, re) = sum(Ntot_allcol(x,:,REGIONS==re),3);
            end %regions
            x=x+1;

        end %sim
    end %ens
end

regionalSum_med=zeros(3,nt,numRegions);
for re=1:numRegions
    regionalSum_med(:,:,re) = quantile(regionalSum(:,:,re),[0.05,0.5,0.95]);
end


%% Plot Nmature for each region 

fig = figure;
set(fig, 'Position', [0, 0, 1800, 900]);

regions_name = {'48.1', '48.5', '48.6', '58.4.2', '58.4.1', '88.1', '88.2', '88.3'};
ylist = [20000 140000 90000 50000 80000 160000 30000 20000];


tiledlayout(2,4)

for re=1:numRegions
    nexttile
    

    Nreg = regionalSum(:,:,re);

    Nreg_med = regionalSum_med(:, :, re);

    % PLOT

    ttt = [timePOP, fliplr(timePOP)];
    inBetween = [Nreg_med(1,:), fliplr(Nreg_med(3,:))];
    f = fill(ttt, inBetween, couleur(1,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot the 50 ens
    %plot(timePOP, Nreg(:,:),'color', couleur(2,:),'linewidth',1)
    hold on

    % Plot median and quantiles
    plot(timePOP, Nreg_med(2,:),'color',couleur(1,:),'linewidth',5)
    hold on
    y = plot(timePOP, Nreg_med(1,:),'color',couleur(1,:),'linewidth',2);
    hold on
    plot(timePOP, Nreg_med(3,:),'color',couleur(1,:),'linewidth',2)
    hold on


    fontsi = 24;
    % FIGURE PARAMETERS
    xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'Bold');
   
    set(gca, 'FontSize', fontsi, 'Fontweight', 'Bold');
    titre=regions_name{re};
    title(titre, 'FontSize', fontsi, 'Fontweight', 'Bold');
    ylim([0 ylist(re)])

     ylabel('Nb of mature individuals', 'FontSize', 20, 'Fontweight', 'Bold');



    %save as high resolution
    %fig_name = sprintf("Fig_regCCMALR_%s.png", regions_name{re});
    %print(fig_name, '-dpng', '-r300')
    %close
end % Region



% -------------------------------------------------------------------------
%% REGIONS GENETIC

load Ntot_allcol_july.mat

% Calculate Region sizes

% Atlantic ocean
%1 48.1
%2 48.5
%3 48.6
% Indian Ocean
%4 58.4.2
%5 58.4.1
% Pacific Ocean
%6 88.1
%7 88.2
%8 88.3

REGIONS=[1 1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7];
numRegions=max(REGIONS);

regionalSum = zeros(3*nens*nsim, nt, numRegions);

x=1;
for scenEXT = [1 2 4]
    for ens=1:nens
        for sim=1:nsim

            for re=1:numRegions
                regionalSum(x, :, re) = sum(Ntot_allcol(x,:,REGIONS==re),3);
            end %regions
            x=x+1;

        end %sim
    end %ens
end

regionalSum_med=zeros(3,nt,numRegions);
for re=1:numRegions
    regionalSum_med(:,:,re) = quantile(regionalSum(:,:,re),[0.05,0.5,0.95]);
end


%% Plot Nmature for each region and save as pdf


regions_name = {'StoS' 'WEDD' 'StoKP' 'MAWS' 'AMPG' 'ROSS' 'AtoBe'};
ylist = [70000 80000 130000 30000 90000 200000 50000 ];

fig = figure;
set(fig, 'Position', [0, 0, 1800, 900]);


tiledlayout(2,4)

for re=1:numRegions
    nexttile

    %fig=figure;

    Nreg = regionalSum(:,:,re);

    Nreg_med = regionalSum_med(:, :, re);

    % PLOT

    ttt = [timePOP, fliplr(timePOP)];
    inBetween = [Nreg_med(1,:), fliplr(Nreg_med(3,:))];
    f = fill(ttt, inBetween, couleur(1,:));
    set(f,'EdgeColor','none','FaceAlpha', 0.1)
    hold on

    % Plot the 50 ens
    %plot(timePOP, Nreg(:,:),'color', couleur(2,:),'linewidth',1)
    hold on

    % Plot median and quantiles
    plot(timePOP, Nreg_med(2,:),'color',couleur(1,:),'linewidth',5)
    hold on
    y = plot(timePOP, Nreg_med(1,:),'color',couleur(1,:),'linewidth',2);
    hold on
    plot(timePOP, Nreg_med(3,:),'color',couleur(1,:),'linewidth',2)
    hold on


    fontsi = 24;
    % FIGURE PARAMETERS
    xlabel('Years', 'FontSize', fontsi, 'Fontweight', 'Bold');
    ylabel('Abundance', 'FontSize', fontsi, 'Fontweight', 'Bold');
    set(gca, 'FontSize', fontsi, 'Fontweight', 'Bold');
    titre=regions_name{re};
    title(titre, 'FontSize', fontsi, 'Fontweight', 'Bold');
    ylim([0 ylist(re)])

    
    %pause
    %close

    %save as high resolution
    fig_name = sprintf("Fig_reg_genet_%d.png", re);
    %print(fig_name, '-dpng', '-r300')
    %close
    
end % Region

