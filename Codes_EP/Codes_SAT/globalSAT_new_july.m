%% Global Satellite Model 
% Code for main projections, for projections with mean of ensembles

% DIRECTORY to the folder Codes_EP
ordi = '/Users/aliceeparvier/Desktop'

addpath(sprintf('%s/Codes_EP/Codes_SAT/', ordi))


%% Parameters

% Observed satellite data 2009-2018
dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));

% Import projected environmental data 1909-2100
% transformed data
load('SIC_proj_trans_total.mat'); %obtained with obtain_SIC_proj_data.m

mod="CESM2";

nens=50; %ensemble members
nsim=100; % to draw parameters in posterior distribution


old_col =[1:18, 20:26, 28, 31:34, 38,39, 41, 43:46, 48:52, 54,55, 57:59, 61, 62, 64]; %50 colonies
ncol=66;

yearstart = 1909;
yearstart = 2009;
yearend = 2100;
timePOP=yearstart:yearend;
nt=length(timePOP);
Ntot = zeros(nt, nsim, nens, ncol); %t, col, ens, simu

% obtain alphabetical order of 50 colonies
ind_letters = [45 27 19 43 23 31 25 16 46 20 37 4 41 3 29 36 24 49 00 2 28 21 47 5 15 1 00 6 00 00 26 42 9 35 00 00 00 18 34 00 32 00 17 38 13 50 00 22 7 14 12 40 00 30 48 00 8 11 33 00 10 44 00 39 00 00];

%% Simulations

idx = randi(5000, nsim, 1); %random simulations among the 5000 posteriors

% Parameter chains from Bilgecan
param_bilgecan = readmatrix('param_chains.csv'); %alpha beta sigma_time
param_bilgecan(1,:)=[];
param_bilgecan(:,1)=[];

beta0 = param_bilgecan(idx,1); %alpha
beta1 = param_bilgecan(idx,2); %beta
sigma_time = param_bilgecan(idx,4:53); %sigma_time 50 colonies: colonies in alphabetical order!
mu_sigma_time = param_bilgecan(idx,3); %mean of 50 colonies

eps = readmatrix('eps_chains.csv'); %epsilon
eps(1,:)=[];
eps(:,1)=[];
eps = eps(idx,:);


for col=1:66 %colonies
    display(col)

    for ens=1:nens %ensemble members

        Rtot=zeros(nt, nsim); %growth rate
        Rtot(1,:) = 1;

        for t=1:nt

            % SIC projections
            if yearstart==2009
                SIC=SIC_proj_total(t+100, col, ens); %select years
            else
                SIC=SIC_proj_total(t, col, ens);
            end

            if SIC==0
                SIC=0.0001;
            end

            SIC = log(SIC);

            for sim=1:nsim

                if ismember(col, old_col)==1 %if col is one of the 50 old colonies
                    col_bis = ind_letters(col); %obtain index in alphabetical order because sigma_time ordered alphabetically
                    s = sigma_time(sim, col_bis); %standard deviation
                else %for the 16 new colonies
                    s = mu_sigma_time(sim); %mean of 50 colonies
                end

                c = normrnd(0, s);
                Ntot(t, sim, ens, col) = exp(beta0(sim) + beta1(sim)*SIC + eps(sim, col) + c);

                % Calculate Growth rates
                %if t>1
                    %Rtot(t,sim) = Ntot(t, sim, ens, col) ./ Ntot(t-1, sim, ens, col);
                %end

            end %sim

        end %time

        % Save one txt file for each col, each ens : (nt, nsimu)

        % Save growth rates
        % if yearstart==1909
        %     filename=sprintf('%s/Codes_EP/GRdisp_mature_results_tot/GR_results_Sat/LE_CESM2/GR_Sat_%s_ens%d_col%d.txt', ordi, mod, ens, col);
        % else
        %     filename=sprintf('%s/Codes_EP/GRdisp_mature_results_tot/2009/GR_results_Sat/LE_CESM2/GR_Sat_%s_ens%d_col%d.txt', ordi, mod, ens, col);
        % end
        %writematrix(Rtot, filename);

        % Save the colony size data for each col
        if yearstart==1909
            filename=sprintf('%s/Codes_EP/N_mature_results_tot/N_results_Sat/LE_CESM2/N_Sat_LE_CESM2_ens%d_col%d.txt', ordi, ens, col);
        else
            filename=sprintf('%s/Codes_EP/N_mature_results_tot/2009/N_results_Sat/LE_CESM2/N_Sat_LE_CESM2_ens%d_col%d.txt', ordi, ens, col);
        end
        writematrix(Ntot(:, :, ens, col), filename);

    end %50 ens
end %col


% Save files for the total population size

Nntot=sum(Ntot(:,:,:,:), 4);

for ens=1:nens

    % GRtot = ones(nt, nsim); %TOTAL POPULATION
    % for t=2:nt
    %     GRtot(t,:)=Nntot(t,:,ens)./Nntot(t-1,:,ens);
    % end

    % SAVE growth rates
    % if yearstart==1909
    %     filename=sprintf('%s/Codes_EP/GRdisp_mature_results_tot/GR_results_Sat/GRtot_Sat/LE_CESM2/GRtot_Sat_%s_ens%d.txt', ordi, mod, ens);
    % else
    %     filename=sprintf('%s/Codes_EP/GRdisp_mature_results_tot/2009/GR_results_Sat/GRtot_Sat/LE_CESM2/GRtot_Sat_%s_ens%d.txt', ordi, mod, ens);
    % end
    % writematrix(GRtot, filename);


    if yearstart==1909
        filename=sprintf('%s/Codes_EP/N_mature_results_tot/N_results_Sat/Ntot_Sat/LE_CESM2/Ntot_Sat_LE_CESM2_ens%d.txt', ordi, ens);
    else
        filename=sprintf('%s/Codes_EP/N_mature_results_tot/2009/N_results_Sat/Ntot_Sat/LE_CESM2/Ntot_Sat_LE_CESM2_ens%d.txt', ordi, ens);
    end
    writematrix(Nntot(:,:, ens), filename);

end



%% FIGURE

figure
yearstart=1909;
yearend=2100;
%yearstart=2009;

timePOP=yearstart:yearend;
nt=length(timePOP);

Nntot = zeros(nt, nsim*nens);

x=1;
for ens=1:nens %calculate mean of simulations for each ens
    filename=sprintf('%s/Codes_EP/N_mature_results_tot/N_results_Sat/Ntot_Sat/LE_CESM2/Ntot_Sat_LE_CESM2_ens%d.txt', ordi, ens);
    N = readmatrix(filename);
    Nntot(:,x:(x+nsim-1))=N;
    plot(timePOP, mean(N(:,:), 2), 'linewidth', 0.5)
    hold on
    x=x+nsim;
end

Nntot=permute(Nntot, [2 1]);

hold on
%plot(2009:2018, dataOBSSAT_bilgecan(:,2),'-k','linewidth',5);

% Plot 95% interval

Ntot_med = quantile(Nntot, [0.025 0.5 0.975]);
ttt = [timePOP, fliplr(timePOP)];
inBetween = [Ntot_med(1,:), fliplr(Ntot_med(3,:))];
couleur=[0.5 0.5 0.5];
f = fill(ttt, inBetween, couleur);
set(f,'EdgeColor','none','FaceAlpha', 0.3)

% Plot min and max
plot(timePOP, min(Nntot), 'k')
hold on 
plot(timePOP, max(Nntot), 'k')


%% ---------------------------------------------------------------------------
% SIMULATIONS WITH THE MEAN OF ENSEMBLE MEMBERS - INTERNAL VARIABILITY

choix = 2 %1 : starts in 1909, 2 starts in 2009

% PARAMETERS

% Import projected environmental data 1909-2100
%load('SIC_proj_trans_total.mat'); %from obtain_SIC_proj_data

mod="CESM2";

nens=50;
nsim=100; % to draw parameters in posterior distribution

old_col =[1:18, 20:26, 28, 31:34, 38,39, 41, 43:46, 48:52, 54,55, 57:59, 61, 62, 64]; %50 colonies
ncol=66;

if choix==1
    yearstart = 1909;
else
    yearstart = 2009;
end

yearend = 2100;
timePOP=yearstart:yearend;
nt=length(timePOP);


%% Simulations
Ntot = zeros(nt, ncol, nsim); %t, col, ens, simu
mean_tot=zeros(nt, 2, ncol); %mean and std 

idx = randi(5000, nsim, 1); %random simulations

% Parameter chains from Bilgecan
param_bilgecan = readmatrix('param_chains.csv'); %alpha beta sigma_time
param_bilgecan(1,:)=[];
param_bilgecan(:,1)=[];

beta0 = param_bilgecan(idx,1); %alpha
beta1 = param_bilgecan(idx,2); %beta
sigma_time = param_bilgecan(idx,4:53); %sigma_time 50 colonies: colonies in alphabetical order!
mu_sigma_time = param_bilgecan(idx,3); %mean of 50 colonies

eps = readmatrix('eps_chains.csv'); %epsilon
eps(1,:)=[];
eps(:,1)=[];
eps = eps(idx,:);


for col=1:66
    display(col)

    Rtot=zeros(nt, nsim);
    Rtot(1,:) = 1;

    for t=1:nt

        % SIC projections, standardized
        if choix==1
            SIC = mean(SIC_proj_total(t, col, :),3); %mean of ensembles
        else
            SIC = mean(SIC_proj_total(t+100, col, :), 3); %for yearstart = 2009
        end

        if SIC==0
            SIC=0.0001;
        end

        %SIC = log(SIC_proj_total(t+100, col, ens)); %2009
        SIC = log(SIC);

        % Posterior / Residual variance across space (process variance)
        parfor sim=1:nsim

            if ismember(col, old_col)==1 %if col is one of the 50 old colonies
                col_bis = ind_letters(col); %obtain index in alphabetical order because sigma_time ordered alphabetically
                s = sigma_time(sim, col_bis);
            else %for the 16 new colonies
                s = mu_sigma_time(sim); %mean of 50 colonies
            end
            c = normrnd(0, s);

            Ntot(t, col, sim) = exp(beta0(sim) + beta1(sim)*SIC + eps(sim, col) + c);

        end %sim param

    end %time

    % Save the colony size data for each col
    if choix==1
        filename=sprintf('%s/Codes_EP/N_mature_results_tot/mean_ensemble/N_results_Sat/LE_CESM2/Nmean_Sat_LE_CESM2_col%d.txt', ordi, col);
    else
        filename=sprintf('%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_Sat/LE_CESM2/Nmean_Sat_LE_CESM2_col%d.txt', ordi, col);
    end
    writematrix(Ntot(:, col, :), filename);


end %col

Ntot = permute(Ntot, [1 3 2]); %t, sim, col
Nntot=sum(Ntot(:,:,:), 3);

% Save files for total population

filename=sprintf('%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_Sat/Ntot_Sat/LE_CESM2/Ntot_Sat_LE_CESM2.txt', ordi);
writematrix(Nntot(:,:), filename);
