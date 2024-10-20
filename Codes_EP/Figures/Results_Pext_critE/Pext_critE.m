%% CRITERIA E - TOTAL POPULATION

% Calculate probability of extinction for each colony for each ecological
% model, under criteria E
% Proba to be under a threshold or 10 or 100 mature individuals in 2073

% 10 and 100 mature individuals / col


% INITIALISATION

load(sprintf("%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat", ordi))
dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));
sites = readmatrix(sprintf('%s/Codes_EP/Codes_CMR/empe_sitesNewNB.csv', ordi));

e_mod={'CMR', 'IPM', 'Sat'};
yearstart=2009;
yearstop=2100;
timePOP=yearstart:yearstop;
nt=length(timePOP);
ncol = 66;
nsim = 100;
nens=50;

time_present= 2024;
GL = 16.4; % Generation length = average age of parents in the current cohort
year_eval_list = [round(time_present + 3*GL) 2100];

Pext_col = zeros(66, 7); %for the csv file : columns = indice_col, lat, long, 10_2073, 100_2073, 10_2073, 100_2100


mod='LE_CESM2';


% SIMULATIONS

for m=[1:3] %ecological models

    eco_mod=e_mod{m};
    if m ==1
        scenEXT =[1 2 4];
    elseif m==2
        scenEXT =[1];
    else
        scenEXT =[1];
    end

    Ntot = zeros(nt, length(scenEXT)*nens*nsim);
    Pext_col = zeros(66, 7);
    
    for col=1:66 % colonies

        sad_tot=zeros(50*length(scenEXT),1); %sad in 2024 mean of 100 sim and 50 ens

        display(col)

        x=1;
        xi=1;
        for scen=scenEXT % Extreme event scenario
            for ens = 1:nens
                if m ==1
                    file_name = sprintf("%s/Codes_EP/N_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_scen%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, scen,col);
                elseif m==2 %IPM and SAT
                    file_name = sprintf("%s/Codes_EP/N_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens,col);
                else
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens,col);
                end
                N = readmatrix(file_name); %(nt, nsim)

                if m ==1
                    file_name2 = sprintf("%s/Codes_EP/MAT_results_tot/2009/MAT_results_Stef/%s/MAT_Stef_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scen,col);
                    sad = readmatrix(file_name2); %(nt, nsim)
                    sad_tot(xi) = mean(sad(16,:),2);
                elseif m==2
                    file_name2 = sprintf("%s/Codes_EP/MAT_results_tot/2009/MAT_results_%s/CESM2/MAT_%s_CESM2_ens%d_col%d.txt", ordi, eco_mod, eco_mod, ens,col);
                    sad = readmatrix(file_name2); %(nt, nsim)
                    sad_tot(ens) = mean(sad(16,:),2);
                end

                Ntot(:,x:x+nsim-1) = N(:,:);

                x=x+nsim;
                xi=xi+1;
                
            end %ens
        end %scenEXT

        if m<3 %CMR IPM
            sad_2024=mean(sad_tot); %sad 2024 for each col
        else
            sad_2024=1; % SAT we don't have
        end

        Nall = permute(Ntot, [2 1]); %(nb_scenEXT*sim*ens, nt)
        N_all_med=quantile(Nall, [0.025 0.5 0.975]);

        % CALCULATE PEXT

        %pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)

        pop_present = N_all_med(2,1);

        % 10 or 100 mature individuals -> need to know SAD
        thr_10 = 10/sad_2024;
        thr_100 = 100/sad_2024;

        for year_eval=year_eval_list

            y_eval = year_eval - yearstart;
            if m==3 % SAT 2099
                y_eval=y_eval-1;
            end

            % Probability to go under the threshold
            proba_tot=[0 0]; %10 100
            j=1;

            for thr = [thr_10 thr_100]
                proba = 0;
                under_thr = 0;
                for i=1:length(Nall(:,1)); %total simulations
                    if Nall(i, y_eval) < thr;
                        under_thr = under_thr + 1;
                        proba = under_thr / length(Nall(:,1));
                    end
                end
                proba_tot(j)=proba;
                j=j+1;
            end

            if year_eval==2073
                %for the csv file : columns = indice_col, lat, long, Pext_vul, Pext_endan, Pext_crit, Pext_extinct
                Pext_col(col, 4) = proba_tot(1);
                Pext_col(col, 5) = proba_tot(2);

            else %2100
                Pext_col(col, 6) = proba_tot(1);
                Pext_col(col, 7) = proba_tot(2);

            end

        end %year eval (2068, 2100)

    end %col

    Pext_col(:, 1) = 1:66; % col#
    Pext_col(:, 2) = sites(:,5); %latitude
    Pext_col(:, 3) = sites(:,6); % longitude

    if m==1
        Pext_col_CMR = Pext_col;
    elseif m==2
        Pext_col_IPM = Pext_col;
    else
        Pext_col_SAT = Pext_col;
    end

end %eco mod


%save('Pext_col_CMR_critE', 'Pext_col_CMR')
%save('Pext_col_IPM_critE', 'Pext_col_IPM')
%save('Pext_col_SAT_critE', 'Pext_col_SAT')


%% Transform in csv file

% load(sprintf('%s/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_ecomod/Pext_col_CMR_critE.mat', ordi))
% load(sprintf('%s/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_ecomod/Pext_col_IPM_critE.mat', ordi))
% load(sprintf('%s/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_ecomod/Pext_col_SAT_critE.mat', ordi))

%colony size
N_avg_66=[3324	1788	3771	3198	4424	4614	12337	4307	6012	7365	4970	6831	4416	3252	217	4626	2789	112	70	69	1289	287	798	4035	2664	4267	70	359	70	70	3289	4083	1320	1553	70	70	70	9711	1578	70	5292	70	2212	7983	26477	14040	70	6348	1416	1071	22266	1753	70	616	3403	70	3684	3872	3046	70	1036	4314	70	797 70 70];

for m=[1:3] %ecological model

    if m==1
        Pext_col = Pext_col_CMR;
    elseif m==2
        Pext_col = Pext_col_IPM;
    else
        Pext_col = Pext_col_SAT;
    end

    Pext_col(:,8) = N_avg_66; %satellite mean 2009-2018

    csv_file = array2table(Pext_col);
    csv_file.Properties.VariableNames = {'Colony', 'Lat', 'Long', 'Pext_10_2073', 'Pext_100_2073', 'Pext_10_2100', 'Pext_100_2100', 'N'};

    if m==1
        writetable(csv_file, 'Pext_col_CMR_SSP370_critE.csv')
    elseif m==2
        writetable(csv_file, 'Pext_col_IPM_SSP370_critE.csv')
    else
        writetable(csv_file, 'Pext_col_SAT_SSP370_critE.csv')
    end

end %eco mod