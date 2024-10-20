%% CALCULATE PEXT FOR EACH COLONY AND FOR TOTAL POPULATION
% FOR EACH ECOLOGICAL MODELS (separate results)

% Criteria A
% We take the mature population
% Calculate proba to be under 3 IUCN thresholds (Vulnerable 30%, Endangered
% 50%, Crit endangered 80%)


% Directory to Codes_EP
ordi = 'D:/Documents/alice';

%% INITIALISATION

load(sprintf("%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat", ordi))
dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_CMR/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));
sites = readmatrix(sprintf('%s/Codes_EP/Codes_CMR/empe_sitesNewNB.csv', ordi));

e_mod={'CMR', 'IPM', 'Sat'};
yearstart=2009;
ncol = 66;
nsim = 100;
nens=50;

time_present= 2024;
GL = 16.4; % Generation length = average age of parents in the current cohort
year_eval_list = [round(time_present + 3*GL) 2100];

mod='LE_CESM2';


% SIMULATIONS

for m=[1:3] %ecological models

    eco_mod=e_mod{m};
    if m ==1
        timePOP=2009:2100;
        nt=length(timePOP);
        scenEXT =[1 2 4];
    elseif m==2
        timePOP=2009:2100;
        nt=length(timePOP);
        scenEXT =[1];
    else
        timePOP=2009:2099;
        nt=length(timePOP);
        scenEXT =[1];
    end

    Ntot = zeros(nt, length(scenEXT)*nens*nsim);
    Pext_col = zeros(66, 11);

    for col=1:66 % colonies

        display(col)

        x=1;
        for scen=scenEXT % Extreme event scenario
            for ens = 1:nens
                if m ==1
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_scen%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, scen,col);
                else %IPM and SAT
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens,col);
                end

                N = readmatrix(file_name); %(nt, nsim)

                Ntot(:,x:x+nsim-1) = N(:,:);

                x=x+nsim;
            end %ens
        end %scenEXT

        Nall = permute(Ntot, [2 1]); %(nb_scenEXT*sim*ens, nt)
        N_all_med=quantile(Nall, [0.025 0.5 0.975]);

        % CALCULATE PEXT

        %pop_present = mean(dataOBSSAT_bilgecan(:,2)); %observed current population (mean 2009-2018)

        pop_present = N_all_med(2,1);
        threshold_vul = 0.7 * pop_present; % Vulnerable criteria
        threshold_endan = 0.5 * pop_present; % Endangered
        threshold_critend = 0.2 * pop_present; % Critically endangered
        thr_ext = 0.05 * pop_present; % Quasi-extinct


        for year_eval=year_eval_list

            y_eval = year_eval - yearstart;

            % Probability to go under the threshold
            proba_tot=[0 0 0 0];
            j=1;

            for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
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
            Pext_col(col, 6) = proba_tot(3);
            Pext_col(col, 7) = proba_tot(4);
        else %2100
            Pext_col(col, 8) = proba_tot(1);
            Pext_col(col, 9) = proba_tot(2);
            Pext_col(col, 10) = proba_tot(3);
            Pext_col(col, 11) = proba_tot(4);
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


%save('Pext_col_CMR_mature', 'Pext_col_CMR')
%save('Pext_col_IPM_mature', 'Pext_col_IPM')
%save('Pext_col_SAT_mature', 'Pext_col_SAT')

%% Transform in csv file

% load(sprintf('%s/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_colonies/Pext_col_CMR_mature.mat', ordi))
% load(sprintf('%s/Codes_EP/Figures/Fig_Pext_ecomod/Pext_col_IPM_mature.mat', ordi))
% load(sprintf('%s/Codes_EP/Figures/Fig_Pext_ecomod/Pext_col_SAT_mature.mat', ordi))

%colony size
N_avg_66=[3324	1788	3771	3198	4424	4614	12337	4307	6012	7365	4970	6831	4416	3252	217	4626	2789	112	70	69	1289	287	798	4035	2664	4267	70	359	70	70	3289	4083	1320	1553	70	70	70	9711	1578	70	5292	70	2212	7983	26477	14040	70	6348	1416	1071	22266	1753	70	616	3403	70	3684	3872	3046	70	1036	4314	70	797 70 70];

for m=[1 2 3] %ecological model

    % Create file with status for each col
    % Threshold over which col is likely to be at this level : 50-59%
    % Status: 1=OK, 2=Vulnerable, 3=Endangered, 4=Critically Endangered,
    % 5=Quasi extinct
    stat = 2:5;
    status=ones(66,1);

    if m==1
        Pext_col = Pext_col_CMR;
    elseif m==2
        Pext_col = Pext_col_IPM;
    else
        Pext_col = Pext_col_SAT;
    end


    % 2068
    for col=1:66
        for thr=[1:4]
            if Pext_col(col, thr+3) > 0.40 %cf risk tolerance, RedListGuidelines
                status(col) = stat(thr);
            end
        end
    end

    Pext_col(:,12) = status; %add to previous file

    % 2100
    for col=1:66
        for thr=[1:4]
            if Pext_col(col, thr+7) > 0.40 %cf risk tolerance, RedListGuidelines
                status(col) = stat(thr);
            end
        end
    end
    Pext_col(:,13) = status; %status in 2100


    Pext_col(:,14) = N_avg_66; %satellite mean 2009-2018

    csv_file = array2table(Pext_col);
    csv_file.Properties.VariableNames = {'Colony', 'Lat', 'Long', 'Pext_vul_2068', 'Pext_endan_2068', 'Pext_crit_2068', 'Pext_ext_2068', 'Pext_vul_2100', 'Pext_endan_2100', 'Pext_crit_2100', 'Pext_ext_2100','status_2068', 'status_2100', 'N'};

    if m==1
        writetable(csv_file, 'Pext_col_CMR_SSP370_mature.csv')
    elseif m==2
        writetable(csv_file, 'Pext_col_IPM_SSP370_mature.csv')
    else
        writetable(csv_file, 'Pext_col_SAT_SSP370_mature.csv')
    end

end %eco mod



