%% Pext regions 

% Calculate Pext for each region and each ecological model
% Criteria A : three IUCN thresholds

% DIRECTORY
ordi='/Users/aliceeparvier/Desktop'
addpath(sprintf('%s/Codes_EP/Figures/Fig_july/N_mature/Fig_Pext_regions', ordi))

%% Calculate colony sizes
ncol=66;
nens=50;
nsim=100;
yearstart=2009;
yearstop=2100;
e_mod={'CMR', 'IPM', 'Sat'};

for m=3
    eco_mod=e_mod{m};

    timePOP=yearstart:yearstop;
    nt=length(timePOP);
    scenEXT=[1 2 4];

    mod="LE_CESM2";

    Ntot_allcol = zeros(length(scenEXT)*nens*nsim, nt, ncol);

    for col=1:ncol

        display(col)

        Ncol = zeros(nt, length(scenEXT)*nens*nsim);

        x=1;
        for scen=scenEXT % Extreme event scenario

            for ens = 1:nens

                %load file with new K for each col
                if m==1
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_scen%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, scen,col);
                else
                    file_name = sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens,col);
                end
                    N = readmatrix(file_name); %(nt, nsim)

                Ncol(:,x:x+nsim-1) = Ncol(:,x:x+nsim-1) + N(:,:);
                x=x+nsim;
            end %ens

        end %scenEXT

        Ncol = permute(Ncol, [2 1]); %(sim*ens, nt)

        Ntot_allcol(:,:,col) = Ncol;

    end %col

    if m==1
        % save('Ntot_allcol_CMR.mat', 'Ntot_allcol')
    elseif m==2
        % save('Ntot_allcol_IPM.mat', 'Ntot_allcol')
    else
        % save('Ntot_allcol_SAT.mat', 'Ntot_allcol')
    end

end


%% Calculate Pext of regions

e_mod={'CMR', 'IPM', 'Sat'};
yearstart=2009;
timePOP=2009:2100;
nt=length(timePOP);
ncol = 66;
nsim = 100;
nens=50;
mod='LE_CESM2';

time_present= 2024;
GL = 16.4; % Generation length = average age of parents in the current cohort
year_eval_list = [round(time_present + 3*GL) 2100];

Pext_reg = zeros(15, 10); % 15 regions and for the csv file : columns = region type(1 or 2), indice_reg, Pext_vul, Pext_endan, Pext_crit, Pext_extinct

for m=3
    eco_mod=e_mod{m};

    if m==1
        load('Ntot_allcol_CMR.mat', 'Ntot_allcol')
        scenEXT =[1 2 4];
    elseif m==2
        load('Ntot_allcol_IPM.mat', 'Ntot_allcol')
        scenEXT =[1];
    else
        load('Ntot_allcol_SAT.mat', 'Ntot_allcol')
        scenEXT =[1];
    end


    % Calculate region sizes

    % REGIONS CCAMLR

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

    REGIONS1=[1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 1 2];
    numRegions1=max(REGIONS1);

    regionalSum1 = zeros(length(scenEXT)*nens*nsim, nt, numRegions1);

    x=1;
    for scenEXT = [1 2 4]
        for ens=1:nens
            for sim=1:nsim

                for re=1:numRegions1
                    regionalSum1(x, :, re) = sum(Ntot_allcol(x,:,REGIONS1==re),3);
                end %regions
                x=x+1;

            end %sim
        end %ens
    end

    regionalSum_med1=zeros(3,nt,numRegions1);
    for re=1:numRegions1
        regionalSum_med1(:,:,re) = quantile(regionalSum1(:,:,re),[0.05,0.5,0.95]);
    end


    % Genetic regions
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

    REGIONS2=[1 1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7];
    numRegions2=max(REGIONS2);

    regionalSum2 = zeros(length(scenEXT)*nens*nsim, nt, numRegions2);

    x=1;
    for scenEXT = [1 2 4]
        for ens=1:nens
            for sim=1:nsim

                for re=1:numRegions2
                    regionalSum2(x, :, re) = sum(Ntot_allcol(x,:,REGIONS2==re),3);
                end %regions
                x=x+1;

            end %sim
        end %ens
    end

    regionalSum_med2=zeros(3,nt,numRegions2);
    for re=1:numRegions2
        regionalSum_med2(:,:,re) = quantile(regionalSum2(:,:,re),[0.05,0.5,0.95]);
    end

    for i_re=1:2 %region type
        if i_re==1
            numRegions=numRegions1;
            regionalSum=regionalSum1;
            regionalSum_med=regionalSum_med1;
        else
            numRegions=numRegions2;
            regionalSum=regionalSum2;
            regionalSum_med=regionalSum_med2;
        end

        for re=1:numRegions
            display(re)

            pop_present = regionalSum_med(2,1,re);
            threshold_vul = 0.7 * pop_present; % Vulnerable criteria
            threshold_endan = 0.5 * pop_present; % Endangered
            threshold_critend = 0.2 * pop_present; % Critically endangered
            thr_ext = 0.05 * pop_present; % Quasi-extinct

            for year_eval=year_eval_list

                y_eval = year_eval - yearstart;
                if m==3
                    y_eval=y_eval-1;
                end

                % Probability to go under the threshold
                proba_tot=[0 0 0 0];
                j=1;

                for thr = [threshold_vul, threshold_endan, threshold_critend, thr_ext]
                    proba = 0;
                    under_thr = 0;
                    for i=1:length(regionalSum(:,1, re)); %total simulations
                        if regionalSum(i, y_eval,re) < thr;
                            under_thr = under_thr + 1;
                            proba = under_thr / length(regionalSum(:,1,re));
                        end
                    end
                    proba_tot(j)=proba;
                    j=j+1;
                end

                if i_re==2
                    reg=re+8;
                else
                    reg=re;
                end

                if year_eval==2073
                    %for the csv file : columns = indice_col, lat, long, Pext_vul, Pext_endan, Pext_crit, Pext_extinct
                    Pext_reg(reg, 3) = proba_tot(1);
                    Pext_reg(reg, 4) = proba_tot(2);
                    Pext_reg(reg, 5) = proba_tot(3);
                    Pext_reg(reg, 6) = proba_tot(4);
                else %2100
                    Pext_reg(reg, 7) = proba_tot(1);
                    Pext_reg(reg, 8) = proba_tot(2);
                    Pext_reg(reg, 9) = proba_tot(3);
                    Pext_reg(reg, 10) = proba_tot(4);
                end

            end %year eval (2068, 2100)

        end %region

    end % region type

    Pext_reg(1:8, 1) = 1; % type 1
    Pext_reg(9:15, 1) =2;
    Pext_reg(1:8, 2) = 1:8;
    Pext_reg(9:15, 2) = 1:7;

    if m==1
        Pext_reg_CMR = Pext_reg;
    elseif m==2
        Pext_reg_IPM = Pext_reg;
    else
        Pext_reg_SAT = Pext_reg;
    end

end %ecological model


%save('Pext_reg_CMR_mature', 'Pext_reg_CMR')
%save('Pext_reg_IPM_mature', 'Pext_reg_IPM')
%save('Pext_reg_SAT_mature', 'Pext_reg_SAT')

%% Transform in csv file

% load(sprintf('%s/Codes_EP/Figures/Fig_Pext_ecomod/Pext_col_CMR_mature.mat', ordi))
% load(sprintf('%s/Codes_EP/Figures/Fig_Pext_ecomod/Pext_col_IPM_mature.mat', ordi))
% load(sprintf('%s/Codes_EP/Figures/Fig_Pext_ecomod/Pext_col_SAT_mature.mat', ordi))

%colony size
N_avg_66=[3324	1788	3771	3198	4424	4614	12337	4307	6012	7365	4970	6831	4416	3252	217	4626	2789	112	70	69	1289	287	798	4035	2664	4267	70	359	70	70	3289	4083	1320	1553	70	70	70	9711	1578	70	5292	70	2212	7983	26477	14040	70	6348	1416	1071	22266	1753	70	616	3403	70	3684	3872	3046	70	1036	4314	70	797 70 70];

for m=[3] %ecological model

    % Create file with status for each col
    % Threshold over which col is likely to be at this level : 50-59%
    % Status: 1=OK, 2=Vulnerable, 3=Endangered, 4=Critically Endangered,
    % 5=Quasi extinct
    stat = 2:5;
    status=ones(15,1);

    if m==1
        Pext_reg = Pext_reg_CMR;
    elseif m==2
        Pext_reg = Pext_reg_IPM;
    else
        Pext_reg = Pext_reg_SAT;
    end


    % 2068
    for re=1:15
        for thr=[1:4]
            if Pext_reg(re, thr+2) > 0.40 %cf risk tolerance, RedListGuidelines
                status(re) = stat(thr);
            end
        end
    end

    Pext_reg(:,11) = status; %add to previous file

    % 2100
    for re=1:15
        for thr=[1:4]
            if Pext_reg(re, thr+6) > 0.40 %cf risk tolerance, RedListGuidelines
                status(re) = stat(thr);
            end
        end
    end
    Pext_reg(:,12) = status; %status in 2100


    csv_file = array2table(Pext_reg);
    csv_file.Properties.VariableNames = {'Reg_type', 'reg', 'Pext_vul_2068', 'Pext_endan_2073', 'Pext_crit_2073', 'Pext_ext_2073', 'Pext_vul_2100', 'Pext_endan_2100', 'Pext_crit_2100', 'Pext_ext_2100','status_2078', 'status_2100'};

    if m==1
        writetable(csv_file, 'Pext_reg_CMR_mature.csv')
    elseif m==2
        writetable(csv_file, 'Pext_reg_IPM_mature.csv')
    else
        writetable(csv_file, 'Pext_reg_SAT_mature.csv')
    end

end %eco mod
