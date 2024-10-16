% Meta-pop model for 66 colonies and new SIC from CMIP6 - LE2 CESM2
% from 1900 to 2100 which start with 66 colonies!

clear all; clc; %close all;

rand('state',sum(100.*clock)) % reset of the rand function

ordi = "/Users/%s/Desktop";
ordi = '/OneDrive/Bureau';
ordi = '/Users/cepar/Desktop';
ordi = 'D:/Documents/alice';

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% Initialization
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

load(sprintf('%s/Codes_EP/Codes_Stef/mat_data/infos_models.mat', ordi))

% number of simulation
nsimulation=1000; %8000 %should be smaller than 8000 because 8000 posterior for dispersive data

dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_Stef/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));

% SAT POP average *2 from BILGECAN
%BESATav=[5506	2956	6495	5138	7039	7604	19697	6203	9581	12568	8377	11430	6706	4974	355	7102	3902	247	180	113	1914	394	1264	5641	3375	7467	180	569	180	180	5642	5133	1993	2454	180	180	180	14963	2829	180	8201	180	3480	13859	46539	23393	180	10833	2370	1835	36018	2578	180	1141	4746	180	5344	5251	4525	180	1560	6528	180	1141	180	180]';
BESATav=2*[3324	1788	3771	3198	4424	4614	12337	4307	6012	7365	4970	6831	4416	3252	217	4626	2789	112	70	69	1289	287	798	4035	2664	4267	70	359	70	70	3289	4083	1320	1553	70	70	70	9711	1578	70	5292	70	2212	7983	26477	14040	70	6348	1416	1071	22266	1753	70	616	3403	70	3684	3872	3046	70	1036	4314	70	797 70 70]';

% N SAD calculated from previous simulation
NSAD=[5161	7579	13474	4667	6048	7157	16964	7339	9425	10798	8597	11823	6955	5124	351	6334	2676	136	114	74	1420	888	2826	2779	3580	7808	294	407	227	405	6527	3582	2084	2338	388	175	179	14986	2896	405	7548	405	3489	21999	31240	9456	700	11122	5372	4051	4996	2772	607	1407	3379	593	2056	2040	2123	555	1609	2273	654	2086	610	556]';
SAD=NSAD./sum(NSAD);
BE=SAD.*dataOBSSAT_bilgecan(1,2);

% BE=1 for 2009 proj; BE=1.3 for IPM 1920; BE=0.85 for CMR 1900
BE=0.85*BE;
BE=2.2*BE % for CMR mature
%BE=1.3*BE; %IPM

%BE=1.5*BE; %IPM mature

% Carrying capacity
%load('/Users/fledgelab/Desktop/Codes_EP/Codes_Stef/Kmax_50ens.mat') %Kmax

Kmax=4.5;
K  = Kmax*BE; % carrying capacity

%Load estimated dispersal parameters
% D  = D_post;   % Mean distance dispersal
% Pm = Pm_post; % Emigration
load(sprintf('%s/Codes_EP/Codes_Stef/post_proc_EP_7pm_informed_jd0.mat', ordi),'D_post','Pm_post');

% Dispersal parameters
choice = 2;
% choice 1) Run with median of mean dispersal and emigration
% choice 2) Run with each mean dispersal and emigration (WARNING a lot of computation)

if (choice == 1)
    D_post = median(D_post);
    Pm_post = median(Pm_post);
end

nD = length(D_post); % there are 8000 posteriors in post_proc_EP_7pm_informed_jd0.mat

% so we reduce to nsimREDUIT to reduce computing time
D_post=D_post(randi(nD,nsimulation,1));
Pm_post=Pm_post(randi(nD,nsimulation,1),:);
nD = length(D_post);

% file with information for blinking colonies
load(sprintf('%s/Codes_EP/Codes_Stef/Blinking.mat', ordi))

% Regions for dispersal
i_subgroup = [1,5,8,22,26,44,52,67];
% 1-4 Snowhill to Smith
% 5-7 Gould Bay to Halley Bay
% 8-21 Dawson to Kloa Point
% 22-25 Fold Island to Cape Darnley
% 26-43 Amanda Bay Point Geologie Davis Bay
% 44-51 Ross Sea
% 52-66 Amundsen Bellington
n_subgroup = length(i_subgroup)-1;

%%%% calcule la distance cotiere entre chaque colonnie en km
lat  = xlsread(sprintf('%s/Codes_EP/Codes_Stef/empe_sitesNewNB.xlsx', ordi),'C2:D67'); %'D2:D67'
long = xlsread(sprintf('%s/Codes_EP/Codes_Stef/empe_sitesNewNB.xlsx', ordi),'D2:E67'); %'E2:D67'
ncol = length(lat);  % Nb of colonies

d    = zeros(1,ncol);
lat  = lat*pi/180;
long = long*pi/180;
for c = 1:ncol-1
    d(c) = 6374.8925* acos(...
        sin(lat(c))*sin(lat(c+1))+...
        cos(lat(c))*cos(lat(c+1))*cos(long(c+1)-long(c)));
end
d(end) = 6374.8925* acos(sin(lat(end))*sin(lat(1))+...
    cos(lat(end))*cos(lat(1))*cos(long(1)-long(end)));
dm = dispersion(d,ncol);   % Distance between colonies

% Growth Function

% With Density dependance K
f=@(N,r) N.*exp(log(1+r).*(1-N./K)).*(r>0) + (1+r).*N.*(r<=0);
% g: Individual growth rate f(N)/N
g = @(N,r) exp(log(1+r).*(1-N./K)).*(r>0) + (1+r).*(r<=0);


%% SIMULATIONS FOR 7 MODELS

% Load Growth rate from one of the 6 models

mod_pop= ["Stef", "IPM", "Sat"];
e_mod=2;
eco_mod=mod_pop(e_mod);

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-Paris2"];

if e_mod==1 %model CMR
    scen_list = [1 2 4];
    mod_list=[2:3 7];
else % other models don't have scenEXT
    scen_list = [1];
    mod_list=6;
end


for m=mod_list %climate model
    mod=model(m)

    if e_mod==1
        yearstart = info_models(m,3) ; yearstop = info_models(m,4);
    elseif e_mod==2
        %yearstart = 2021; yearstop=2100;
        yearstart = 2015; yearstop=2100;
        yearstart=1921; yearstop=2100;
        yearstart=2009;
    else
        yearstart = 1909; yearstop=2099;
    end

    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    ncol=66;
    nens= info_models(m,5);
    nsim=100;

    for scenEXT=scen_list
        display(scenEXT)

        % SIMULATIONS

        jd = 0; % Choose 1) jd=0 Random search strategy; 2) jd=1 Oriented search strategy

        for ens=1:nens

            display(ens)
            R=zeros(ncol, nt, nsim);

            % Import data of all the colonies
            parfor col=1:ncol
                if e_mod==1
                    file_name=sprintf('%s/Codes_EP/GR_results_tot/GR_results_%s/%s/GR_%s_%s_ens%d_scen%d_col%d.txt', ordi, eco_mod, mod, eco_mod, mod, ens, scenEXT, col);

                else
                    file_name=sprintf('%s/Codes_EP/GR_results_tot/GR_results_%s/LE_CESM2/GR_%s_CESM2_ens%d_col%d.txt', ordi, eco_mod, eco_mod, ens, col);
                end

                R(col, :, :) = readmatrix(file_name); %(nt, nsim)
            end

            N_allsim=zeros(ncol, nt, nsim);

            for i=1:nsim

                r = R(:,:,i); %(ncol, nt)
                r  = log(r); %(nt, ncol)

                A_col = [];

                % Matrix of connexion %
                %%%%%%%%%%%%%%%%%%%%%%%
                D = D_post(i); % Mean distance dispersal
                % Dispersal kernel        % Rem: k(0)=0, the probability of staying where we are, is zero
                k = @(x) (x<D)/D .*(x>0); % Uniform kernel with D = mean distance dispersal

                % Matrix of connexion without renormalisation
                dmat = k(dm);    % Unif kernel
                % Renormalization of the connexion matrix %
                den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
                Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
                Dmat = dmat*diag(1./Den);

                pm = [];
                parfor j = 1:n_subgroup
                    Is = i_subgroup(j+1)-i_subgroup(j);
                    pm = [pm;ones(Is,1)*Pm_post(i,j)];
                end

                % Movement scenario %
                %%%%%%%%%%%%%%%%%%%%%
                rM = 0.25;      % maximal death rate
                rm = rM*(1-pm);
                if (pm==1)
                    movement=@(r)  (r<0); %HIGH
                else
                    movement=@(r)  (r<-rm) - (r./rm).*(-rm<=r).*(r<0); %HIGH
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Initialization of size of colonies %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Nntot = zeros(nt, nsim);
                Nncol = zeros(ncol,nt); % n(i,t,j)=  # individuals in colony j at time t in simu i
                Nncol(:,1) = BE; %no uncertainties in BE
                Nnew = Nncol(:,1);

                % Time loop forward
                % Inidividuals
                t = yearstart;
                a_col = [];
                it = 1;
                while (t<yearstop)

                    Nold = Nnew;

                    %%% Reproduction at time t+1
                    %%% ###################################################################################################################################
                    fN    = f(Nold,r(:,it));
                    rstar = g(Nold,r(:,it))-1; % Effective growth rate at time t when blicking happens
                    %%% ###################################################################################################################################


                    %%%% Construction of informed dispersal matrix
                    mat_r = meshgrid(rstar)';
                    [rmax_vec, ind_col_sel] = max(mat_r.*(dm<D)+((dm<D)-1)); % ((dm<D)-1) to eliminate zeros and the max provieds the largest negative growth rate
                    Mat_select = zeros(ncol,ncol);
                    for ii=1:ncol
                        Mat_select(ind_col_sel(ii),ii) = 1;
                    end
                    rmax_mat = meshgrid(rmax_vec);

                    if (jd==1)
                        % Oriented dispersal
                        disp_sel = Mat_select; display(scenEXT)
                    else
                        % Random Dispersal
                        disp_sel = Dmat;
                    end

                    %%% Exploration
                    lambda = movement(rstar); % mean proportion of individual that leaves colonies

                    % Blinking of colonies (everybody dipserse and nobody arrive directly) ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    %%% ###################################################################################################################################

                    iblink = Blinking(Blinking(:,2)==(t+1),1); % Give the blinking colonies at time t
                    s_disp_blink = 0.5*(0.9438+0.9253);   % survival rate after living a blinking colony / ADD a 10% decrease from what is optimal**
                    pm_blink     = 0.9;                       % emigration rate of blinking colonies ****************************

                    lambda(iblink) = s_disp_blink*pm_blink;    % proportion of emigrants that survive
                    disp_sel(iblink,:) = 0;                    % no emigration to the blinking colony


                    den_disp = sum(disp_sel,1)';
                    Disp_sel = disp_sel-diag(den_disp);        % pour assurer que la somme des colonnes vaut 0

                    NI = (Disp_sel*(lambda.*fN));              % emigrants
                    NO = fN;                                   % newborn and alive individuals
                    N_E = NI + NO;

                    % blinking ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    N_E(iblink) = (1-pm_blink)*fN(iblink);     % survivers in blinking colonies x

                    %%% ###################################################################################################################################


                    %           Random settlment
                    full    = (N_E>K);                     % localization of overcrowded colonies
                    Overpgn = (N_E-K).*full;               % surplus de manchots
                    Overpgn_mov = Dmat*Overpgn;            % surplus re-distribue en fonction du kernel des distance
                    N_RS    = N_E - Overpgn + Overpgn_mov; % variable intermediare: locaux+ surplus dans chaque col
                    Nnew = N_RS;

                    %           Nnew = N_E; % informed dispersal
                    %%% ###################################################################################################################################

                    Nncol(:,it+1) = Nnew;
                    it = it+1;
                    a_col = [a_col,lambda];
                    t = t+1;

                end % year

                N_allsim(:, :, i) = Nncol; % all sim for one colony

            end % sim

            % SAVE one file txt per model, scen, colony, ens
            
            N_allsim_mature=zeros(ncol,nt,nsim);
            
            for col=1:ncol
                if e_mod==1
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/MAT_results_Stef/%s/MAT_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);
                else
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/MAT_results_IPM/%s/MAT_IPM_CESM2_ens%d_col%d.txt', ordi, mod, ens, col);
                end

                MAT=readmatrix(filename2);

                N=permute(N_allsim(col,:,:), [2 3 1]); % nt, nsim TOTAL POPULATION
                GR = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
                GR_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER

                for t=1:length(N(:,1))
                    for s=1:length(N(1,:))
                        N_allsim_mature(col,t,s) = N(t,s)*(MAT(t,s));
                    end
                    if t>1
                        GR(t,:)=N(t,:)./N(t-1,:);
                        GR_mature(t,:)=N_allsim_mature(col,t,:)./N_allsim_mature(col,t-1,:);
                    end
                end

                N_mature=permute(N_allsim(col, :,:), [2 3 1]);
               
                % Files for BE=1.1 BE and Kmax=3.15
                if e_mod==1 %CMR, we precise the scenEXT in the name of the file
                    file_name =sprintf("%s/Codes_EP/N_results_tot/N_results_CMR/%s/N_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/N_results_CMR/%s/N_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/GR_results_CMR/%s/GR_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/GR_results_CMR/%s/GR_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                else %other models
                    file_name =sprintf("%s/Codes_EP/N_results_tot/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/GR_results_%s/%s/GR_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/GR_results_%s/%s/GR_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                end

                writematrix(N, file_name) %(nt, nsim)
                writematrix(N_mature, file_name2) %(nt, nsim)
                writematrix(GR, file_name3) %(nt, nsim)
                writematrix(GR_mature, file_name4) %(nt, nsim)
            end %col

            % SAVE total population for each model, scen, ens

            Nntot=sum(permute(N_allsim, [2 3 1]), 3); % nt, nsim
            Nntot_mature=sum(permute(N_allsim_mature, [2 3 1]), 3); % nt, nsim
            GRtot = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
            GRtot_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER
            
            for t=2:nt
                GRtot(t,:)=Nntot(t,:)./Nntot(t-1,:);
                GRtot_mature(t,:)=Nntot_mature(t,:)./Nntot_mature(t-1,:);
            end

            if e_mod==1
                file_name =sprintf("%s/Codes_EP/N_results_tot/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);

            else
                file_name =sprintf("%s/Codes_EP/N_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/GR_results_%s/GRtot_%s/%s/GRtot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/GR_results_%s/GRtot_%s/%s/GRtot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
            end

            writematrix(Nntot, file_name)
            writematrix(Nntot_mature, file_name2)
            writematrix(GRtot, file_name3)
            writematrix(GRtot_mature, file_name4)

        end %ens

    end %scenEXT

end %climate model



%% MEAN OF ENSEMBLES FOR CMR MODEL

% Load Growth rate from one of the 6 models

mod_pop= ["Stef", "IPM", "Sat"];
e_mod=1;
eco_mod=mod_pop(e_mod);

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-Paris2"];

if e_mod==1 %model CMR
    scen_list = [1 2 4];
    mod_list=[1:3 6 7];
else % other models don't have scenEXT
    scen_list = [1];
    mod_list=6;
end


for m=mod_list %climate model
    mod=model(m)

    if e_mod==1
        yearstart = info_models(m,3); yearstop = info_models(m,4);
    elseif e_mod==2
        %yearstart = 2021; yearstop=2100;
        yearstart = 2015; yearstop=2100;
        yearstart = 1920;
    else
        yearstart = 1909; yearstop=2099;
    end

    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    ncol=66;
    nens= info_models(m,5);
    nsim=100;

    for scenEXT=scen_list
        display(scenEXT)

        % SIMULATIONS

        jd = 0; % Choose 1) jd=0 Random search strategy; 2) jd=1 Oriented search strategy

        %for ens=1:nens

        R=zeros(ncol, nt, nsim);

            % Import data of all the colonies
            parfor col=1:ncol
                if e_mod==1
                    file_name=sprintf('%s/Codes_EP/GR_results_tot/mean_ensemble/GR_results_%s/%s/GRmean_%s_%s_scen%d_col%d.txt', ordi, eco_mod, mod, eco_mod, mod, scenEXT, col);
                else
                   file_name=sprintf('%s/Codes_EP/GR_results_tot/mean_ensemble/GR_results_%s/CESM2/GR_%s_CESM2_ens%d_col%d.txt', ordi, eco_mod, eco_mod, ens, col);
                end

                R(col, :, :) = readmatrix(file_name); %(nt, nsim)
            end

            N_allsim=zeros(ncol, nt, nsim);

            for i=1:nsim

                r = R(:,:,i); %(ncol, nt)
                r  = log(r); %(nt, ncol)

                A_col = [];

                % Matrix of connexion %
                %%%%%%%%%%%%%%%%%%%%%%%
                D = D_post(i); % Mean distance dispersal
                % Dispersal kernel        % Rem: k(0)=0, the probability of staying where we are, is zero
                k = @(x) (x<D)/D .*(x>0); % Uniform kernel with D = mean distance dispersal

                % Matrix of connexion without renormalisation
                dmat = k(dm);    % Unif kernel
                % Renormalization of the connexion matrix %
                den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
                Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
                Dmat = dmat*diag(1./Den);

                pm = [];
                parfor j = 1:n_subgroup
                    Is = i_subgroup(j+1)-i_subgroup(j);
                    pm = [pm;ones(Is,1)*Pm_post(i,j)];
                end

                % Movement scenario %
                %%%%%%%%%%%%%%%%%%%%%
                rM = 0.25;      % maximal death rate
                rm = rM*(1-pm);
                if (pm==1)
                    movement=@(r)  (r<0); %HIGH
                else
                    movement=@(r)  (r<-rm) - (r./rm).*(-rm<=r).*(r<0); %HIGH
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Initialization of size of colonies %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Nntot = zeros(nt, nsim);
                Nncol = zeros(ncol,nt); % n(i,t,j)=  # individuals in colony j at time t in simu i
                Nncol(:,1) = BE; %no uncertainties in BE
                Nnew = Nncol(:,1);

                % Time loop forward
                % Inidividuals
                t = yearstart;
                a_col = [];
                it = 1;
                while (t<yearstop)

                    Nold = Nnew;

                    %%% Reproduction at time t+1
                    %%% ###################################################################################################################################
                    fN    = f(Nold,r(:,it));
                    rstar = g(Nold,r(:,it))-1; % Effective growth rate at time t when blicking happens
                    %%% ###################################################################################################################################


                    %%%% Construction of informed dispersal matrix
                    mat_r = meshgrid(rstar)';
                    [rmax_vec, ind_col_sel] = max(mat_r.*(dm<D)+((dm<D)-1)); % ((dm<D)-1) to eliminate zeros and the max provieds the largest negative growth rate
                    Mat_select = zeros(ncol,ncol);
                    for ii=1:ncol
                        Mat_select(ind_col_sel(ii),ii) = 1;
                    end
                    rmax_mat = meshgrid(rmax_vec);

                    if (jd==1)
                        % Oriented dispersal
                        disp_sel = Mat_select; display(scenEXT)
                    else
                        % Random Dispersal
                        disp_sel = Dmat;
                    end

                    %%% Exploration
                    lambda = movement(rstar); % mean proportion of individual that leaves colonies

                    % Blinking of colonies (everybody dipserse and nobody arrive directly) ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    %%% ###################################################################################################################################

                    iblink = Blinking(Blinking(:,2)==(t+1),1); % Give the blinking colonies at time t
                    s_disp_blink = 0.5*(0.9438+0.9253);   % survival rate after living a blinking colony / ADD a 10% decrease from what is optimal**
                    pm_blink     = 0.9;                       % emigration rate of blinking colonies ****************************

                    lambda(iblink) = s_disp_blink*pm_blink;    % proportion of emigrants that survive
                    disp_sel(iblink,:) = 0;                    % no emigration to the blinking colony


                    den_disp = sum(disp_sel,1)';
                    Disp_sel = disp_sel-diag(den_disp);        % pour assurer que la somme des colonnes vaut 0

                    NI = (Disp_sel*(lambda.*fN));              % emigrants
                    NO = fN;                                   % newborn and alive individuals
                    N_E = NI + NO;

                    % blinking ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    N_E(iblink) = (1-pm_blink)*fN(iblink);     % survivers in blinking colonies x

                    %%% ###################################################################################################################################


                    %           Random settlment
                    full    = (N_E>K);                     % localization of overcrowded colonies
                    Overpgn = (N_E-K).*full;               % surplus de manchots
                    Overpgn_mov = Dmat*Overpgn;            % surplus re-distribue en fonction du kernel des distance
                    N_RS    = N_E - Overpgn + Overpgn_mov; % variable intermediare: locaux+ surplus dans chaque col
                    Nnew = N_RS;

                    %           Nnew = N_E; % informed dispersal
                    %%% ###################################################################################################################################

                    Nncol(:,it+1) = Nnew;
                    it = it+1;
                    a_col = [a_col,lambda];
                    t = t+1;

                end % year

                N_allsim(:, :, i) = Nncol; % all sim for one colony

            end % sim

            % SAVE one file txt per model, scen, colony, ens

            N_allsim_mature=zeros(ncol,nt,nsim);

            for col=1:ncol
                if e_mod==1
                filename2=sprintf('%s/Codes_EP/MAT_results_tot/mean_ensemble/MAT_results_Stef/%s/MAT_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);
                else
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/mean_ensemble/MAT_results_Stef/%s/MAT_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);
                end
                MAT=readmatrix(filename2);

                N=permute(N_allsim(col,:,:), [2 3 1]); % nt, nsim TOTAL POPULATION
                GR = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
                GR_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER

                for t=1:length(N(:,1))
                    for s=1:length(N(1,:))
                        N_allsim_mature(col,t,s) = N(t,s)*(MAT(t,s));
                    end
                    if t>1
                        GR(t,:)=N(t,:)./N(t-1,:);
                        GR_mature(t,:)=N_allsim_mature(col,t,:)./N_allsim_mature(col,t-1,:);
                    end
                end

                N_mature=permute(N_allsim(col, :,:), [2 3 1]);

                N=permute(N_allsim(col,:,:), [2 3 1]); % nt, nsim

                % if e_mod==1 %CMR, we precise the scenEXT in the name of the file
                %     file_name =sprintf("/Users/fledgelab/Desktop/Codes_EP/N_results_tot/N_results_CMR/%s/N_CMR_%s_ens%d_scen%d_col%d.txt", mod, mod, ens, scenEXT, col);
                % else %other models
                %     file_name =sprintf("/Users/fledgelab/Desktop/Codes_EP/N_results_tot/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", eco_mod, mod, eco_mod, mod, ens, col);
                % end

                % Files for BE=1.1 BE and Kmax=3.15
                if e_mod==1 %CMR, we precise the scenEXT in the name of the file
                    file_name =sprintf("%s/Codes_EP/N_results_tot/mean_ensemble/N_results_CMR/%s/Nmean_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/mean_ensemble/N_results_CMR/%s/N_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/mean_ensemble/GR_results_CMR/%s/GR_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/mean_ensemble/GR_results_CMR/%s/GR_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    
                end
 
                writematrix(N, file_name) %(nt, nsim)
                writematrix(N_mature, file_name2) %(nt, nsim)
                writematrix(GR, file_name3) %(nt, nsim)
                writematrix(GR_mature, file_name4) %(nt, nsim)

            end %col

            % SAVE total population for each model, scen, ens

            Nntot=sum(permute(N_allsim, [2 3 1]), 3); % nt, nsim
            Nntot_mature=sum(permute(N_allsim_mature, [2 3 1]), 3); % nt, nsim
            GRtot = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
            GRtot_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER
            
            for t=2:nt
                GRtot(t,:)=Nntot(t,:)./Nntot(t-1,:);
                GRtot_mature(t,:)=Nntot_mature(t,:)./Nntot_mature(t-1,:);
            end

            if e_mod==1
                %file_name =sprintf("/Users/fledgelab/Desktop/Codes_EP/N_results_tot/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_ens%d_scen%d.txt", mod, mod, ens, scenEXT);
                file_name =sprintf("%s/Codes_EP/N_results_tot/mean_ensemble/N_results_CMR/Ntot_CMR/%s/Ntotmean_CMR_%s_scen%d.txt", ordi, mod, mod, scenEXT);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/mean_ensemble/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/mean_ensemble/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/mean_ensemble/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);

            end

            writematrix(Nntot, file_name)
            writematrix(Nntot_mature, file_name2)
            writematrix(GRtot, file_name3)
            writematrix(GRtot_mature, file_name4)

        %end %ens

    end %scenEXT

end %climate model




%% 2009

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% Initialization
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

load(sprintf('%s/Codes_EP/Codes_Stef/mat_data/infos_models.mat', ordi))

% number of simulation
nsimulation=1000; %8000 %should be smaller than 8000 because 8000 posterior for dispersive data

dataOBSSAT_bilgecan=readmatrix(sprintf('%s/Codes_EP/Codes_Stef/DataOBS_SAT_BILGECAN/N_glb.csv', ordi));

% SAT POP average *2 from BILGECAN
%BESATav=[5506	2956	6495	5138	7039	7604	19697	6203	9581	12568	8377	11430	6706	4974	355	7102	3902	247	180	113	1914	394	1264	5641	3375	7467	180	569	180	180	5642	5133	1993	2454	180	180	180	14963	2829	180	8201	180	3480	13859	46539	23393	180	10833	2370	1835	36018	2578	180	1141	4746	180	5344	5251	4525	180	1560	6528	180	1141	180	180]';
BESATav=2*[3324	1788	3771	3198	4424	4614	12337	4307	6012	7365	4970	6831	4416	3252	217	4626	2789	112	70	69	1289	287	798	4035	2664	4267	70	359	70	70	3289	4083	1320	1553	70	70	70	9711	1578	70	5292	70	2212	7983	26477	14040	70	6348	1416	1071	22266	1753	70	616	3403	70	3684	3872	3046	70	1036	4314	70	797 70 70]';

% N SAD calculated from previous simulation
NSAD=[5161	7579	13474	4667	6048	7157	16964	7339	9425	10798	8597	11823	6955	5124	351	6334	2676	136	114	74	1420	888	2826	2779	3580	7808	294	407	227	405	6527	3582	2084	2338	388	175	179	14986	2896	405	7548	405	3489	21999	31240	9456	700	11122	5372	4051	4996	2772	607	1407	3379	593	2056	2040	2123	555	1609	2273	654	2086	610	556]';
SAD=NSAD./sum(NSAD);
BE=SAD.*dataOBSSAT_bilgecan(1,2);


BE=1.7*BE; %so that mature fits obsCMR

BE=1.2*BE %IPM

% Carrying capacity
%load('/Users/fledgelab/Desktop/Codes_EP/Codes_Stef/Kmax_50ens.mat') %Kmax

Kmax=4.5;
K  = Kmax*BE; % carrying capacity

%Load estimated dispersal parameters
% D  = D_post;   % Mean distance dispersal
% Pm = Pm_post; % Emigration
load(sprintf('%s/Codes_EP/Codes_Stef/post_proc_EP_7pm_informed_jd0.mat', ordi),'D_post','Pm_post');

% Dispersal parameters
choice = 2;
% choice 1) Run with median of mean dispersal and emigration
% choice 2) Run with each mean dispersal and emigration (WARNING a lot of computation)

if (choice == 1)
    D_post = median(D_post);
    Pm_post = median(Pm_post);
end

nD = length(D_post); % there are 8000 posteriors in post_proc_EP_7pm_informed_jd0.mat

% so we reduce to nsimREDUIT to reduce computing time
D_post=D_post(randi(nD,nsimulation,1));
Pm_post=Pm_post(randi(nD,nsimulation,1),:);
nD = length(D_post);

% file with information for blinking colonies
load(sprintf('%s/Codes_EP/Codes_Stef/Blinking.mat', ordi))

% Regions for dispersal
i_subgroup = [1,5,8,22,26,44,52,67];
% 1-4 Snowhill to Smith
% 5-7 Gould Bay to Halley Bay
% 8-21 Dawson to Kloa Point
% 22-25 Fold Island to Cape Darnley
% 26-43 Amanda Bay Point Geologie Davis Bay
% 44-51 Ross Sea
% 52-66 Amundsen Bellington
n_subgroup = length(i_subgroup)-1;

%%%% calcule la distance cotiere entre chaque colonnie en km
lat  = xlsread(sprintf('%s/Codes_EP/Codes_Stef/empe_sitesNewNB.xlsx', ordi),'C2:D67'); %'D2:D67'
long = xlsread(sprintf('%s/Codes_EP/Codes_Stef/empe_sitesNewNB.xlsx', ordi),'D2:E67'); %'E2:D67'
ncol = length(lat);  % Nb of colonies

d    = zeros(1,ncol);
lat  = lat*pi/180;
long = long*pi/180;
for c = 1:ncol-1
    d(c) = 6374.8925* acos(...
        sin(lat(c))*sin(lat(c+1))+...
        cos(lat(c))*cos(lat(c+1))*cos(long(c+1)-long(c)));
end
d(end) = 6374.8925* acos(sin(lat(end))*sin(lat(1))+...
    cos(lat(end))*cos(lat(1))*cos(long(1)-long(end)));
dm = dispersion(d,ncol);   % Distance between colonies

% Growth Function

% With Density dependance K
f=@(N,r) N.*exp(log(1+r).*(1-N./K)).*(r>0) + (1+r).*N.*(r<=0);
% g: Individual growth rate f(N)/N
g = @(N,r) exp(log(1+r).*(1-N./K)).*(r>0) + (1+r).*(r<=0);


%% SIMULATIONS FOR 7 MODELS

% Load Growth rate from one of the 6 models

mod_pop= ["Stef", "IPM", "Sat"];
e_mod=2;
eco_mod=mod_pop(e_mod);

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-Paris2"];

if e_mod==1 %model CMR
    scen_list = [1 2 4];
    mod_list=[1:3 6 7];
else % other models don't have scenEXT
    scen_list = [1];
    mod_list=6;
end


for m=mod_list %climate model
    mod=model(m)

    if e_mod==1
        yearstart = 2009 ; yearstop = info_models(m,4);
    elseif e_mod==2
        yearstart = 2009; yearstop=2100;
    else
        yearstart = 1909; yearstop=2099;
    end

    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    ncol=66;
    nens= info_models(m,5);
    nsim=100;

    for scenEXT=scen_list
        display(scenEXT)

        % SIMULATIONS

        jd = 0; % Choose 1) jd=0 Random search strategy; 2) jd=1 Oriented search strategy

        for ens=1:nens

            display(ens)
            R=zeros(ncol, nt, nsim);

            % Import data of all the colonies
            parfor col=1:ncol
                if e_mod==1
                    file_name=sprintf('%s/Codes_EP/GR_results_tot/2009/GR_results_%s/%s/GR_%s_%s_ens%d_scen%d_col%d.txt', ordi, eco_mod, mod, eco_mod, mod, ens, scenEXT, col);

                else
                    file_name=sprintf('%s/Codes_EP/GR_results_tot/2009/GR_results_%s/CESM2/GR_%s_CESM2_ens%d_col%d.txt', ordi, eco_mod, eco_mod, ens, col);
                end

                R(col, :, :) = readmatrix(file_name); %(nt, nsim)
            end

            N_allsim=zeros(ncol, nt, nsim);

            for i=1:nsim

                r = R(:,:,i); %(ncol, nt)
                r  = log(r); %(nt, ncol)

                A_col = [];

                % Matrix of connexion %
                %%%%%%%%%%%%%%%%%%%%%%%
                D = D_post(i); % Mean distance dispersal
                % Dispersal kernel        % Rem: k(0)=0, the probability of staying where we are, is zero
                k = @(x) (x<D)/D .*(x>0); % Uniform kernel with D = mean distance dispersal

                % Matrix of connexion without renormalisation
                dmat = k(dm);    % Unif kernel
                % Renormalization of the connexion matrix %
                den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
                Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
                Dmat = dmat*diag(1./Den);

                pm = [];
                parfor j = 1:n_subgroup
                    Is = i_subgroup(j+1)-i_subgroup(j);
                    pm = [pm;ones(Is,1)*Pm_post(i,j)];
                end

                % Movement scenario %
                %%%%%%%%%%%%%%%%%%%%%
                rM = 0.25;      % maximal death rate
                rm = rM*(1-pm);
                if (pm==1)
                    movement=@(r)  (r<0); %HIGH
                else
                    movement=@(r)  (r<-rm) - (r./rm).*(-rm<=r).*(r<0); %HIGH
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Initialization of size of colonies %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Nntot = zeros(nt, nsim);
                Nncol = zeros(ncol,nt); % n(i,t,j)=  # individuals in colony j at time t in simu i
                Nncol(:,1) = BE; %no uncertainties in BE
                Nnew = Nncol(:,1);

                % Time loop forward
                % Inidividuals
                t = yearstart;
                a_col = [];
                it = 1;
                while (t<yearstop)

                    Nold = Nnew;

                    %%% Reproduction at time t+1
                    %%% ###################################################################################################################################
                    fN    = f(Nold,r(:,it));
                    rstar = g(Nold,r(:,it))-1; % Effective growth rate at time t when blicking happens
                    %%% ###################################################################################################################################


                    %%%% Construction of informed dispersal matrix
                    mat_r = meshgrid(rstar)';
                    [rmax_vec, ind_col_sel] = max(mat_r.*(dm<D)+((dm<D)-1)); % ((dm<D)-1) to eliminate zeros and the max provieds the largest negative growth rate
                    Mat_select = zeros(ncol,ncol);
                    for ii=1:ncol
                        Mat_select(ind_col_sel(ii),ii) = 1;
                    end
                    rmax_mat = meshgrid(rmax_vec);

                    if (jd==1)
                        % Oriented dispersal
                        disp_sel = Mat_select; display(scenEXT)
                    else
                        % Random Dispersal
                        disp_sel = Dmat;
                    end

                    %%% Exploration
                    lambda = movement(rstar); % mean proportion of individual that leaves colonies

                    % Blinking of colonies (everybody dipserse and nobody arrive directly) ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    %%% ###################################################################################################################################

                    iblink = Blinking(Blinking(:,2)==(t+1),1); % Give the blinking colonies at time t
                    s_disp_blink = 0.5*(0.9438+0.9253);   % survival rate after living a blinking colony / ADD a 10% decrease from what is optimal**
                    pm_blink     = 0.9;                       % emigration rate of blinking colonies ****************************

                    lambda(iblink) = s_disp_blink*pm_blink;    % proportion of emigrants that survive
                    disp_sel(iblink,:) = 0;                    % no emigration to the blinking colony


                    den_disp = sum(disp_sel,1)';
                    Disp_sel = disp_sel-diag(den_disp);        % pour assurer que la somme des colonnes vaut 0

                    NI = (Disp_sel*(lambda.*fN));              % emigrants
                    NO = fN;                                   % newborn and alive individuals
                    N_E = NI + NO;

                    % blinking ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    N_E(iblink) = (1-pm_blink)*fN(iblink);     % survivers in blinking colonies x

                    %%% ###################################################################################################################################


                    %           Random settlment
                    full    = (N_E>K);                     % localization of overcrowded colonies
                    Overpgn = (N_E-K).*full;               % surplus de manchots
                    Overpgn_mov = Dmat*Overpgn;            % surplus re-distribue en fonction du kernel des distance
                    N_RS    = N_E - Overpgn + Overpgn_mov; % variable intermediare: locaux+ surplus dans chaque col
                    Nnew = N_RS;

                    %           Nnew = N_E; % informed dispersal
                    %%% ###################################################################################################################################

                    Nncol(:,it+1) = Nnew;
                    it = it+1;
                    a_col = [a_col,lambda];
                    t = t+1;

                end % year

                N_allsim(:, :, i) = Nncol; % all sim for one colony

            end % sim

            % SAVE one file txt per model, scen, colony, ens
            
            N_allsim_mature=zeros(ncol,nt,nsim);
            
            for col=1:ncol

                if e_mod==1
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/MAT_results_Stef/%s/MAT_Stef_%s_ens%d_scen%d_col%d.txt', ordi, mod, mod, ens, scenEXT, col);
                else
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/MAT_results_IPM/CESM2/MAT_IPM_CESM2_ens%d_col%d.txt', ordi, ens, col);
                end
                MAT=readmatrix(filename2);

                N=permute(N_allsim(col,:,:), [2 3 1]); % nt, nsim TOTAL POPULATION
                GR = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
                GR_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER

                for t=1:length(N(:,1))
                    for s=1:length(N(1,:))
                        N_allsim_mature(col,t,s) = N(t,s)*(MAT(t,s));
                    end
                    if t>1
                        GR(t,:)=N(t,:)./N(t-1,:);
                        GR_mature(t,:)=N_allsim_mature(col,t,:)./N_allsim_mature(col,t-1,:);
                    end
                end

                N_mature=permute(N_allsim(col, :,:), [2 3 1]);
               
                % Files for BE=1.1 BE and Kmax=3.15
                if e_mod==1 %CMR, we precise the scenEXT in the name of the file
                    file_name =sprintf("%s/Codes_EP/N_results_tot/2009/N_results_CMR/%s/N_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_CMR/%s/N_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/GR_results_CMR/%s/GR_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/GR_results_CMR/%s/GR_CMR_%s_ens%d_scen%d_col%d.txt", ordi, mod, mod, ens, scenEXT, col);
                else %other models
                    file_name =sprintf("%s/Codes_EP/N_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/%s/N_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/GR_results_%s/%s/GR_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/GR_results_%s/%s/GR_%s_%s_ens%d_col%d.txt", ordi, eco_mod, mod, eco_mod, mod, ens, col);
                end

                writematrix(N, file_name) %(nt, nsim)
                writematrix(N_mature, file_name2) %(nt, nsim)
                writematrix(GR, file_name3) %(nt, nsim)
                writematrix(GR_mature, file_name4) %(nt, nsim)
            end %col

            % SAVE total population for each model, scen, ens

            Nntot=sum(permute(N_allsim, [2 3 1]), 3); % nt, nsim
            Nntot_mature=sum(permute(N_allsim_mature, [2 3 1]), 3); % nt, nsim
            GRtot = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
            GRtot_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER
            
            for t=2:nt
                GRtot(t,:)=Nntot(t,:)./Nntot(t-1,:);
                GRtot_mature(t,:)=Nntot_mature(t,:)./Nntot_mature(t-1,:);
            end

            if e_mod==1
                file_name =sprintf("%s/Codes_EP/N_results_tot/2009/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_ens%d_scen%d.txt", ordi, mod, mod, ens, scenEXT);
            
                else
                file_name =sprintf("%s/Codes_EP/N_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/N_results_%s/Ntot_%s/%s/Ntot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/GR_results_%s/GRtot_%s/%s/GRtot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/GR_results_%s/GRtot_%s/%s/GRtot_%s_%s_ens%d.txt", ordi, eco_mod, eco_mod, mod, eco_mod, mod, ens);
            end

            writematrix(Nntot, file_name)
            writematrix(Nntot_mature, file_name2)
            writematrix(GRtot, file_name3)
            writematrix(GRtot_mature, file_name4)

        end %ens

    end %scenEXT

end %climate model


%% MEAN OF ENSEMBLES 2009

% Load Growth rate from one of the 6 models

mod_pop= ["Stef", "IPM", "Sat"];
e_mod=2;
eco_mod=mod_pop(e_mod);

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3","LE_MPI-ESM", "LE_CESM2", "LE_CESM1-CAM5-Paris2"];

if e_mod==1 %model CMR
    scen_list = [1 2 4];
    mod_list=[1:3 6 7];
else % other models don't have scenEXT
    scen_list = [1];
    mod_list=6;
end


for m=mod_list %climate model
    mod=model(m)

    if e_mod==1 %CMR
        yearstart = 2009; yearstop = info_models(m,4);
    elseif e_mod==2 %IPM
        yearstart = 2009; yearstop=2100;
    end

    timePOP   = yearstart:yearstop;
    nt        = length(timePOP);
    ncol=66;
    nens= info_models(m,5);
    nsim=100;

    for scenEXT=scen_list
        display(scenEXT)

        % SIMULATIONS

        jd = 0; % Choose 1) jd=0 Random search strategy; 2) jd=1 Oriented search strategy

        %for ens=1:nens

        R=zeros(ncol, nt, nsim);

            % Import data of all the colonies
            parfor col=1:ncol
                if e_mod==1
                    file_name=sprintf('%s/Codes_EP/GR_results_tot/2009/mean_ensemble/GR_results_%s/%s/GRmean_%s_%s_scen%d_col%d.txt', ordi, eco_mod, mod, eco_mod, mod, scenEXT, col);
                else
                   file_name=sprintf('%s/Codes_EP/GR_results_tot/2009/mean_ensemble/GR_results_%s/CESM2/GRmean_%s_CESM2_col%d.txt', ordi, eco_mod, eco_mod, col);
                end

                R(col, :, :) = readmatrix(file_name); %(nt, nsim)
            end

            N_allsim=zeros(ncol, nt, nsim);

            for i=1:nsim

                r = R(:,:,i); %(ncol, nt)
                r  = log(r); %(nt, ncol)

                A_col = [];

                % Matrix of connexion %
                %%%%%%%%%%%%%%%%%%%%%%%
                D = D_post(i); % Mean distance dispersal
                % Dispersal kernel        % Rem: k(0)=0, the probability of staying where we are, is zero
                k = @(x) (x<D)/D .*(x>0); % Uniform kernel with D = mean distance dispersal

                % Matrix of connexion without renormalisation
                dmat = k(dm);    % Unif kernel
                % Renormalization of the connexion matrix %
                den  = sum(dmat)';         % sum_{l\neq i}(d_{i->l})
                Den  = den+(den<=0);       % delete zeros in Den et les remplace par des 1
                Dmat = dmat*diag(1./Den);

                pm = [];
                parfor j = 1:n_subgroup
                    Is = i_subgroup(j+1)-i_subgroup(j);
                    pm = [pm;ones(Is,1)*Pm_post(i,j)];
                end

                % Movement scenario %
                %%%%%%%%%%%%%%%%%%%%%
                rM = 0.25;      % maximal death rate
                rm = rM*(1-pm);
                if (pm==1)
                    movement=@(r)  (r<0); %HIGH
                else
                    movement=@(r)  (r<-rm) - (r./rm).*(-rm<=r).*(r<0); %HIGH
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Initialization of size of colonies %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Nntot = zeros(nt, nsim);
                Nncol = zeros(ncol,nt); % n(i,t,j)=  # individuals in colony j at time t in simu i
                Nncol(:,1) = BE; %no uncertainties in BE
                Nnew = Nncol(:,1);

                % Time loop forward
                % Inidividuals
                t = yearstart;
                a_col = [];
                it = 1;
                while (t<yearstop)

                    Nold = Nnew;

                    %%% Reproduction at time t+1
                    %%% ###################################################################################################################################
                    fN    = f(Nold,r(:,it));
                    rstar = g(Nold,r(:,it))-1; % Effective growth rate at time t when blicking happens
                    %%% ###################################################################################################################################


                    %%%% Construction of informed dispersal matrix
                    mat_r = meshgrid(rstar)';
                    [rmax_vec, ind_col_sel] = max(mat_r.*(dm<D)+((dm<D)-1)); % ((dm<D)-1) to eliminate zeros and the max provieds the largest negative growth rate
                    Mat_select = zeros(ncol,ncol);
                    for ii=1:ncol
                        Mat_select(ind_col_sel(ii),ii) = 1;
                    end
                    rmax_mat = meshgrid(rmax_vec);

                    if (jd==1)
                        % Oriented dispersal
                        disp_sel = Mat_select; display(scenEXT)
                    else
                        % Random Dispersal
                        disp_sel = Dmat;
                    end

                    %%% Exploration
                    lambda = movement(rstar); % mean proportion of individual that leaves colonies

                    % Blinking of colonies (everybody dipserse and nobody arrive directly) ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    %%% ###################################################################################################################################

                    iblink = Blinking(Blinking(:,2)==(t+1),1); % Give the blinking colonies at time t
                    s_disp_blink = 0.5*(0.9438+0.9253);   % survival rate after living a blinking colony / ADD a 10% decrease from what is optimal**
                    pm_blink     = 0.9;                       % emigration rate of blinking colonies ****************************

                    lambda(iblink) = s_disp_blink*pm_blink;    % proportion of emigrants that survive
                    disp_sel(iblink,:) = 0;                    % no emigration to the blinking colony


                    den_disp = sum(disp_sel,1)';
                    Disp_sel = disp_sel-diag(den_disp);        % pour assurer que la somme des colonnes vaut 0

                    NI = (Disp_sel*(lambda.*fN));              % emigrants
                    NO = fN;                                   % newborn and alive individuals
                    N_E = NI + NO;

                    % blinking ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    N_E(iblink) = (1-pm_blink)*fN(iblink);     % survivers in blinking colonies x

                    %%% ###################################################################################################################################


                    %           Random settlment
                    full    = (N_E>K);                     % localization of overcrowded colonies
                    Overpgn = (N_E-K).*full;               % surplus de manchots
                    Overpgn_mov = Dmat*Overpgn;            % surplus re-distribue en fonction du kernel des distance
                    N_RS    = N_E - Overpgn + Overpgn_mov; % variable intermediare: locaux+ surplus dans chaque col
                    Nnew = N_RS;

                    %           Nnew = N_E; % informed dispersal
                    %%% ###################################################################################################################################

                    Nncol(:,it+1) = Nnew;
                    it = it+1;
                    a_col = [a_col,lambda];
                    t = t+1;

                end % year

                N_allsim(:, :, i) = Nncol; % all sim for one colony

            end % sim

            % SAVE one file txt per model, scen, colony, ens

            N_allsim_mature=zeros(ncol,nt,nsim);

            for col=1:ncol

                if e_mod==1
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/mean_ensemble/MAT_results_Stef/%s/MAT_Stef_%s_scen%d_col%d.txt', ordi, mod, mod, scenEXT, col);
                else
                    filename2=sprintf('%s/Codes_EP/MAT_results_tot/2009/mean_ensemble/MAT_results_IPM/CESM2/MAT_IPM_CESM2_col%d.txt', ordi, col);
                end
                    
                MAT=readmatrix(filename2);

                N=permute(N_allsim(col,:,:), [2 3 1]); % nt, nsim TOTAL POPULATION
                GR = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
                GR_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER

                for t=1:length(N(:,1))
                    for s=1:length(N(1,:))
                        N_allsim_mature(col,t,s) = N(t,s)*(MAT(t,s));
                    end
                    if t>1
                        GR(t,:)=N(t,:)./N(t-1,:);
                        GR_mature(t,:)=N_allsim_mature(col,t,:)./N_allsim_mature(col,t-1,:);
                    end
                end

                N_mature=permute(N_allsim(col, :,:), [2 3 1]);

                N=permute(N_allsim(col,:,:), [2 3 1]); % nt, nsim

                if e_mod==1 %CMR, we precise the scenEXT in the name of the file
                    file_name =sprintf("%s/Codes_EP/N_results_tot/2009/mean_ensemble/N_results_CMR/%s/Nmean_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_CMR/%s/N_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/mean_ensemble/GR_results_CMR/%s/GR_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/mean_ensemble/GR_results_CMR/%s/GR_CMR_%s_scen%d_col%d.txt", ordi, mod, mod, scenEXT, col);
                    
                else
                    file_name =sprintf("%s/Codes_EP/N_results_tot/2009/mean_ensemble/N_results_IPM/%s/Nmean_IPM_%s_col%d.txt", ordi, mod, mod, col);
                    file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_IPM/%s/N_IPM_%s_col%d.txt", ordi, mod, mod, col);
                    file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/mean_ensemble/GR_results_IPM/%s/GR_IPM_%s_col%d.txt", ordi, mod, mod, col);
                    file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/mean_ensemble/GR_results_IPM/%s/GR_IPM_%s_col%d.txt", ordi, mod, mod, col);

                end
 
                writematrix(N, file_name) %(nt, nsim)
                writematrix(N_mature, file_name2) %(nt, nsim)
                writematrix(GR, file_name3) %(nt, nsim)
                writematrix(GR_mature, file_name4) %(nt, nsim)

            end %col

            % SAVE total population for each model, scen, ens

            Nntot=sum(permute(N_allsim, [2 3 1]), 3); % nt, nsim
            Nntot_mature=sum(permute(N_allsim_mature, [2 3 1]), 3); % nt, nsim
            GRtot = ones(length(N(:,1)), length(N(1,:))); %TOTAL POPULATION
            GRtot_mature = ones(length(N(:,1)), length(N(1,:))); %BREEDER AND NON-BREEDER
            
            for t=2:nt
                GRtot(t,:)=Nntot(t,:)./Nntot(t-1,:);
                GRtot_mature(t,:)=Nntot_mature(t,:)./Nntot_mature(t-1,:);
            end

            if e_mod==1
                %file_name =sprintf("/Users/fledgelab/Desktop/Codes_EP/N_results_tot/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_ens%d_scen%d.txt", mod, mod, ens, scenEXT);
                file_name =sprintf("%s/Codes_EP/N_results_tot/2009/mean_ensemble/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_scen%d.txt", ordi, mod, mod, scenEXT);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_CMR/Ntot_CMR/%s/Ntot_CMR_%s_scen%d.txt", ordi, mod, mod, scenEXT);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/mean_ensemble/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_scen%d.txt", ordi, mod, mod, scenEXT);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/mean_ensemble/GR_results_CMR/GRtot_CMR/%s/GRtot_CMR_%s_scen%d.txt", ordi, mod, mod, scenEXT);

            else %IPM
                file_name =sprintf("%s/Codes_EP/N_results_tot/2009/mean_ensemble/N_results_IPM/Ntot_IPM/%s/Ntot_IPM_%s.txt", ordi, mod, mod);
                file_name2 =sprintf("%s/Codes_EP/N_mature_results_tot/2009/mean_ensemble/N_results_IPM/Ntot_IPM/%s/Ntot_IPM_%s.txt", ordi, mod, mod);
                file_name3 =sprintf("%s/Codes_EP/GRdisp_results_tot/2009/mean_ensemble/GR_results_IPM/GRtot_IPM/%s/GRtot_IPM_%s.txt", ordi, mod, mod);
                file_name4 =sprintf("%s/Codes_EP/GRdisp_mature_results_tot/2009/mean_ensemble/GR_results_IPM/GRtot_IPM/%s/GRtot_IPM_%s.txt", ordi, mod, mod);

            end

            writematrix(Nntot, file_name)
            writematrix(Nntot_mature, file_name2)
            writematrix(GRtot, file_name3)
            writematrix(GRtot_mature, file_name4)

        %end %ens

    end %scenEXT

end %climate model