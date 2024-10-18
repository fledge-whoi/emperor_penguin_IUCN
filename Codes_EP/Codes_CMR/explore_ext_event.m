% Explore extreme event scenarios

% approximately 42% of emperor penguin colony locations (i.e., 28 of the 66)
% experienced early fast-ice break-up in at least one year out of five (ATCMXLV, WP52), 
% possibly resulting in breeding failures. Moreover, four of five breeding 
% colonies in the central and eastern Bellingshausen Sea experienced total 
% breeding failure during 2022, after sea-ice break-up prior to fledging 
% (Fretwell et al., 2023). This latter event is the first recorded incident 
% of a regionally widespread breeding failure of emperor penguins that is 
% clearly linked with large-scale contractions in sea-ice extent.

clear; close all; clc;

% DIRECTORY to Codes_EP
ordi = "/Users/aliceeparvier/Desktop";

addpath(sprintf('%s/Codes_EP/Codes_CMR/', ordi))

%load informations on climate models (years, number of ens)
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/infos_models.mat', ordi));

model= ["LE_CanESM2", "LE_CESM1-CAM5", "LE_CSIRO-Mk3-6-0", "LE_GFDL-CM3", "LE_MPI-ESM", "LE_CESM2", "Paris2"] % "CanESM2","CESM1-CAM5", "CSIRO-Mk3-6-0", "GFDL-CM3","MPI-ESM"]
ncol=66


for i=[1:3 6:7]% Climate model

    mod=model(i) %name
    display(mod)
    file_name=sprintf('%s/Codes_EP/Codes_CMR/mat_data/%s_seaice_std.mat', ordi, mod)

    load(file_name)

    load NewNameColony.mat

    % colony
    yearstart=info_models(i,3);
    yearend=info_models(i,4);
    timePOP=yearstart:yearend;
    nt=length(timePOP);
    nc=66;

    % climate
    ntime=yearstart:yearend;
    nens=info_models(i,5); % climate ensemble, nens minimum is 5

    %Matrix containing when extreme events happen
    nsim=300;
    extreme_event=zeros(2, ncol, nt, nens, nsim); % scenEXT, colony, time, simu with scenEXT=1 or 2 (historical level vs increase)
    yearstart=2018;
    yearend=2022;
    ystart=find(ntime==yearstart);yend=find(ntime==yearend);
    nb_col_ext = zeros(1,nens*nsim);


    % perturbation - extreme events
    percentSI=13/355; % =3.66% of extreme events counted on historical period

    y2009=find(ntime==2009);y2018=find(ntime==2018);

    % MAIN LOOP

        for c=1:66 % colony
            display(c)

            for ens=1:nens % climate ensemble
                
                % SEA ICE FORECAST ENS (climatic ens 50) matrix SICa_500
                % time*season*colony

                % Extreme sea ice loss threshold
                threshold_SIC2=prctile(ENS(ens).SICa(y2009:y2018,2,c),percentSI);
                threshold_SIC4=prctile(ENS(ens).SICa(y2009:y2018,4,c),percentSI);
                EXT=find(ENS(ens).SICa(:,2,c)<threshold_SIC2 | ENS(ens).SICa(:,4,c)<threshold_SIC4);
                

                for s=1:nsim
                    % Depending on sea ice levels
                    iEXT=EXT(binornd(1,0.47,length(EXT),1)==1); % likely as not to observe an extreme events
                    extreme_event(2, c, iEXT,ens, s) = 1;
                    
                    % Historical level
                    PERTclimat=binornd(1, percentSI, 1, nt);
                    extreme_event(1, c, PERTclimat==1,ens, s) = 1;
                end % sim
                
            end  % ensemble

        end % colony

end %model


% Explore compatibility
yearstart=2018;
yearend=2022;
ystart=find(ntime==yearstart);yend=find(ntime==yearend);

nb_col_ext_hist = zeros(1,nens*nsim);
nb_col_ext_increase = zeros(1,nens*nsim);

x=1;
for ens=1:nens
    for s=1:nsim
        for c=1:ncol
            %historical level
            if sum(extreme_event(1, c, ystart:yend, ens, s)) > 0 %if col experiences at least one complete breeding failure
                nb_col_ext_hist(x) = nb_col_ext_hist(x) + 1;
            end
            %binomial rand in iEXT
            if sum(extreme_event(2, c, ystart:yend, ens, s)) > 0 %if col experiences at least one complete breeding failure
                nb_col_ext_increase(x) = nb_col_ext_increase(x) + 1;
            end
        end
        x=x+1;
    end %sim
end %ens

% figure
histogram(nb_col_ext_increase)
hold on
xline(median(nb_col_ext_increase), 'color', 'r')
hold on 
xline(28)
% figure
% histogram(nb_col_ext_hist)

median(nb_col_ext_increase)
median(nb_col_ext_hist)


%% Results
% Models 1-5 and 6=CEMS2
% The proportion of extreme event that lead to breeding failure to have the
% same mean number of colonies experiencing total breeding failure betweem
% 2018-2022 than the observation.

%proportion_ext_event = [0.66 0.82 0.66 0.76 0.64 0.60]
proportion_ext_event = [0.67 0.83 0.76 0.76 0.64 0.61 0.47]
%save('proportion_ext_event', 'proportion_ext_event')