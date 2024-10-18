%% Environmental projections 

% Import SIC_trans data from Bilgecan and create data mat file
% These data are modified to have the same mean and std than the
% observation 2009-2018
% File created SIC_proj_trans_total = (nt, col, simu_env)

ncol=66; %colonies
nt=192; %time
nens=50; %ensemble

SIC_proj_total = zeros(nens, ncol, nt); % SIC Projections for 1909-2100 all colonies, 50 ensembles

for ens=1:nens

    file_name=sprintf('%s/Codes_EP/Codes_SAT/SIC_proj_trans/SIC_proj_trans_%d.csv', ordi, ens);
    env_proj=readmatrix(file_name);
    env_proj(1,:)=[]; env_proj(:,1)=[];

    for col=1:ncol
        for t=1:nt

        
            SIC_proj_total(ens, col, t) = env_proj(col, t);

        end %t
    end %col 
end %50 simu env

SIC_proj_total=permute(SIC_proj_total, [3 2 1]); % time, col, sim_env

% save('SIC_proj_trans_total', 'SIC_proj_total') 


%% Visualize environmental data

%load SIC_proj_trans_total.mat

ncol=66;
nt=201;
nens=50;

figure
for ens=1:50
    plot(1909:2100, mean(SIC_proj_total(:,:,ens), 2)); %mean on colonies
    hold on
end





