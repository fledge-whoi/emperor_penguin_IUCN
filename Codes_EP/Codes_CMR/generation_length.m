%Calculate Generation Length
% Visser et al 2021

% DIRECTORY to Codes_EP
ordi = 

%% l, V, W, A stables

%load Astable.mat %"Atot", "Vtot", "ltot","Wtot", "uftot", "umtot", "thetatot"
%load Astable1.mat 
% 1 1900-1950  2 1950-2000  3 2000-2050 present  4 2050-2100

% Most recent version, with standardized sea ice
% one value per colony and simulation
%load("Astable_SICstd.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot");

% Version with observed data
load("Astable_SICobs.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot");

ncol = 66;
nsim = 100;

% Select stable population lambda
lambda=ltot(300,:,:); %(nsim, ncol)
lambda=median(median(lambda)); %lambda equilibre mean on colonies and simulations

nsimulation = ncol*nsim;

% Calculate median of simulations and colonies

A=zeros(5);
for indice=1:25
    med=zeros(nsimulation,1);
    x=1;
    for simu=1:nsim
        for col=1:ncol
            med(x)=Atot{simu, col}(indice);
            x=x+1;
        end
    end
A(indice)=median(med);
end

V=zeros(1,5);
W=zeros(5,1);
for indice=1:5
    med=zeros(nsimulation,1);
    med2=zeros(nsimulation,1);
    x=1;
    for simu=1:nsim
        for col=1:ncol
            med(x)=Vtot{simu,col}(indice);
            med2(x)=Wtot{simu,col}(indice);
            x=x+1;
        end
    end
V(indice)=median(med);
W(indice)=median(med2);
end

[w d v]=eig(A);
d=diag(d);
ivp=find(d==max(d)); %indice de la vp dominante
V=v(:,ivp);
W=w(:,ivp);

% Elasticities
% Parameters to calculate fertilities
% Values at equilibrium from proj_pop_glb_nomov2
uf=median(median(uftot)); %median on 66 colonies and sim
um=median(median(umtot));
thetatot=median(thetatot(:,:,:), 3); %median on colonies
thetatot=permute(thetatot, [2 1]);
theta=median(thetatot(:,:)); %median on simulations

m=(theta(2))/((theta(8)^(8/12))*(theta(9)^(8/12))); %proba that a breeding pair raises offspring
f1=(1-theta(1))*m;
f3=theta(1)*m;

%Fertilities
f15=A(1,5);
f35=A(3,5);
f55=uf*theta(3)*theta(4)*f1 + um*theta(3)*theta(4)*f3;

ef15=(f15*V(1)*W(5))/(lambda*V'*W);
ef35=(f35*V(3)*W(5))/(lambda*V'*W);
ef55=(f55*V(5)*W(5))/(lambda*V'*W);

T=1/(ef15+ef35+ef55)

% Formula 13 Bienvenu and Legendre

% Decomposing A=U+F
F= zeros(5); %Fertilities
F(1,5)=f15;
F(3,5)=f35;
F(5,5)=f55;

U=A-F;

% Generation length
GT1= 1 + ((V'*U*W)/(V'*F*W))



%% Calculate one GL per colony

% Most recent version, with standardized sea ice, scenEXT=[1 2 4] and
% the number of total breeding failures adapted to observations
% one value per colony and simulation

clear A; clear V; clear W
load("Astable_SICobs.mat", "Atot", "Vtot", "ltot","Wtot", "uftot", "umtot","thetatot");

ncol = 66;
nsim = 100;

lambda=ltot(300,:,:); %(nsim, ncol)
lambda=median(lambda(:,:)); %lambda equilibre mean simulations


% Calculate median of simulations for each colony

A(ncol,1)={zeros(5)};
A2=zeros(5);
for col=1:ncol
    for indice=1:25
        med=zeros(nsim,1);
        for simu=1:nsim
            med(simu)=Atot{simu, col}(indice);
        end
        A2(indice)=median(med);

    end %indice matrix
    A{col,1} = A2;
end %col


V(ncol,1)={zeros(1,5)};
W(ncol,1)={zeros(5,1)};

for col=1:ncol
    V2=zeros(1,5);
    W2=zeros(5,1);

    for indice=1:5
        med=zeros(nsim,1);
        med2=zeros(nsim,1);

        for simu=1:nsim
            med(simu)=Vtot{simu,col}(indice);
            med2(simu)=Wtot{simu,col}(indice);

        end
        V2(indice)=median(med);
        W2(indice)=median(med2);
    end
    V{col,1}=V2;
    W{col,1}=W2;
end %col


GL_tot=zeros(ncol,1);

for col=1:ncol
    
    A2 = A{col};
    uftot2=uftot(:,col);
    umtot2=umtot(:,col);
    thetatot2 = thetatot(:,:,col);

    [w d v]=eig(A2);
    d=diag(d);
    ivp=find(d==max(d)); %indice de la vp dominante
    V=v(:,ivp);
    W=w(:,ivp);

    % Elasticities
    %Parameters to calculate fertilities
    % Values at equilibrium from proj_pop_glb_nomov2
    uf=median(uftot2); %median of simulations
    um=median(umtot2);
    thetatot2=permute(thetatot2, [2 1]); %(simulations, 9 param)
    theta=median(thetatot2(:,:)); %median of simulations

    m=(theta(2))/((theta(8)^(8/12))*(theta(9)^(8/12))); %proba that a breeding pair raises offspring
    f1=(1-theta(1))*m;
    f3=theta(1)*m;

    %Fertilities
    f15=A2(1,5);
    f35=A2(3,5);
    f55=uf*theta(3)*theta(4)*f1 + um*theta(3)*theta(4)*f3;

    ef15=(f15*V(1)*W(5))/(lambda*V'*W);
    ef35=(f35*V(3)*W(5))/(lambda*V'*W);
    ef55=(f55*V(5)*W(5))/(lambda*V'*W);

    T=1/(ef15+ef35+ef55);

    % Formule 13 Bienvenu et Legendre

    % Decomposing A=U+F
    F= zeros(5); %Fertilities
    F(1,5)=f15;
    F(3,5)=f35;
    F(5,5)=f55;

    U=A2-F;

    GT1= 1 + ((V'*U*W)/(V'*F*W));
    GL_tot(col,1) = GT1;

end %col


%% Figure GL for all colonies

figure
scatter([1:66], GL_tot, 50, 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);
xlabel('Colonies', 'Fontsize', 20)
ylabel('Generation length (Years)', 'Fontsize', 20)
set(gca, 'FontSize', 18);
title('Generation length','Fontsize', 25)


%% Figure GL as function of SIC
% Calculate mean SIC for the 50 ens for each season

%load seaiceLE.mat %LE2 CESM CMIP6 SSP370
load(sprintf('%s/Codes_EP/Codes_CMR/mat_data/LE_CESM2_seaice_std.mat', ordi));

ncol=66;
SIC=zeros(ncol,4);
nens=50;

% Choose period 
tstart=1; % 0 =1900
tend=51; %50 =1950
nt=tend-tstart+1;

for c=1:ncol
    for si=1:4 %season
        SICmoy=zeros(nens*nt,1);
        x=1;
        for ens=1:nens
            for t=1:50
                SICmoy(x,1)=ENS(ens).SICa(t, si, c);
                x=x+1;
            end
        end
        SIC(c,si)=median(SICmoy(:,1));
    end
end

SIC = mean(SIC(:,:), 2); %mean on season

scatter(SIC, GL_tot, 50, 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);
xlabel('SICa mean')
ylabel('Generation length')

%% Observed data

ncol = 66;
nt = 40;
SIC_tot = zeros(ncol, 4, nt);
for col=1:66
    SIC_col = readmatrix(sprintf("%s/Codes_EP/Codes_CMR/SIC_obs_allcol/SIC_col%d.csv", ordi, col));
    for season=1:4
        SIC_tot(col, season, :)= SIC_col(:,season);
        SIC_mean = mean(SIC_tot(col,season,1:30));
        SIC_tot(col,season,:) = (SIC_tot(col,season,:) - SIC_mean) ./ SIC_mean;
    end
end

% Calculate anomalies
% for s=1:4
%     SIC_mean = mean(mean(SIC_tot(:,s,:), 3));
%     SIC_tot(:,s,:) = (SIC_tot(:,s,:) - SIC_mean) ./ SIC_mean;
% end

SIC=mean(mean(SIC_tot(:,:,:), 3),2); %mean on the 40 years

scatter(SIC, GL_tot, 50, 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);
xlabel('SICa mean')
ylabel('Generation length')