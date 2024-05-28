% Program to calculate the carrying capacity
% We use projections of the population 
% NRBE contains intial population sizes BE obtained from projection of the
% population without movement without density dependance with constant SICa
% (fixed at the predicted mean betweem 1900-1950)
% Rcol and Ncol were obtained 

close all; clear all

%% Initialisation


load NRBE_new.mat %BE, Rcol, Ncol

Rcol=Rcol_mean; %(sim, ens, nt, ncol)
Ncol=Ncol_mean; %(class, nt, ncol, ens, sim)

ncol=66;
nt=51;
nens=50;
time=(1:nt);
nsim=length(Rcol(:,1,1,1));

% variables: BE, r, MED

%Median on simulations

%MED=median(Ncol(5,:,:,:,:), 5); %Breeding pairs   (class 5, t, c, ens)

%r=median(Rcol(:,:,:,:)); %(1, nens, nt, ncol)

MED=Ncol(5,:,:,:,:); %Breeding pairs   (class 5, t, c, ens, sim)
r=Rcol(:,:,:,:); %(nsim, nens, nt, ncol)

BE(46)=0;              % Colony 46 is extinct
%r(:,:,:,end)=-0.25; %r       % growth rate unkown

lat  = xlsread('empe_sitesNewNB.xlsx','D2:D67');
long = xlsread('empe_sitesNewNB.xlsx','E2:E67');
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
dm=dispersion(d,ncol); % Distance between colonies


IK=[]; % Initialisation iK

% Without density dependance
f=@(N,r) N.*(1+r);
% g: Individual growth rate f(N)/N
g=@(N,r) (1+r);


% % With Density dependance K
% 
% f=@(N,r,K) N.*exp(r.*(1-N./K)).*(r>0) + (1+r).*N.*(r<=0);
% % g: Individual growth rate f(N)/N
% g=@(N,r) exp(r.*(1-N./K)).*(r>0) + (1+r).*(r<=0);


% Fonction K en fonction de N et r 
fK=@(N1,N2,r) log(1+r).*N1./(log(1+r)-log(N2)+log(N1)).*(r>0);


%% Calcul des K

%Ktot=zeros(nens, ncol);

Ktot=zeros(nsim, nens, ncol);
r=permute(r, [1 3 4 2]); %(nsim, nt, ncol, nens)

for sim=1:nsim
for ens=1:nens
    display(ens)
    K=[];k=[];Ik=[];
    for t=2:nt-1
        kinter=zeros(ncol,1);
        for c=1:ncol
            %kinter(c)=fK(MED(1, t,c,ens), MED(1, t+1,c,ens), r(1, t, c,ens));
            kinter(c)=fK(MED(1, t,c,ens,sim), MED(1, t+1,c,ens,sim), r(sim, t, c,ens));

            %     ki=log(1+r(:,t)).*MED(:,t)./log((1+r(:,t))./(1+R(:,t)));
        end
        ki=kinter(:,1);
        k=[k,ki];
    end
    ikp=(k>0); ikm=(k<0);

    for i=1:ncol
        if((sum(ikp(i,:)+ikm(i,:))>0))
            kip=k(i,ikp(i,:)); %les k>0 dans col i
            kimax=max(kip);
            nkim=sum(ikm(i,:));
            kim=kimax*ones(1,nkim);
            ki=max([median([kip,kim]),BE(1,i)]);
            %         ki=max([kip,BE(i)]);
            %k(i,ik(i,:))/BE(i);
        else
            ki=BE(1,i);
        end
        K=[K;ki];
    end %col
    K(end)=200;
    KK=K*ones(1,91);

    Ktot(sim, ens,:)=K;

end %ens

end %sim
%save('K_sim.mat', 'K')

r=permute(r, [1 4 2 3]);


%% Knorm for each ensemble

%K_norm=zeros(nens, ncol); %K=1*BE, 2*BE, etc
K_norm=zeros(nsim, nens, ncol); %K=1*BE, 2*BE, etc
for c=1:ncol
    K_norm(:,:,c) = Ktot(:,:,c)./BE(1,c);
end

%K_norm=permute(K_norm, [ 2 1 ]);
%K_norm_median=median(K_norm);

% figure
% for ens=1:50
%     clf
%     histogram(K_norm(:,ens),7)
%     xline(K_norm_median(1,ens), 'Color', 'r', 'Linewidth', 2)
%     axis([0 5 0 inf])
%     pause
% end

%m=median(K_norm);
%Kmoytot=median(m);
%K_norm=permute(K_norm, [ 2 1 ]);

Kmax = max(max(K_norm(K_norm<10))) 
%save('Kmax_50ens', "Kmax")


%% Normalized K - Figure for each colony
       
K_norm=zeros(nens, ncol); %K=1*BE, 2*BE, etc
for c=1:ncol
    K_norm(:,c) = Ktot(:,c)/BE(1,c);
end
  
K_norm_median=median(K_norm);

figure
for c=[1:45 47:ncol]
    clf
    histogram(K_norm(:,c),7)
    xline(K_norm_median(1,c), 'Color', 'r', 'Linewidth', 2)
    axis([0 5 0 50])
    pause
end

m=median(K_norm);
Kmoytot=median(m([1:45 47:end])); %K median des colonies
   

%% Figure K, BE
Kmean=mean(Ktot(:,:))
figure

plot(Kmean(1,:),'k*-', 'linewidth',2) 
hold on
plot(BE,'r*-', 'linewidth',2)
hold off
xlabel('Colony')
legend('K', 'BE')
pause


%% For 5 scenario
load NRBE_temp2.mat
load NRBE_5scen.mat
ncol=66;
nt=51;
time=(1:nt);
nsim=length(Rcol(:,1,1));

lat  = xlsread('empe_sitesNewNB.xlsx','D2:D67');
long = xlsread('empe_sitesNewNB.xlsx','E2:E67');
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
dm=dispersion(d,ncol); % Distance between colonies


% Without density dependance
f=@(N,r) N.*(1+r);
% g: Individual growth rate f(N)/N
g=@(N,r) (1+r);

% With Density dependance K
% f=@(N,r,K) N.*exp(r.*(1-N./K)).*(r>0) + (1+r).*N.*(r<=0);
% % g: Individual growth rate f(N)/N
% g=@(N,r) exp(r.*(1-N./K)).*(r>0) + (1+r).*(r<=0);

% Fonction K en fonction de N et r
fK=@(N1,N2,r) log(1+r).*N1./(log(1+r)-log(N2)+log(N1)).*(r>0);

Kmoytot=zeros(5,1)

for scenEXT=1:5
    Ncol=SCEN2(scenEXT).Ncol;
    Rcol=SCEN2(scenEXT).GRcol;

    % variables: BE, r, MED
    MED=Ncol(5,:,:,:); % (class, t, c, x) %only Breeding pairs

    r=Rcol(:,:,:); %nsim, nt, nc

    BE(46)=0;              % Colony 46 is extinct
    r(:,:,end)=-0.25; %r       % growth rate unkown

    IK=[]; % Initialisation iK


    % Calcul des K
    Ktot=zeros(nsim, ncol);
    r=permute(r, [2 3 1]);
    for s=1:nsim
        display(s)
        K=[];k=[];Ik=[];
        for t=2:nt-1
            kinter=zeros(ncol,1);
            for c=1:ncol
                kinter(c)=fK(MED(1,t,c,s),MED(1,t+1,c,s),r(t,c,s));
                %     ki=log(1+r(:,t)).*MED(:,t)./log((1+r(:,t))./(1+R(:,t)));
            end
            ki=kinter(:,1);
            k=[k,ki];
        end
        ikp=(k>0);ikm=(k<0);

        for i=1:ncol
            if((sum(ikp(i,:)+ikm(i,:))>0))
                kip=k(i,ikp(i,:)); %les k>0 dans col i
                kimax=max(kip);
                nkim=sum(ikm(i,:));
                kim=kimax*ones(1,nkim);
                ki=max([median([kip,kim]),BE(1,i)]);
                %         ki=max([kip,BE(i)]);
                %k(i,ik(i,:))/BE(i);
            else
                ki=BE(1,i);
            end
            K=[K;ki];
        end
        K(end)=200;
        KK=K*ones(1,91);

        Ktot(s,:)=K;

    end
    %save('K_sim.mat', 'K')
    r=permute(r, [3 1 2]);


    % Normalized K

    K_norm=zeros(nsim, ncol); %K=1*BE, 2*BE, etc
    for c=1:ncol
        K_norm(:,c) = Ktot(:,c)/BE(1,c);
    end

    K_norm_median=median(K_norm);

    figure
    for c=[1:45 47:ncol]
        clf
        histogram(K_norm(:,c),20)
        xline(K_norm_median(1,c), 'Color', 'r', 'Linewidth', 2)
        axis([0 6 0 150])
        pause
    end

    m=mean(K_norm);
    Kmoytot(scenEXT)=mean(m([1:45 47:end])); %K moyen des colonies

end