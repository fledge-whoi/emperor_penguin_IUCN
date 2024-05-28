function [theta]=parameter_RAND3(xA,SICa,bR,bS,SICf,SICm,W,tlimit)

theta=zeros(9,tlimit);
theta(1,:)=0.5; %sex-ratio
theta(3,:)=0.0561/(1-xA);% Proportion of 1st return
theta(4,:)=0.7715/(1-xA);% 1st year at sea survival

% stochasticity due to interannual variability
theta(2,:)=invlogit(-0.0059-1.7535*SICa(:,4)+(randn(tlimit,1)*1.1646)); %BS 
theta(6,:)=invlogit(bR(ceil(17.*rand(tlimit,1)))); % Proportion of return NB 
theta(7,:)=invlogit(bR(17+ceil(13.*rand(tlimit,1)))); % Proportion of return B 

% % no ENV stochasticity 
% theta(2,:)=invlogit(-0.0059-1.7535*SICa(:,4)); %BS 
% theta(6,:)=0.0446; % Proportion of return NB 
% theta(7,:)=0.8619; % Proportion of return B 


% survie uncertainties in model selection - survival is a weigthed average
% wAIC
sf=zeros(1,tlimit); %female survival
for j=1:length(bS(:,1))
     sf(1,:)=sf(1,:)+invlogit(bS(j,1)+(bS(j,2)*SICa(:,SICf(j)))+(bS(j,3)*SICa(:,SICf(j)).^2))'.*repmat(W(j,1),1,tlimit);     
end  
theta(8,:)=sf;

sm=zeros(1,tlimit); %male survival

for j=1:length(bS(:,1))
    sm(1,:)=sm(1,:)+invlogit(bS(j,4)+(bS(j,5)*SICa(:,SICm(j)))+(bS(j,6)*SICa(:,SICm(j)).^2))'.*repmat(W(j,1),1,tlimit);
end
theta(9,:)=sm;

theta(5,:)=(theta(9,:)+theta(8,:))/2;% pre-breeders survival, inconnue, on fait moy de FS et MS

return