function [bBS bPB bR bS]=betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)

%breeding success
bBS=mvnrnd(betaBS,covBS,1);

% Proportion of 1st return - 1st year survival ***********************************
bPB=mvnrnd(betaPB,covPB,1);

% return probabilities
bR=mvnrnd(betaR,covR,1); % time varying estimates

% survival probabilities
% survival probabilities
bS=zeros(length(betaS(:,1)),6);
for j=1:length(betaS(:,1)) % weigth
    bS(j,:)=mvnrnd(betaS(j,:),covS(:,:,j),1);
    if bS(j,3)>0, bS(j,3)=0;end;
    if bS(j,6)>0, bS(j,6)=0;end;
end


return
