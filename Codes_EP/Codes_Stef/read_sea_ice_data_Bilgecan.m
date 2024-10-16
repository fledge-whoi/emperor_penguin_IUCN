clear; close all; clc;
timeSICcesm=1900:2100;ntcesm=length(timeSICcesm);
nc=66; nENS=50;

Indexa=find(timeSICcesm==2007);

for ens=1:nENS % ensemble
    ENS(ens).SIC=nan(ntcesm,4,nc);
    ENS(ens).SICa=nan(ntcesm,4,nc);

    for c=1:nc % colony

        filename=['/Users/sjenouvrier/Dropbox/_SJ_work/_projet/Emperor_penguin/NewforecastEP2023_LE2/SIC_CESMLE2/' ...
            num2str(ens),'/col',num2str(c),'.csv'];
        ENS(ens).SIC(:,:,c)=readmatrix(filename,'Range','B2:E202');

        for s=1:4 %seasons
            ENS(ens).SICa(:,s,c)=...
            (ENS(ens).SIC(:,s,c)-(mean(ENS(ens).SIC(1:Indexa,s,c))))./mean(ENS(ens).SIC(1:Indexa,s,c));
        end  %seasons

    end % colony
end% ensemble

%%
% POINTE GEOLOGIE
load SICaOBS.mat
SICa_stef=SICa;
clear SICa

figure
for s=1:4

    subplot(2,2,s)
    hold on
    for ens=1:nENS
        plot(1979:2018,ENS(ens).SICa(find(timeSICcesm==1979):find(timeSICcesm==2018),s,39),'k','linewidth', 1)
    end

    %plot(1979:2018,SICa(:,s,39),'r','linewidth', 3)

    plot(1979:2010,SICa_stef(:,s,30),'b','linewidth', 2)

end

% Kloa point
figure
for s=1:4

    subplot(2,2,s)
    hold on
    for ens=1:nENS
        plot(1979:2018,ENS(ens).SICa(find(timeSICcesm==1979):find(timeSICcesm==2018),s,21),'k','linewidth', 1)
    end

    %plot(1979:2018,SICa(:,s,39),'r','linewidth', 3)

    plot(1979:2010,SICa_stef(:,s,19),'r','linewidth', 2)

end

% crozet
figure
for s=1:4

    subplot(2,2,s)
    hold on
    for ens=1:nENS
        plot(1979:2018,ENS(ens).SICa(find(timeSICcesm==1979):find(timeSICcesm==2018),s,50),'k','linewidth', 1)
    end

    %plot(1979:2018,SICa(:,s,39),'r','linewidth', 3)

    plot(1979:2010,SICa_stef(:,s,38),'g','linewidth', 2)

end

% halley
figure
for s=1:4

    subplot(2,2,s)
    hold on
    for ens=1:nENS
        plot(1979:2018,ENS(ens).SICa(find(timeSICcesm==1979):find(timeSICcesm==2018),s,7),'k','linewidth', 1)
    end

    %plot(1979:2018,SICa(:,s,39),'r','linewidth', 3)

    plot(1979:2010,SICa_stef(:,s,7),'c','linewidth', 2)

end