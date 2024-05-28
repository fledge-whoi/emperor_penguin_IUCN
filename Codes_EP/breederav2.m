function [Fav Mav]=breederav2(THETA,nt)
% calculate the number of males Mav and Females Fav avalaible for mating
% at time t+1
% ASR: proportion of males in the breeders class at time t
%          1 2  3   4   5   6  7   8  9
%   THETA=[s;BS;Ppb;Pnb;Pb;Spb;Sf;Sm;S0];
    Fav=zeros(5,nt);
    Mav=zeros(5,nt);
    Fav(1,:) =   THETA(3,:).*THETA(5,:);
    Fav(2,:) =   THETA(6,:).*THETA(8,:);
    Fav(5,:) = 0.5*THETA(7,:).*THETA(8,:);
    Mav(3,:)=    THETA(3,:).*THETA(5,:);
    Mav(4,:)=   THETA(6,:).*THETA(9,:);
    Mav(5,:)= 0.5*THETA(7,:).*THETA(9,:);
return
