function [A]=popmat(theta,uf,um) 
% define population matrix 
% theta
%   1   2   3    4   5    6    7  8   9
%[rho, BS, Rpb, S0, Spb, Rnb, Rb, Sf, Sm] 

%m=(0.5*theta(2))/((theta(8)^(8/12))*(theta(9)^(8/12)));
m=(theta(2))/((theta(8)^(8/12))*(theta(9)^(8/12))); %proba that a breeding pair raises offspring

f1=(1-theta(1))*m;
f3=theta(1)*m;

A=zeros(5,5);
A(1,1)=  (1-uf*theta(3))*theta(5);
A(1,2)=0; A(1,3)=0; A(1,4)=0;
A(1,5)=(1-uf*theta(3))*theta(4)*f1;

A(2,1)=0; A(2,3)=0; A(2,4)=0;
A(2,2)=(1-uf*theta(6))*theta(8);
A(2,5)=0.5*(1-uf*theta(7))*theta(8);

A(3,1)=0; A(3,2)=0; A(3,4)=0;
A(3,3)=(1-um*theta(3))*theta(5);
A(3,5)=(1-um*theta(3))*theta(4)*f3;

A(4,1)=0; A(4,2)=0; A(4,3)=0;
A(4,4)=(1-um*theta(6))*theta(9);
A(4,5)=0.5*(1-um*theta(7))*theta(9);

A(5,1)=uf*theta(3)*theta(5);
A(5,2)=uf*theta(6)*theta(8);
A(5,3)=um*theta(3)*theta(5);
A(5,4)=um*theta(6)*theta(9);
A(5,5)=0.5*uf*theta(7)*theta(8)+0.5*um*theta(7)*theta(9)+...
uf*theta(3)*theta(4)*f1+...
um*theta(3)*theta(4)*f3;
