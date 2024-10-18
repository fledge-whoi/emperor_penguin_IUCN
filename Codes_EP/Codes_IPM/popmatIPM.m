function [A]=popmatIPM(theta) 
% define population matrix 
% theta
%   1     2     3      4      5      6    
% [fec, phiA, phiJ, Pbreed, recr5, recr6]


A=zeros(7,7);
A(1,1) = theta(2)*theta(4);
A(1,2) = theta(2)*theta(4);
A(1,6) = theta(2)*theta(5);
A(1,7) = theta(2)*theta(6);

A(2,1) = theta(2)*(1-theta(4));
A(2,2) = theta(2)*(1-theta(4));

A(3,1) = 0.5 * theta(1) * theta(3);

A(4,3) = theta(2);
A(5,4) = theta(2);

A(6,5) = theta(2);

A(7,6) = theta(2) * (1- theta(5));
A(7,7) = theta(2) * (1- theta(6));