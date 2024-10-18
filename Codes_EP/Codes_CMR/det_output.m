function [lambda1, SAD, RVw, Sensitivity, Elasticity, t20, dr, GT]= det_output(Umat, Fmat, theta) % S, GT
%function [lambda1, SAD, RVw,t20]= det_output(mat,lv)

%**********************************************************************
%INPUT
%mat: population matrix
%OUTPUT
%Lambda, stable age distribution and reproductive value, 
%Sensitivity and Elasticity matrices
%sensitivity lower parameter:S: 1) sensitivity, 2) elasticity, 3) sinus
%sensitivities (Doherty et al)

%*************************************************************************
%Lambda, stable age distribution and reproductive value routine;
mat=Umat+Fmat;
% Calculate the eigenvalues (dmat) and eigenvectors (wmat, vmat) of the matrix;
[wmat, dmat,vmat ]=eig(mat);

% Calculate the rate of population change (lambda1 = the dominant eigenvalue), ...
% the subdominant eigenvalue (lambda2), the damping ratio (rho) and...
% the time to convergence (t20);
lambda=diag(dmat);
lambvec=sort(abs(lambda));
dim=size(lambvec,1);
lambda1=lambvec(dim);
lambda2=lambvec(dim-1);
rho=lambda1/lambda2;
dr=lambda1/abs(lambda2);
t20=log(20)/log(rho);


% Calculate the stable age distribution (SAD) and the weighted reproductive value (RVw);
imax=find(lambda==max(lambda));
v = vmat(:,imax); % vecteur propre a gauche
w = wmat(:,imax); % vecteur propre a droite


RawAge=wmat(:,imax);
SAD=RawAge./sum(RawAge);
RawRV=real(vmat(imax,:))';
RVw=RawRV/RawRV(1,1);

%**********************************************************************
%  Sensitivity and Elasticity subroutine.

%Calculate a new matrix;
num=RVw*SAD';
%Calculate product of stable age and reproductive value;
den=sum(RVw.*SAD);
%Calculate sensitivity;
Sensitivity=num./den;

%Divide original matrix by lambda.
temp=mat./lambda1;
%Calculate elasticity.
Elasticity=Sensitivity.*temp;

%calculate sensitivities for lower level
%S: 1) sensitivity, 2) elasticity, 3) sinus sensitivities (Doherty et al)

%**********************************************************************************************
% [PDlv]= PD_lowerlevel(lv);
% S(:,:)=zeros(length(lv),3);
% for j=1:length(lv)
%  S(j,1)=sum(sum(PDlv(:,:,j).*Sensitivity));
%  S(j,2)=lv(j)/lambda1*S(j,1);
%  %S(j,3)=(((lv(j)*(1-lv(j)))^0.5)/lambda1)*S(j,1);
% end
% 
% GT=1/S(1,2);
% %GT2=log(rho)/log(lambda1);

 GT= 1 + ((v'*Umat*w)/(v'*Fmat*w));
 
return