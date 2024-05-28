function dm=dispersion(d,ncol)

dm=zeros(ncol,ncol);
% calcul coastal distance between colony i and j
for i=1:ncol
    for j=i+1:ncol
% if j<i
%    travel1=[i-1:-1:j];
%    travel2=[i:1:ncol 1:1:j-1];
%  elseif j>i
   travel1=[i:1:j-1];
   travel2=[i-1:-1:1 46:-1:j];
%end
    
d1=sum(d(travel1));
d2=sum(d(travel2));
dm(i,j)=min(d1,d2);
dm(j,i)=min(d1,d2);
    end
end

return

