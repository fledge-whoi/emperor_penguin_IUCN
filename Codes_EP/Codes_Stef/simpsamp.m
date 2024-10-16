function c=simpsamp(k)
% H Caswell
% function c=simpsamp(k)
% returns a set of k  weights uniformly sampled from the k-simplex
% sampling by broken stick method

breaks=[0 rand(1,k-1) 1];

breaks=sort(breaks);

intervals=breaks(2:k+1)-breaks(1:k);

c=intervals;

end
