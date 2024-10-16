function r=rmse(F,A)

    n=length(F);
    
    S=0;
    for i=1:n
        S = S + ((F(i)-A)^2);
    end

    r=sqrt(S./n);
