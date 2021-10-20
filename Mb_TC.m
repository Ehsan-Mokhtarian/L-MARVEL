function [Mb,nTests,sep_sets] = Mb_TC(D, alpha)
% Total Conditioning, see Pellet and Elisseeff

    nVars = size(D,2);  
    nTests = 0;
    Mb = zeros(nVars);
    sep_sets = false(nVars,nVars,nVars);
    for X=1:(nVars-1)
        for Y=(X+1):nVars
            S = [1:(X-1),(X+1):(Y-1),(Y+1):nVars];
            CI = CI_Test(X,Y,S,D,alpha);
            nTests = nTests+1;
            if ~CI
                Mb(X,Y)=1;
                Mb(Y,X)=1;
            else
                sep_sets(X,Y,S) = true;
                sep_sets(Y,X,S) = true;
            end
        end
    end
end
