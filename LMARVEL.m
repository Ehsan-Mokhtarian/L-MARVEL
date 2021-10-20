function [Adj, sep_sets, nTests, condSize] = ...
    LMARVEL(D, Mb, sep_sets, alpha, alphaMb)

% input:
% D: Data matrix with size (number of samples)*n
% Mb: n*n Markov boundary matrix. Entries can be 0 or 1
% sep_sets: separatin sets (n*n*n) found by Mb with true false entries
% alpha: significance level for CI tests. for example 0.01
% alpha_Mb: significance level for updating Mb. for example 2/n^2
%
% output:
% Adj: The matrix of learned skeleton. Entries can be 0 or 1
% sep_sets: all the separatin sets (n*n*n) 
% nTests: the number of performed CI tests
% condSize: Size of Conditioning sets. It is a (n+1)*1 vector such that
% i-th entry indicates the number of conditioning sets with size i-1

%**************************************************************************
%********************** Start of the Algorithm ****************************
%**************************************************************************


nVar = size(Mb, 1);             % number of variables
V = ones(1,nVar);               % Remaining Variables
Adj = -ones(nVar);              % Learned skeleton
nTests = 0;                     % Number of CI tests performed
condSize = zeros(nVar+1,1);     % Average conditioning set size
remFlag = true(nVar);           % Already checked for removability?

% Remove Edges between vertices that do not share Mb
Adj(~Mb)=0;

% Recursive part
for iter=nVar:-1:1
    Mb_size = sum(Mb);
    [~,ind] = sort(Mb_size);
    ind = ind(V(ind)>0 & remFlag(ind));
    Xrem = 0;             % removable variable to be identified
    % Identify Xrem:
    for X = ind
        remFlag(X) = false;
        Mb_X = find(Mb(X,:));
        [Adj, sep_sets, nTests, condSize] = findAdjacents(D,X,Mb_X,Adj,...
            nTests,condSize,sep_sets,alpha);
        [isR, Adj, nTests, condSize, sep_sets] = ...
            isRemovable(D,X,Mb_X,Adj,nTests,condSize,sep_sets,alpha);
        if isR
            Xrem = X;
            break
        end
    end
    
    % No removable variables found. Faulty CI tests!
    if Xrem==0
        break
    end
    % Remove Xrem and Update Mb info:
    [Mb, Adj, V, nTests, condSize, sep_sets, remFlag] = UpdateProcess...
        (D,Xrem,Mb,Adj,V,nTests,condSize,sep_sets,alphaMb,remFlag);
end
end

%**************************************************************************
%************************* Main Functions *********************************
%**************************************************************************

function [Adj, sep_sets, nTests, condSize] = findAdjacents(D,X,Mb_X,Adj,...
    nTests,condSize,sep_sets,alpha)
%findAdjacencies updates the list of variables adjacent to X given the
%Markov blanket of X, i.e., Mb_X.

for Y = Mb_X
    if Adj(X, Y) == -1
        [is_adj, sep_sets, nTests, condSize] = ...
            isAdjacent(D,X,Y,Mb_X,nTests,condSize,sep_sets,alpha);
        Adj(X, Y) = is_adj;
        Adj(Y, X) = is_adj;
    end
end
end


function [is_adj, sep_sets, nTests, condSize] = ...
    isAdjacent(D,X,Y,Mb_X,nTests,condSize, sep_sets,alpha)
%isAdjacent decides whether X and Y are adjacent, given the Markov blanket
%of the variable X, i.e., Mb_X. This is done through a search for
%separating set for X, Y in Mb_X. If such a set does not exist, X,Y are
%adjacent. Otherwise, the sep set is added to sep_sets.

is_adj = 1;
search_set = mysetdiff(Mb_X, Y);
for sep_size = 0:(length(search_set)-1)
    mb_subsets = subsets1(search_set,sep_size);
    for t=1:length(mb_subsets)
        S = mb_subsets{t};
        CI = CI_Test(X,Y,S,D,alpha);
        nTests = nTests+1;
        condSize(1+length(S)) = condSize(1+length(S)) + 1;
        if CI
            n = size(sep_sets,1);
            sep = false(1,n);
            sep(S) = true;
            sep_sets(X,Y,:) = sep;
            sep_sets(Y,X,:) = sep;
            is_adj = 0;
            return;
        end
    end
end
end


function [isR, Adj, nTests, condSize, sep_sets] = ...
    isRemovable(D,X,Mb_X,Adj,nTests,condSize,sep_sets,alpha)

size_sep_set = size(sep_sets,1);

% Adj(X):
adj_X = myintersect(find(Adj(X,:)), Mb_X);

% Testing removability
isR = true;   % first assume X is removable

% A set of variable pairs that need not be checked
skip_check = Adj>0;   % if Y,Z are connected, no need to check them.
for Y = adj_X
    for Z = mysetdiff(Mb_X, Y)
        if Adj(Y, Z) == 0 && skip_check(Y, Z) == 0
            if ~sep_sets(Y,Z,X)
                skip_check(Y, Z) = 1;
                skip_check(Z, Y) = 1;
                continue
            end
        end
    end
end

for Yind = 1:length(adj_X)
    Y = adj_X(Yind);
    for Z = mysetdiff(Mb_X, adj_X(1:Yind))
        if skip_check(Y,Z)
            continue
        end
        search_set = mysetdiff(Mb_X, [Y Z]);
        for sep_size = 0:length(Mb_X)-2
            mb_subsets = subsets1(search_set,sep_size);
            CI = false;
            for t=1:length(mb_subsets)
                cond_set = mb_subsets{t};
                
                CI = CI_Test(Y,Z,cond_set,D,alpha);
                nTests = nTests+1;
                condSize(1+length(cond_set)) = condSize(1+length(cond_set)) + 1;
                if CI
                    Adj(Y, Z) = 0;
                    Adj(Z, Y) = 0;
                    sepYZ = false(1,size_sep_set);
                    sepYZ(cond_set) = true;
                    sep_sets(Z,Y,:) = sepYZ;
                    sep_sets(Y,Z,:) = sepYZ;
                    skip_check(Y, Z) = 1;
                    skip_check(Z, Y) = 1;
                    break
                end
            end
            if CI
                break
            end
        end
        if skip_check(Y,Z)
            continue
        end
        for sep_size = 0:length(Mb_X)-2
            mb_subsets = subsets1(search_set,sep_size);
            CI = false;
            for t=1:length(mb_subsets)
                cond_set = mb_subsets{t};
                
                CI = CI_Test(Y,Z,[X cond_set],D,0.8);
                nTests = nTests+1;
                condSize(2+length(cond_set)) = condSize(2+length(cond_set)) + 1;
                if CI
                    Adj(Y, Z) = 0;
                    Adj(Z, Y) = 0;
                    sepYZ = false(1,size_sep_set);
                    sepYZ([X cond_set]) = true;
                    sep_sets(Z,Y,:) = sepYZ;
                    sep_sets(Y,Z,:) = sepYZ;
                    isR = false;
                    return
                end
            end
            if CI
                break
            end
        end
    end
end

end



function [Mb,Adj,V,nTests,condSize,sep_sets,remFlag] = ...
    UpdateProcess(D,X,Mb,Adj,V,nTests,condSize,sep_sets,alpha,remFlag)

Mb_X = find(Mb(X,:));

% Update V
V(X) = 0;

% Update Markov boundaries
for Y = Mb_X
    Mb(X,Y) = 0;
    Mb(Y,X) = 0;
    remFlag(Y) = true;
end

mbs = length(Mb_X);
for Yind=1:(mbs-1)
    Y = Mb_X(Yind);
    for Wind=(Yind+1):mbs
        W = Mb_X(Wind);
        if nnz(Mb(Y,:))>nnz(Mb(W,:))
            S = mysetdiff(find(Mb(W,:)),Y);
        else
            S = mysetdiff(find(Mb(Y,:)),W);
        end
        CI = CI_Test(Y,W,S,D,alpha);
        nTests = nTests+1;
        condSize(1+length(S)) = condSize(1+length(S)) + 1;
        if CI
            Mb(W,Y) = 0;
            Mb(Y,W) = 0;
            Adj(W,Y) = 0;
            Adj(Y,W) = 0;
            sep = false(1,size(sep_sets,1));
            sep(S) = true;
            sep_sets(W,Y,:) = sep;
            sep_sets(Y,W,:) = sep;
        end
    end
end
end

%**************************************************************************
%************************* Minor Functions ********************************
%**************************************************************************

function C = mysetdiff(A,B)
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }
% Original by Kevin Murphy, modified by Leon Peshkin

if isempty(A)
    C = [];
    return;
elseif isempty(B)
    C = A;
    return;
else % both non-empty
    bits = zeros(1, max(max(A), max(B)));
    bits(A) = 1;
    bits(B) = 0;
    C = A(logical(bits(A)));
end
end

function C = myintersect(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = myintersect(A,B)

A = A(:)'; B = B(:)';

if isempty(A)
    ma = 0;
else
    ma = max(A);
end

if isempty(B)
    mb = 0;
else
    mb = max(B);
end

if ma==0 || mb==0
    C = [];
else
    %bits = sparse(1, max(ma,mb));
    bits = zeros(1, max(ma,mb));
    bits(A) = 1;
    C = B(logical(bits(B)));
end
end

function sub_s = subsets1(s,k)
% SUBSETS1 creates sub-sets of a specific from a given set
% SS = subsets1(S, k)
%
% S is the given set
% k is the required sub-sets size
%
% Example:
%
% >> ss=subsets1([1:4],3);
% >> ss{:}
% ans =
%      1     2     3
% ans =
%      1     2     4
% ans =
%      1     3     4
% ans =
%      2     3     4
%
% Written by Raanan Yehezkel, 2004

if k<0 % special case
    error('subset size must be positive');
elseif k==0 % special case
    sub_s={[]};
else
    l=length(s);
    ss={};
    if l>=k
        if k==1 % Exit condition
            for I=1:l
                ss{I}=s(I);
            end
        else
            for I=1:l
                ss1=subsets1(s([(I+1):l]),k-1);
                for J=1:length(ss1)
                    ss{end+1}=[s(I),ss1{J}];
                end
            end
        end
    end
    sub_s=ss;
end
end