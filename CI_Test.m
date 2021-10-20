function [CI, p] = CI_Test(X,Y,S,D,alpha)

% Input arguments:
% X,Y: the variables to check the conditional independence of
% S: conditioning set
% D: Data matrix (number_of_Samples * n)
% alpha: the parameter of Fischer's Z transform
%
% Output arguments:
% CI: the conditional indpendence relation between X,Y, true if
% independent, false if dependent
% p: p-value (The null hypothesis is for independence)
%--------------------------------------------------------------------------

n = size(D, 1);
DD = D(:,[X,Y,S]);
R = corrcoef(DD);
P = inv(R);
ro = abs(P(1,2)/sqrt(P(1,1)*P(2,2)));
zro = 0.5*log((1+ro)/(1-ro));
df = n-size(S,2)-3;
W = abs(zro)*sqrt(df);
p = 2*tcdf(-abs(W),df);
CI = p>alpha;
end
