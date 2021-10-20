path(pathdef)
close all
clear
clc


% load data
mag_skelet = load("data.mat").mag;
data = load("data.mat").data;

% parameters
alpha = 0.01;
n = size(mag_skelet,1);
alpha_mb = 2/(n^2);

% Markov boundary discovery
[Mb,tests_mb,sep_sets] = Mb_TC(data, alpha_mb);

% Call L-MARVEL function
[Adj,sep_sets,nTests,condSize] = LMARVEL(data,Mb,sep_sets,alpha,alpha_mb);

% number of CI tests and accuracy of L-MARVEL
avg_size_cond_set = ((0:n)*condSize)/nTests;
[extra_edges,missing_edges,precision,recall,F1_score]=...
    learning_errors(mag_skelet,Adj);

fprintf("number of CI tests: %d", nTests);
fprintf("\nprecision: %.2f", precision);
fprintf("\nrecall: %.2f", recall);
fprintf("\nF1 score: %.2f\n", F1_score);
