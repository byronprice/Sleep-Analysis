function [ CDF,alpha ] = Power_Law_Estimation(X,xmax)
%Power_Law_Estimation.m
%   Fits a discrete power-law distribution, cumulative distribution, 
%    exponent equal to alpha.
% INPUT: X - x values at which to calculate maximum likelihood estimator
%       for alpha
%        xmax - max value at which to calculate the CDF
% 
% OUTPUT - CDF - cumulative distribution function using the calculated
%           maximum likelihood estimator
%          alpha - maximum likelihood estimator for the exponent alpha
%

% CALCULATE MAXIMUM LIKELIHOOD ESTIMATOR FOR ALPHA
xmin = min(X);
n = length(X);
W = 5000; % for eval of zeta function
alphas = 1:0.01:5;
L = zeros(length(alphas),1);
count = 1;
for alpha = 1:0.01:5
    zeta = zeros(W,1);
    zeta(1) = (0+xmin)^(-alpha);
    for j=1:W
        zeta(j+1) = zeta(j)+(j+xmin)^(-alpha);
    end
    L(count) = -n*log(zeta(W))-alpha*sum(log(X));
    count = count+1;
end
[~,I] = max(L);
alpha = alphas(I);

% CALCULATE THEORETICAL CDF
zeta_bottom = zeros(W,1);
zeta_bottom(1) = (0+xmin)^(-alpha);
for i=1:W
    zeta_bottom(i+1) = zeta_bottom(i)+(i+xmin)^(-alpha);
end
zeta_top = zeros(xmax,1);
for xval = 1:xmax
    zeta = zeros(W,1);
    zeta(1) = (0+xval)^(-alpha);
    for i=1:W
        zeta(i+1) = zeta(i) + (i+xval)^(-alpha);
    end
    zeta_top(xval) = zeta(W);
end
cCDF = zeta_top./zeta_bottom(W);

CDF = 1-cCDF;


end
