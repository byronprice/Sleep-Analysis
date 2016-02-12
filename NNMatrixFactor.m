function [B,C] = NNMatrixFactor(Matrix,K,Tolerance,maxITER)
%NNMatrixFactor.m
%   Perform non-negative matrix factorization according to 
%    Lee & Seung, Algorithms for Non-negative Matrix Factorization
%    Matrix = B*C
%     if Matrix has size I-by-J and K < I && K < J
%     then B will be I-by-K and C K-by-J
%
%INPUT: Matrix - the matrix to be factorized
%       K - size of second dimension of B and first of C
%       Tolerance - difference in the cost from one step to the next such
%        that the algorithm will stop
%       maxITER - algorithm will not do more than maxITER steps
%        if desired tolerance is not reached, a warning will display
%
%OUTPUT: B & C - factorized matrices
%
%Created: 2016/02/12, 24 Cummington, Boston
% Byron Price
%Updated: 2016/02/12
% By: Byron Price

I = size(Matrix,1);
J = size(Matrix,2);
minVal = min(min(Matrix));
maxVal = max(max(Matrix));
if K > I || K > J
    display('K must be less than I and J')
end
B = randi([minVal,maxVal],[I,K]);
C = randi([minVal,maxVal],[K,J]);

cost = zeros(maxITER,1);
for iter = 2:maxITER
    C = C.*((B'*Matrix)./(B'*B*C));
    B = B.*((Matrix*C')./(B*C*C'));
    B(B<0) = 0;
    C(C<0) = 0;
    cost(iter) = norm(Matrix-B*C,2);
    if abs(cost(iter)-cost(iter-1)) < Tolerance
        display(sprintf('Reached tolerance after %i iterations',iter'))
        return;
    end
end
if iter == maxITER
    display('Desired tolerance not achieved')
end
end

