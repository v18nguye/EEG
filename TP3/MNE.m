function s_h = MNE(x,A,lambda)

% Estimate the source signal by using MNE

[n,~] = size(A);
s_h = A.'*inv((A*A.' + lambda*eye(n)))*x;

end