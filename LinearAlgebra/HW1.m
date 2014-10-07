%% Adaptive Filters Homework 1
% Kelvin Lin
clc; clear all; close all;
%% Problem 2
%b) P_u(x) = x;
%d)
M = 10;
n = 5;
a = randi([0 10],n,M) + 1j*randi([0 10],n,M);
b = randi([0 10],n,M) + 1j*randi([0 10],n,M);
spans = [chkSpan(a,b);chkSpan(a,a+b);chkSpan(b,a+b);chkSpan(2*a+b,2*a-b)];
%e) it is the set that contains only the zero vector;
%f)
A = randi([0 10],2,4);
A = [A;sum(A)];
rankc = rank(A);
%% Problem 2 part g
[U,S,V] = svd(A);
dS = (diag(S)>eps(max(max((S))))).*diag(S);
dS = dS(dS ~= 0);
dS = 1./dS;
dS = diag(dS);
if size(dS,1) < size(S,1);
    dS = [dS;zeros(size(S,1)-size(dS,1),size(dS,2))];
end
if size(dS,2) < size(S,2)
    dS = [dS,zeros(size(dS,1),size(S,2)-size(dS,2))];
end
dS = dS.';
Ainv = V*dS*U';
Pinv = pinv(A);
Compareinv = max(max(abs(Ainv - Pinv))) < 1e-6; % they match to within a tolerance
%2 - it is possible because it's over the complex field
rankA = rank(A); % However, A is not full rank, so it is not possible
%3
A1 = A*A';
A2 = A'*A;
eA1 = eig(A1); eA1 = eA1(eA1>eps('double')*max([size(eA1,1),size(eA1,2)])*max(max(eA1)));
eA2 = eig(A2); eA2 = eA2(eA2>eps('double')*max([size(eA2,1),size(eA2,2)])*max(max(eA2)));
sA1 = flipud(svd(A1)); sA1 = sA1(sA1>eps('double')*max([size(sA1,1),size(sA1,2)])*max(max(sA1)));
sA2 = flipud(svd(A2)); sA2 = sA2(sA2>eps('double')*max([size(sA2,1),size(sA2,2)])*max(max(sA2)));
% the corresponding values match
% 4)
RA = orthProj(A);
RAp = orthProj(null(A'));
NA = orthProj(null(A));
NAp = orthProj(A');
% 5)
PRAc = RA*RAp;
PRAc = PRAc.*(PRAc > eps('double')*max(size(PRAc)));
PNAc = NA*NAp;
PNAc = PNAc.*(PNAc > eps('double')*max(size(PNAc)));
% both PRAc and PNAc are zero matrices.

%% Problem 5
aoa = [40,0;80,50];
L = 100;
alpha = [0 5]';
sn = 20;
[x,y] = meshgrid(-2:2,-2:2);
r = [x(:),y(:),zeros(size(x(:)))];
d = 1;
lambda = d*3;

[A,U,S,V,s] = arraysvd(aoa,alpha,r,L,sn,lambda);

% calculate correlation matrices
covs = s'*s;
covsU = s'*U;
covsU = abs(covsU.*(covsU > eps('double')*max(size(covsU))));

% plot singular values
sv = diag(S);
figure
subplot(2,1,1)
stem(20*log10(sv));
title('Singular Values')
ylabel('Singular Values (dB)');

%estimate noise power
N = numel(alpha);
dS = diag(S); dS = dS(N+1:end);
nvar = sqrt(sum(dS.^2));
dS = diag(S); dS = [dS(1:N);nvar];
dS = dS.^2;
subplot(2,1,2)
stem(10*log10(dS/max(dS)));
title('Approximate Signal Power Relative to Max Power');
ylabel('Signal Power (dB)');

