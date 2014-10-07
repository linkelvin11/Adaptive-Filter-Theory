function [ A,U,S,V,s ] = arraysvd( aoa,alpha,r,L,sn,lambda )
%ARRAYRCV Generate the Left-Singular Vector of an Antenna Array
%   N = number of sources (depth of matrix)
%   M = number of sensors (# of rows)
%   L = number of time samples (# of columns)
%   alpha = relative amplitude of signals
%   npwr = noise power relative to adb(1)
%   aoa = column vector of incidence angles in degrees

% initialize dimension variables
M = size(r,1);
N = size(aoa,1);

% convert aoa to rectangular
aoa = aoa/180*pi;
theta = aoa(:,1);
phi = aoa(:,2);
aoa = [sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta)];
k = 2*pi/lambda*aoa;

% generate random phase offset 
psi = 2*pi*rand(1,L,N);
psi = repmat(psi,M,1,1);
psi(:,1,:) = zeros(size(psi(:,1,:)));
phase = exp(1j*psi);

% combine amplitude and phase
alpha = 10.^(-alpha/20);
alpha = reshape(alpha,1,1,numel(alpha));
alpha = repmat(alpha,M,L,1);

% calculate k dot r
dir = reshape(r*k.',M,1,N);
dir = repmat(dir,1,L,1);
dir = exp(-1j*dir);
s = 1/sqrt(M)*dir(:,1,:);
s = reshape(s,M,N,1);

% combine everything
A = alpha.*phase.*dir;
A = sum(A,3);

% add noise
sn = 10^(-sn/10);
v = sqrt(sn/2)*randn(size(A)) + 1j*sqrt(sn/2)*randn(size(A));
A = A + v;

% svd
[U,S,V] = svd(A);

end