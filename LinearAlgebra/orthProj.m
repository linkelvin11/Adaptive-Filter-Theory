function [ T ] = orthProj( A )
%ORTHPROJ generates an orthogona projection onto the rangespace of A
T = orth(A)*orth(A).';
end

