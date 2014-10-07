function [ B ] = chkSpan( A, B )
%CHKSPAN Checks if two sets of vectors span the same space
%   Check if the projection onto the rangespace of both matrices is the
%   same. Some numerical tolerance is needed to to quantization error
if sum(sum(orthProj(A)-orthProj(B))) < 1e-3
   B = 1;
   return;
else
B = 0;  
end
end

