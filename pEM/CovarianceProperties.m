function [diagonals correlations C] = CovarianceProperties(diffX)
%-------------------------------------------------------------------------- 
% Summary: This function calculates the empirical covariance properties for
% each particle track, assuming normal diffusion with a constant
% diffusivity and static localization noise throughout each particle track.
% Each dimension is calculated separately.
% 
% Input:
%       deltaX = cell with the particle track displacements
%
% Output:
%       diagonals = average diagonal elements of empirical covariance matrix
%       correlations = average nearest-neighbor correlations of empirical
%       covariance matrix for each particle track and each dimension
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% calculate MSD and correlation for all particle tracks
numTracks = length(diffX);
dim = size(diffX{1},2);
diagonals = zeros(numTracks,dim);
correlations = zeros(numTracks,dim);
C = cell(numTracks,dim);
for i = 1:numTracks
    for j = 1:dim
        deltaX = diffX{i}(:,j);
        C{i,j} = deltaX*deltaX';
        diagonals(i,j) = mean(deltaX.^2);
        correlations(i,j) = mean(deltaX(1:end-1).*deltaX(2:end));
    end
end 

