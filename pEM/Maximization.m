function [D p S2] = Maximization(gamma,diagonals,correlations,dt,R)
%-------------------------------------------------------------------------- 
% Summary: This function maximizes the log-likelihood for the diffusivities, 
% static localization noise variance, and population fractions for given a 
% posterior probability, gamma.
% 
% Input:
%       deltaX = cell with the particle track displacements
%       gamma = posterior probability otherwise known as responsibilities
%       trackInfo = track information
%               trackInfo.numberOfTracks = number of particle tracks
%               trackInfo.dimensions = number of dimensions
%               trackInfo.dt = frame duration
%               trackInfo.R = motion blur coefficient
%               trackInfo.trackLength = vector of particle track lengths 
%               trackInfo.uniqueLength = vector of unique track lengths
%               trackInfo.diagonals = vector of empirical covariance matrix diagonals
%               trackInfo.correlations = vector of empirical covariance matrix correlations
%
% Output:
%       D = maximized estimates of diffusivities
%       P = maximized estimates of population fractions
%       S2 = maximized estimates of static localization noise VARIANCE
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% get track population parameters
[numTracks numStates dim] = size(gamma);

% population fraction maximization
Nk = zeros(dim,numStates);
for i = 1:dim
    Nk(i,:) = sum(gamma(:,:,i));
end
meanNk = mean(Nk,1);
p = meanNk/numTracks;

% diffusion coefficient maximization
diagterms = 0; corrterms = 0;
for j = 1:dim
    diagterms = diagterms + sum(gamma(:,:,j).*(diagonals(:,j)*ones(1,numStates)))./Nk(j,:);
    corrterms = corrterms + sum(gamma(:,:,j).*(correlations(:,j)*ones(1,numStates)))./Nk(j,:);
end
diagterms = diagterms/dim;
corrterms = corrterms/dim;
D = max((diagterms + 2*corrterms)/2/dt,1e-6);

% localization noise variance maximization
S2 = max((diagterms - 2*D*dt*(1-2*R))/2,1e-8);

