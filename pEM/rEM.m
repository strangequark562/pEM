function [D S P Lmax] = rEM(deltaX,D0,P0,S0,params,trackInfo,numTrial)
%-------------------------------------------------------------------------- 
% Summary: This function employs reinitialization expectation maximization 
% for a Gaussian Mixture Model specific for single particle tracking data.
% The EM is remployed and the trial which yields the highest likelihood
% value is saved.
% 
% Input:
%       deltaX = cell with the particle track displacements
%       D0 = vector of initial diffusivity values for each diffusive state
%       S0 = vector of initial population fractions for each diffusive state
%       P0 = vecotr of initial static localization noise for each diffusive state
%       params = parameters for pEM
%               params.maxiter = maximum number of iterations for EM
%               params.converged = convergence condition for change in log-likelihood
%               params.numPerturbation = number of perturbations trials
%               params.showplot = displays progress of parameter estimates (0,1)
%               params.verbose = display progress on command window (0,1)
%       trackInfo = track information
%               trackInfo.numberOfTracks = number of particle tracks
%               trackInfo.dimensions = number of dimensions
%               trackInfo.dt = frame duration
%               trackInfo.R = motion blur coefficient
%   
% Output:
%        D = vector of diffusivity estimates
%        S = vector of static localization noise estimates
%        P = vector of population fraction estimates
%        Lmax = maximum likelihood found with optimal parameter estimates
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% apply reinitialization EM to protein trajectories
numStates = length(D0);
baseD = zeros(numTrial,numStates);
baseP = zeros(numTrial,numStates);
baseS = zeros(numTrial,numStates);
Lmax = zeros(numTrial,1);
for i = 1:numTrial
    disp('-------------------------------------------------------');
    disp(['EM trial #' num2str(i)]);

    % run em algorithm
    [D_est S_est P_est L] = EM(deltaX,D0,P0,S0,params,trackInfo);
    baseD(i,:) = D_est(end,:);
    baseS(i,:) = S_est(end,:);
    baseP(i,:) = P_est(end,:);
    Lmax(i) = L(end);
    
    disp(['D_k = ' num2str(baseD(i,:)) ' um^2/s']);
    disp(['sigma_k = ' num2str(baseS(i,:)) ' um']);
    disp(['pi_k = ' num2str(baseP(i,:)) ]);
    disp(['L = ' num2str(Lmax(i))]);
    disp('-------------------------------------------------------');

    % reinitialize EM parameters
    [D0 P0 S0] = RandomInitialization(length(D0),trackInfo.D_cve,trackInfo.sigma_cve);
end

% find maximum likelihood across trials
[MAX index] = max(Lmax);
D  = baseD(index,:);
P  = baseP(index,:);
S  = baseS(index,:);
Lmax = MAX;