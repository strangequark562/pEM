function [baseD baseS baseP Lmax posteriorProb] = pEM(deltaX,baseD,baseP,baseS,Lmax,params,trackInfo)
%-------------------------------------------------------------------------- 
% Summary: This function employs perturbation expectation maximization 
% for a Gaussian Mixture Model specific for single particle tracking data.
% First, the EM is employed.  Then, the likelihood surface is perturbed
% with a bootstrap set of particle tracks.  The EM is remployed on the
% bootstrap set of particle tracks.  If the new parameters from the
% perturbed likelihood landscape yields a higher likelihood on the original
% likelihood landscape, then the parameters are updated.  This is repeated
% for a predetermined number of perturbations.
% 
% Input:
%       X = cell with the particle track displacements
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
%        baseD = vector of diffusivity estimates
%        baseS = vector of static localization noise estimates
%        baseP = vector of population fraction estimates
%        Lmax = maximum likelihood found with optimal parameter estimates
%        posteriorProb = posterior for M particle tracks and K diffusive 
%        states (M x K matrix)
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% begin perturbation EM to protein trajectories
for j = 1:params.numPerturbation
    
    % display current iteration
    disp(['Perturbation trial #' num2str(j)]);
    
    % perturb likelihood surface with bootstrap resampled data
    trackIndex = ceil(rand(trackInfo.numberOfTracks,1)*trackInfo.numberOfTracks);
    bootStrapDeltaX = cell(trackInfo.numberOfTracks,1);
    for i = 1:trackInfo.numberOfTracks
        bootStrapDeltaX{i} = deltaX{trackIndex(i)};        
    end

    % calculate new track properties of resampled data
    bsTrackInfo = trackInfo;
    bsTrackInfo.trackLength = trackInfo.trackLength(trackIndex);
    bsTrackInfo.diagonals = trackInfo.diagonals(trackIndex,:);
    bsTrackInfo.correlations = trackInfo.correlations(trackIndex,:);
    
    % Expectation-maximization with bootstrap data
    params2 = params;
    params2.showplot = 0;
    [D_est S_est P_est] = EM(bootStrapDeltaX,baseD,baseP,baseS,params2,bsTrackInfo);
    D = D_est(end,:);
    P = P_est(end,:);
    S = S_est(end,:);

    % calculate log-likelihood on original data set with bootstrap parameters
    [gamma Lnew] = Expectation(deltaX,D,P,S.^2,trackInfo);

    % if bootstrap parameters yield higher likelihood, then employ a full EM with new parameters
    if Lnew > Lmax
        disp('Higher likelihood found');
        disp('-------------------------------------------------------');

        [D_est S_est P_est L] = EM(deltaX,D,P,S,params,trackInfo);
        baseD = D_est(end,:);
        baseP = P_est(end,:);
        baseS = S_est(end,:);
        Lmax = L(end);
        
        
        disp(['D_k = ' num2str(baseD) ' um^2/s']);
        disp(['sigma_k = ' num2str(baseS) ' um']);
        disp(['pi_k = ' num2str(baseP) ]);
        disp(['L = ' num2str(Lmax)]);
        disp('-------------------------------------------------------');

    end
end

% calculate the posterior probability of each diffusive state
gamma = Expectation(deltaX,baseD,baseP,baseS.^2,trackInfo);
posteriorProb = 0;
for i = 1:trackInfo.dimensions
    posteriorProb = posteriorProb + gamma(:,:,i)/trackInfo.dimensions;
end


