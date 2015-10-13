function  [D_est S_est P_est L] = EM(deltaX,D,P,S,params,trackInfo)
%-------------------------------------------------------------------------- 
% Summary: This function applies the expectation-maximization to find the 
% maximum log-likelihood for single particle tracking data.  Starting from 
% initial parameter values for the diffusivity, static localization noise, 
% and population fraction for each diffusive state, the EM learns the 
% parameters which maximized the log-likelihood of a Gaussian Mixture Mode.
% 
% Input:
%       deltaX = cell with the particle track displacements
%       D = vector of starting guess of diffusivities for each state 
%       P = vector of starting guess of population fractions for each state
%       S = vector of starting guess of static localization noise for each state
%       params = parameters for EM
%               params.converged = log-likelihood convergence criterion
%               params.maxiter = maximum number of iterations for EM
%               params.showplot = plots the parameter estimates (0,1)
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
%        D_est = diffusivities estimates at each iteration 
%        P_est = population fraction estimates at each iteration
%        S_est = static localization noise estimates at each iteration 
%        L = log-likelihood at each iteration for parameters at that iteration
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% set the temperature schedule
showplot = params.showplot;
converged = params.converged;
maxiter = params.maxiter;

diagonals = trackInfo.diagonals;
correlations = trackInfo.correlations;
dt = trackInfo.dt;
R = trackInfo.R;

% initialize variables
numStates = length(D);
D_est = zeros(maxiter,numStates); 
P_est =  zeros(maxiter,numStates); 
S2_est = zeros(maxiter,numStates); 
L = zeros(maxiter,1);

% setup plot (optional)
if showplot == 1
    h = figure(500); clf; box on;
    set(h,'position',[20 20 1200 300]); 
end
  
% expectation step with a initial parameters
[gamma Lnew] = Expectation(deltaX,D,P,S.^2,trackInfo);
 
% expectation-maximization at each annealing temperature
D_est(1,:) = D; P_est(1,:) = P; S2_est(1,:)= S.^2; L(1) = Lnew;
for i = 2:maxiter    
    
    % maximization step
    [D P S2] = Maximization(gamma,diagonals,correlations,dt,R);

    % expectation step    
    [gamma L(i,1)] = Expectation(deltaX,D,P,S2,trackInfo);
        
    % check for convergence of the log-likelihood
    if L(i) - L(i-1) < converged 
        break;
    else
        % if not converged, then store updated parameters
        D_est(i,:) = D;
        P_est(i,:) = P;
        S2_est(i,:) = S2;
    end
    
    % plot parameters
    if showplot == 1
        figure(h);
        subplot(1,3,1); plot(D_est(1:i,:)); title('D_k'); xlabel('iteration'); ylabel('Diffusivity [\mum^2s^{-1}]');
        subplot(1,3,2); plot(P_est(1:i,:)); title('\pi_k'); xlabel('iteration'); ylabel('Population fraction');
        subplot(1,3,3); plot(L(1:i)); title('L'); xlabel('iteration'); ylabel('log-Likelihood');
    end
end 

% clean parameter estimates
index = find(D_est(:,1)==0,1,'first');
D_est = D_est(1:index-1,:);
P_est = P_est(1:index-1,:);
S_est = sqrt(S2_est(1:index-1,:));
L = L(1:index-1);

% sort diffusivities and store parameter values at iterations
[D index] = sort(D_est(end,:));
D_est = D_est(:,index);
S_est = S_est(:,index);
P_est = P_est(:,index); 


