clear all;
clc;
close all;

%% load file

[filename,dirpath] = uigetfile('select file');
load(fullfile(dirpath,filename));

%% user set parameters

% movie parameters
dt = .032;  %frame interval
dE = .032; %exposure time

% rEM parameters
numReinitialize = 3;    % number of reinitializations  (number of initial pints to try to avoid getting stuck in a local maximum)  5 is Plenty

% pEM parameters
minStates = 1;          % minimum number of states to explore
maxStates = 20;         % maximum number of states to explore
numPerturb = 20;         % messing with the likelikhood surface to escpace from local minima...number of perturbation trials (use 100 if you want publish something like 5 to 20 to mess around)
maxiter = 10000;        % maximum number of iterations within EM trial
convergence = 1e-7;     % convergence criteria for change in log-likelihood 
showplot = 0;           % display the progress


%% run pEM

% structure for track info
trackInfo.numberOfTracks = length(X);   % number of tracks
trackInfo.dimensions = size(X{1},2);    % particle track dimensions
trackInfo.dt = dt;                      % frame duration
trackInfo.R = 1/6*dE/dt;                % motion blur coefficient

% structure for pEM
params.numPerturbation = numPerturb;    % number of perturbations trials
params.converged = convergence;         % convergence condition for EM
params.maxiter = maxiter;               % maximum number of iterations for EM
params.showplot = showplot;             % displays progress of parameter estimates (0,1)
params.verbose = 1;                     % display progress on command window (0,1)


% calculate the displacements for each particle track
deltaX = cell(trackInfo.numberOfTracks,1);
for i = 1:trackInfo.numberOfTracks
    deltaX{i} = diff(X{i});
end

% calculate relevant properties to enhance compuatational time
[trackInfo.trackLength trackInfo.uniqueLength] = TrackLengthParameters(deltaX);
[trackInfo.diagonals trackInfo.correlations trackInfo.C] = CovarianceProperties(deltaX);

% diffusivity and static localization estimate from covariance-based estimator
trackInfo.D_cve = mean((trackInfo.diagonals+2*trackInfo.correlations)/(2*trackInfo.dt),2);
trackInfo.sigma_cve = mean(trackInfo.diagonals,2)/2 - trackInfo.D_cve*trackInfo.dt*(1-2*trackInfo.R); 

% BIC Model Selection Loop
results = struct;
BIC = zeros(maxStates,1); 
for numStates = minStates:maxStates
    startTime = tic;
    disp([num2str(numStates) ' state model']);
    
    % random initialization
    [D0 P0 S0] = RandomInitialization(numStates,trackInfo.D_cve,trackInfo.sigma_cve);

    % run rEM
    [baseD baseS baseP Lmax] = rEM(deltaX,D0,P0,S0,params,trackInfo,numReinitialize);
    
    % run pEM
    [baseD baseS baseP Lmax posteriorProb] = pEM(deltaX,baseD,baseP,baseS,Lmax,params,trackInfo);

    % calculate BIC
    BIC(numStates) = Lmax - numStates*log(trackInfo.numberOfTracks);

    % display results    
    disp('-------------------------------------------------------');
    disp([num2str(numStates) ' state model results:']);
    disp(['D_k = ' num2str(baseD) ' um^2/s']);
    disp(['sigma_k = ' num2str(baseS) ' um']);
    disp(['pi_k = ' num2str(baseP) ]);
    disp(['L = ' num2str(Lmax)]);
    disp(['BIC = ' num2str(BIC(numStates))]);
    disp('-------------------------------------------------------');
    
    % store results
    results(numStates).numberOfStates = numStates;
    results(numStates).BIC = BIC(numStates);
    results(numStates).optimalD = baseD;
    results(numStates).optimalS = baseS;
    results(numStates).optimalP = baseP;
    results(numStates).optimalL = Lmax;
    results(numStates).posteriorProb = posteriorProb;
    results(numStates).elapsedTime = toc(startTime);
    
    % check BIC model selection
    if numStates > 1 
        if  BIC(numStates-1) ~= 0
            if BIC(numStates) < BIC(numStates-1)
                display(['Lower BIC found.  Optimal number of State:' num2str(numStates-1)]);
                break;
            end
            if numStates == maxStates
                display(['Optimal number of states not found is larger than ' num2str(maxStates)]);
                break;
            end
        end
    end
end
[MAX numStates] = max(BIC);

% store results
data.X = X;
data.params = params;
data.trackInfo = trackInfo;
data.results = results;
data.BIC = BIC(numStates);
data.optimalD = results(numStates).optimalD;
data.optimalP = results(numStates).optimalP;
data.optimalS = results(numStates).optimalS;
data.optimalL = results(numStates).optimalL;
data.posteriorProb = results(numStates).posteriorProb;

% display results
disp('Finished analysis');
disp('-------------------------------------------------------');
disp(['Optimal size: ' num2str(numStates) ' states']);
disp(['D_k = ' num2str(data.optimalD) ' um^2/s']);
disp(['sigma_k = ' num2str(data.optimalS) ' um']);
disp(['pi_k = ' num2str(data.optimalP) ]);
disp(['L = ' num2str(data.optimalL(end))]);
disp(['BIC = ' num2str(BIC(numStates))]);
disp('-------------------------------------------------------');


%% Display posterior-weighted tracks

DisplayPosteriorTracks(X,data.posteriorProb);

%% Display posterior-weighted MSD

numLags = 10;
DisplayWeightedMSD(X,posteriorProb,numLags,trackInfo.dt);


%%





