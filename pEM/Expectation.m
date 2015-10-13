function  [gamma L] = Expectation(deltaX,D,P,S2,trackInfo)
%-------------------------------------------------------------------------- 
% Summary: This function calculates the posterior probabilities and
% likelihood for a given set of parameters.  
% 
% Input:
%       deltaX = cell with the particle track displacements
%       D = vector of starting value of diffusivities for each state 
%       P = vector of starting value of population fractions for each state
%       S2 = vector of starting value of static localization noise VARIANCE for each state
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
%       gamma = posterior probabilities for each particle track for each 
%       diffusive state.  gamma is a cell for each spatial dimension.
%       L = log-likelihood
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% get track population parameters
trackLength = trackInfo.trackLength;
uniqueLength = trackInfo.uniqueLength;
numTracks = trackInfo.numberOfTracks;

dim = trackInfo.dimensions;
numStates = length(D);

% calculate analytical covariance terms for each unique track length
diagterms = 2*D*trackInfo.dt + 2*S2 - 4*trackInfo.R*D*trackInfo.dt;
corrterms = -S2 + 2*trackInfo.R*D*trackInfo.dt;
invC = cell(uniqueLength(end),numStates);
detC = zeros(uniqueLength(end),numStates);
for i = 1:length(uniqueLength)
	for k = 1:numStates
		C = toeplitz([diagterms(1,k),corrterms(1,k),zeros(1,uniqueLength(i)-2)]);
        [logdetC invS]= LogDeterminant(C);        
        invC{uniqueLength(i),k} = invS;
		detC(uniqueLength(i),k) = logdetC;
	end
end

% calculate normalized posterior probability and log-likelihood
logN = zeros(numTracks,numStates,dim); 
for k = 1:numStates
    for i = 1:numTracks
        logN(i,k,:) = -.5*diag((deltaX{i}'*invC{trackLength(i),k}*deltaX{i}));
    end
end

L = 0; 
gamma = zeros(numTracks,numStates,dim); 
for j = 1:dim
    logpiN = logN(:,:,j) + ones(numTracks,1)*log(P) - trackLength*ones(1,numStates)/2*log(2*pi) - .5*detC(trackLength,:);

    % fix numerical issues associated with log-sum-exp
    MAX = max(logpiN,[],2);
    logpiN_norm = bsxfun(@minus,logpiN,MAX);
    logSumPiN = MAX + log(sum(exp(logpiN_norm),2));
    index = find(~isfinite(MAX));
    if ~isempty(index)
        logSumPiN(index) = MAX(index);
    end
    logr = logpiN - repmat(logSumPiN,1,numStates);
    gammank = exp(logr);
    gammank(find(gammank(:) < 1e-100)) = 1e-100;
    gamma(:,:,j) = gammank; 
    L = L + sum(logSumPiN);
end


