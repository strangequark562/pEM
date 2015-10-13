function [D0 P0 S0] = RandomInitialization(numStates,D_cve,sigma_cve)

% random population fractions greater than 5%
P0 = zeros(1,numStates);
for i = 1:numStates
    P0(i) = rand;
end
P0 = P0/sum(P0);

% CDF of diffusivities to generate a good random initialization
bin = 1000;
MIN = min(D_cve);
MAX = max(D_cve);
edges = MIN:(MAX-MIN)/bin:MAX;
DN = histc(D_cve,edges);
cumprob = cumsum(DN)/sum(DN);

% mean location of each random population fraction along CDF 
Dpoints = cumsum([0 P0(1:end-1)] + diff([0 P0])/2);

% find which diffusivity corresponds to each random population
% fraction
D0 = zeros(1,numStates);
for i = 1:numStates
    D0(i) = edges(find(cumprob > Dpoints(i),1,'first'));
end

S0 = ones(1,numStates)*mean(sigma_cve);