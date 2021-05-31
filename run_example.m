% Run examples
clear
%% Simple EIV linear model
data = datagenSLR(1234,4,250,0.15,0.25,1);
% Compute Estimate
est = estimatorSLR(data.Y, data.W);
% Compute covariance matrix of estimates for the phase function estimates using a modified
% bootstrap procedure
var_est = phase_variance(data.Y,data.W,est.tGrid,est.bPhase,100);

%% Multiple EIV linear model
clear
rng('default')
data = datagenMLR(1234,1,500,0.075,0.15,1,0.5);

est = estimatorMLR(data.Y, [data.W data.Z]);  % This does not include the intercept
