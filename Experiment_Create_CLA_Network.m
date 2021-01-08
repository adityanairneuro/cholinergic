% Create a SOC Claustrum Network

gamma = 9;
rate = 10;
desired_SA = 0.3;

% Create claustrum network with connectivity from literature:
% E-E connectivity: 0.03 (Kim et al., JNeurosci 2014)
% E-I & I-E connectivity: 0.4 (Kim et al., JNeurosci 2014)
% I-I connectivity: 0.5 (Kim et al., JNeurosci 2014, Chia et al., Current Biology 2020)

W0 = create_matrix(270,30,0.03,0.4,0.4,0.5,10,gamma);

% Perform optimization for creating SOC
Wsoc = soc_function(W0, rate, desired_SA, gamma, 270);