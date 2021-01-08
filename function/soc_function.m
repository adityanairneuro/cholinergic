function Wsoc = soc_function(W_initial, rate, desired_SA, gamma, end_exc)
% Code from Vogels et al., 2018, Nature Neuroscience
% By Jake Stroud

%% Create a stability-optimised circuit

Wsoc = W_initial;

% Plot spectrum of W_initial
figure
subplot(2,1,1)
plot(eig(Wsoc),'*');
xlim([-10 10]);
ylim([-10 10]);
xlabel('Real part');
ylabel('Imaginary part');

i = 1;
e = max(real(eig(Wsoc)));   %Initial largest real part in the spcetrum of Wsoc
while e(i) > desired_SA
    
    i = i+1; %Output iteration number
    
    Wsoc = ssaCode(Wsoc,rate,gamma,end_exc);  %Run one iteration of the gradient descent
    e(i) = max(real(eig(Wsoc)));
    
    if mod(i,20) == 0 %Plot spectrum of Wsoc every 20 iterations
        
        subplot(2,1,1)
        plot(eig(Wsoc),'*');
        xlim([-10 10]);
        ylim([-10 10]);
        xlabel('Real part');
        ylabel('Imaginary part');
        
        subplot(2,1,2)
        plot(e)
        ylabel('Spectral abscissa')
        xlabel('Number of iterations')
        spec_abscissa = e(i)    %Print the current maximum real eigenvalue
        pause(0.01)
        
    end
    
end


end

%% Code that produces one gradient descent iteration
function [Wo, Emax] = ssaCode(Wi,rate,gamma, end_exc)

N = length(Wi);         %Number of neurons
start_inh = end_exc + 1;    %Index of first inhibitory neuron

Emax = max(real(eig(Wi)));

% Setup lyapunov equations
s = max(Emax*1.5, Emax + 0.2);

A = Wi - s*eye(size(Wi));
X = 2*eye(size(Wi));

%Solve the lyapunov equations to obtain the observability and
%controllability gramians
Q = lyap(A',X); 
P = lyap(A,X);

%Calculate gradient with which to move the inhibitory weights
grad = Q*P/(trace(Q*P)); 

Wo = Wi;

%Change the inhibitory weights according to the gradient
Wo(:,start_inh:end) = Wi(:,start_inh:end) - (rate*grad(:,start_inh:end));

%% Now perform all constraints

%Set any positive inhibitory weights to 0
I = Wo > 0; 
I(:,1:end_exc) = 0;
Wo(I) = 0;

%Make all inhibitory weights on average gamma times stronger than the
%excitatory ones.

exc_to_exc = Wo(1:end_exc,1:end_exc);
meanEE = mean(exc_to_exc(:));

exc_to_inh = Wo(start_inh:end,1:end_exc);
meanEI = mean(exc_to_inh(:));

inh_to_exc = Wo(1:end_exc,start_inh:end);
meanIE = mean(inh_to_exc(:));

inh_to_inh = Wo(start_inh:end,start_inh:end);
meanII = mean(inh_to_inh(:));

Wo(1:end_exc,start_inh:end) = -gamma*(meanEE/meanIE) * Wo(1:end_exc,start_inh:end);
Wo(start_inh:end,start_inh:end) = -gamma*(meanEI/meanII) * Wo(start_inh:end,start_inh:end);

%Make only 40% of the inh weights non-zero by setting the smallest
%inhitiroy weights to 0. (This step is slightly different to that used in
%Hennequin et al., Neuron, 2014.)

inh_weights = Wo(:,start_inh:end);
inh_weights = inh_weights(:);
[~,I] = sort(inh_weights,'ascend');
thres = round(0.4*length(inh_weights));
inh_weights(I(thres:end)) = 0;
Wo(:,start_inh:end) = reshape(inh_weights,N,N-end_exc);

% Remove any self-loops
Wo = Wo - diag(diag(Wo));

end