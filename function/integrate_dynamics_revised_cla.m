% Integrate Dynamics
% Nair et al., 2021, Aditya Nair
function [output] = integrate_dynamics_revised_cla(W, params, x0,It, input_c)

    %It = linspace(0,params.tfinal, params.n_timepoints);
   
    %input_c = input_C(x0,It);
    
    [output.t, X] = ode45(@rate_dynamics_ode, ...
        linspace(0, params.tfinal, params.n_timepoints), x0, [], W, params, It, input_c);

    % Convert neuronal activities into firing rates
    output.r = mixed_gain_2(X, params.ACh);
    
    output.X = X;

end

function x_dot = rate_dynamics_ode(t, X, W, params, It, input_c)
   
    input_c = interp1(It, input_c, t);
    x_dot = params.over_tau*(-X + input_c+ W*mixed_gain(X,params.ACh));

end

% % Input Current
% function I = input_C(X, t)
% 
% if t < 100
%    I = zeros(length(X),100);  
% elseif (t >= 100) & (t <= 200)
% 
%     I = zeros(length(X),100);  
%     I(1:end,100) = normrnd(400,10,[length(X),100]);
%         
% %elseif (t >= 700) && (t <= 750)
% %    I = zeros(size(X));  
% %    I(1:end,1) =  normrnd(400,10);
% 
% else
%    I=zeros(length(X),1);  
%    
% end
% 
% end


% Gain function
function out_X = mixed_gain(X, ACh)
    CC = X(1:180);
    CS = X(181:270);
    VIP = X(271:280);
    SST = X(281:290);
    PV = X(291:300);
    
    if ACh == 0
        CC_out = gain_fn(CC,20,44.1,0.3243);
        CS_out = gain_fn(CS,20,45.38,0.2092);
        VIP_out = gain_fn(VIP,20,48.21,0.2799);
    else
        CC_out = gain_fn(CC,20,43.53,0.2098);
        CS_out = gain_fn(CS,20,40.68,0.3221);
        VIP_out = gain_fn(VIP,20,48.92,1.117);
    end
    
    SST_out = gain_fn(SST,20, 71.53,0.536);
    PV_out = gain_fnPV(PV,20,100,1.826);
    out_X = [CC_out;CS_out;VIP_out;SST_out;PV_out];
end

function out = gain_fn(x, r0, rmax, g)
    out = zeros(size(x));
    for n = 1 : length(x)
        
        if x(n) < 0
        out(n) = r0*tanh(g*x(n)/r0); 
        else
        out(n) = (rmax-r0)*tanh(g*x(n)/(rmax-r0));    
            
        end
                
    end
end

function out = gain_fnPV(x, r0, rmax, g)
    out = zeros(size(x));
    x1 = x - 140;
    for n = 1 : length(x1)
        
        if x1(n) < 0
        out(n) = r0*tanh(g*x1(n)/r0); 
        else
        out(n) = (rmax-r0)*tanh(g*x1(n)/(rmax-r0));    
            
        end
                
    end
end


function out_X = mixed_gain_2(X, ACh)
     CC = X(:,1:180);
    CS = X(:,181:270);
    VIP = X(:,271:280);
    SST = X(:,281:290);
    PV = X(:,291:300);
    
    if ACh == 0
        CC_out = gain_fn_2(CC,20,44.1,0.3243);
        CS_out = gain_fn_2(CS,20,45.38,0.2092);
        VIP_out = gain_fn_2(VIP,20,48.21,0.2799);
    else
        CC_out = gain_fn_2(CC,20,43.53,0.2098);
        CS_out = gain_fn_2(CS,20,40.68,0.3221);
        VIP_out = gain_fn_2(VIP,20,48.92,1.117);
    end
    
    SST_out = gain_fn_2(SST,20, 71.53,0.536);
    PV_out = gain_fn_2PV(PV,20,100,1.826);
    out_X = [CC_out,CS_out,VIP_out,SST_out,PV_out];
end

function out = gain_fn_2(x, r0, rmax, g)
    out = zeros(size(x));
    for row = 1 : size(x,1)
        for col = 1 : size(x,2)
            if x(row,col) < 0
              out(row,col) =  r0*tanh(g*x(row,col)/r0);
               
            else
              out(row,col) = (rmax-r0)*tanh(g*x(row,col)/(rmax-r0));    
            end
                       
        end
    end
end

function out = gain_fn_2PV(x, r0, rmax, g)
    out = zeros(size(x));
    x1 = x - 140;
    for row = 1 : size(x1,1)
        for col = 1 : size(x1,2)
            if x1(row,col) < 0
              out(row,col) =  r0*tanh(g*x1(row,col)/r0);
               
            else
              out(row,col) = (rmax-r0)*tanh(g*x1(row,col)/(rmax-r0));    
            end
                       
        end
    end
end