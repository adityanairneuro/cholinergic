function [CC_All, CS_All] = simulation_cholinergic(i_curr,curr_length)
% Nair et al., 2020, Aditya Nair
    load('params.mat');
    CC_All = zeros(1,curr_length);
    CS_All = zeros(1,curr_length);
    jValues = 50:50:400;
    for jdx = 1:numel(jValues)
       j_curr = jValues(jdx);
         

    %% Define noise and input
    It = linspace(0,params.tfinal, params.n_timepoints);

    input_signal = (normrnd(20,20,[1,1200])).*((It>=0)&(It<100)) + ...
           (normrnd(i_curr,i_curr/20,[1,1200])).*((It>100)&(It<500)) + ...
           (normrnd(20,20,[1,1200])).*(It>200);

    input_noise = (normrnd(20,20,[1,1200])).*((It>=0)&(It<100)) + ...
           (normrnd(j_curr,j_curr,[1,1200])).*((It>100)&(It<500)) + ...
           (normrnd(20,20,[1,1200])).*(It>200);

%% Ctrl Signal

    params.ACh = 0;
    dynamics = integrate_dynamics_revised_cla(Wsoc, params, X0, It,input_signal);
    f_rates = dynamics.r;
    CC.signal_ctrl_f_t1000_ini_vec = f_rates(:,1:180);
    CS.signal_ctrl_f_t1000_ini_vec = f_rates(:,181:270);

%% Stim Signal
    params.ACh = 1;
    dynamics = integrate_dynamics_revised_cla(Wsoc, params, X0, It,input_signal);
    f_rates = dynamics.r;
    CC.signal_stim_f_t1000_ini_vec = f_rates(:,1:180);
    CS.signal_stim_f_t1000_ini_vec = f_rates(:,181:270);

%% Ctrl Noise

    params.ACh = 0;
    dynamics = integrate_dynamics_revised_cla(Wsoc, params, X0, It,input_noise);
    f_rates = dynamics.r;
    CC.noise_ctrl_f_t1000_ini_vec = f_rates(:,1:180);
    CS.noise_ctrl_f_t1000_ini_vec = f_rates(:,181:270);


%% Stim Noise

    params.ACh = 1;
    dynamics = integrate_dynamics_revised_cla(Wsoc, params, X0, It,input_noise);
    f_rates = dynamics.r;
    CC.noise_stim_f_t1000_ini_vec = f_rates(:,1:180);
    CS.noise_stim_f_t1000_ini_vec = f_rates(:,181:270);

%% CC Signal

    CCCtrl = CC.signal_ctrl_f_t1000_ini_vec(1:end,:);
    CCStim = CC.signal_stim_f_t1000_ini_vec(1:end,:);

    CCCtrlcorr = corr(CCCtrl);
    CCStimcorr = corr(CCStim);

    CCCtrlcorrT1 = triu(CCCtrlcorr);
    CCCtrlCorrDist1 =  single.empty;
    for i = 1:size(CCCtrl,2)

        for j = 1:size(CCCtrl,2)
            if CCCtrlcorrT1(i,j) == 0
            else

                CCCtrlCorrDist1 = [CCCtrlCorrDist1;CCCtrlcorrT1(i,j)];


            end
        end

    end

    CCStimcorrT2 = triu(CCStimcorr);
    CCStimCorrDist2 =  single.empty;
    for i = 1:size(CCCtrl,2)

        for j = 1:size(CCCtrl,2)
            if CCStimcorrT2(i,j) == 0
            else

                CCStimCorrDist2 = [CCStimCorrDist2;CCStimcorrT2(i,j)];


            end
        end

    end


    Signal_Corr_Ctrl_CC = CCCtrlCorrDist1;
    Signal_Corr_Stim_CC = CCStimCorrDist2;

%% CC Noise

    CCCtrl = CC.noise_ctrl_f_t1000_ini_vec(1:end,:);
    CCStim = CC.noise_stim_f_t1000_ini_vec(1:end,:);

    CCCtrlcorr = corr(CCCtrl);
    CCStimcorr = corr(CCStim);

    CCCtrlcorrT1 = triu(CCCtrlcorr);
    CCCtrlCorrDist1 =  single.empty;
    for i = 1:size(CCCtrl,2)

        for j = 1:size(CCCtrl,2)
            if CCCtrlcorrT1(i,j) == 0
            else

                CCCtrlCorrDist1 = [CCCtrlCorrDist1;CCCtrlcorrT1(i,j)];


            end
        end

    end

    CCStimcorrT2 = triu(CCStimcorr);
    CCStimCorrDist2 =  single.empty;
    for i = 1:size(CCCtrl,2)

        for j = 1:size(CCCtrl,2)
            if CCStimcorrT2(i,j) == 0
            else

                CCStimCorrDist2 = [CCStimCorrDist2;CCStimcorrT2(i,j)];


            end
        end

    end


    Noise_Corr_Ctrl_CC = CCCtrlCorrDist1;
    Noise_Corr_Stim_CC = CCStimCorrDist2;



%% CC Calc
    Signal_Noise_Ctrl = [Signal_Corr_Ctrl_CC,Noise_Corr_Ctrl_CC];
    Signal_Noise_Stim = [Signal_Corr_Stim_CC,Noise_Corr_Stim_CC];


    P_Ctrl = polyfit(Signal_Noise_Ctrl(:,1),Signal_Noise_Ctrl(:,2),1);
    CC_slope_ctrl = P_Ctrl(1);

    P_Stim = polyfit(Signal_Noise_Stim(:,1),Signal_Noise_Stim(:,2),1);
    CC_slope_stim = P_Stim(1);

    %% CS Signal

    CSCtrl = CS.signal_ctrl_f_t1000_ini_vec(1:end,:);
    CSStim = CS.signal_stim_f_t1000_ini_vec(1:end,:);

    CSCtrlcorr = corr(CSCtrl);
    CSStimcorr = corr(CSStim);

    CSCtrlcorrT1 = triu(CSCtrlcorr);
    CSCtrlCorrDist1 =  single.empty;
    for i = 1:size(CSCtrl,2)

        for j = 1:size(CSCtrl,2)
            if CSCtrlcorrT1(i,j) == 0
            else

                CSCtrlCorrDist1 = [CSCtrlCorrDist1;CSCtrlcorrT1(i,j)];


            end
        end

    end

    CSStimcorrT2 = triu(CSStimcorr);
    CSStimCorrDist2 =  single.empty;
    for i = 1:size(CSCtrl,2)

        for j = 1:size(CSCtrl,2)
            if CSStimcorrT2(i,j) == 0
            else

                CSStimCorrDist2 = [CSStimCorrDist2;CSStimcorrT2(i,j)];


            end
        end

    end


    Signal_Corr_Ctrl_CS = CSCtrlCorrDist1;
    Signal_Corr_Stim_CS = CSStimCorrDist2;

%% CS Noise


    CSCtrl = CS.noise_ctrl_f_t1000_ini_vec(1:end,:);
    CSStim = CS.noise_stim_f_t1000_ini_vec(1:end,:);

    CSCtrlcorr = corr(CSCtrl);
    CSStimcorr = corr(CSStim);

    CSCtrlcorrT1 = triu(CSCtrlcorr);
    CSCtrlCorrDist1 =  single.empty;
    for i = 1:size(CSCtrl,2)

        for j = 1:size(CSCtrl,2)
            if CSCtrlcorrT1(i,j) == 0
            else

                CSCtrlCorrDist1 = [CSCtrlCorrDist1;CSCtrlcorrT1(i,j)];


            end
        end

    end

    CSStimcorrT2 = triu(CSStimcorr);
    CSStimCorrDist2 =  single.empty;
    for i = 1:size(CSCtrl,2)

        for j = 1:size(CSCtrl,2)
            if CSStimcorrT2(i,j) == 0
            else

                CSStimCorrDist2 = [CSStimCorrDist2;CSStimcorrT2(i,j)];


            end
        end

    end


    Noise_Corr_Ctrl_CS = CSCtrlCorrDist1;
    Noise_Corr_Stim_CS = CSStimCorrDist2;

%% CS Calc

    Signal_Noise_Ctrl = [Signal_Corr_Ctrl_CS,Noise_Corr_Ctrl_CS];
    Signal_Noise_Stim = [Signal_Corr_Stim_CS,Noise_Corr_Stim_CS];


    P_Ctrl = polyfit(Signal_Noise_Ctrl(:,1),Signal_Noise_Ctrl(:,2),1);
    CS_slope_ctrl = P_Ctrl(1);

    P_Stim = polyfit(Signal_Noise_Stim(:,1),Signal_Noise_Stim(:,2),1);
    CS_slope_stim = P_Stim(1);


%%

    CC_All(jdx) = (CC_slope_ctrl - CC_slope_stim)/CC_slope_ctrl;
    CS_All(jdx) = (CS_slope_ctrl - CS_slope_stim)/CS_slope_ctrl;

% disp(j_index)
    end

CC_All;
CS_All;
end
    