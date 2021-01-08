% This script performs the primary analysis in Nair et al., 2020 and
% reproduces Fig 3 E,F  & Fig S2 

% This makes use of the parallel computing toolbox in matlab, please change
% the parfor loop to a for loop if this toolbox is unavailable:

CC_All_parfor = zeros(8,8)
CS_All_parfor = zeros(8,8)


%% Loop!:
% Range of input current values
iValues = 50:50:400;
CC_iter = zeros(8,8,5);
CS_iter = zeros(8,8,5);

for iter = 1:10
    CC_All_parfor = zeros(8,8);
    CS_All_parfor = zeros(8,8);
    tic
parfor idx = 1:numel(iValues)

    i_curr = iValues(idx);
    [CC_All, CS_All] = simulation_cholinergic(i_curr,8)
    CC_All_parfor(idx,:) = CC_All;
    CS_All_parfor(idx,:) = CS_All;
end
    toc
    CC_iter(:,:,iter) = CC_All_parfor;
    CS_iter(:,:,iter) = CS_All_parfor;
    disp_X = ['Iteration ',num2str(iter)];
    disp(disp_X)
end

%% Analyse matrices
figure; subplot(1,2,1)
imagesc(-mean(CC_iter,3)*100);
xlabel('mean Noise')
ylabel('mean Signal')
title('CC')

 subplot(1,2,2);imagesc(-mean(CS_iter,3)*100);
xlabel('mean Noise')
ylabel('mean Signal')
title('CS')