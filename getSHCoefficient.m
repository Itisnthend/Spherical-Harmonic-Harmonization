function [S,C] = getSHCoefficient(GradTable,BValLowTHD,BValHighTHD, dMRI, SPHMaxOrder)

B0_Ind = find(GradTable(:,4)<100);

DTI_Ind = find(GradTable(:,4)>BValLowTHD & GradTable(:,4)<BValHighTHD);
GradVec = GradTable(DTI_Ind,1:3);
BVal = GradTable(DTI_Ind,4);


params.MDirection = GradVec;
params.maxOrder = SPHMaxOrder;
params.bval = BVal;
params.bval = round(params.bval/50)*50; %Assume all bvals are rounded to 50s. For HCP, bvals are around 1000,2000,3000. To speed up the matrixG_Multishell and matrixG_grad_MultiShell calculation. 


%Generate matrices with consisits of samples of spherical harmonics on
%point samples on the unit sphere
%This is the matrix for point samples of gradient directions
params.B = matrixB(params.MDirection,params.maxOrder);
params.s = zeros(length(params.MDirection),1);
%R = (params.maxOrder+1)*(params.maxOrder+2)/2; % number of SH basis functions

% %calculate FOD reconstruction for the selected slice.

S0 = mean(dMRI(1,1,1,B0_Ind),4);       % S0 (121*145) double

params.bval = round(params.bval/50)*50; %Assume all bvals are rounded to 50s. For HCP, bvals are around 1000,2000,3000. To speed up the matrixG_Multishell and matrixG_grad_MultiShell calculation. 
%Recompute B Matrix
params.B = matrixB(params.MDirection,params.maxOrder);

%Data samples on gradient directions
params.s(:) = double(dMRI(1,1,1,DTI_Ind));
params.s = params.s./double(S0);
if max(params.s)>1
    params.s = params.s./max(params.s);
end

coefficient = inv(params.B'*params.B)*params.B'*params.s;
C = coefficient;
S = params.s;
end
