% addpath('/ifs/loni/faculty/shi/spectrum/qwang/util/sh')
% 
% GradTable = '/ifs/loni/faculty/shi/spectrum/yqiao/EEAJ/processed_ACPC/AD_010/Diffusion/GradientTable.txt';
% BValLowTHD = '1000';
% BValHighTHD = '2000';
% Data = '/ifs/loni/faculty/shi/spectrum/yqiao/EEAJ/processed_ACPC/AD_010/Diffusion/data.nii.gz';
% SPHMaxOrder = '8';
% NumOptiSteps = '100';

coefficient = ones(45,1);
indx = 1;
res = zeros(5,1);
for order = 0 : 4
    for count = 1:order*2+1
        res(order+1,1) = res(order+1,1) + coefficient(indx);
        indx = indx+1;
    end
end
