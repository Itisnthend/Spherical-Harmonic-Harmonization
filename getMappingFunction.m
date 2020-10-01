%addpath('/Users/itisnthend/INI/spectrum/qwang/util/sh')
addpath('/ifs/loni/faculty/shi/spectrum/qwang/util/sh')
addpath('/ifs/loni/faculty/shi/spectrum/qwang/matlab/NIfTI_20140122')
%DataDir = '/Users/itisnthend/INI/spectrum/qwang/HarmExperiment';
DataDir = '/ifs/loni/faculty/shi/spectrum/qwang/HarmExperiment';
%saveDir = '/Users/itisnthend/INI/spectrum/qwang/HarmExperiment';
saveDir = '/ifs/loni/faculty/shi/spectrum/qwang/HarmExperiment';

BValLowTHD = 2500;
BValHighTHD = 3500;
SPHMaxOrder = 8;
basic = (SPHMaxOrder+1)*(SPHMaxOrder+2)/2;

targetSiteFolder = 'HCPDevelopment';
referenceSiteFolder = 'HCP';
targetSite = {dir(fullfile(DataDir, targetSiteFolder)).name};
targetSite = targetSite(5:27);
referenceSite = {dir(fullfile(DataDir, referenceSiteFolder)).name};
referenceSite = referenceSite(3:27);

subGradTableTarget = cell(1, length(targetSite));
subGradTableReference = cell(1, length(referenceSite));
subMaskTarget = [];
subMaskReference = [];
subCoefTarget = [];
%subCoefReference = [];

for index = 1:length(targetSite)
    dMRIDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/crgdata.nii.gz');
    %fprintf(dMRIDir);
    GradTableDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/subGradientTable.txt');
    MaskDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/crgnodif_brain.nii.gz');
    %subCoefTarget = [subCoefTarget load_nii(dMRIDir)];
    subGradTableTarget(1,index) = {load(GradTableDir)};
    subMaskTarget = [subMaskTarget load_untouch_nii(MaskDir)];
end
dMRIDir = fullfile(DataDir, targetSiteFolder, targetSite{1},'Diffusion/crgdata.nii.gz');
subCoefTarget = load_nii(dMRIDir);


for index = 1:length(referenceSite)
    dMRIDir = fullfile(DataDir, referenceSiteFolder, referenceSite{index},'Diffusion/crgdata.nii.gz');
    GradTableDir = fullfile(DataDir, referenceSiteFolder, referenceSite{index},'Diffusion/subGradientTable.txt');
    subGradTableReference(1,index) = {load(GradTableDir)};
    %subCoefReference = [subCoefReference load_untouch_nii(MaskDir)];
end


NumOfSlices = subCoefTarget.hdr.dime.dim(4);
sz = subCoefTarget.hdr.dime.dim(2:5);
%[182,218,182,53]
mapFunction = zeros(sz(1),sz(2),sz(3),SPHMaxOrder/2+1);

% for index = 1: length(targetSite)
%     subCoefTarget(index).hdr.dime.dim(5) = basic;
%     subCoefTarget(index).img = zeros(sz(1),sz(2),sz(3),basic);
% end

for k = 1:NumOfSlices
    fprintf(int2str(k)+" ");
    flag = false;
    for index = 1:length(targetSite)
        nonZeroIndex = length(find(subMaskTarget(index).img(:,:,k)>0));
        if nonZeroIndex>0
            flag = true;
            break;
        end
    end
    if flag==true
        subDataTarget = [];
        for index = 1:length(targetSite)
            dMRIDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/crgdata.nii.gz');
            subDataTarget = [subDataTarget load_untouch_nii(dMRIDir,[],[],[],[],[],k)];
        end
        
        subDataReference = [];
        for index = 1:length(referenceSite)
            dMRIDir = fullfile(DataDir, referenceSiteFolder, referenceSite{index},'Diffusion/crgdata.nii.gz');
            subDataReference = [subDataReference load_untouch_nii(dMRIDir,[],[],[],[],[],k)];
        end


        for i = 1:sz(1)
            for j = 1:sz(2)
                if subDataTarget(1).img(i,j,1) > 0 || subDataTarget(2).img(i,j,1) > 0 || subDataReference(1).img(i,j,1) > 0 || subDataReference(2).img(i,j,1) > 0
                    tCoef = zeros(basic, length(targetSite));
                    tRISH = zeros(SPHMaxOrder/2+1,length(targetSite));
                    for index = 1: length(targetSite)
                        
                        %only load the specific slice
                        % nii = load_untouch_nii(dMRIDir,[],[],[],[],[],k);  % size: [121,145,1,99]
                        nii = subDataTarget(index).img(i,j,1,:);
                        % [tS_Org(:, index), tCoef(:,index)] = getSHCoefficient(subGradTable{1,index}, BValLowTHD,BValHighTHD, nii, i, j, SPHMaxOrder);
                        [tS_Org, tCoef(:,index)] = getSHCoefficient(subGradTableTarget{1,index}, BValLowTHD,BValHighTHD, nii, SPHMaxOrder);
                        tRISH(:,index) = getRISH(tCoef(:,index), SPHMaxOrder);

                    end
                    tRISH_bar = mean(tRISH, 2);
                    
                    rCoef = zeros(basic, length(referenceSite));
                    rRISH = zeros(SPHMaxOrder/2+1,length(referenceSite));
                    %rS_Org = zeros(46,2);
                    for index = 1: length(referenceSite)

                        %only load the specific slice
                        nii = subDataReference(index).img(i,j,1,:);
                        [rS_Org, rCoef(:,index)] = getSHCoefficient(subGradTableReference{1,index}, BValLowTHD,BValHighTHD, nii, SPHMaxOrder);
                        rRISH(:,index) = getRISH(rCoef(:,index), SPHMaxOrder);

                    end
                    rRISH_bar = mean(rRISH, 2);
                    temppie = rRISH_bar ./ tRISH_bar;
%                     [f, ~] = size(temppie);
%                     pie = zeros(f,length(targetSite));
%                     for num = 1: length(targetSite)
%                         pie(:, num) = temppie;
%                     end
%                     
%                     pie = pie.^(1/2);
                    mapFunction(i,j,k,:) = temppie.^(1/2);
%                     tCoef = getMappedCoef(tCoef, pie, SPHMaxOrder);
                    
                    %S_hat = [];
%                     for index = 1: length(targetSite)
%                         subCoefTarget(index).img(i,j,k,:) = tCoef(:,index);
%                         %S_hat = HarmMapping(tCoef(:,index), subGradTableTarget{1,index}, BValLowTHD, BValHighTHD, SPHMaxOrder);
%                         
% %                         [numOfDirections, numOfSubjects] = size(S_hat);
% %                         for indx = 1:numOfDirections
% %                             if S_hat(indx,1) < 0
% %                                 S_hat(indx,1) = 0;
% %                             end
% %                         end
% %                         DTI_Ind = find(subGradTableTarget{1,index}(:,4)>BValLowTHD & subGradTableTarget{1,index}(:,4)<BValHighTHD);
% %                         B0_Ind = find(subGradTableTarget{1,index}(:,4)<100);
% %                         S0 = mean(subData(1,index).img(i,j,k,B0_Ind),4);
% %                         subDataTarget(index).img(i,j,k,DTI_Ind) = S_hat*S0;
%                     end
                    
                    
                end
            end
        end
    end
end


% for index = 1: length(targetSite)
% %     filename = fullfile(saveDir, targetSiteFolder, targetSite{index}, 'Diffusion/harmdata.nii.gz');
% %     save_untouch_nii(subDataTarget(index), filename);  
%     filename = fullfile(saveDir, targetSiteFolder, targetSite{index}, 'Diffusion/mappedtCoef.nii.gz');
%     save_nii(subCoefTarget(index), filename); 
% end

filename = fullfile(saveDir, targetSiteFolder, 'mapFunction.mat');
save(filename,'mapFunction');






