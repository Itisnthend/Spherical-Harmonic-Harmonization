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
targetSite = targetSite(5:7);
referenceSite = {dir(fullfile(DataDir, referenceSiteFolder)).name};
referenceSite = referenceSite(3:5);


subDataTarget = [];
subGradTableTarget = cell(1, length(targetSite));
subGradTableRef = cell(1, length(referenceSite));
subMaskTarget = [];
subCoefTarget = [];


for index = 1:length(targetSite)
    dMRIDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/crgdata.nii.gz');
    %fprintf(dMRIDir);
    GradTableDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/subGradientTable.txt');
    MaskDir = fullfile(DataDir, targetSiteFolder, targetSite{index},'Diffusion/crgnodif_brain_mask.nii.gz');
    subDataTarget =  [subDataTarget load_untouch_nii(dMRIDir)];
    subCoefTarget = [subCoefTarget load_nii(dMRIDir)];
    subGradTableTarget(1,index) = {load(GradTableDir)};
    subMaskTarget = [subMaskTarget load_untouch_nii(MaskDir)];
end

subDataRef = [];
subCoefRef = [];
subMaskRef = [];
for index = 1:length(referenceSite)
    dMRIDir = fullfile(DataDir, referenceSiteFolder, referenceSite{index},'Diffusion/crgdata.nii.gz');
    GradTableDir = fullfile(DataDir, referenceSiteFolder, referenceSite{index},'Diffusion/subGradientTable.txt');
    MaskDir = fullfile(DataDir, referenceSiteFolder, referenceSite{index},'Diffusion/crgnodif_brain_mask.nii.gz');
    subDataRef =  [subDataRef load_untouch_nii(dMRIDir)];
    subCoefRef =  [subCoefRef load_nii(dMRIDir)];
    subGradTableRef(1,index) = {load(GradTableDir)};
    subMaskRef = [subMaskRef load_untouch_nii(MaskDir)];
end

NumOfSlices = subCoefTarget(1).hdr.dime.dim(4);
sz = subCoefTarget(1).hdr.dime.dim(2:5);
%[182,218,182,53]
%mapFunction = zeros(sz(1),sz(2),sz(3),SPHMaxOrder/2+1);

for indx = 1: length(targetSite)
    subCoefTarget(indx).hdr.dime.dim(2:5) = [sz(1) sz(2) sz(3) basic];
    subCoefRef(indx).hdr.dime.dim(2:5) = [sz(1) sz(2) sz(3) basic]; 
end

% signalO = [];
% signalh = [];
% signalind = find(subGradTable{1,1}(:,4)>BValLowTHD & subGradTable{1,index}(:,4)<BValHighTHD);
% signalO = subData(1,1).img(77,92,53,signalind);
for indx = 1: length(targetSite)
    subCoefTarget(indx).img = zeros(182, 218, 182,basic); 
end
for indx = 1: length(referenceSite)
    subCoefRef(indx).img = zeros(182, 218, 182,basic); 
end

mapFunction = zeros(sz(1),sz(2),sz(3),SPHMaxOrder/2+1);

for i = 1:sz(1)
    fprintf(int2str(i)+" ");
    for j = 1:sz(2)
        for k = 1:sz(3)
            if subMaskTarget(1).img(i,j,k) > 0 || subMaskTarget(2).img(i,j,k) > 0 || subMaskRef(1).img(i,j,k) > 0 || subMaskRef(2).img(i,j,k) > 0
            %if subMaskTarget(1).img(i,j,k) > 0 || subMaskRef(1).img(i,j,k) > 0
                tCoef = zeros(basic, length(targetSite));
                tRISH = zeros(SPHMaxOrder/2+1,length(targetSite));
                %tS_Org = zeros(46,2);
                for index = 1: length(targetSite)

                    %only load the specific slice
                    % nii = load_untouch_nii(dMRIDir,[],[],[],[],[],k);  % size: [121,145,1,99]
                    nii = subDataTarget(index).img(i,j,k,:);
                    % [tS_Org(:, index), tCoef(:,index)] = getSHCoefficient(subGradTable{1,index}, BValLowTHD,BValHighTHD, nii, i, j, SPHMaxOrder);
                    [tS_Org, tCoef(:,index)] = getSHCoefficient(subGradTableTarget{1,index}, BValLowTHD,BValHighTHD, nii, SPHMaxOrder);
                    %subCoefTarget(index).img(i,j,k,:) = tCoef(:,index);
                    tRISH(:,index) = getRISH(tCoef(:,index), SPHMaxOrder);

                end

                tRISH_bar = mean(tRISH, 2);

                rCoef = zeros(basic, length(referenceSite));
                rRISH = zeros(SPHMaxOrder/2+1,length(referenceSite));
                %rS_Org = zeros(46,2);
                for index = 1: length(referenceSite)

                    %only load the specific slice
                    nii = subDataRef(index).img(i,j,k,:);
                    [rS_Org, rCoef(:,index)] = getSHCoefficient(subGradTableRef{1,index}, BValLowTHD,BValHighTHD, nii, SPHMaxOrder);
                    subCoefRef(index).img(i,j,k,:) = rCoef(:,index);
                    rRISH(:,index) = getRISH(rCoef(:,index), SPHMaxOrder);

                end

                rRISH_bar = mean(rRISH, 2);

                temppie = rRISH_bar ./ tRISH_bar;
                [f, ~] = size(temppie);
                pie = zeros(f,length(targetSite));
                for num = 1: length(targetSite)
                    pie(:, num) = temppie;
                end
                pie = pie.^(1/2);
                mapFunction(i,j,k,:) = temppie.^(1/2);
                
                
                tCoef = getMappedCoef(tCoef, pie, SPHMaxOrder);
                
                S_hat = [];
                for index = 1: length(targetSite)

                    %only load the specific slice
                    %nii = load_untouch_nii(dMRIDir,[],[],[],[],[],k);
                    subCoefTarget(index).img(i,j,k,:) = tCoef(:,index);
                    S_hat = HarmMapping(tCoef(:,index), subGradTableTarget{1,index}, BValLowTHD, BValHighTHD, SPHMaxOrder);
                    % tRISH(:,index) = getRISH(tCoef(:,index), SPHMaxOrder); 
                    [numOfDirections, numOfSubjects] = size(S_hat);
                    for indx = 1:numOfDirections
                        if S_hat(indx,1) < 0
                            S_hat(indx,1) = 0;
                        end
                    end
                    DTI_Ind = find(subGradTableTarget{1,index}(:,4)>BValLowTHD & subGradTableTarget{1,index}(:,4)<BValHighTHD);
                    subDataTarget(index).img(i,j,k,DTI_Ind) = S_hat;

                end

                
            end
            
        end
        
    end
    
end

filename = fullfile(saveDir, targetSiteFolder, 'mapFunction.mat');
save(filename,'mapFunction');

for index = 1: length(targetSite)
    filename = fullfile(saveDir, targetSiteFolder, targetSite{index},'harmdata.nii.gz');
    save_untouch_nii(subDataTarget(index), filename);  
    filename = fullfile(saveDir, targetSiteFolder, targetSite{index},'tCoef.nii.gz');
    save_nii(subCoefTarget(index), filename); 
end

for index = 1: length(referenceSite) 
    filename = fullfile(saveDir, referenceSiteFolder, referenceSite{index},'rCoef.nii.gz');
    save_nii(subCoefRef(index), filename); 
end











