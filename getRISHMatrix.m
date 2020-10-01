function [tCoef,tRISH] = getRISHMatrix(DataDir, Site, BValLowTHD,BValHighTHD, i, j, k, SPHMaxOrder)
    tCoef = zeros(basic, length(Site));
    tRISH = zeros(SPHMaxOrder/2+1,length(Site));
    for index = 1: length(Site)
        dMRIDir = fullfile(DataDir, Site{index},'Diffusion/data.nii.gz');
        GradTableDir = fullfile(DataDir, Site{index},'Diffusion/GradientTable.txt');
        MaskDir = fullfile(DataDir, Site{index},'Diffusion/nodif_brain.nii.gz');

        i = 77;
        j = 92;
        k = 53;

        GradTable = load(GradTableDir);
        Mask = load_untouch_nii(MaskDir);

        if Mask.img(i,j,k) > 0
            %only load the specific slice
            nii = load_untouch_nii(dMRIDir,[],[],[],[],[],k);
            tCoef(:,index) = getSHCoefficient(GradTable, BValLowTHD,BValHighTHD, nii, i, j, SPHMaxOrder);
            tRISH(:,index) = getRISH(tCoef(:,index), SPHMaxOrder);

        else
            fprintf(join(['Index of subject ', Site{1}, ' exceed brain region.','\n']));
        end
    end
end

