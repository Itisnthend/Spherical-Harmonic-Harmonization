function res = HarmMapping(coefficient, GradTable, BValLowTHD, BValHighTHD, SPHMaxOrder)

    DTI_Ind = find(GradTable(:,4)>BValLowTHD & GradTable(:,4)<BValHighTHD);
    GradVec = GradTable(DTI_Ind,1:3);

    params.MDirection = GradVec;
    params.maxOrder = SPHMaxOrder;
    params.B = matrixB(params.MDirection,params.maxOrder);
    res = params.B*coefficient;
end

