function res = getMappedCoef(coefficient, pie, SPHMaxOrder)
    
    [~, numSubjects] = size(coefficient);
    orderMatrix = [0,2,4,6,8];
    indx = 1;
    res = zeros((SPHMaxOrder+1)*(SPHMaxOrder+2)/2, numSubjects);
    for order = 1:SPHMaxOrder/2+1
        for count = 1:orderMatrix(order)*2+1
            for num = 1: numSubjects
                res(indx,num) = coefficient(indx,num)*pie(order,num);
            end
            indx = indx+1;
        end
    end
end

