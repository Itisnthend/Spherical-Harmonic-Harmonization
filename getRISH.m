function C = getRISH(coefficient, SPHMaxOrder)


coefficient = coefficient.^2;
orderMatrix = [0,2,4,6,8];
indx = 1;
res = zeros(SPHMaxOrder/2+1,1);
for order = 1:SPHMaxOrder/2+1
    for count = 1:orderMatrix(order)*2+1
        res(order,1) = res(order,1) + coefficient(indx);
        indx = indx+1;
    end
end
   
C = res;
end

