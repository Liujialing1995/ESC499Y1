function y = evaluatePopulationFitness(P, MASK_LOWER, MASK_UPPER, PHASE, BASIS)
% Performs a fitness evaluation on each member of the entire population
    size = length(P(1,:));
    y = zeros(size, 1);
    for i = 1:size
        y(i) = evaluateFitness(P(:,i),MASK_LOWER,MASK_UPPER,PHASE,BASIS);
    end
end

