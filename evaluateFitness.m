function fitness = evaluateFitness(chromosome, MASK_LOWER, MASK_UPPER, PHASE, BASIS)
% Evaluates the fitness for a single chromosome given the mask definition
% of error
    coefficient = exp(1j*arrayfun(@decode_gene,chromosome));
    % Coefficients need to be passed in as a row vector
    error =  multibeam_error_sumsqr_points_outside_mask(MASK_LOWER, MASK_UPPER, coefficient.', PHASE, BASIS);
    fitness = 1/(1 + sqrt(error));
end

