function total_error = multibeam_error_sumsqr_points_outside_mask(MASK_LOWER, MASK_UPPER, X_, PHASE, BASIS)
% Calculates the error as the sum of squares of difference of squares of
% array factor points outside of either mask
%     Array_factor = Series_coefficients.' * Basis ./ length(Series_coefficients);
%     total_error = (sumsqr(abs(Array_factor(abs(Array_factor) < Mask_lower)).^2 - abs(Mask_lower(abs(Array_factor) < Mask_lower)).^2) + ...
%             sumsqr(abs(Array_factor(abs(Array_factor) > Mask_upper)).^2 - abs(Mask_upper(abs(Array_factor) > Mask_upper)).^2));

%     AF = abs(X_.' * PHASE_BASIS);
    AF = ones(length(X_),1)*X_.*PHASE*BASIS;
    UPPER_EXCEED = AF > MASK_UPPER;
    LOWER_EXCEED = AF < MASK_LOWER;
    total_error = sumsqr(AF(UPPER_EXCEED) - MASK_UPPER(UPPER_EXCEED)) + ...
                  sumsqr(MASK_LOWER(LOWER_EXCEED) - AF(LOWER_EXCEED));
end

