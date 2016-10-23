function total_error = multibeam_error_sumsqr_points_outside_mask(Series_coefficients, Basis, Mask_upper, Mask_lower)
% Calculates the error as the sum of squares of difference of squares of
% array factor points outside of either mask
    Array_factor = Series_coefficients.' * Basis;
    total_error = (sumsqr(abs(Array_factor(abs(Array_factor) < Mask_lower)).^2 - abs(Mask_lower(abs(Array_factor) < Mask_lower)).^2) + ...
            sumsqr(abs(Array_factor(abs(Array_factor) > Mask_upper)).^2 - abs(Mask_upper(abs(Array_factor) > Mask_upper)).^2));
end

