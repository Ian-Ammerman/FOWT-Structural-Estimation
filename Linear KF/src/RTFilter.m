function [a_filtered, filter_state] = RTFilter(a_raw, filter_coeffs, filter_state)
    % a_raw:       Current raw acceleration measurement
    % band_range:  [low_cutoff, high_cutoff] in Hz
    % dt:          Time step (s)
    % filter_coeffs: Structure containing filter coefficients
    % filter_state: Previous state of the filter

    % Apply the filter to the current sample
    [a_filtered, filter_state] = filter(filter_coeffs.b, filter_coeffs.a, a_raw, filter_state);
end
