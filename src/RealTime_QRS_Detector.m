
function [idx] = RealTime_QRS_Detector(file, sampling_rate, plots_enabled)
    S = load(file);
    original_signal = S.val(1, :);
    
    % RESAMPLE AND REMOVE DC
    MLII_raw = resample_signal(original_signal, 200, sampling_rate);
    MLII_no_DC = MLII_raw - mean(MLII_raw);
    original_sr = sampling_rate;
    sampling_rate = 200;
    delay_compensation = 0; 

    % BANDPASS
    lowpass = lowpass(MLII_no_DC);
    delay_compensation += 6;
    bandpass = highpass(lowpass);
    delay_compensation += 16;

    % DERIVATIVE
    d = derivative(bandpass);
    delay_compensation += 2;

    % SQUARE
    sq = d .^ 2;

    % INTEGRATION
    qrs_time_window_ms = 150;
    N = round(qrs_time_window_ms / 1000 * sampling_rate);
    mwi = integrator(sq, N);
    
    % QRS DETECTION
    RP_ms = 200;
    Tw_ms = 360;
    r_peaks = QRS_PT(bandpass, mwi, sampling_rate, RP_ms, Tw_ms);
    delay_compensation += 2;

    % QRS ALIGNMENT
    r_peaks_time_adjusted = r_peaks - delay_compensation - 2;
    r_peaks_old_frequency = round( r_peaks_time_adjusted * (original_sr / sampling_rate) );
    

    idx = r_peaks_old_frequency;

    if plots_enabled
        start_ms = 0;
        end_ms = 7000;
        plot_signal_time_domain(original_signal, original_sr, start_ms, end_ms, "0_raw.png");
        plot_signal_time_domain(MLII_no_DC, sampling_rate, start_ms, end_ms, "1_no_DC_resample.png");
        plot_signal_time_domain(lowpass, sampling_rate, start_ms, end_ms, "2a_lowpass.png");
        plot_signal_time_domain(bandpass, sampling_rate, start_ms, end_ms, "2b_bandpass.png");
        plot_signal_time_domain(d, sampling_rate, start_ms, end_ms, "3_derivative.png");
        plot_signal_time_domain(sq, sampling_rate, start_ms, end_ms, "4_squared.png");
        plot_signal_time_domain(mwi, sampling_rate, start_ms, end_ms, "5_mwa.png");
        plot_r_peaks_on_raw(MLII_raw, r_peaks_time_adjusted, sampling_rate, start_ms, end_ms, "6_peaks.png");
        plot_r_peaks_on_raw(original_signal, r_peaks_old_frequency, original_sr, start_ms, end_ms, "7_final.png");
    end

end



function y = resample_signal(x, new_fs, old_fs)
    [p, q] = rat(new_fs / old_fs, 1e-12);
    x = x(:)';
    fc = 0.9 * min(1/p, 1/q) / 2;
    L = 10 * max(p, q);
    n = -L : L;
    h = 2 * fc * sinc(2 * fc * n);
    w = blackman(length(h))';
    h = h .* w;
    h = p * h / sum(h);
    xu = zeros(1, length(x) * p);
    xu(1:p:end) = x;
    xf = conv(xu, h, "same");
    y = xf(1:q:end);
end


function y = lowpass(x)
    n = length(x);
    y = zeros(1,n);

    for i = 13:n
        y(i) = (2*y(i-1) - y(i-2) + x(i) - 2*x(i-6) + x(i-12)) / 36;
    end
end

function y = highpass(x)
    n = length(x);
    y = zeros(1,n);

    for i = 33:n
        y(i) = (32*x(i-16) - ( y(i-1) + x(i) - x(i-32))) / 32;
    end
end

function y = derivative(x)
    n = length(x);
    y = zeros(1,n);
    for i = 3:(n-2)
        y(i) = (-x(i-2) - 2*x(i-1) + 2*x(i+1) + x(i+2)) / 8;
    end
end

function y = integrator(x, N)
    n = length(x);
    y = zeros(1,n);

    buffer = zeros(1, N);
    buf_sum = 0;
    idx = 1;

    for i = 1:n
        buf_sum = buf_sum - buffer(idx);
        buffer(idx) = x(i);
        buf_sum = buf_sum + x(i);

        idx = idx + 1;
        if idx > N
            idx = 1;
        end

        if i < N
            y(i) = buf_sum / i;
        else
            y(i) = buf_sum / N;
        end
    end
end




function [RR_latest, RR_selected, RR_average1, RR_average2, RR_missed] = rr_update(RR_latest, RR_selected, new_RR)
    RR_average1 = [];
    RR_average2 = [];
    RR_missed   = [];

    % update RR_latest
    if ~isempty(new_RR)
        RR_latest(end+1) = new_RR;
        if length(RR_latest) > 8
            RR_latest = RR_latest(end-7:end);
        end
    end

    % calculate mean of RR_latest
    if length(RR_latest) == 8
        RR_average1 = mean(RR_latest);
    end


    if isempty(RR_selected)
        if length(RR_latest) >= 1
            RR_selected(end+1) = new_RR;
            if length(RR_selected) > 8
                RR_selected = RR_selected(end-7:end);
            end
        end
    else
        % USING AVERAGE 2 TO INSERT BASED OFF LIMITS
        if length(RR_selected) == 8
            RR_average2 = mean(RR_selected);
            RR_LOW  = 0.92 * RR_average2;
            RR_HIGH = 1.16 * RR_average2;
        else
            RR_average2 = [];
        end

        if ~isempty(RR_average2)
            if new_RR > RR_LOW && new_RR < RR_HIGH
                RR_selected(end+1) = new_RR;
                if length(RR_selected) > 8
                    RR_selected = RR_selected(end-7:end);
                end
            end
        else
            RR_selected(end+1) = new_RR;
            if length(RR_selected) > 8
                RR_selected = RR_selected(end-7:end);
            end
        end
    end


    if length(RR_selected) == 8
        RR_average2 = mean(RR_selected);
    end


    if ~isempty(RR_average2)
        RR_missed = round(1.66 * RR_average2);
    end
end



function [THR_F1, THR_F2, THR_I1, THR_I2] = update_thresholds(SPKF, NPKF, SPKI, NPKI)
    THR_F1 = NPKF + 0.25*(SPKF - NPKF);
    THR_F2 = 0.5 * THR_F1;
    THR_I1 = NPKI + 0.25*(SPKI - NPKI);
    THR_I2 = 0.5 * THR_I1;
end

function QRS_locations = QRS_PT(filtered, integrated, sampling_rate, RP_ms, Tw_ms)
    % constants
    refractory_samples = round(RP_ms / 1000 * sampling_rate);
    T_wave_period_samples = round(Tw_ms / 1000 * sampling_rate);
    T_wave_slope_samples = round(75 / 1000 * sampling_rate);

    % global algorithm variables
    QRS_locations = [];
    previous_QRS_sample_ind = -1;

    SPKI = 0; NPKI = 0;
    THR_I1 = 0; THR_I2 = 0;

    SPKF = 0; NPKF = 0;
    THR_F1 = 0; THR_F2 = 0;

    RR_latest = [];
    RR_selected = [];

    RR_average1 = []; 
    RR_average2 = []; 
    RR_missed = [];

    % implementation specific variables
    rising_f = false;
    rising_i = false;
    candidate_peak_f = 0;
    candidate_peak_i = 0;

    for n = 3:length(filtered)
        
        % 1) PEAK DETECTION
        if filtered(n) > filtered(n-1)
            rising_f = true;
            candidate_peak_f = n;
            is_peak_f = false;
        else
            is_peak_f = rising_f && candidate_peak_f == n-1;
            rising_f = false;
        end

        if integrated(n) > integrated(n-1)
            rising_i = true;
            candidate_peak_i = n;
            is_peak_i = false;
        else
            is_peak_i = rising_i && candidate_peak_i == n-1;
            rising_i = false;
        end

        % NO PEAK
        if ~is_peak_f && ~is_peak_i
            continue;
        end

        % PEAK ALWAYS AT SAME POS IF DETECTED
        peak_sample_index = n-1;
        pF = filtered(peak_sample_index);
        pI = integrated(peak_sample_index);


        % 2) SIGNAL CLASSIFICATION AND UPDATING ESTIMATES
        is_signal_F = false;
        is_signal_I = false;

        if pF > THR_F1
            is_signal_F = true;
            SPKF = 0.125*pF + 0.875*SPKF;
        else
            NPKF = 0.125*pF + 0.875*NPKF;
        end

        if pI > THR_I1
            is_signal_I = true;
            SPKI = 0.125*pI + 0.875*SPKI;
        else
            NPKI = 0.125*pI + 0.875*NPKI;
        end


        % 3) VALIDATE QRS CANDIDATE
        RR = peak_sample_index - previous_QRS_sample_ind;
        is_valid_QRS_candidate = is_signal_F && is_signal_I && ((RR >= refractory_samples) || previous_QRS_sample_ind < 0);

        if is_valid_QRS_candidate

            % 3a) T-WAVE DISCRIMINATION
            possible_t_wave = RR < T_wave_period_samples && previous_QRS_sample_ind > 0;

            if possible_t_wave
                prev_QRS = previous_QRS_sample_ind;
                prev_start = max(prev_QRS - T_wave_slope_samples, 1);
                prev_end   = prev_QRS;

                curr_peak = peak_sample_index;
                curr_start = max(curr_peak - T_wave_slope_samples, 1);
                curr_end   = curr_peak;

                prev_slope = max(abs(diff(filtered(prev_start:prev_end))));
                new_slope  = max(abs(diff(filtered(curr_start:curr_end))));

                if new_slope < 0.5 * prev_slope
                    % T-WAVE IS NOISE
                    NPKI = 0.125*pI + 0.875*NPKI;
                    NPKF = 0.125*pF + 0.875*NPKF;
                    [THR_F1, THR_F2, THR_I1, THR_I2] = update_thresholds(SPKF, NPKF, SPKI, NPKI);
                    continue;
                end
            end

            % 3b) SURVIVING CANDIDATE IS A QRS PEAK - update RR intervals and averages

            if previous_QRS_sample_ind > 0
                new_RR = peak_sample_index - previous_QRS_sample_ind;
                [RR_latest, RR_selected, RR_average1, RR_average2, RR_missed] = rr_update(RR_latest, RR_selected, new_RR);
            end

            QRS_locations(end+1) = peak_sample_index;
            previous_QRS_sample_ind = peak_sample_index;

        end

        % 4) UPDATE THRESHOLDS
        [THR_F1, THR_F2, THR_I1, THR_I2] = update_thresholds(SPKF, NPKF, SPKI, NPKI);


        % 5) SEARCHBACK
        if (n - previous_QRS_sample_ind) > RR_missed

            window_start = previous_QRS_sample_ind + refractory_samples;
            window_end   = n - 1;

            if window_end > window_start
                [max_val, rel_idx] = max(integrated(window_start:window_end));
                max_idx = window_start + rel_idx - 1;

                if max_val > THR_I2
                    QRS_locations(end+1) = max_idx;
                    previous_QRS_sample_ind = max_idx;

                    SPKI = 0.25 * pI + 0.75 * SPKI;
                    SPKF = 0.25 * pF + 0.75 * SPKF;

                    [THR_F1, THR_F2, THR_I1, THR_I2] = update_thresholds(SPKF, NPKF, SPKI, NPKI);
                end
            end
        end

    end
end





function plot_signal_time_domain(signal, fs, ms_start, ms_end, save_path)
    start_idx = max(1, round((ms_start / 1000) * fs) + 1);
    end_idx   = min(length(signal), round((ms_end / 1000) * fs));
    y = signal(start_idx:end_idx);
    t = (0:length(y)-1) / fs;
    figure('visible', 'off'); plot(t, y, 'r'); print(save_path, '-dpng'); close all;
end

function plot_r_peaks_on_raw(raw, R_locs, fs, ms_start, ms_end, save_path)
    start_idx = max(1, round((ms_start / 1000) * fs) + 1);
    end_idx   = min(length(raw), round((ms_end / 1000) * fs));
    sig = raw(start_idx:end_idx); t = (0:length(sig)-1) / fs;
    valid = R_locs(R_locs>=start_idx & R_locs<=end_idx);
    R_times = (valid - start_idx) / fs;
    figure('visible','off'); plot(t,sig,'b'); hold on;
    plot(R_times, raw(valid),'rv','MarkerFaceColor','r'); print(save_path,'-dpng'); close all;
end
