
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
    
    % PEAK DETECTION
    RP_ms = 200;
    Tw_ms = 360;
    r_peaks = peak_detection(bandpass, mwi, sampling_rate, RP_ms, Tw_ms);
    r_peaks_time_adjusted = r_peaks - 16;
    r_peaks_old_frequency = round( r_peaks_time_adjusted * (original_sr / sampling_rate) );

    idx = r_peaks_old_frequency;

    if plots_enabled
        start_ms = 0;
        end_ms = 5000;
        plot_signal_time_domain(original_signal, original_sr, start_ms, end_ms, "0_raw.png");
        plot_signal_time_domain(MLII_no_DC, sampling_rate, start_ms, end_ms, "1_no_DC_resample.png");
        plot_signal_time_domain(lowpass, sampling_rate, start_ms, end_ms, "2a_lowpass.png");
        plot_signal_time_domain(bandpass, sampling_rate, start_ms, end_ms, "2b_bandpass.png");
        plot_signal_time_domain(d, sampling_rate, start_ms, end_ms, "3_derivative.png");
        plot_signal_time_domain(sq, sampling_rate, start_ms, end_ms, "4_squared.png");
        plot_signal_time_domain(mwi, sampling_rate, start_ms, end_ms, "5_mwa.png");
        plot_r_peaks_on_raw(MLII_raw, r_peaks_time_adjusted, sampling_rate, start_ms, end_ms, "6_final_peaks.png");
        plot_r_peaks_on_raw(original_signal, r_peaks_old_frequency, original_sr, start_ms, end_ms, "7_test.png");
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
    c = cumsum(x);

    for i = 1:n
        if i < N
            y(i) = c(i) / i;
        else
            y(i) = (c(i) - c(i-N+1)) / N;
        end
    end
end






function [THF1,THF2,THI1,THI2] = update_thresholds(SPKF,NPKF,SPKI,NPKI)
    THF1 = NPKF + 0.25*(SPKF - NPKF);
    THF2 = 0.5 * THF1;
    THI1 = NPKI + 0.25*(SPKI - NPKI);
    THI2 = 0.5 * THI1;
end

function QRS_locs = peak_detection(filt, mwa, fs, RP_ms, Tw_ms)

    refractory_samples = round(RP_ms / 1000 * fs);
    T_wave_period_samples = round(Tw_ms / 1000 * fs);

    last_QRS = -inf;
    QRS_locs = [];

    SPKF = 0; NPKF = 0;
    SPKI = 0; NPKI = 0;

    THF1 = 0; THF2 = 0;
    THI1 = 0; THI2 = 0;

    RR_intervals = [];
    RR_missed = inf;


    rising_f = false;
    rising_i = false;

    candidate_peak_f = 0;
    candidate_peak_i = 0;

    for n = 3:length(mwa)

        % 1) FIDUCIAL MARK DETECTION

        if filt(n) > filt(n-1)
            rising_f = true;
            candidate_peak_f = n;
            is_peak_f = false;
        else
            is_peak_f = rising_f && candidate_peak_f == n-1;
            rising_f = false;
        end

        if mwa(n) > mwa(n-1)
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
        peak_index = n-1;
        peakF = filt(peak_index);
        peakI = mwa(peak_index);



        is_QRS_candidate = peakF > THF1 && peakI > THI1 && (peak_index - last_QRS) > refractory;

        if is_QRS_candidate

            % ======================================
            % 3. T-wave suppression (causal version)
            % ======================================
            dt = peak_index - last_QRS;

            if dt < Twave_limit && last_QRS > 10
                % slope of last QRS
                window_prev_start = max(last_QRS-10,1);
                window_prev_end   = last_QRS;

                prev_slope = max(abs(diff( filt(window_prev_start:window_prev_end) )));

                % slope of current region — ONLY past portion
                window_new_start = max(peak_index-10,1);
                window_new_end   = peak_index;

                new_slope = max(abs(diff( filt(window_new_start:window_new_end) )));

                % slope rule
                if new_slope < 0.5 * prev_slope
                    NPKI = 0.125*peakI + 0.875*NPKI;
                    NPKF = 0.125*peakF + 0.875*NPKF;
                    [THF1,THF2,THI1,THI2] = update_thresholds(SPKF,NPKF,SPKI,NPKI);
                    continue;
                end
            end

            % ======================================
            % 4. Valid QRS detected
            % ======================================
            QRS_locs(end+1) = peak_index;
            RR_intervals(end+1) = peak_index - last_QRS;
            last_QRS = peak_index;

            SPKF = 0.125*peakF + 0.875*SPKF;
            SPKI = 0.125*peakI + 0.875*SPKI;

            % RR average for search-back
            if length(RR_intervals) > 8
                RR_last = RR_intervals(end-7:end);
                RR_avg = mean(RR_last);
                RR_missed = round(1.66 * RR_avg);
            end

        else
            % noise update
            NPKF = 0.125*peakF + 0.875*NPKF;
            NPKI = 0.125*peakI + 0.875*NPKI;
        end

        % update thresholds
        [THF1,THF2,THI1,THI2] = update_thresholds(SPKF,NPKF,SPKI,NPKI);

        % ================================
        % 5. Search-Back (fully causal)
        % ================================
        if (n - last_QRS) > RR_missed

            window_start = last_QRS + refractory;
            window_end   = n - 1;

            if window_end > window_start
                [max_val, rel_idx] = max(mwa(window_start:window_end));
                max_idx = window_start + rel_idx - 1;

                if max_val > THI2
                    QRS_locs(end+1) = max_idx;
                    last_QRS = max_idx;

                    SPKI = 0.25*mwa(max_idx)  + 0.75*SPKI;
                    SPKF = 0.25*filt(max_idx) + 0.75*SPKF;

                    [THF1,THF2,THI1,THI2] = update_thresholds(SPKF,NPKF,SPKI,NPKI);
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
