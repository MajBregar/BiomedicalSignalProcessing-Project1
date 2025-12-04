
function [idx] = RealTime_QRS_Detector(fileName, sampling_rate, plots_enabled)

    S = load(fileName);
    original_signal = S.val(1, :);

    MLII_raw = resample_signal(original_signal, 200, sampling_rate);
    MLII_no_DC = MLII_raw - mean(MLII_raw);
    original_sr = sampling_rate;
    sampling_rate = 200;
    delay_compensation = 0; 

    % Bandpass filtering (paper)
    lowpass = lowpass(MLII_no_DC);
    delay_compensation += 6;
    bandpass = highpass(lowpass);
    delay_compensation += 16;

    % Derivative (paper)
    d = derivative(bandpass);
    delay_compensation += 2;


    % Squaring (paper)
    sq = d .^ 2;

    % Moving-window integration (paper)
    N = round(0.150 * sampling_rate);
    mwi = integrator(sq, N);
    
    r_peaks = peak_detect_PT(bandpass, mwi, sampling_rate);
    idx_corrected = r_peaks - 16;
    idx_360 = round( idx_corrected * (360 / 200) );


    idx = idx_360;

    if plots_enabled
        start_ms = 0;
        end_ms = 20000;

        % *** RAW SIGNAL ***
        plot_signal_time_domain(original_signal, original_sr, start_ms, end_ms, "0_raw.png");
        plot_signal_time_domain(MLII_no_DC, sampling_rate, start_ms, end_ms, "1_no_DC_resample.png");


        % *** BANDPASS OUTPUT ***
        plot_signal_time_domain(lowpass, sampling_rate, start_ms, end_ms, "2a_lowpass.png");
        plot_signal_time_domain(bandpass, sampling_rate, start_ms, end_ms, "2b_bandpass.png");


        % *** DERIVATIVE OUTPUT ***
        plot_signal_time_domain(d, sampling_rate, start_ms, end_ms, "3_derivative.png");

        % *** SQUARING OUTPUT ***
        plot_signal_time_domain(sq, sampling_rate, start_ms, end_ms, "4_squared.png");

        % *** MOVING WINDOW INTEGRATION OUTPUT ***
        plot_signal_time_domain(mwi, sampling_rate, start_ms, end_ms, "5_mwa.png");

        % *** RAW ECG WITH FINAL DETECTED QRS PEAKS ***
        plot_r_peaks_on_raw(MLII_raw, idx_corrected, sampling_rate, start_ms, end_ms, "6_final_peaks.png");
        plot_r_peaks_on_raw(original_signal, idx_360, original_sr, start_ms, end_ms, "7_test.png");
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

function QRS_locs = peak_detect_PT(filt, mwa, fs)

    refractory = round(0.200 * fs);
    Twave_limit = round(0.360 * fs);

    % peak buffers
    last_QRS = -inf;
    QRS_locs = [];

    % running estimates (paper eq. 12–21)
    SPKF = 0; NPKF = 0;
    SPKI = 0; NPKI = 0;

    THF1 = 0; THF2 = 0;
    THI1 = 0; THI2 = 0;

    RR_intervals = [];
    RR_missed = inf;

    % ============================================================
    % MAIN LOOP
    % ============================================================
    for n = 2:length(mwa)-1

        % --------------------------------------------------------
        % 1) locate peaks in EACH signal separately
        % --------------------------------------------------------
        is_peak_f = (filt(n-1) < filt(n) && filt(n) > filt(n+1));
        is_peak_i = (mwa(n-1)  < mwa(n)  && mwa(n)  > mwa(n+1));

        if !is_peak_f && !is_peak_i
            continue;
        end

        % note: peak candidate = current sample n
        peakF = filt(n);
        peakI = mwa(n);

        % --------------------------------------------------------
        % 2) classify as signal or noise based on thresholds
        % --------------------------------------------------------
        is_QRS_candidate = ...
            peakF > THF1 && peakI > THI1 && ...
            (n - last_QRS) > refractory;

        if is_QRS_candidate

            % ----------------------------------------------------
            % 3) T-wave discrimination (paper)
            % ----------------------------------------------------
            if (n - last_QRS) < Twave_limit && last_QRS > 5 && n < length(filt)-5
                prev_slope = max(abs(diff(filt(last_QRS-5:last_QRS+5))));
                new_slope  = max(abs(diff(filt(n-5:n+5))));
                if new_slope < 0.5 * prev_slope
                    % classify as noise
                    NPKI = 0.125*peakI + 0.875*NPKI;
                    NPKF = 0.125*peakF + 0.875*NPKF;
                    [THF1,THF2,THI1,THI2] = update_thresholds(SPKF,NPKF,SPKI,NPKI);
                    continue;
                end
            end

            % ----------------------------------------------------
            % 4) VALID QRS DETECTED
            % ----------------------------------------------------
            QRS_locs(end+1) = n;
            RR_intervals(end+1) = n - last_QRS;
            last_QRS = n;

            % update signal peak estimates
            SPKF = 0.125*peakF + 0.875*SPKF;
            SPKI = 0.125*peakI + 0.875*SPKI;

            % update RR averages
            if length(RR_intervals) > 8
                RR_last = RR_intervals(end-7:end);
                RR_avg = mean(RR_last);
                RR_missed = round(1.66 * RR_avg);
            end

        else
            % ----------------------------------------------------
            % noise peak
            % ----------------------------------------------------
            NPKF = 0.125*peakF + 0.875*NPKF;
            NPKI = 0.125*peakI + 0.875*NPKI;
        end

        % --------------------------------------------------------
        % update thresholds (paper equations)
        % --------------------------------------------------------
        [THF1,THF2,THI1,THI2] = update_thresholds(SPKF,NPKF,SPKI,NPKI);

        % --------------------------------------------------------
        % 5) SEARCH-BACK (paper)
        % --------------------------------------------------------
        if (n - last_QRS) > RR_missed

            % search max peak in interval (paper: no buffer needed here)
            window_start = last_QRS + refractory;
            window_end   = n - 1;

            if window_end > window_start
                [max_val, max_idx_rel] = max(mwa(window_start:window_end));
                max_idx = window_start + max_idx_rel - 1;

                if max_val > THI2
                    QRS_locs(end+1) = max_idx;
                    last_QRS = max_idx;

                    % update signal peak (paper: faster update)
                    SPKI = 0.25*mwa(max_idx)  + 0.75*SPKI;
                    SPKF = 0.25*filt(max_idx) + 0.75*SPKF;

                    % refresh thresholds
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
