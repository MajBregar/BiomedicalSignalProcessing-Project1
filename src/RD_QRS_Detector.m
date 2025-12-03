function [idx] = RD_QRS_Detector(fileName, sampling_rate, plots_enabled)
    S = load(fileName);
    MLII_raw = S.val(1, :);
    
    % 1) NOISE FILTER
    K = 5;
    L = 200;
    filtered = noise_filter(MLII_raw, K, L);

    % 2) DIFFERENTIATOR
    differentiated = differentiator(filtered);

    % 3) ENERGY COLLECTOR
    squared_energy = multiplicator(differentiated);

    N = round(0.03 * sampling_rate);
    integrated = integrator(squared_energy, N);

    % 4) ADAPTIVE MIN DIST CLASSIFIER
    W = 15;
    [peak_values, peak_indices, classes, mu_QRS_hist, mu_nonQRS_hist] = min_dist_classifier(integrated, W, sampling_rate);

    % 5) MINMAX SEARCHER
    QRS_peaks = peak_indices(classes == 1);
    [R_locs, S_locs] = minimax_searcher(QRS_peaks, differentiated, MLII_raw);

    % DEBUG
    if plots_enabled
        plot_ms_start = 0;
        plot_ms_end = 300;
        plot_signal_time_domain(MLII_raw,       sampling_rate, plot_ms_start, plot_ms_end, '0_MLII_raw.png');
        plot_signal_time_domain(filtered,       sampling_rate, plot_ms_start, plot_ms_end, '1_filtered_signal.png');
        plot_signal_time_domain(differentiated, sampling_rate, plot_ms_start, plot_ms_end, '2_differentiator.png');
        plot_signal_time_domain(squared_energy, sampling_rate, plot_ms_start, plot_ms_end, '3_squared.png');
        plot_signal_time_domain(integrated,     sampling_rate, plot_ms_start, plot_ms_end, '4_integrated.png');
        plot_signal_time_domain(peak_values,    sampling_rate, plot_ms_start, plot_ms_end, '5_peak_value.png');
        plot_signal_time_domain(classes,        sampling_rate, plot_ms_start, plot_ms_end, '6_classes.png');
        plot_signal_time_domain(mu_QRS_hist,    sampling_rate, plot_ms_start, plot_ms_end, '7_mu_QRS_hist.png');
        plot_signal_time_domain(mu_nonQRS_hist, sampling_rate, plot_ms_start, plot_ms_end, '8_mu_nonQRS_hist.png');
        plot_r_peaks_on_raw(MLII_raw, R_locs,   sampling_rate, plot_ms_start, plot_ms_end, '9_detected.png');
    end

    % OUTPUT
    idx = R_locs;
end



function y = noise_filter(x, K, L)
    hK = ones(1, K) / K;
    hL = ones(1, L) / L;

    yK = conv(conv(x, hK, 'same'), hK, 'same');
    yL = conv(conv(x, hL, 'same'), hL, 'same');

    y = yK - yL;
end

    pri

function y = differentiator(x)
    h = [-1, -2, 0, 2, 1] / 3;
    y = conv(x, h, 'same');
end

function y = multiplicator(x)
    y = x .^ 2;
end


function y = integrator(x, N)
    h = ones(1, N);
    y = conv(x, h, 'same');
end


function [peak_values, peak_indices, classes, mu_QRS_hist, mu_nonQRS_hist] = min_dist_classifier(yb, W, fs)

    half = floor(W/2);
    peak_values  = [];
    peak_indices = [];
    classes      = [];

    init_seconds = 2;
    init_seg = yb(1:min(length(yb), init_seconds*fs));
    mu_QRS     = max(init_seg);
    mu_nonQRS  = min(init_seg);

    mu_QRS_hist    = mu_QRS;
    mu_nonQRS_hist = mu_nonQRS;

    k = 1;
    N = length(yb);

    while k <= N-W+1

        % extract W-sample window
        win = yb(k : k+W-1);

        % find its single global max
        [pk, idx] = max(win);
        peak_idx = k + idx - 1;

        % append only once per window
        peak_values(end+1)  = pk;
        peak_indices(end+1) = peak_idx;

        % ---- classification ----
        d_QRS    = abs(pk - mu_QRS);
        d_nonQRS = abs(pk - mu_nonQRS);

        if d_QRS < d_nonQRS
            classes(end+1) = 1;
            mu_QRS = 0.9*mu_QRS + 0.1*pk;
        else
            classes(end+1) = 0;
            mu_nonQRS = 0.9*mu_nonQRS + 0.1*pk;
        end

        mu_QRS_hist(end+1)    = mu_QRS;
        mu_nonQRS_hist(end+1) = mu_nonQRS;

        % important: advance by W, not by 1
        k = k + W;
    end
end



function [R_locs, S_locs] = minimax_searcher(peak_indices, derivative, raw)
    R_locs = [];
    S_locs = [];

    N = length(raw);

    for i = 1:length(peak_indices)
        p = peak_indices(i);

        % -----------------------------
        % 1. FIND LEFT ZERO CROSSING
        % -----------------------------
        k = p;
        while k > 2 && sign(derivative(k)) == sign(derivative(k-1))
            k = k - 1;
        end
        left_zc = k;

        % -----------------------------
        % 2. FIND RIGHT ZERO CROSSING
        % -----------------------------
        k = p;
        while k < N-1 && sign(derivative(k)) == sign(derivative(k+1))
            k = k + 1;
        end
        right_zc = k;

        if right_zc <= left_zc
            continue;
        end

        % -----------------------------
        % 3. MINIMAX SEARCH
        % -----------------------------
        segment = raw(left_zc:right_zc);

        [~, local_max] = max(segment);
        [~, local_min] = min(segment);

        R_locs(end+1) = left_zc + local_max - 1;
        S_locs(end+1) = left_zc + local_min - 1;
    end
end








% VISUALIZATION HELPERS

function plot_signal_time_domain(signal, fs, ms_start, ms_end, save_path)
    start_idx = max(1, round((ms_start / 1000) * fs) + 1);
    end_idx   = min(length(signal), round((ms_end / 1000) * fs));

    y = signal(start_idx:end_idx);
    t = (0:length(y)-1) / fs;

    figure('visible', 'off');
    plot(t, y, 'r', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Signal Slice (%d ms → %d ms)', ms_start, ms_end));
    grid on;

    print(save_path, '-dpng');
    close all;
end


function plot_r_peaks_on_raw(raw, R_locs, fs, ms_start, ms_end, save_path)
    start_idx = max(1, round((ms_start / 1000) * fs) + 1);
    end_idx   = min(length(raw), round((ms_end / 1000) * fs));

    sig = raw(start_idx:end_idx);

    t = (0:length(sig)-1) / fs;

    valid_mask = (R_locs >= start_idx) & (R_locs <= end_idx);
    valid_R = R_locs(valid_mask);

    R_times = (valid_R - start_idx) / fs;

    figure('visible','off');
    plot(t, sig, 'b', 'LineWidth', 1); hold on;

    plot(R_times, raw(valid_R), 'rv', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'r');

    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    title(sprintf('ECG Slice %d–%d ms with %d QRS markers', ...
                  ms_start, ms_end, length(valid_R)));
    grid on;


    print(save_path, '-dpng');
    close all;
end
