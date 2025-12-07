function [idx] = QRS_Detector_IMPROVED(fileName, sampling_rate, plots_enabled)
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

    N = 20;
    integrated = integrator(squared_energy, N);
    
    % 4) ADAPTIVE MIN DIST CLASSIFIER
    W = 15;
    [peak_values, peak_indices, classes, mu_QRS_hist, mu_nonQRS_hist] = min_dist_classifier(integrated, W, sampling_rate);

    % 5) MINMAX SEARCHER
    QRS_peaks = peak_indices(classes == 1);
    [R_locs, S_locs] = minimax_searcher(QRS_peaks, differentiated, MLII_raw);

    % DEBUG
    if plots_enabled
        plot_ms_start = 0000;
        plot_ms_end = 1000;
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


function y = differentiator(x)
    h = [-1, 8, 0, -8, 1] / 12;
    y = conv(x, h, 'same');
end

function y = multiplicator(x)
    y = x .^ 2;
end


function y = integrator(x, N)
    h = ones(1, N);
    y = conv(x, h, 'same');
end


function [peaks, values] = global_maxima(yb, W)
    half = floor(W/2);
    N = length(yb);

    peaks = [];
    values = [];

    for k = half+1 : N-half
        win = yb(k-half : k+half);
        
        if yb(k) == max(win) && sum(win == yb(k)) == 1
            peaks(end+1) = k;
            values(end+1) = yb(k);
        end
    end

end


function [peak_values, peak_indices, classes, SPK_hist, NPK_hist, TH1_hist] = ...
         min_dist_classifier(yb, W, fs)

    refractory_samples = floor(0.2 * fs);

    [peak_indices, peak_values] = global_maxima(yb, W);

    init_seconds = 2;
    init_seg = yb(1 : min(length(yb), init_seconds*fs));

    SPK = max(init_seg);
    NPK = mean(init_seg(init_seg < mean(init_seg)));
    if isnan(NPK), NPK = 0; end
    TH1 = NPK + 0.25 * (SPK - NPK);

    SPK_hist = SPK;
    NPK_hist = NPK;
    TH1_hist = TH1;

    classes = zeros(1, length(peak_indices));
    last_QRS = -inf;

    for k = 1:length(peak_indices)

        ind = peak_indices(k);
        pk  = peak_values(k);

        is_QRS_candidate = pk > TH1 && (ind - last_QRS) > refractory_samples;

        if is_QRS_candidate
            classes(k) = 1;
            last_QRS = ind;
            SPK = 0.125 * pk + 0.875 * SPK;
        else
            classes(k) = 0;
            NPK = 0.125 * pk + 0.875 * NPK;
        end

        TH1 = NPK + 0.25 * (SPK - NPK);

        SPK_hist(end+1) = SPK;
        NPK_hist(end+1) = NPK;
        TH1_hist(end+1) = TH1;
    end
end






function [R_locs, S_locs] = minimax_searcher(peak_indices, derivative, raw)
    R_locs = [];
    S_locs = [];

    N = length(raw);

    for i = 1:length(peak_indices)
        p = peak_indices(i);

        % nearest left zero crossing
        k = p;
        while k > 2 && sign(derivative(k)) == sign(derivative(k-1))
            k = k - 1;
        end
        left_zc = k;

        % nearest right zero crossing
        k = p;
        while k < N-1 && sign(derivative(k)) == sign(derivative(k+1))
            k = k + 1;
        end
        right_zc = k;

        % reject faulty pairs
        if right_zc <= left_zc
            continue;
        end

        % minmax
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
