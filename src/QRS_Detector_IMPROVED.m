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
    QRS_peaks= min_dist_classifier(integrated, W, sampling_rate);

    % DEBUG
    if plots_enabled
        plot_ms_start = 0000;
        plot_ms_end = 60000;
        plot_signal_time_domain(MLII_raw,       sampling_rate, plot_ms_start, plot_ms_end, '0_MLII_raw.png');
        plot_signal_time_domain(filtered,       sampling_rate, plot_ms_start, plot_ms_end, '1_filtered_signal.png');
        plot_signal_time_domain(differentiated, sampling_rate, plot_ms_start, plot_ms_end, '2_differentiator.png');
        plot_signal_time_domain(squared_energy, sampling_rate, plot_ms_start, plot_ms_end, '3_squared.png');
        plot_signal_time_domain(integrated,     sampling_rate, plot_ms_start, plot_ms_end, '4_integrated.png');
        plot_r_peaks_on_raw(MLII_raw, QRS_peaks,   sampling_rate, plot_ms_start, plot_ms_end, '5_detected.png');
    end

    % OUTPUT
    idx = QRS_peaks;
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


function peaks = global_maxima(yb, W)
    half = floor(W/2);
    N = length(yb);

    peaks = [];
    for k = half+1 : N-half
        win = yb(k-half : k+half);
        
        if yb(k) == max(win) && sum(win == yb(k)) == 1
            peaks(end+1) = k;
        end
    end

end


function QRS_peaks = min_dist_classifier(yb, W, fs)

    init_seconds = 2;
    refractory_period_samples = floor(0.2 * fs);
    class_separation_coef = 0.15;

    peak_indices = global_maxima(yb, W);

    init_seg = yb(1 : min(length(yb), init_seconds*fs));

    QRS = max(init_seg);
    nQRS = mean(init_seg(init_seg < mean(init_seg)));
    THR = nQRS + class_separation_coef * (QRS - nQRS);

    last_QRS = -inf;
    QRS_peaks = [];

    for k = 1:length(peak_indices)

        ind = peak_indices(k);
        val  = yb(ind);

        is_QRS = val > THR && (ind - last_QRS) > refractory_period_samples;

        if is_QRS
            last_QRS = ind;
            QRS = 0.125 * val + 0.875 * QRS;
            QRS_peaks(end+1) = ind;
        else
            nQRS = 0.125 * val + 0.875 * nQRS;
        end

        THR = nQRS + class_separation_coef * (QRS - nQRS);

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
