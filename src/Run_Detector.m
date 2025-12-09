function Run_Detector( record_path, sampling_rate, output_plots )

    [record_folder, record_id, ~] = fileparts(record_path);

    fprintf('Record folder: %s\n', record_folder);
    fprintf('Record ID: %s\n', record_id);

    record_matlab_file = fullfile(record_folder, sprintf('%sm.mat', record_id));
    fprintf('Loaded File: %s with sampling freq %f\n', record_matlab_file, sampling_rate);

    t = cputime();

    #idx = QRS_Detector_IMPROVED(record_matlab_file, sampling_rate, output_plots);
    idx = QRS_Detector_ORIGINAL(record_matlab_file, sampling_rate, output_plots);


    ascii_label_file = fullfile(record_folder, sprintf('%s.asc', record_id));
    fid = fopen(ascii_label_file, 'wt');

    for i = 1:size(idx, 2)
        fprintf(fid, '0:00:00.00 %d N 0 0 0\n', idx(1, i));
    end

    fclose(fid);

    fprintf('FINISHED - Ran for time: %f\n', cputime() - t);
end
