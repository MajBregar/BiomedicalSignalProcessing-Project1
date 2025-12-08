function Run_Detector( record_path, sampling_rate, output_plots )

    record_matlab_file = sprintf('%sm.mat', record_path);
    fprintf('Loaded File: %s with sampling freq %f:\n', record_matlab_file, sampling_rate)

    t=cputime();

    % run the detector
    idx = QRS_Detector_IMPROVED(record_matlab_file, sampling_rate, output_plots);

    [~, base, ~] = fileparts(record_path);
    ascii_label_file = sprintf('../converted_ltst/%s.asc', base);
    fid = fopen(ascii_label_file, 'wt');
    for i=1:size(idx,2)
        fprintf(fid,'0:00:00.00 %d N 0 0 0\n', idx(1,i) );
    end
    fclose(fid);

    fprintf('FINISHED - Ran for time: %f\n', cputime() - t);

end