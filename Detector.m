function Detector(record)
  % Summary: This function performs QRS detection on an ECG record.
  % It converts the record into a MATLAB-compatible format, detects QRS indices,
  % saves them to an ASCII file, and prepares annotations for evaluation.

  % Step 1: Convert the ECG record into a MATLAB-compatible format.
  % This step assumes that the WFDB toolkit is installed and the wfdb2mat command is available.
  % wfdb2mat -r record

  % Construct the filename for the MATLAB-compatible file (recordm.mat)
  fileName = sprintf('./path/%sm.mat', record);
  disp(['Processing file: ', fileName]); % Display the file being processed

  % Step 2: Start timing the QRS detection process
  t = cputime();

  % Call the QRS detection function (QRSDetect) and get indices of detected QRS complexes
  idx = QRSDetect(fileName);

  % Print the time taken for QRS detection
  fprintf('Running time: %f seconds\n', cputime() - t);

  % Step 3: Create an ASCII file to save the QRS detection results
  % The format follows WFDB annotation format conventions
  asciName = sprintf('%s.asc', record); % Define the output ASCII file name
  fid = fopen(asciName, 'wt');          % Open the file in write-text mode

  % Write QRS indices to the ASCII file in the annotation format
  for i = 1:size(idx, 2)
      fprintf(fid, '0:00:00.00 %d N 0 0 0\n', idx(1, i));
  end

  % Close the file after writing
  fclose(fid);

  % Step 4: Convert the ASCII output to binary WFDB annotation format
  % This requires the WFDB toolkit and the `wrann` command.
  % Command: wrann -r record -a qrs <record.asc

  % Step 5: Evaluate the annotations against reference annotations (atr) using bxb
  % This step compares the detected annotations with the reference annotations.
  % Command: bxb -r record -a atr qrs

end
