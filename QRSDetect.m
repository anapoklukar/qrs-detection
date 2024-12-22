function idx = QRSDetect(fileName)
    % Load the ECG signal
    % This function detects QRS complexes in an ECG signal stored in a .mat file.
    data = load(fileName);  % Load the .mat file containing ECG data
    
    x = data.val;  % Extract ECG data (assuming it's stored in 'val')
    
    % Extract individual channels (assuming a two-channel ECG signal)
    sig1 = x(1, :);  % Channel 1
    sig2 = x(2, :);  % Channel 2

    fs = 250;  % Sampling frequency in Hz

    %coef = 360/250; % The article wants the sampling frequency 250Hz
    coef = 1;

    % Step 1: Noise Filtering
    fprintf('Applying noise filter to Channel 1...\n');
    filtered = NoiseFilter(sig1, coef);

    % Step 2: Band-Limited Differentiation
    fprintf('Applying band-limited differentiator to filtered signal...\n');
    differentiated = BandLimitedDifferentiator(filtered);

    % Step 3: Energy Collector
    fprintf('Applying energy collector...\n');
    energy = EnergyCollector(differentiated, coef);

    % Step 4: Adaptive Minimum Distance Classifier
    fprintf('Classifying QRS complexes...\n');
    qrs_indices = AdaptiveMinimumDistanceClassifier(energy, fs);

    % Step 5: Minimax Searcher
    fprintf('Refining QRS complex locations using minimax search...\n');
    %refined_qrs_indices = MiniMaxSearcher(sig1, differentiated, qrs_indices);
    refined_qrs_indices = MiniMaxSearcherImproved(sig1, qrs_indices); % Improved minimax searcher

    % Improvement: Adaptive Thresholding to Refine QRS Indices
    fprintf('Applying adaptive thresholding to refine QRS indices...\n');
    % Define the threshold as a fraction of the maximum signal amplitude
    threshold_fraction = 0.05; % Adjust this value to tune the sensitivity
    adaptive_threshold = threshold_fraction * max(abs(sig1));
    % Retain indices where the signal exceeds the adaptive threshold
    refined_qrs_indices = refined_qrs_indices(sig1(refined_qrs_indices) >= adaptive_threshold);

    % % Plot Original ECG Signal with Detected R Peaks
    % figure;
    % plot(sig1, 'k', 'LineWidth', 1.2);  % Plot the raw ECG signal
    % hold on;
    % % Plot detected R-peaks
    % plot(refined_qrs_indices, sig1(refined_qrs_indices), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
    % % Add labels and legend
    % title('ECG Signal with Detected R Peaks');
    % xlabel('Sample Number');
    % ylabel('Amplitude');
    % legend('ECG Signal', 'R Peaks');
    % grid on;

    % Return only the unique detected QRS indices as output
    idx = unique(refined_qrs_indices);

end

% NoiseFilter: Removes low-frequency drifts and high-frequency noise using
% two moving average filters.
% NoiseFilter: Removes low-frequency drifts and high-frequency noise using
% two moving average filters.
function y = NoiseFilter(ecg_signal, coef)
    % Define window sizes
    K = round(5 * coef);   % Short-term moving average window
    L = round(200 * coef); % Long-term moving average window

    % Ensure the input is a row vector
    ecg_signal = ecg_signal(:)';

    % Initialize moving averages and output
    N = length(ecg_signal);
    short_term_avg = zeros(1, N);
    long_term_avg = zeros(1, N);
    y = zeros(1, N);

    % Compute short-term moving average
    short_sum = sum(ecg_signal(1:min(K, N)));
    for i = 1:N
        short_term_avg(i) = short_sum / K;
        if i + 1 <= N
            short_sum = short_sum - ecg_signal(max(1, i - K + 1)) + ecg_signal(min(N, i + 1));
        end
    end

    % Compute long-term moving average
    long_sum = sum(ecg_signal(1:min(L, N)));
    for i = 1:N
        long_term_avg(i) = long_sum / L;
        if i + 1 <= N
            long_sum = long_sum - ecg_signal(max(1, i - L + 1)) + ecg_signal(min(N, i + 1));
        end
    end

    % Subtract long-term average from short-term average
    y = short_term_avg - long_term_avg;
end

% BandLimitedDifferentiator: Highlights the steep slopes of the QRS complex
% while filtering out high-frequency noise.
function differentiated_signal = BandLimitedDifferentiator(signal)
    % Impulse response coefficients for the band-limited differentiator
    hd = (1/3) * [-1, -2, 0, 2, 1];  % Derived from the paper

    % Apply the differentiator using convolution
    differentiated_signal = conv(signal, hd, 'same');  % Convolve and retain original length
end

% EnergyCollector: Amplifies QRS complexes by squaring the signal and
% smoothing it with a moving average filter.
function y_b = EnergyCollector(x_k, N)

    % Step 1: Nonlinear squaring (Equation [4])
    y_a = x_k .^ 2;

    % Step 2: Moving average integrator (Equation [5])
    % Initialize the output
    y_b = zeros(size(y_a));

    % Compute the moving average using manual summation
    for k = 1:length(y_a)
        if k >= N
            y_b(k) = sum(y_a(k-N+1:k)) / N;  % Full window
        else
            y_b(k) = sum(y_a(1:k)) / k;  % Partial window for boundary
        end
    end
end

% AdaptiveMinimumDistanceClassifier: Classifies energy maxima as QRS or non-QRS
% and updates the cluster centers adaptively.
function qrs_indices = AdaptiveMinimumDistanceClassifier(energy, fs)
    % Initialization using the first 2 seconds of the signal
    initial_window = 2 * fs;  % Number of samples in the first 2 seconds
    mu_qrs = max(energy(1:initial_window));  % Initial estimate for QRS cluster mean
    mu_nonqrs = min(energy(1:initial_window));  % Initial estimate for non-QRS mean

    % Parameters
    window_size = 15;  % Window size for finding maxima
    qrs_indices = [];  % Initialize array for QRS indices

    % Iterate through the energy signal in windows
    for k = window_size:window_size:length(energy)-window_size
        % Find the maximum value in the current window
        [yb, max_idx] = max(energy(k-window_size+1:k+window_size));
        max_pos = k - window_size + max_idx;  % Adjust position relative to full signal

        % Calculate distances to the QRS and non-QRS cluster means
        d_qrs = abs(yb - mu_qrs);
        d_nonqrs = abs(yb - mu_nonqrs);

        % Classify the maximum based on the minimum distance
        if d_qrs < d_nonqrs
            qrs_indices = [qrs_indices, max_pos];  % Add index to QRS list
            mu_qrs = 0.9 * mu_qrs + 0.1 * yb;  % Update QRS cluster mean
        else
            mu_nonqrs = 0.9 * mu_nonqrs + 0.1 * yb;  % Update non-QRS cluster mean
        end
    end
end

% MiniMaxSearcher Improvement: Refines the detected QRS locations by finding exact R peaks
% in the ECG signal within a specified window around each QRS index.
function [R_peaks] = MiniMaxSearcherImproved(ecg_signal, qrs_indices)
    % Convert window from milliseconds to samples
    window_ms = 250;
    window_samples = round(window_ms / 1000 * 250);

    R_peaks = [];  % Initialize array for R peak locations

    i = 1;
    while i <= length(qrs_indices)
        qrs_index = qrs_indices(i);

        % Define the search window around the QRS index
        segment_start = max(1, qrs_index - window_samples);
        segment_end = min(length(ecg_signal), qrs_index + window_samples);

        % Find the maximum value in the window
        segment = ecg_signal(segment_start:segment_end);
        [R_value, R_index] = max(segment);

        % Compute the absolute index of the maximum point
        refined_index = segment_start + R_index - 1;

        % Append the refined R peak index
        R_peaks = [R_peaks, refined_index];

        % Move to the next window, skipping over the current window
        i = find(qrs_indices > segment_end, 1);
        if isempty(i)
            break;
        end
    end

    % Ensure unique and sorted R peaks
    R_peaks = unique(sort(R_peaks));
end

% MiniMaxSearcher: Refines the detected QRS locations by finding exact R and S
% peaks in the ECG signal using zero crossings in the derivative.
function [R_peaks] = MiniMaxSearcher(ecg_signal, derivative_signal, qrs_indices)
    R_peaks = [];  % Initialize array for R peak locations

    % Iterate over each detected QRS complex
    for i = 1:length(qrs_indices)
        qrs_index = qrs_indices(i);

        % Find the nearest zero crossings in the derivative
        zero_crossings = find_zero_crossings(derivative_signal, qrs_index);

        if length(zero_crossings) < 2
            continue;  % Skip if not enough zero crossings are found
        end

        % Define the segment between the zero crossings
        segment_start = zero_crossings(1);
        segment_end = zero_crossings(2);

        % Find the R peak (maximum) in this segment
        [R_value, R_index] = max(ecg_signal(segment_start:segment_end));
        R_peaks = [R_peaks, segment_start + R_index - 1];
    end
end

% Helper function to find zero crossings in the derivative signal
function zero_crossings = find_zero_crossings(derivative_signal, qrs_index)
    N = length(derivative_signal);

    % Search backward for the previous zero crossing
    prev_zero = qrs_index;
    while prev_zero > 1 && derivative_signal(prev_zero) * derivative_signal(prev_zero - 1) >= 0
        prev_zero = prev_zero - 1;
    end

    % Search forward for the next zero crossing
    next_zero = qrs_index;
    while next_zero < N && derivative_signal(next_zero) * derivative_signal(next_zero + 1) >= 0
        next_zero = next_zero + 1;
    end

    % Collect zero crossings
    zero_crossings = [prev_zero, next_zero];
end
