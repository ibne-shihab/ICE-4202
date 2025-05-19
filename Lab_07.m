clc;
clear all;
close all;


% Parameters
f = 1000;            % Frequency of audio signal (Hz)
Fs = 4000;           % Sampling frequency (samples/sec)
t = 1/Fs : 1/Fs : 1; % Time vector for 1 second
Am = 1.0;            % Amplitude

% Generate sinusoidal signal
signal = Am * sin(2 * pi * f * t);

% Plot original signal segment
figure(1);
plot(t(1:200), signal(1:200));
set(gca, 'ytick', [-1 0 1]);
title('Original Audio Signal Segment');
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
grid on;

% Quantization parameters
maximumvalue = max(signal);
minimumvalue = min(signal);
interval = (maximumvalue - minimumvalue) / 255;
partition = minimumvalue : interval : maximumvalue;
codebook = (minimumvalue - interval) : interval : maximumvalue;

% Quantize the signal
[index, quants, distor] = quantiz(signal, partition, codebook);

% Convert to 8-bit binary representation
indxtrn = index';
matrix = zeros(length(indxtrn), 8);
for i = 1:length(indxtrn)
    matrix(i, :) = bitget(uint8(indxtrn(i)), 1:8);
end

% Create baseband bitstream
matrixtps = matrix';
baseband = reshape(matrixtps, [], 1);
Tb = 1/32000;
time = 0 : Tb : (length(baseband)-1)*Tb;

% Plot baseband signal
figure(2);
stairs(time(1:500), baseband(1:500));
title('Baseband Signal Segment');
xlabel('Time (sec)');
ylabel('Binary value');
set(gca, 'ytick', [0 1]);
axis([0 time(500) 0 1]);

% Convolutional Encoding
trellis = poly2trellis(7, [171 133]);
code = convenc(baseband', trellis);

% Interleaving
interleaver_depth = 4831;
data_interleave = randintrlv(code, interleaver_depth);

% QAM Modulation
M = 4;
k = log2(M);
% Pad with zeros if needed
if mod(length(data_interleave), k) ~= 0
    data_interleave = [data_interleave zeros(1, k-mod(length(data_interleave),k))];
end
symbol = bi2de(reshape(data_interleave, k, []).', 'left-msb');
qam_modulated_data = qammod(symbol, M);

% QAM Demodulation
qam_demodulated_symbol = qamdemod(qam_modulated_data, M);

% Symbol error rate
[number, ratio] = symerr(symbol, qam_demodulated_symbol);
fprintf('Symbol errors: %d, Symbol error rate: %f\n', number, ratio);

% Convert back to bits
retrieved_bits = de2bi(qam_demodulated_symbol, k, 'left-msb')';
retrieved_bits = reshape(retrieved_bits, [], 1);

% Deinterleaving
data_deinterleave = randdeintrlv(retrieved_bits, interleaver_depth);

% Viterbi Decoding
tblen = 3;
decoded_bits = vitdec(data_deinterleave, trellis, tblen, 'cont', 'hard');

% Trim and align vectors
decoded_bits_trimmed = decoded_bits(tblen+1:end)';
baseband_trimmed = baseband(tblen+1:end)';

% Ensure equal length
min_len = min(length(decoded_bits_trimmed), length(baseband_trimmed));
decoded_bits_trimmed = decoded_bits_trimmed(1:min_len);
baseband_trimmed = baseband_trimmed(1:min_len);

% Calculate BER
[~, ber1] = biterr(decoded_bits_trimmed, baseband_trimmed);
fprintf('Bit Error Rate after decoding: %f\n', ber1);

% Reshape to bytes
if mod(length(decoded_bits_trimmed), 8) ~= 0
    padding = zeros(1, 8-mod(length(decoded_bits_trimmed),8));
    decoded_bits_trimmed = [decoded_bits_trimmed padding];
end
decoded_matrix = reshape(decoded_bits_trimmed, 8, [])';

% Compare with original
if size(decoded_matrix,1) == size(matrixtps,2)
    [~, ber2] = biterr(decoded_matrix', matrixtps);
    fprintf('Bit Error Rate vs original bit matrix: %f\n', ber2);
else
    fprintf('Size mismatch in bit matrix comparison\n');
end

% Convert to samples
int_decoded = bi2de(decoded_matrix(1:length(index),:));
sample_value = minimumvalue + int_decoded * interval;

% Plot comparison
figure(3);
subplot(2,1,1);
plot(t(1:100), signal(1:100));
title('Original Signal');
xlabel('Time (sec)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t(1:100), sample_value(1:100));
title('Decoded Signal');
xlabel('Time (sec)');
ylabel('Amplitude');
grid on;