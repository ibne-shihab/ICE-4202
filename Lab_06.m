clear all;
close all;
clc;

% Parameters
f = 1000;         % Frequency of the audio signal
Fs = 4000;        % Sampling rate
t = 1/Fs : 1/Fs : 1; % Time vector for 1 second
Am = 1.0;         % Amplitude

% Generate original sinusoidal signal
signal = Am * sin(2*pi*f*t);
figure(1);
plot(t(1:200), signal(1:200))
set(gca,'ytick',[-1 0 1])
title('A segment of synthetically generated sinusoidal waveform')
grid on
xlabel('Time (sec)')
ylabel('Amplitude (volt)')

% Quantization parameters
maximumvalue = max(signal);
minimumvalue = min(signal);
interval = (maximumvalue - minimumvalue) / 255; % 256 levels
partition = minimumvalue:interval:maximumvalue;
codebook = (minimumvalue - interval):interval:maximumvalue;

% Quantize signal
[index, ~, ~] = quantiz(signal, partition, codebook);

% Convert quantized decimal values to binary bits (8 bits per sample)
indxtrn = index';
matrix = zeros(length(indxtrn),8);
for i = 1:length(indxtrn)
    matrix(i, :) = bitget(uint8(indxtrn(i)), 1:8);
end

% Transpose matrix to get bits in columns
matrixtps = matrix'; % 8 rows x 4000 columns

% Convert bit matrix into a serial baseband bitstream
baseband = reshape(matrixtps, [], 1); % 32000 bits (8*4000)
Tb = 1/32000;    % Bit duration
time = 0 : Tb : (length(baseband)-1)*Tb;

% Plot baseband signal segment
figure(2);
stairs(time(1:500), baseband(1:500))
title('A segment of baseband signal')
xlabel('Time (sec)')
ylabel('Binary value')
set(gca,'ytick',[0 1])
axis([0 time(500) 0 1])

% Convolutional Encoder Setup
input_to_Convolutional_encoder = baseband'; % Row vector
trellis = poly2trellis(7, [171 133]);
code = convenc(input_to_Convolutional_encoder, trellis); % Encoded bits

% Interleaving
interleaver_depth = 4831;
data_interleave = randintrlv(code, interleaver_depth);

% QPSK Modulation
M = 4;
k = log2(M);
% Ensure length is multiple of k
if mod(length(data_interleave), k) ~= 0
    data_interleave = [data_interleave zeros(1, k - mod(length(data_interleave), k))];
end
symbol = bi2de(reshape(data_interleave, k, length(data_interleave)/k).', 'left-msb');
qpsk_modulated_data = pskmod(symbol, M);

% QPSK Demodulation
qpsk_demodulated_symbol = pskdemod(qpsk_modulated_data, M);

% Symbol error calculation
[number, ratio] = symerr(symbol, qpsk_demodulated_symbol);
fprintf('Symbol errors: %d, Symbol error rate: %f\n', number, ratio);

% Convert demodulated symbols back to bits
retrieved_bits = de2bi(qpsk_demodulated_symbol, k, 'left-msb')';
retrieved_bits = reshape(retrieved_bits, [], 1);

% Deinterleaving
errors = zeros(size(retrieved_bits));
interleaved_with_errors = bitxor(retrieved_bits, errors);
data_deinterleave = randdeintrlv(interleaved_with_errors, interleaver_depth);

% Viterbi Decoding
tblen = 3;
decoded_bits = vitdec(data_deinterleave, trellis, tblen, 'cont', 'hard');

% Trim first tblen bits (due to decoding delay)
decoded_bits_trimmed = decoded_bits(tblen+1:end)';
baseband_trimmed = baseband(tblen+1:end)';

% Ensure vectors have same length for BER calculation
min_length = min(length(decoded_bits_trimmed), length(baseband_trimmed));
decoded_bits_trimmed = decoded_bits_trimmed(1:min_length);
baseband_trimmed = baseband_trimmed(1:min_length);

% Bit error rate after decoding
[~, ber1] = biterr(decoded_bits_trimmed, baseband_trimmed);
fprintf('Bit Error Rate after decoding: %f\n', ber1);

% Reshape trimmed decoded bits to 8-bit rows (bytes)
if mod(length(decoded_bits_trimmed), 8) ~= 0
    padding = zeros(1, 8 - mod(length(decoded_bits_trimmed), 8));
    decoded_bits_trimmed = [decoded_bits_trimmed padding];
end
decoded_matrix = reshape(decoded_bits_trimmed, 8, [])';

% Check if sizes match before BER check
if size(decoded_matrix,1) == size(matrixtps,2)
    [~, ber2] = biterr(decoded_matrix', matrixtps);
    fprintf('Bit Error Rate between decoded bits and original bit matrix: %f\n', ber2);
else
    fprintf('Size mismatch between decoded matrix (%d) and original bit matrix (%d).\n', ...
            size(decoded_matrix,1), size(matrixtps,2));
end

% Convert decoded bits to decimal values
int_decoded = bi2de(decoded_matrix(1:length(index),:));

% BER between decoded integers and quantization indices
[~, ber3] = biterr(int_decoded, index');
fprintf('Bit Error Rate between decoded integers and original quantized indices: %f\n', ber3);

% Reconstruct signal from decoded quantized levels
sample_value = minimumvalue + int_decoded * interval;

% Plot original and retrieved audio segments
figure(3);
subplot(2,1,1)
plot(t(1:100), signal(1:100));
set(gca, 'ytick', [-1 0 1])
axis([0 t(100) -1 1])
title('Original audio signal segment')
xlabel('Time (sec)')
ylabel('Amplitude')
grid on

subplot(2,1,2)
plot(t(1:100), sample_value(1:100));
set(gca, 'ytick', [-1 0 1])
axis([0 t(100) -1 1])
title('Retrieved audio signal segment after decoding')
xlabel('Time (sec)')
ylabel('Amplitude')
grid on