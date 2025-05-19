clear all;
close all;
clc;
%% Parameters for Synthetic Signal
f = 1000;       % Frequency of audio signal (Hz)
Fs = 4000;      % Sampling frequency (samples/sec)
t = (1/Fs:1/Fs:1); % Time vector (1 sec duration)
Am = 1.0;       % Amplitude of signal

% Generate original sinusoidal signal
signal = Am * sin(2 * pi * f * t);

% Plot original sinusoidal waveform (Segment)
figure(1);
plot(t(1:200), signal(1:200));
set(gca, 'ytick', [-1 0 1]);
title('Segment of Original Sinusoidal Waveform');
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
grid on;

%% Quantization (8-bit uniform)
maxVal = max(signal);
minVal = min(signal);
interval = (maxVal - minVal) / 255;
partition = minVal:interval:maxVal;
codebook = (minVal - interval):interval:maxVal;

[index, quants, distor] = quantiz(signal, partition, codebook);

% Convert decimal indices to binary (8 bits per sample)
indxtrn = index';
matrix = zeros(length(indxtrn), 8);

for i = 1:length(indxtrn)
    matrix(i, :) = bitget(uint8(indxtrn(i)), 1:8);
end

matrixtps = matrix'; % 8 x N matrix

% Create baseband bit stream (column vector)
baseband = reshape(matrixtps, [], 1);

Tb = 1 / (Fs * 8);  % Bit duration (32 kbps bit rate)
time = 0:Tb:Tb * (length(baseband) - 1);

% Plot baseband binary signal segment
figure(2);
stairs(time(1:500), baseband(1:500));
title('Segment of Baseband Binary Signal');
xlabel('Time (sec)');
ylabel('Binary value');
set(gca, 'ytick', [0 1]);
axis([0 time(500) 0 1]);
grid on;

%% Convolutional Encoding
input_to_Convolutional_encoder = baseband';  % 1 x total bits

trellis = poly2trellis(7, [171 133]);  % Constraint length 7, polynomials in octal
code = convenc(input_to_Convolutional_encoder, trellis);

%% Interleaving
interleaver_depth = 4831;
data_interleave = randintrlv(code, interleaver_depth);

%% BPSK Modulation
M = 2;
k = log2(M);

% Convert bits to symbols
symbol = bi2de(reshape(data_interleave, k, []).', 'left-msb');
symbol = double(symbol);

bpsk_modulated_data = pskmod(symbol, M);

%% BPSK Demodulation
bpsk_demodulated = pskdemod(bpsk_modulated_data, M);

% Symbol error count
[num_err, ratio_err] = symerr(symbol, bpsk_demodulated);

% Convert symbols back to bits
retrieved_bits = de2bi(bpsk_demodulated, k, 'left-msb');

%% Deinterleaving and Decoding
% Include burst error simulation (zero errors here)
errors = zeros(size(retrieved_bits));
inter_err = bitxor(retrieved_bits, errors);

data_deinterleave = randdeintrlv(inter_err, interleaver_depth);

tblen = 3; % Traceback length for Viterbi decoding
decodx = vitdec(data_deinterleave, trellis, tblen, 'cont', 'hard');

% Fix indexing for decoded bits
N3 = length(decodx);
decod2 = zeros(1, N3);
decod2(1:(N3 - tblen)) = decodx(tblen+1:end);
decod2(N3) = decodx(1);

decod2 = decod2';  % Column vector

%% Error Checking
[number_bits_err, bit_error_ratio] = biterr(decod2, baseband);

convert = reshape(decod2, 8, length(signal));  % 8 x 4000 matrix
convert = convert';

[number_conv_err, conv_error_ratio] = biterr(convert, matrixtps');

% Convert binary back to decimal
intconv = bi2de(convert);

[number_final_err, final_error_ratio] = biterr(intconv, index');

% Reconstruct signal from quantized levels
reconstructed_signal = minVal + intconv * interval;

%% Plot original vs reconstructed signals
figure(3);
subplot(2,1,1);
plot(time(1:100), signal(1:100));
set(gca, 'ytick', [-1 0 1]);
axis([0 time(100) -1 1]);
title('Segment of Original Audio Signal');
xlabel('Time (sec)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(time(1:100), reconstructed_signal(1:100));
set(gca, 'ytick', [-1 0 1]);
axis([0 time(100) -1 1]);
title('Segment of Retrieved Audio Signal after Decoding');
xlabel('Time (sec)');
ylabel('Amplitude');
grid on;
