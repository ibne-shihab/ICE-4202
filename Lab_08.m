clc;
clear all;
close all;

%% 1. Generate Synthetic Audio Signal
Fs = 4000;              
f = 1000;               
t = 0:1/Fs:1-1/Fs;      
Am = 1.0;
signal = Am * sin(2*pi*f*t);  

% Plot 1: Original waveform
figure(1);
plot(t(1:200), signal(1:200));
title('Original Audio Signal Segment');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 2. Quantization
max_val = max(signal);
min_val = min(signal);
interval = (max_val - min_val) / 255;
partition = min_val:interval:max_val;
codebook = [min_val-interval:interval:max_val];
[index, quants, distor] = quantiz(signal, partition, codebook);

%% 3. Convert to Binary Bit Stream
index = index';
binary_matrix = zeros(length(index), 8);
for i = 1:length(index)
    binary_matrix(i,:) = bitget(uint8(index(i)), 8:-1:1);  
end
binary_matrix = binary_matrix';
bit_stream = reshape(binary_matrix, [], 1);  

% Plot 2: Baseband signal (binary)
Tb = 1 / length(bit_stream);
time_bit = 0:Tb:(length(bit_stream)-1)*Tb;
figure(2);
stairs(time_bit(1:500), bit_stream(1:500));
title('Baseband Binary Signal Segment');
xlabel('Time (s)');
ylabel('Binary Value');
axis([0, time_bit(500), 0, 1]);
grid on;

%% 4. FEC Encoding (Convolutional)
trellis = poly2trellis(7, [171 133]);
coded_bits = convenc(bit_stream', trellis);  

%% 5. Interleaving
seed = 4831;
interleaved_bits = randintrlv(coded_bits, seed);

%% 6. 16-QAM Modulation
M = 16;
k = log2(M);  
symbols = bi2de(reshape(interleaved_bits, k, []).', 'left-msb');
modulated_signal = qammod(symbols, M, 'UnitAveragePower', true);

%% 7. 16-QAM Demodulation
demodulated_symbols = qamdemod(modulated_signal, M, 'UnitAveragePower', true);
retrieved_bits = de2bi(demodulated_symbols, k, 'left-msb')';
retrieved_bits = reshape(retrieved_bits, [], 1);

%% 8. Deinterleaving
deinterleaved_bits = randdeintrlv(retrieved_bits', seed);

%% 9. Viterbi Decoding
tblen = 35;  
decoded_bits = vitdec(deinterleaved_bits, trellis, tblen, 'cont', 'hard');
decoded_bits = decoded_bits(tblen+1:end);  
decoded_bits = decoded_bits(:);

%% 10. Bit Stream Trimming and Reshape
min_len = min(length(bit_stream), length(decoded_bits));
bit_stream_trimmed = bit_stream(1:min_len);
decoded_bits_trimmed = decoded_bits(1:min_len);

% Make length divisible by 8
L = floor(min_len / 8) * 8;
bit_stream_trimmed = bit_stream_trimmed(1:L);
decoded_bits_trimmed = decoded_bits_trimmed(1:L);

% Reconstruct quantized values
reshape_bits = reshape(decoded_bits_trimmed, 8, []).'; 
int_values = bi2de(reshape_bits, 'left-msb');
reconstructed_signal = min_val + int_values * interval;

% Plot 3: Original vs Reconstructed Signal
figure(3);
subplot(2,1,1);
plot(t(1:100), signal(1:100));
title('Original Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t(1:100), reconstructed_signal(1:100));
title('Retrieved Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 11. Bit Error Rate
[bit_errors, bit_error_rate] = biterr(bit_stream_trimmed, decoded_bits_trimmed);
fprintf('\nTotal Bit Errors: %d\nBit Error Rate: %.5f\n', bit_errors, bit_error_rate);
