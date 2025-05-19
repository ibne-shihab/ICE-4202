%% DS-CDMA with 1/2-rate Convolutional Encoding in AWGN Channel
% Author: [Your Name]
% Description: BER performance evaluation of a convolutionally encoded DS-CDMA system in AWGN.

clear; clc; close all;

%% Parameters
msgLength     = 1000;                   % Number of bits
bitRate       = 1e3;                    % Bit rate (Hz)
chipRate      = 10e3;                   % Chip rate (Hz)
fc            = 5e3;                    % Carrier frequency (Hz)
Eb            = 0.5;                    % Energy per bit (Joules)
tb            = 1/bitRate;              % Bit duration (s)
tc            = 1/chipRate;             % Chip duration (s)
chipsPerBit   = chipRate / bitRate;
samplesPerBit = chipsPerBit;            % 1 chip = 1 sample for simplicity
snr_dB        = 0:1:10;                 % SNR in dB

%% Generate Random Message
msg = randi([0 1], 1, msgLength);

%% Convolutional Encoding (1/2 rate)
trellis = poly2trellis(3, [6 7]);
encodedBits = convenc(msg, trellis);
dataSymbols = 2 * encodedBits - 1;  % BPSK Mapping: 0 → -1, 1 → +1

%% Time Vector
t = 0:tc:(length(dataSymbols)*tb - tc);

%% Generate Baseband Signal
basebandSig = repelem(dataSymbols, samplesPerBit);

%% BPSK Modulation
carrier = cos(2 * pi * fc * t);
bpskMod = sqrt(2 * Eb) * basebandSig .* carrier;

%% PN Sequence Generator (Simple LFSR)
pnSeed = [1 -1 1 -1];
pnSeq = generatePN(pnSeed, length(dataSymbols), chipsPerBit);
pnUpsampled = repelem(pnSeq, 1);

%% Spread Signal
spreadSig = bpskMod .* pnUpsampled;

%% Plot Original Baseband Signal
figure('Name', 'Baseband Signal', 'NumberTitle', 'off');
plot(t(1:1000), basebandSig(1:1000), 'b', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Amplitude');
title('Baseband Signal (First 1000 Samples)');
grid on; axis tight;

%% Frequency Spectrum
figure('Name', 'BPSK Spectrum', 'NumberTitle', 'off');
N = length(t);
f = (0:N-1)*(2*fc/N);
plot(f, abs(fft(bpskMod)), 'r', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of BPSK Modulated Signal');
grid on; axis tight;

%% Transmitted Signal Plot
figure('Name', 'Transmitted DS-CDMA Signal', 'NumberTitle', 'off');
plot(t(1:1000), spreadSig(1:1000), 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Segment of Transmitted DS-CDMA Signal');
grid on; axis tight;

%% BER Simulation over AWGN Channel
ber = zeros(1, length(snr_dB));

for idx = 1:length(snr_dB)
    % Add AWGN
    rxSig = awgn(spreadSig, snr_dB(idx), 'measured');
    
    % Despread
    despreadSig = rxSig .* pnUpsampled;
    
    % Coherent Demodulation
    demodCarrier = sqrt(2 * Eb) * cos(2 * pi * fc * t);
    demodulated = despreadSig .* demodCarrier;
    
    % Integrate over bit intervals
    rxBits = [];
    for i = 1:length(dataSymbols)
        bitSamples = demodulated((i-1)*samplesPerBit + 1 : i*samplesPerBit);
        rxBits(i) = sum(bitSamples) > 0;
    end
    
    % Viterbi Decoder (hard decision)
    traceback = 3;
    decodedBits = vitdec(rxBits, trellis, traceback, 'cont', 'hard');
    
    % Calculate BER (ignoring traceback delay)
    delay = traceback;
    [~, ber(idx)] = biterr(decodedBits(delay+1:end), msg(1:end-delay));
end

%% Plot BER Performance
figure('Name', 'BER Performance', 'NumberTitle', 'off');
semilogy(snr_dB, ber, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR for Convolutionally Coded DS-CDMA in AWGN');
legend('Simulation', 'Location', 'southwest');
grid on; axis tight;



function pnSeq = generatePN(seed, numBits, chipsPerBit)
    % Generates a PN sequence for a given seed and length
    pnSeq = [];
    reg = seed;
    for i = 1:numBits
        for j = 1:chipsPerBit
            pnSeq(end+1) = reg(4);
            % LFSR tap: XOR(4,3) mapped to ±1
            temp = xor(reg(4)==1, reg(3)==1) * 2 - 1;
            reg = [temp reg(1:3)];
        end
    end
end

