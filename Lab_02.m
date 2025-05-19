% Convolutionally Encoded DS-CDMA System in AWGN and Rayleigh Fading
clear; clc; close all;

%% Parameters
msg = randi([0 1], 1, 1000);                   % Random binary message
trellis = poly2trellis(3, [6 7]);              % 1/2-rate convolutional encoder
encodedBits = convenc(msg, trellis);           % Convolutional encoding
encodedBits(encodedBits==0) = -1;              % BPSK Mapping: 0 -> -1, 1 -> +1

fc = 5e3;                                      % Carrier frequency (Hz)
Eb = 0.5;                                      % Energy per bit
bitrate = 1e3;                                 % Bitrate (bps)
tb = 1/bitrate;                                % Time per bit
chiprate = 10e3;                               % Chiprate (10x bit rate)
tc = 1/chiprate;                               % Time per chip
chipsPerBit = chiprate / bitrate;

t = 0:tc:tb*(length(encodedBits))-tc;          % Time vector

%% Baseband Signal
basebandsig = repelem(encodedBits, chipsPerBit);

figure;
stairs(t(1:800), basebandsig(1:800), 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Binary Value');
title('Baseband Signal for Encoded Data');
ylim([-1.5, 1.5]); grid on;

%% BPSK Modulation
bpskmod = sqrt(2*Eb) * basebandsig .* cos(2*pi*fc*t);

%% Spectrum Analysis
L = length(bpskmod);
spectrum = abs(fft(bpskmod));
f = (0:L-1) * (chiprate / L);                  % Frequency vector

figure;
plot(f, spectrum, 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of BPSK Modulated Signal');
grid on;

%% PN Sequence Generator (for 1 user)
pn_seed = [1 -1 1 -1];                         % Bipolar NRZ format
pn = zeros(1, length(encodedBits) * chipsPerBit);
register = pn_seed;

for i = 1:length(encodedBits)
    for j = 1:chipsPerBit
        pn((i-1)*chipsPerBit + j) = register(4);
        feedback = xor(register(4) == 1, register(3) == 1);
        register = [1 -2*feedback+1 register(1:3)];
    end
end

pn(pn == 0) = -1;

%% Spreading (DS-CDMA Transmit Signal)
sigtx = bpskmod .* pn;

figure;
plot(t(1:200), sigtx(1:200), 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Segment of DS-CDMA Transmitted Signal');
grid on;

%% Channel: Rayleigh Fading + AWGN
rayleighChan = comm.RayleighChannel(...
    'SampleRate', chiprate, ...
    'PathDelays', 0, ...
    'AveragePathGains', 0, ...
    'MaximumDopplerShift', 100);

%% BER Simulation over SNR
snr_dB = 0:1:10;
ber = zeros(size(snr_dB));

for k = 1:length(snr_dB)
    reset(rayleighChan);                               % Reset channel
    fadedSig = rayleighChan(sigtx.');                  % Apply fading
    noisySig = awgn(fadedSig, snr_dB(k), 'measured');  % Add AWGN
    
    %% Receiver: Despreading
    rx = (noisySig.').* pn;                            % Despread
    
    %% BPSK Demodulation
    demod = rx .* cos(2*pi*fc*t);                      % Coherent demod
    demod_reshaped = reshape(demod, chipsPerBit, []);
    sum_samples = sum(demod_reshaped, 1);
    rxbits = double(sum_samples > 0);
    
    %% Viterbi Decoding
    traceback = 3;
    decoded = vitdec(rxbits, trellis, traceback, 'cont', 'hard');
    
    %% BER Calculation (delay compensation)
    delay = traceback;
    [~, ber(k)] = biterr(decoded(delay+1:end), msg(1:end-delay));
end

%% BER Plot
figure;
semilogy(snr_dB, ber, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR under Rayleigh Fading with AWGN');
legend('1/2-Rate Coded DS-CDMA');
grid on;
