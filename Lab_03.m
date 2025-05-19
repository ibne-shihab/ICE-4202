% 1/2-Rated Convolutionally Encoded DS-CDMA System BER Simulation
% under AWGN and Rician fading channel

clear; close all; clc;

% Parameters
msgLen = 1000;             % Number of bits in message
fc = 5000;                 % Carrier frequency (Hz)
Eb = 0.5;                  % Energy per bit
bitrate = 1000;            % Bitrate (bits/sec)
tb = 1/bitrate;            % Bit duration (sec)
chiprate = 10000;          % Chip rate (chips/sec)
tc = 1/chiprate;           % Chip duration (sec)

% Generate random message bits
msg = randi([0 1],1,msgLen);

% --- Figure 1: Original Message bits (first 100 bits) ---
figure(1);
stem(msg(1:100), 'filled');
title('Original Message Bits (Segment)');
xlabel('Bit Index');
ylabel('Bit Value');
ylim([-0.1 1.1]);
grid on;

% Convolutional Encoder: rate 1/2, constraint length 3, generators [6 7] octal
trellis = poly2trellis(3,[6 7]);
coded = convenc(msg, trellis);

% Map 0->-1 and 1->1 (bipolar NRZ)
codedBipolar = 2*coded - 1;

% Length of coded sequence
N = length(codedBipolar);

% Time vector for chips
t = tc:tc:tc*N;

% BPSK Modulation
bpskMod = sqrt(2*Eb) * codedBipolar .* cos(2*pi*fc*t);

% --- Figure 2: BPSK Modulated Signal (segment) ---
figure(2);
plot(t(1:500), bpskMod(1:500));
title('BPSK Modulated Signal (Segment)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% PN Sequence Generation (length = number of chips)
seed = [1 -1 1 -1]; % Bipolar initial seed
pn = zeros(1,N);
state = seed;

for i=1:N
    pn(i) = state(end);
    newBit = -state(end)*state(end-1); % XOR in bipolar: product
    state = [newBit state(1:end-1)];
end

% Spread the BPSK modulated signal with PN sequence
spreadSig = bpskMod .* pn;

% --- Figure 3: Spread DS-CDMA Signal (segment) ---
figure(3);
plot(t(1:500), spreadSig(1:500));
title('Spread DS-CDMA Signal (Segment)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Rician channel object using comm.RicianChannel
ricianChan = comm.RicianChannel(...
    'SampleRate', 1/tc, ...
    'PathDelays', 0, ...
    'KFactor', 15, ...
    'DirectPathDopplerShift', 0, ...
    'MaximumDopplerShift', 0, ...
    'DirectPathInitialPhase', 0);

% SNR range (dB)
snr_dB = 0:2:10;
ber = zeros(size(snr_dB));

% Traceback length for Viterbi decoder
% Correct way to get constraint length:
constraintLength = log2(trellis.numStates) + 1;
tblen = 5 * constraintLength;

for idx = 1:length(snr_dB)
    % Reset channel for each SNR to get independent fading realizations
    reset(ricianChan);
    
    % Pass signal through Rician channel
    fadedSig = ricianChan(spreadSig')';  % transpose for column vector, then back
    
    % Add AWGN noise
    rxSig = awgn(fadedSig, snr_dB(idx), 'measured');
    
    % Despread (multiply by PN sequence)
    despread = rxSig .* pn;
    
    % Demodulate BPSK by correlating with carrier
    carrier = sqrt(2*Eb) * cos(2*pi*fc*t);
    demod = despread .* carrier;
    
    % Sum demodulated chips over each bit interval
    chipsPerBit = chiprate / bitrate; % 10 chips per bit
    numBitsCoded = length(demod) / chipsPerBit;
    decisionVar = zeros(1, numBitsCoded);
    for k = 1:numBitsCoded
        decisionVar(k) = sum(demod((k-1)*chipsPerBit + 1 : k*chipsPerBit));
    end
    
    % Decide bits from sign of decision variable
    rxBits = decisionVar > 0;
    
    % Viterbi decoding (hard decision)
    decodedBits = vitdec(rxBits, trellis, tblen, 'cont', 'hard');
    
    % Remove delay from decoding
    delay = tblen;
    decodedBits = decodedBits(delay+1:end);
    msgTrimmed = msg(1:length(decodedBits));
    
    % Calculate Bit Error Rate (BER)
    [~, ber(idx)] = biterr(decodedBits, msgTrimmed);
end

% --- Figure 4: BER vs SNR ---
figure(4);
semilogy(snr_dB, ber, '-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance of 1/2-rate Convolutionally Encoded DS-CDMA System');
legend('Rician Fading + AWGN');