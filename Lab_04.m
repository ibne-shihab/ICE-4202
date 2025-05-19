clear all;
close all;
clc;

% Input bit stream
xbit = [1 0 1 1 0 1 0 0 0 1 1 0];

% Differential encoding
difencod(1) = ~(1 - xbit(1));
for i = 2:length(xbit)
    difencod(i) = ~(difencod(i-1) - xbit(i));
end

% Differential decoding (to re-derive xbit)
xbit(1) = 1 - ~(difencod(1));
for i = 2:length(xbit)
    xbit(i) = difencod(i-1) - ~(difencod(i));
    if xbit(i) == -1
        xbit(i) = 1;
    end
end

% Prepare inphase (I) and quadrature (Q) streams (unipolar)
inp = zeros(1,length(difencod));
qp = zeros(1,length(difencod));
for i = 1:2:(length(difencod)-1)
    inp(i) = difencod(i);
    inp(i+1) = inp(i);
end
for i = 2:2:length(difencod)
    qp(i) = difencod(i);
    qp(i-1) = qp(i);
end

% Convert unipolar to bipolar NRZ
it = 2*inp - 1;  % 1 -> 1, 0 -> -1
qt = 2*qp - 1;

% Raised cosine filter design using rcosdesign
filtorder = 40;
nsamp = 4;
rolloff = 0.5;
span = filtorder / nsamp; % Number of symbol periods the filter spans
rrcfilter = rcosdesign(rolloff, span, nsamp, 'normal');

% Plot impulse response
figure(1);
impz(rrcfilter,1);
grid on;
title('Impulse Response of Raised Cosine Filter');

% Transmit filtering (upsample and filter)
itx = upfirdn(it, rrcfilter, nsamp, 1);
qtx = upfirdn(qt, rrcfilter, nsamp, 1);

% Time vector for plots
Drate = 64000; % bit rate
T = 1 / Drate;
Ts = T / nsamp;
time_itx = (0:length(itx)-1) * Ts;
time_qtx = (0:length(qtx)-1) * Ts;

% Plot filtered I and Q signals
figure(2);
plot(time_itx, itx);
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
title('Low pass filtered InPhase Component');
grid on;

figure(3);
plot(time_qtx, qtx);
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
title('Low pass filtered Quadrature Component');
grid on;

% Carrier frequency and modulation
fc = 900e6; % 900 MHz carrier frequency
dd = 2*pi*fc*time_itx';
ddd = 2*pi*fc*time_qtx';

% Offset quadrature by nsamp samples (half symbol delay)
delay = zeros(size(qtx));
delay(nsamp+1:end) = qtx(1:end-nsamp);

% OQPSK modulated signal
mt = cos(dd).*itx + sin(ddd).*delay';

figure(4);
plot(time_itx, mt);
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
title('Differentially encoded OQPSK modulated signal');
grid on;

% Add AWGN noise with SNR 10 dB
snr = 10;
madd = awgn(mt, snr, 'measured');

figure(5);
plot(time_itx, madd);
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
title('OQPSK Signal with AWGN (SNR=10 dB)');
grid on;

% Demodulation: multiply with carrier
cscomp = madd .* cos(dd);
sincomp = madd .* sin(ddd);

% Low pass filtering of demodulated signals
lpfin = upfirdn(cscomp, rrcfilter, 1, nsamp);
lpfqu = upfirdn(sincomp, rrcfilter, 1, nsamp);

% Time vectors for filtered demod signals
tmx = (0:length(lpfin)-1) * Ts * nsamp;
tmy = (0:length(lpfqu)-1) * Ts * nsamp;

figure(6);
plot(tmx, lpfin);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Low pass filtered I-channel after demodulation');
grid on;

figure(7);
plot(tmy, lpfqu);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Low pass filtered Q-channel after demodulation');
grid on;

% Sampling points for decision (accounting for filter delay)
len_lpfin = length(lpfin);
samplePoints = round(span*nsamp + nsamp*(0:length(xbit)-1));
samplePoints(samplePoints > len_lpfin) = [];  % Prevent out-of-range

% Initialize decoded bit arrays
chk1 = zeros(1,length(samplePoints));
chk2 = zeros(1,length(samplePoints));

for i = 1:length(samplePoints)
    chk1(i) = lpfin(samplePoints(i)) > 0;
    chk2(i) = lpfqu(samplePoints(i)) > 0;
end

% Convert logicals to bipolar NRZ (1 -> 1, 0 -> -1)
chk1 = 2*chk1 - 1;
chk2 = 2*chk2 - 1;

disp('I channel bit stream checking');
distortionI = sum((it(1:length(chk1)) - chk1).^2) / length(chk1);
disp(distortionI);

disp('Q channel bit stream checking');
distortionQ = sum((qt(1:length(chk2)) - chk2).^2) / length(chk2);
disp(distortionQ);

% Differential decoding from I and Q
dfd = zeros(1,length(xbit));
dfd(1:2:end) = chk1(1:2:end);
dfd(2:2:end) = chk2(2:2:end);

dfdecod = zeros(1,length(dfd));
dfdecod(dfd==1) = 1;
dfdecod(dfd==-1) = 0;

detected(1) = 1 - ~dfdecod(1);
for i = 2:length(dfd)
    detected(i) = dfdecod(i-1) - ~dfdecod(i);
    if detected(i) == -1
        detected(i) = 1;
    end
end

disp('Distortion between transmitted and received NRZ bit stream');
distortionTotal = sum((xbit - detected).^2) / length(detected);
disp(distortionTotal);

% Plot transmitted and detected bit streams
tmx_bits = (0:length(xbit)-1) / Drate;

figure(8);
subplot(2,1,1);
stairs(tmx_bits, xbit, 'LineWidth', 1.5);
ylim([-0.2 1.2]);
xlabel('Time (sec)');
ylabel('Bit Value');
title('Transmitted Bit Stream');
grid on;

subplot(2,1,2);
stairs(tmx_bits, detected, 'LineWidth', 1.5);
ylim([-0.2 1.2]);
xlabel('Time (sec)');
ylabel('Bit Value');
title('Received Bit Stream');
grid on;
