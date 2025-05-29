%% clear
clc
clear all
close all

%% receiver

receive = comm.SDRuReceiver(...
              'Platform','B210', ...
              'SerialNum','327AB44', ...
              'ChannelMapping',1, ...
              'Gain', 40, ...
              'CenterFrequency', 1200e6, ...
              'DecimationFactor', 64, ...
              'SamplesPerFrame', 15000);
disp(receive)

rx = [];
try
    for i = 1:20
        [receivedSignal, len] = step(receive);
        if len > 0 && any(receivedSignal(:) ~=0)
            rx = [rx; receivedSignal];
        end
    end
catch exception
    disp('error: ');
    disp(exception.message);
end
release(receive);
%% plot signal
rx = double(rx);
rx = rx'./max(rx);
plot(abs(rx))
%% load data and preamble

load("data.mat")
load("preamble.mat")
load("preamble_cp.mat")
%%

PL = 8; %payload length
N = 2; %bits per sample (qpsk)
FFT_N = 64; %fft  size
CP_size = FFT_N/4; %cyclic prefix size = 16
sym_len = FFT_N*N; % symbol len = 128
L = FFT_N / 2; %length of a half preamble
frame_len = (PL + 2)*(FFT_N + CP_size);
xq = [1:1:64];
pilot_index = [11 20 45 54];
null_index = [1 2 29:36 63 64]; 
%% Coarse Time Syncronization

%Schmidl and Cox metric is used for coarse time sycronization

P = zeros(length(rx), 0);
R = zeros(length(rx), 0);
for d = 1:(length(rx) - 2*L + 1)
    P(d) = 0;
    R(d) = 0;
    for m = 0:(L-1)
        P(d) = P(d) + conj(rx(d + m)).*rx(d + m + L);
        R(d) = R(d) + abs(rx(d + m + L)).^2;
    end
end

M = abs(P).^2./R.^2;
figure
plot(abs(M))
hold on
plot(abs(rx))
hold off

rectpulse = ones(CP_size, 1)';
a_pick = conv(rectpulse, M);
figure
plot(a_pick)

[val_max_cf, max_cf] = max(a_pick); %maximum value of Schmidl and Cox metric
if (a_pick(max_cf - 80) >= 0.5*a_pick(max_cf))
    max_cf = max_cf - 80;
end
max_cf
%% Coarse Frequency Syncronization
% Schmidl and Cox metric is used for coarse frequency syncronization
df = angle(P(max_cf))/pi;
rx = rx.*exp(-1i*2.*pi.*df.*[1:length(rx)]/64);
%% Fine Time Syncronization
% Fine time syncronization is done using cross-correlation of received
% signal and known preamble 
[c, lags] = xcorr(rx, preamble_cp(1:40));
c = abs(c).^2;
padding = length(rx) - 40 + 1;
c_pad = c(padding:end);
figure
plot(c_pad)
[val_max, max_f] = max(c_pad);
th = val_max * 0.65; %threshold value for peak detection
if (0.75*c_pad(max_f) <= c_pad(max_f + 32))
    max_f = max_f + 32;
end
max_f

%% Frame Extraction
frames = [];
for i = 161:length(rx) - frame_len
    if c_pad(i) >= th && max(c_pad(i + 30:i + 35)) < th && c_pad(i) >= max(c_pad(i - 3:i + 3))
        peak = c_pad(i);
        frame = rx(i - FFT_N + 1:i + frame_len - 2*FFT_N - 2*CP_size);
        frames = [frames; frame];
    end
end

%% Fine Frequency Syncronization
for i = 1:height(frames)
    pre_detected = frames(i, 1:FFT_N);
    sum = 0;
    for j = 1:32
        sum = sum + conj(pre_detected(j)).*pre_detected(32 + j);
    end
    df_2 = angle(sum)/pi;
    frames(i, :) = frames(i, :).*exp(-1i*2.*pi.*df_2.*[1:length(frames(i, :))]/64);
end

%% Channel Estimation
H = zeros(height(frames) , FFT_N);
for i = 1:height(frames)
    pre_detected = frames(i, 1:FFT_N);
    pre_y = fft(pre_detected, 64);
    H(i, :) = pre_y./preamble;
    H_q(i, :) = H(i, 1:2:end);
    H(i, :) = interp1([1:2:length(H(i, :))], H_q(i, :), xq, 'spline','extrap'); 
    %H(i, :) = smoothdata(H(i, :));
    %H(i, 2:2:end - 2) = (H(i, 1:2:end - 2) + H(i, 3:2:end)) / 2;
    H(i, end) = (H(i, end - 1) + H(i, 1))/2;
   % H(i, :) = smoothdata(H(i, :));

end
scatterplot(H(2, :))
figure
plot(abs(H(2, :)))
figure
plot(angle(H(2, :)))
%% Demodulation and Channel Equalization
sym_all = [];
ber = 0;
active_carriers = [];
for i = 1:height(frames)
    frame = frames(i, FFT_N + 1:end);
    sym = [];
    rx_data = [];
    for j = 1:PL
        sym_time = frame(j * CP_size + (j - 1)* FFT_N + 1: j*FFT_N + j*CP_size);
        sym_freq = fft(sym_time, FFT_N);
        sym_freq = sym_freq.*conj(H(i,:))./abs(H(i,:)).^2; %channel equalization
        %sym_freq = sym_freq.*exp(1i.*angle(conj(H(i,:)))); %channel equalization

        %pilot_equ = (sym_freq(11)*sym_freq(33)*sym_freq(55))^(1/3)/abs(sym_freq(11)*sym_freq(33)*sym_freq(55));
        pilots = interp1(pilot_index, sym_freq(pilot_index), xq, 'linear','extrap');
        pilots = smoothdata(pilots);
        %pilots = sym_freq(13);
        sym_freq = sym_freq./pilots;

        %pilots = sym_freq([13 22 43 56]);
        sym_freq([pilot_index null_index]) = [];
        for k = 1:4:48
            sub_block = sym_freq(k:k+3);
            [sub_data, active_carrier] = demapper(sub_block);
            rx_data = [rx_data sub_data];
            active_carriers = [active_carriers active_carrier];
        end
        sym = [sym sym_freq];
    end
    %scatterplot(sym)
    sym_all = [sym_all sym];
    %rx_data = zeros(1, PL*FFT_N*2);
    %rx_data(2:2:end) = real(sym) < 0;
    %rx_data(1:2:end) = imag(sym) < 0;
    result = bin2text(rx_data)
    ber = ber + biterr(data, rx_data) / length(data);
end
ber = ber / height(frames)
scatterplot(sym_all)
grid on
scatterplot(active_carriers)
grid on
%%
function [d, p] = demapper(sub_block)
    d = zeros(1, 6);
    [p, x] = maxk(sub_block, 2);
    x = sort(x);
    i1 = x(1);
    i2 = x(2);
    d(6) = real(sub_block(i2)) < 0;
    d(5) = imag(sub_block(i2)) < 0;
    d(4) = real(sub_block(i1)) < 0;
    d(3) = imag(sub_block(i1)) < 0;
    if [i1 i2] == [1 2]
        d(1:2) = [0 0];
    elseif [i1 i2] == [2 3]
        d(1:2) = [0 1];
    elseif [i1 i2] == [3 4]
        d(1:2) = [1 0];
    elseif [i1 i2] == [1 4]
        d(1:2) = [1 1];
    end
end

%%
function text = bin2text(binVS)
    btxt  = reshape(binVS,[8, length(binVS)/8])';
    if length(class(btxt))== 6
        text  = char(bin2dec(char(btxt+48)))';
    else
        text  = char(bin2dec(btxt))';
    end
end