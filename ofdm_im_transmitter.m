%Omer Lutfu Karakelle

clc
clear
close all
%%
[data, data2] = text2bin('Lorem ipsum dolor sit amet, consectetur adipiscing elit. Ut porta augue non magna aliquam.')
%% tx options
transmit = comm.SDRuTransmitter(...
              'Platform','B200', ...
              'SerialNum','31FD9D5', ...
              'ChannelMapping',1, ...
              'Gain', 50, ...
              'CenterFrequency', 1200e6, ...
              'InterpolationFactor', 64);
disp(transmit)
%% load data and preamble

load("preamble.mat")
load("preamble_cp.mat")

%% Transmitter
N = 2; %QPSK
FFT_N = 64; %FFT length
sym_len = 90;
CP_size = 16; %cycle prefix length
PL = 8; %payload length
pilot = 1;
%data = randi([0 1], sym_len*PL, 1)'; %random data sequence
tx = [];
for i = 1:PL

    sym = [];
    for j = 1:6:sym_len
        m = mapper(data((i - 1)*sym_len +j:(i - 1)*sym_len + j + 5));
        sym = [sym m];
        if ismember(length(sym), [12 21 42 55])
            sym = [sym pilot];
        end
    end
    time = ifft(sym, FFT_N);
    tx = [tx time(FFT_N - CP_size + 1:end) time];
end

%preamble = 1 - 2*randi([0 1], FFT_N/2, 1)' + 1i - 1i*2*randi([0 1], FFT_N/2, 1)';
%preamble_time = ifft(preamble, FFT_N/2);
%preamble_time = [preamble_time preamble_time];
%preamble_cp = [preamble_time(FFT_N - CP_size + 1:end) preamble_time];
tx = [preamble_cp/2 tx];
%preamble = fft(preamble_time, FFT_N);

%L = length(preamble)/2;
preamble_2 = 1 - 2*randi([0 1], FFT_N, 1)' + 1i - 1i*2*randi([0 1], FFT_N, 1)';
preamble_2([2:4:end]) = 0 + 1i*0;
preamble_2([3:4:end]) = 0 + 1i*0;
preamble_2([4:4:end]) = 0 + 1i*0;
preamble_time_2 = ifft(preamble_2, FFT_N);
preamble_cp_2 = [preamble_time_2(FFT_N - CP_size + 1:end) preamble_time_2];
tx = [preamble_cp_2*sqrt(2) tx];


%% transmit
    while (true)
      transmit(tx');
    end
    release(transmit);

%%
function m = mapper(sub_data)
    sig = 1 - 2*sub_data(3:6);
    s1 = sig(2) + 1i*sig(1);
    s2 = sig(4) + 1i*sig(3);
    if sub_data(1:2) == [0 0]
        m = [s1 s2 0 0];
    elseif sub_data(1:2) == [0 1]
        m = [0 s1 s2 0];
    elseif sub_data(1:2) == [1 0]
        m = [0 0 s1 s2];
    elseif sub_data(1:2) == [1 1]
        m = [s1 0 0 s2];
    end
end
%%
function [binV, binS] = text2bin(text)
    binS = dec2bin(text,8);
    binS = binS';
    binS = binS(:)';
    binV = binS-48;
end