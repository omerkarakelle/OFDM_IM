%Omer Lutfu Karakelle
%Istanbul Technical University
clc
clear
close all

SNR_dB = [10:3:40];
df = 0.05; %frequency shift
%% Transmitter
N = 2; %QPSK
FFT_N = 64; %FFT length
CP_size = 16; %cycle prefix length
PL = 1000; %payload length
pilot = 1;
pilot_index = [];
null_index = []; %subcarrier indexes for guard band
sym_len = FFT_N /4*6;
sym_len_ofdm = FFT_N*2;
data = randi([0 1], sym_len*PL, 1)'; %random data sequence
data_ofdm = randi([0 1], sym_len_ofdm*PL, 1)'; %random data sequence

tx = [];
for i = 1:PL
    sym = [];
    for j = 1:6:sym_len
        m = mapper(data((i - 1)*sym_len +j:(i - 1)*sym_len + j + 5));
        sym = [sym m];
    end
    time = ifft(sym, FFT_N);500
    tx = [tx time(FFT_N - CP_size + 1:end) time];
end
%% OFDM Transmitter
tx_ofdm = [];
for i = 1:PL
    sym_ofdm = [];
    for j = 1:2:sym_len_ofdm
        m = mapper_ofdm(data_ofdm((i - 1)*sym_len_ofdm +j:(i - 1)*sym_len_ofdm + j + 1));
        sym_ofdm = [sym_ofdm m];
    end
    time = ifft(sym_ofdm, FFT_N);
    tx_ofdm = [tx_ofdm time(FFT_N - CP_size + 1:end) time];
end
%% Receiver OFDM-IM
sym_all = [];
active_carriers = [];
ber_ofdm_im = [];
Es = sum(tx.*conj(tx))/PL;
for snr = SNR_dB
    ber = 0;
    SNR_linear = power(10, snr/10);
    N0 = Es/SNR_linear;
    n = sqrt(N0)*(randn(1, length(tx)) + 1i*randn(1, length(tx)));
    rx = tx + n;
    
    for i = 1:1
        frame = rx;
        sym = [];
        rx_data = [];
        for j = 1:PL
            sym_time = frame(j * CP_size + (j - 1)* FFT_N + 1: j*FFT_N + j*CP_size);
            sym_time = sym_time.*exp(1i*2.*pi.*df.*[1:length(sym_time)]/64);
            sym_freq = fft(sym_time, FFT_N);
            for k = 1:4:64
                sub_block = sym_freq(k:k+3);
                [sub_data, active_carrier] = demapper(sub_block);
                rx_data = [rx_data sub_data];
                active_carriers = [active_carriers active_carrier];
            end
            sym = [sym sym_freq];
        end
        
        sym_all = [sym_all sym];
        ber = ber + biterr(data, rx_data) / length(data);
    end
    ber_ofdm_im = [ber_ofdm_im ber];
    ber
end

%% Receiver OFDM
sym_all = [];
active_carriers = [];
ber_ofdm = [];
Es = sum(tx_ofdm.*conj(tx_ofdm))/PL;
for snr = SNR_dB
    ber = 0;
    SNR_linear = power(10, snr/10);
    N0 = Es/SNR_linear;
    n = sqrt(N0)*(randn(1, length(tx_ofdm)) + 1i*randn(1, length(tx_ofdm)));
    rx_ofdm = tx_ofdm + n;
    
    for i = 1:1
        frame = rx_ofdm;
        sym = [];
        rx_data = [];
        for j = 1:PL
            sym_time = frame(j * CP_size + (j - 1)* FFT_N + 1: j*FFT_N + j*CP_size);
            sym_time = sym_time.*exp(1i*2.*pi.*df.*[1:length(sym_time)]/64);
            sym_freq = fft(sym_time, FFT_N);
            for k = 1:1:64
            sub_block = sym_freq(k);
            [sub_data, active_carrier] = demapper_ofdm(sub_block);
            rx_data = [rx_data sub_data];
            active_carriers = [active_carriers active_carrier];
            end
            sym = [sym sym_freq];
        end
        
        sym_all = [sym_all sym];
        ber = ber + biterr(data_ofdm, rx_data) / length(data_ofdm);
    end
    ber_ofdm = [ber_ofdm ber];
    ber
end
%% Figures

figure
semilogy(SNR_dB, ber_ofdm_im, 'LineWidth', 2)
xlabel('E_s/N_0')
ylabel('BER')
grid on
hold on
semilogy(SNR_dB, ber_ofdm, 'LineWidth', 2)
legend('OFDM-IM', 'OFDM')

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

function m = mapper_ofdm(sub_data)
    sig = 1 - 2*sub_data;
    s1 = sig(2) + 1i*sig(1);
    m = [s1/sqrt(2)];
end

function [d, p] = demapper(sub_block)
    lookup = [1 2; 2 3; 3 4; 1 4];
    powers = abs(sub_block).^2;
    pair_scores = powers(lookup(:,1)) + powers(lookup(:,2));

    d = zeros(1, 6);
    [p, x] = max(pair_scores); %maximum likelihood decoding
    x = lookup(x, :);
    i1 = x(1);
    i2 = x(2);
    p = sub_block([i1 i2]);
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

function [d, p] = demapper_ofdm(sub_block)
    d = zeros(1, 2);
    p = sub_block;
    d(2) = real(sub_block) < 0;
    d(1) = imag(sub_block) < 0;
end