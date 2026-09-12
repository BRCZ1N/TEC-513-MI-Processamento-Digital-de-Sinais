clear;
clc;
close all;

pkg load signal;

f_analog = 1e5;            % Frequência analógica simulada (100 kHz)
Ta = 1/f_analog;           % Período de amostragem analógico
T_total = 0.01;             % Tempo total de simulação (1s)0.1
t = 0:Ta:T_total-Ta;        % Vetor de tempo "contínuo"
N_analog = length(t);       % Número total de pontos analógicos

% Frequências do sinal analógico de entrada
f1 = 1e3;                   % Senoide desejada (Banda do sinal: 1KHz)
f2 = 10e3;                  % Segunda senoide (Banda do sinal: 10KHz)
% Amplitudes das senoides
A1 = 1.0;
A2 = 1.2;

senoideA = A1*sin(2*pi*f1*t);

%figure;
subplot(3,1,1);
plot(t, senoideA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal analógico A');

senoideB = A2*sin(2*pi*f2*t);

%figure;
subplot(3,1,2);
plot(t, senoideB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal analógico B');

senoideComposta = senoideA + senoideB;

%figure;
subplot(3,1,3);
plot(t, senoideComposta);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal analógico composto');

fs = 150;                   % Frequência de amostragem ajustável (Hz)
Ts = 1/fs;                  % Período de amostragem digital (s)
k_passo = round(Ts/Ta);

trem_impulsos = zeros(1, N_analog);
k_passo = round(Ts/Ta);
trem_impulsos(1:k_passo:end) = 1;




