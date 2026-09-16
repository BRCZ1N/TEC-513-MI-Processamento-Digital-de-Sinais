clear;
clc;
close all;

pkg load signal;


% Frequências do sinal analógico de entrada
f1 = 1e3;                   % Senoide desejada (Banda do sinal: 1KHz)
f2 = 5e3;                   % Segunda senoide (Banda do sinal: 5KHz)
% Amplitudes das senoides
A1 = 1.0;
A2 = 1.2;

fAnalog = 1e5;             % Frequência analógica simulada (100 kHz)
Ta = 1/fAnalog;            % Período de amostragem analógico
tempoTotal = 1;            % Tempo total de simulação (1s)
t = 0:Ta:tempoTotal-Ta;    % Vetor de tempo "contínuo"
nPAnalog = length(t);      % Número total de pontos analógicos
fs = 10e3;                 % Frequencia de amostragem(Trem de impulsos)
Ts = 1/fs;                 % Periodo de amostragem(Trem de impulsos)
nPlots = 8;                % Número total de plots
ordemFiltro = 100;
% fs >= 2*f1
fNyquist = fs/2;
fCorte = 1e3;
fs = 10e3;
Ts = 1/fs;
kPasso = round(Ts/Ta);

sinalA = A1*sin(2*pi*f1*t);

%figure;
subplot(nPlots,1,1);
plot(t, sinalA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.01]);

sinalB = A2*sin(2*pi*f2*t);

%figure;
subplot(nPlots,1,2);
plot(t, sinalB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.01]);

sinalComposto = sinalA + sinalB;

%figure;
subplot(nPlots,1,3);
plot(t,sinalComposto);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal composto');
xlim([0 0.01]);

filtro = fir1(ordemFiltro, fCorte/fNyquist);
sinalFiltrado = filter(filtro, 1, sinalComposto);

tremImpulsos = zeros(1, nPAnalog);
tremImpulsos(1:kPasso:end) = 1;

subplot(nPlots,1,4);
stem(t, tremImpulsos);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de impulsos');
xlim([0 0.01]);

sinalAmostrado = sinalFiltrado .* tremImpulsos;

subplot(nPlots,1,5);
stem(t, sinalAmostrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal amostrado');
xlim([0 0.01]);

f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog)

sinalCompostoFTT = abs(fftshift(fft(sinalComposto)))/nPAnalog;

subplot(nPlots,1,6);
stem(f, sinalCompostoFTT);
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Sinal Composto');
xlim([-10000 10000]);
xticks(-10000 : 1000 : 10000);

tremImpulsosFTT = abs(fftshift(fft(tremImpulsos)))/nPAnalog;

subplot(nPlots,1,7);
stem(f, tremImpulsosFTT);
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Trem de Impulsos');
xlim([-20000 20000]);
xticks(-200000 : 10000 : 200000);

espectroSinalAmostrado = conv(sinalCompostoFTT,tremImpulsosFTT);

nP_conv = length(espectroSinalAmostrado);
f_conv = linspace(-fAnalog, fAnalog, nP_conv);

subplot(nPlots,1,8);
stem(f_conv, espectroSinalAmostrado);
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Convolução dos sinais');
xlim([-40000 40000]);
xticks(-40000 : 10000 : 40000);








