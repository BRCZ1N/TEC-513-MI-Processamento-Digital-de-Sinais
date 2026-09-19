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
fNyquist = fs/2;
fCorte = 2e3;
ordemFiltro = 100;
% fs >= 2*f1
nPlots = 5;                % Número total de plots

sinalA = A1*sin(2*pi*f1*t);

figure;
subplot(nPlots,3,1);
plot(t, sinalA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.01]);

sinalB = A2*sin(2*pi*f2*t);

##figure;
subplot(nPlots,3,2);
plot(t, sinalB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.01]);

sinalComposto = sinalA + sinalB;

##figure;
subplot(nPlots,3,3);
plot(t,sinalComposto);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal composto');
xlim([0 0.01]);

filtroIdeal = zeros(1, nPAnalog);

for i = 1:nPAnalog
    fAtual = (i - 1) * (fAnalog / nPAnalog);
    % Permite frequências abaixo do corte ou frequências espelhadas
    if fAtual <= fCorte || fAtual >= (fAnalog - fCorte)
        filtroIdeal(i) = 1;
    else
        filtroIdeal(i) = 0;
    end
end

sinalCompostoFFT = fft(sinalComposto);
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal;
sinalFiltrado = real(ifft(sinalFiltradoFFT));

##figure;
subplot(nPlots,3,4);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

kPasso = round(Ts/Ta);
tremImpulsos = zeros(1, nPAnalog);
tremImpulsos(1:kPasso:end) = 1/Ta;

##figure;
subplot(nPlots,3,5);
stem(t, tremImpulsos*Ta);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de impulsos');
xlim([0 0.01]);

sinalAmostrado = sinalFiltrado .* (tremImpulsos*Ta);

##figure;
subplot(nPlots,3,6);
stem(t, sinalAmostrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal amostrado');
xlim([0 0.01]);

f = (-nPAnalog/2 : nPAnalog/2 - 1)*(fAnalog/nPAnalog);

sinalFiltradoFFT = fftshift(fft(sinalFiltrado))/nPAnalog;

##figure;
subplot(nPlots,3,7);
stem(f, abs(sinalFiltradoFFT));
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Sinal Filtrado');
xlim([-5000 5000]);
xticks(-5000 : 1000 : 5000);

tremImpulsosFFT = fftshift(fft(tremImpulsos))/nPAnalog;

##figure;
subplot(nPlots,3,8);
stem(f, abs(tremImpulsosFFT));
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Trem de Impulsos');
xlim([-20000 20000]);
xticks(-200000 : 10000 : 200000);

sinalAmostradoFFT = fftshift(fft(sinalAmostrado))/nPAnalog;

##figure;
subplot(nPlots,3,9);
stem(f, imag(sinalAmostradoFFT));
grid on;
xlabel('Frequência (Hz)');
title('Espectro do Sinal Amostrado');
xlim([-35000 35000]);
xticks(-35000 : 10000 : 35000);

##figure;
subplot(nPlots,3,10);
stem(f, abs(sinalAmostradoFFT));
grid on;
xlabel('Frequência (Hz)');
title('Espectro de Magnitude do Sinal Amostrado');
xlim([-35000 35000]);
xticks(-35000 : 10000 : 35000);










