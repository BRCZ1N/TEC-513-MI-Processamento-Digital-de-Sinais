clear;
clc;
close all;

pkg load signal;

f1 = 2e3;        % Frequencia da senoide a ser amostrada
f2 = 5e3;        % Frequencia da segunda senoide
A1 = 1.0;        % Amplitude da senoide a ser amostrada
A2 = 1.5;        % Amplitude da segunda senoide

fs = 3e3;        % Frequencia de amostragem
Ts = 1/fs;       % Periodo de amostragem
fNyquist = fs/2; % Frequencia de Nyquist
fAnalog = 360e3; % Frequencia da grade computacional
Ta = 1/fAnalog;  % Periodo da grade computacional
tempoTotal = 1;  % Tempo total de simulacao
t = 0:Ta:tempoTotal-Ta; % Grade computacional
nPAnalog = length(t);   % Total de pontos da grade
fCorte = 3e3;    % Frequencia de corte do filtro passa baixas

sinalA = A1*sin(2*pi*f1*t);
sinalB = A2*sin(2*pi*f2*t);
sinalComposto = sinalA + sinalB;

filtroIdeal = zeros(1, nPAnalog);

for i = 1:nPAnalog
    fAtual = (i-1) * (fAnalog / nPAnalog);

    if fAtual <= fCorte || fAtual >= (fAnalog - fCorte)
        filtroIdeal(i) = 1;
    else
        filtroIdeal(i) = 0;
    end
end

sinalCompostoFFT = fft(sinalComposto);
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal;
sinalFiltrado = real(ifft(sinalFiltradoFFT));

% Figura 1: Sinais no Tempo
figure;

subplot(4,1,1);
plot(t,sinalA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.05]);

subplot(4,1,2);
plot(t,sinalB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.05]);

subplot(4,1,3);
plot(t,sinalComposto);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Composto');
xlim([0 0.05]);

subplot(4,1,4);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.05]);

saveas(gcf, sprintf('FlatTop_%gHz_01_Sinais_No_Tempo.png', fs));

kPasso = round(Ts/Ta);

d = Ts/3; % Tempo de abertura do pulso
numPontosPulso = round(d/Ta);
dReal = numPontosPulso * Ta;

amostras = sinalFiltrado(1:kPasso:end); % Obtendo as amostras do sinal filtrado
pulsoFlatTop = [ones(1, numPontosPulso), zeros(1, kPasso - numPontosPulso)]; % Pulso retangular Flat-Top

sinalAmostradoFlatTop = kron(amostras, pulsoFlatTop); % Construindo o sinal Flat-Top
sinalAmostradoFlatTop = sinalAmostradoFlatTop(1:nPAnalog);

% Figura 2: Processo de amostragem no tempo (3 subplots)
figure;

subplot(3,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.05]);

subplot(3,1,2);
stairs(t,sinalAmostradoFlatTop);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Amostrado Flat-Top');
xlim([0 0.05]);

subplot(3,1,3);
plot(t,sinalFiltrado,'--');
hold on;
stairs(t,sinalAmostradoFlatTop);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado e Sinal Flat-Top');
legend('Sinal Filtrado','Flat-Top');
xlim([0 0.05]);
hold off;

saveas(gcf, sprintf('FlatTop_%gHz_02_Processo_De_Amostragem.png', fs));

f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog); % Grade computacional da frequencia

sinalFiltradoFFT_plot = fftshift(fft(sinalFiltrado)) / nPAnalog;
sinalAmostradoFlatTopFFT_plot = fftshift(fft(sinalAmostradoFlatTop)) / nPAnalog;

% Figura 3: Espectros de Magnitude
figure;

subplot(2,1,1);
stem(f, abs(sinalFiltradoFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);
xtickangle(45);

subplot(2,1,2);
stem(f, abs(sinalAmostradoFlatTopFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Flat-Top');
xlim([-15000 15000]);
xticks(-15000:3000:15000);

saveas(gcf, sprintf('FlatTop_%gHz_03_Espectros_De_Magnitude.png', fs));

espectroMag = abs(sinalAmostradoFlatTopFFT_plot);

fSinc = f;
fSinc(fSinc == 0) = 1e-10;

envolventeSinc = abs(sin(pi * fSinc * dReal) ./ (pi * fSinc * dReal)); % Envolvente Sinc
amplitudeMaxima = max(espectroMag);
envolventeSinc = envolventeSinc * amplitudeMaxima;

% Figura 4: Espectro Flat-Top e Envolvente Sinc
figure;

stem(f, espectroMag, '.');
hold on;
plot(f, envolventeSinc, 'k--');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro Flat-Top e Envolvente Sinc');
xlim([-15000 15000]);
xticks(-15000:3000:15000);
ylim([0 amplitudeMaxima * 1.2]);
legend('Espectro Flat-Top', 'Envolvente Sinc', 'Location', 'northeast');
hold off;

saveas(gcf, sprintf('FlatTop_%gHz_04_Envolvente_Sinc.png', fs));

fCorteReconstrucao = fs/2;
filtroReconstrucao = zeros(1, nPAnalog);

for i = 1:nPAnalog
    fAtual = (i-1) * (fAnalog / nPAnalog);

    if fAtual <= fCorteReconstrucao || fAtual >= (fAnalog - fCorteReconstrucao)
        filtroReconstrucao(i) = 1;
    else
        filtroReconstrucao(i) = 0;
    end
end

fatorCompensacao = Ts / dReal;

sinalAmostradoFlatTopFFT = fft(sinalAmostradoFlatTop);
sinalReconstruidoFFT = sinalAmostradoFlatTopFFT .* filtroReconstrucao * fatorCompensacao;
sinalReconstruido = real(ifft(sinalReconstruidoFFT));

filtroReconstrucaoPlot = fftshift(filtroReconstrucao);
sinalReconstruidoFFT_plot = fftshift(sinalReconstruidoFFT) / nPAnalog;

% Figura 5: Reconstrução e Comparação (3 subplots)
figure;

subplot(3,1,1);
stem(f, abs(sinalAmostradoFlatTopFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Sinal Flat-Top');
xlim([-5*fs 5*fs]);
xticks(-5*fs:fs:5*fs);

subplot(3,1,2);
plot(f, filtroReconstrucaoPlot);
grid on;
xlabel('Frequência (Hz)');
ylabel('Ganho');
title('Filtro Ideal de Reconstrução');
xlim([-5*fs 5*fs]);
xticks(-5*fs:fs:5*fs);
ylim([-0.2 1.2]);

subplot(3,1,3);
plot(t, sinalFiltrado);
hold on;
plot(t, sinalReconstruido, '--');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado vs Sinal Reconstruído');
legend('Sinal Filtrado', 'Sinal Reconstruído');
xlim([0 0.05]);
hold off;

saveas(gcf, sprintf('FlatTop_%gHz_05_Reconstrucao_E_Comparacao.png', fs));

% Figura 6: Comparação final detalhada no tempo
figure;

plot(t, sinalFiltrado);
hold on;
plot(t, sinalReconstruido, '--');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Comparação entre Sinal Original e Reconstruído');
legend('Sinal Filtrado', 'Sinal Reconstruído');
xlim([0 0.05]);
hold off;

saveas(gcf, sprintf('FlatTop_%gHz_06_Comparacao_Final.png', fs));
