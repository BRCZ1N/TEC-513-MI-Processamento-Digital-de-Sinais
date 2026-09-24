clear;
clc;
close all;

%% ===================== PARAMETROS =====================
f1 = 2e3;        % Frequencia do Sinal A (util)
f2 = 5e3;        % Frequencia do Sinal B (fora da banda permitida)
A1 = 1.0;        % Amplitude do Sinal A
A2 = 1.5;        % Amplitude do Sinal B

fs = 3e3;        % Frequencia de amostragem (use 8e3 para o caso sem aliasing)
Ts = 1/fs;       % Periodo de amostragem
fAnalog = 360e3; % Frequencia da grade computacional (emula tempo continuo)
Ta = 1/fAnalog;  % Passo da grade computacional

tempoTotal = 1;  % Tempo total (s). Resolucao espectral = 1/tempoTotal (1 Hz)
fatorZP = 8;     % Fator de zero-padding, usado nas Figuras 6 e 7 (zoom nas amostras)
larguraZoom = 5; % Meia-largura (Hz) da janela de zoom em torno de cada linha
freqsZoom = [-2000 -1000 1000 2000]; % Frequencias (Hz) onde dar zoom
parteZoom = @abs;    % Pode ser @abs, @imag ou @real

fCorte = 3e3;    % Corte do filtro anti-aliasing (use fs/2 para o caso correto)
fCorteReconstrucao = fs/2; % Corte do filtro de reconstrucao

duty = 1/3;      % Duty cycle do pulso

% Faixa de frequencias exibida nos espectros
limEsp1 = 10000;      % Figura 3 (sinal filtrado) e Figura 4 (superior)
limEsp2 = 2.5*fs;     % Demais espectros

% Requisitos para evitar vazamento espectral
if abs(fAnalog/fs - round(fAnalog/fs)) > 1e-9
    error('fAnalog/fs deve ser inteiro.');
end
if abs(fs*tempoTotal - round(fs*tempoTotal)) > 1e-9
    error('fs*tempoTotal deve ser inteiro.');
end

%% ===================== GRADE E SINAIS =====================
t = 0:Ta:tempoTotal-Ta;   % Grade computacional
nPAnalog = length(t);     % Total de pontos da grade

sinalA = A1*sin(2*pi*f1*t);   % Sinal a ser amostrado
sinalB = A2*sin(2*pi*f2*t);   % Sinal fora da banda
sinalComposto = sinalA + sinalB;

%% ===================== FILTRO ANTI-ALIASING IDEAL =====================
fVet = (0:nPAnalog-1) * (fAnalog/nPAnalog);   % Eixo de frequencia da FFT (0 a fAnalog)
filtroIdeal = double(fVet <= fCorte | fVet >= (fAnalog - fCorte));

sinalCompostoFFT = fft(sinalComposto);
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal;
sinalFiltrado = real(ifft(sinalFiltradoFFT));

%% ===================== AMOSTRAGEM NATURAL =====================
kPasso = round(fAnalog / fs);
kPulso = round(duty * kPasso);

tremPulsos = zeros(1,nPAnalog);
tremPulsos(mod(0:nPAnalog-1, kPasso) < kPulso) = 1;
dutyReal = mean(tremPulsos);   % c0 medido do trem de pulsos

sinalAmostrado = sinalFiltrado .* tremPulsos;   % x(t) * p(t)

%% ===================== FIGURA 1: SINAIS NO TEMPO =====================
figure;

subplot(4,1,1);
plot(t,sinalA,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinal A'); xlim([0 0.01]);

subplot(4,1,2);
plot(t,sinalB,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinal B'); xlim([0 0.01]);

subplot(4,1,3);
plot(t,sinalComposto,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinal Composto'); xlim([0 0.01]);

subplot(4,1,4);
plot(t,sinalFiltrado,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinal Filtrado'); xlim([0 0.01]);

%% ===================== FIGURA 2: PROCESSO DE AMOSTRAGEM =====================
figure;

subplot(4,1,1);
plot(t,sinalFiltrado,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinal Filtrado'); xlim([0 0.01]);

subplot(4,1,2);
plot(t,tremPulsos,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Trem de Pulsos'); xlim([0 0.01]); ylim([-0.1 1.1]);

subplot(4,1,3);
plot(t,sinalAmostrado,'LineWidth',1.2);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Amostragem (Sinal Amostrado)'); xlim([0 0.01]);

sinalSegmentos = sinalFiltrado;
sinalSegmentos(tremPulsos == 0) = NaN;   % Oculta os trechos fora do pulso

subplot(4,1,4);
plot(t, sinalFiltrado, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;
plot(t, sinalSegmentos, 'r', 'LineWidth', 1.5);
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Segmentos da Amostragem Natural sobre a Senoide');
xlim([0 0.01]);
legend('Sinal Continuo', 'Trechos Amostrados', 'Location', 'northeast');
hold off;

%% ===================== ESPECTROS (VISAO GERAL) =====================
Nfft = nPAnalog;   % Visao geral: sem zero-padding (linhas exatas, 1 Hz de resolucao)
f = (-Nfft/2:Nfft/2-1) * (fAnalog/Nfft);   % Eixo de frequencia centrado

% Normalizacao por nPAnalog mantem a amplitude das linhas espectrais
sinalFiltradoFFT_plot  = fftshift(fft(sinalFiltrado,  Nfft))/nPAnalog;
tremPulsosFFT_plot     = fftshift(fft(tremPulsos,     Nfft))/nPAnalog;
sinalAmostradoFFT_plot = fftshift(fft(sinalAmostrado, Nfft))/nPAnalog;

% Mascaras para plotar apenas a faixa de interesse (evita stem de milhoes de pontos)
idx1 = (f >= -limEsp1) & (f <= limEsp1);
idx2 = (f >= -limEsp2) & (f <= limEsp2);

xt1 = -limEsp1:2000:limEsp1;

% Plota apenas as linhas com amplitude acima de tol (evita o 'borrao' dos zeros)
tol = 1e-4;
stemLinhas = @(X, idx, parte) stem(f(idx & abs(X) > tol), parte(X(idx & abs(X) > tol)), 'filled', 'MarkerSize', 5);
passoTick = 1000;
limTick = floor(limEsp2/passoTick)*passoTick;
xt2 = -limTick:passoTick:limTick;   % ticks em multiplos de 1 kHz

%% ===================== FIGURA 3: MAGNITUDE =====================
figure;

subplot(3,1,1);
stemLinhas(sinalFiltradoFFT_plot, idx1, @abs);
grid on; xlabel('Frequencia (Hz)'); ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Filtrado');
xlim([-limEsp1 limEsp1]); set(gca,'XTick',xt1);

subplot(3,1,2);
stemLinhas(tremPulsosFFT_plot, idx2, @abs);
grid on; xlabel('Frequencia (Hz)'); ylabel('Magnitude');
title('Espectro de Magnitude do Trem de Pulsos');
xlim([-limEsp2 limEsp2]); set(gca,'XTick',xt2);

subplot(3,1,3);
stemLinhas(sinalAmostradoFFT_plot, idx2, @abs);
grid on; xlabel('Frequencia (Hz)'); ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Amostrado');
xlim([-limEsp2 limEsp2]); set(gca,'XTick',xt2);

%% ===================== FIGURA 4: PARTE IMAGINARIA =====================
figure;

subplot(2,1,1);
stemLinhas(sinalFiltradoFFT_plot, idx1, @imag);
grid on; xlabel('Frequencia (Hz)'); ylabel('Amplitude');
title('Espectro do Sinal Filtrado (parte imaginaria)');
xlim([-limEsp1 limEsp1]); set(gca,'XTick',xt1);

subplot(2,1,2);
stemLinhas(sinalAmostradoFFT_plot, idx2, @imag);
grid on; xlabel('Frequencia (Hz)'); ylabel('Amplitude');
title('Espectro do Sinal Amostrado (parte imaginaria)');
xlim([-limEsp2 limEsp2]); set(gca,'XTick',xt2);

%% ===================== RECONSTRUCAO =====================
filtroReconstrucao = double(fVet <= fCorteReconstrucao | fVet >= (fAnalog - fCorteReconstrucao));

% Compensacao da atenuacao do pulso dividindo pelo duty cycle real (c0)
sinalReconstruidoFFT = fft(sinalAmostrado) .* filtroReconstrucao / dutyReal;
sinalReconstruido = real(ifft(sinalReconstruidoFFT));

% Espectro do sinal reconstruido na mesma grade (com zero-padding)
sinalReconstruidoFFT_plot = fftshift(fft(sinalReconstruido, Nfft))/nPAnalog;

% Filtro de reconstrucao desenhado direto no eixo f centrado
filtroReconstrucaoPlot = double(abs(f) <= fCorteReconstrucao);

%% ===================== FIGURA 5: RECONSTRUCAO =====================
figure;

subplot(4,1,1);
stemLinhas(sinalAmostradoFFT_plot, idx2, @imag);
grid on; xlabel('Frequencia (Hz)'); ylabel('Amplitude');
title('Espectro do Sinal Amostrado (parte imaginaria)');
xlim([-limEsp2 limEsp2]); set(gca,'XTick',xt2);

subplot(4,1,2);
plot(f(idx2), filtroReconstrucaoPlot(idx2), 'LineWidth', 2, 'Color', 'r');
grid on; xlabel('Frequencia (Hz)'); ylabel('Ganho');
title('Filtro Ideal de Reconstrucao');
xlim([-limEsp2 limEsp2]); set(gca,'XTick',xt2); ylim([-0.2 1.2]);

subplot(4,1,3);
stemLinhas(sinalReconstruidoFFT_plot, idx2, @imag);
grid on; xlabel('Frequencia (Hz)'); ylabel('Amplitude');
title('Espectro do Sinal Reconstruido (parte imaginaria)');
xlim([-limEsp2 limEsp2]); set(gca,'XTick',xt2);

subplot(4,1,4);
plot(t, sinalFiltrado, 'LineWidth', 2, 'DisplayName', 'Sinal Filtrado');
hold on;
plot(t, sinalReconstruido, '--', 'LineWidth', 2, 'DisplayName', 'Sinal Reconstruido');
grid on; xlabel('Tempo (s)'); ylabel('Amplitude');
title('Comparacao: Sinal Filtrado vs Sinal Reconstruido');
xlim([0 0.01]);
legend('Location', 'northeast');
hold off;

%% ===================== ESPECTROS DENSOS (ZERO-PADDING) PARA ZOOM =====================
NfftZ = fatorZP * nPAnalog;
fZ = (-NfftZ/2:NfftZ/2-1) * (fAnalog/NfftZ);   % Resolucao = 1/(fatorZP*tempoTotal) Hz
XamostZ = fftshift(fft(sinalAmostrado,    NfftZ))/nPAnalog;
XreconZ = fftshift(fft(sinalReconstruido, NfftZ))/nPAnalog;

%% ===================== FIGURA 6: ZOOM NO ESPECTRO DO SINAL AMOSTRADO =====================
figure;
nz = length(freqsZoom);
for k = 1:nz
    subplot(ceil(nz/2), 2, k);
    iz = abs(fZ - freqsZoom(k)) <= larguraZoom;
    stem(fZ(iz), parteZoom(XamostZ(iz)), 'filled', 'MarkerSize', 3);
    grid on; xlabel('Frequencia (Hz)'); ylabel('Amplitude');
    title(sprintf('Sinal Amostrado - zoom em %g Hz', freqsZoom(k)));
    xlim([freqsZoom(k)-larguraZoom freqsZoom(k)+larguraZoom]);
end

%% ===================== FIGURA 7: ZOOM NO ESPECTRO DO SINAL RECONSTRUIDO =====================
figure;
for k = 1:nz
    subplot(ceil(nz/2), 2, k);
    iz = abs(fZ - freqsZoom(k)) <= larguraZoom;
    stem(fZ(iz), parteZoom(XreconZ(iz)), 'filled', 'MarkerSize', 3);
    grid on; xlabel('Frequencia (Hz)'); ylabel('Amplitude');
    title(sprintf('Sinal Reconstruido - zoom em %g Hz', freqsZoom(k)));
    xlim([freqsZoom(k)-larguraZoom freqsZoom(k)+larguraZoom]);
end
