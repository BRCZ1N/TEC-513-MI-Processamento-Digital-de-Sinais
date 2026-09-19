clear;
clc;
close all;
pkg load signal; % Remova ou comente esta linha se estiver usando MATLAB

%% 1. PARÂMETROS DO SINAL E DA SIMULAÇÃO
f1 = 1e3;                   % Senoide 1 (1 kHz) - Sinal Desejado
f2 = 5e3;                   % Senoide 2 (5 kHz) - Sinal Indesejado (Ruído/Aliasing)
A1 = 1.0;                   % Amplitude da senoide 1
A2 = 1.2;                   % Amplitude da senoide 2

fAnalog = 2e5;              % Frequência analógica simulada (200 kHz)
Ta = 1/fAnalog;             % Período de amostragem analógico
tempoTotal = 1;             % Tempo total de simulação (1s)
t = 0:Ta:tempoTotal-Ta;     % Vetor de tempo "contínuo"
nPAnalog = length(t);       % Número total de pontos analógicos

fs = 20e3;                  % Frequência de amostragem ajustada (20 kHz)
Ts = 1/fs;                  % Período de amostragem

% Sinais Individuais e Composto
sinalA = A1*sin(2*pi*f1*t);
sinalB = A2*sin(2*pi*f2*t);
sinalComposto = sinalA + sinalB;

%% 2. FILTRAGEM PRÉVIA (EXTRAÇÃO DO SINAL A ANTES DA AMOSTRAGEM)
fCorte = 3e3;               % Frequência de corte (3 kHz)

% FFT do sinal composto
X_comp = fft(sinalComposto);
filter_ideal = zeros(1, nPAnalog);

% Criar filtro ideal
for i = 1:nPAnalog
    f_atual = (i - 1) * (fAnalog / nPAnalog);
    % Permite frequências abaixo do corte ou frequências espelhadas
    if f_atual <= fCorte || f_atual >= (fAnalog - fCorte)
        filter_ideal(i) = 1;
    else
        filter_ideal(i) = 0;
    end
end

% Aplica o filtro e retorna ao tempo (O resultado será o Sinal A isolado)
X_comp_filtered = X_comp .* filter_ideal;
sinalPreFiltrado = real(ifft(X_comp_filtered));

%% 3. GERAÇÃO DO TREM DE PULSOS
kPasso = round(Ts/Ta);
d = (1/3) * Ts;             % Largura do pulso (duty cycle de 1/3)
numPontosPulso = round(d/Ta);

%% 4. AMOSTRAGEM FLAT-TOP (Feita sobre o sinal JÁ filtrado)
sinalAmostrado = zeros(1, nPAnalog);
for i = 0 : floor(nPAnalog/kPasso) - 1
    idx = i*kPasso + 1;
    if idx <= nPAnalog
        val = sinalPreFiltrado(idx); % Faz o 'Sample'
        idx_fim = min(idx + numPontosPulso - 1, nPAnalog);
        sinalAmostrado(idx:idx_fim) = val; % Aplica o 'Hold' (topo plano)
    end
end

%% 5. RECONSTRUÇÃO (Para demonstrar a eficácia da amostragem)
% Usa-se o mesmo filtro ideal para recuperar o sinal amostrado
X_amostrado = fft(sinalAmostrado);
X_reconstruido = X_amostrado .* filter_ideal;
sinalReconstruido = real(ifft(X_reconstruido));

% Compensação de amplitude pelo Duty Cycle REAL simulado
d_real = numPontosPulso * Ta;
ganhoReconstrucao = Ts / d_real;
sinalReconstruidoCompensado = sinalReconstruido * ganhoReconstrucao;

%% 6. ANÁLISE ESPECTRAL DOS SINAIS PARA VISUALIZAÇÃO
f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog);

sinalCompostoFFT = fftshift(fft(sinalComposto)) / nPAnalog;
sinalPreFiltradoFFT = fftshift(fft(sinalPreFiltrado)) / nPAnalog;
sinalAmostradoFFT = fftshift(fft(sinalAmostrado)) / nPAnalog;
sinalReconstruidoFFT = fftshift(fft(sinalReconstruidoCompensado)) / nPAnalog;

%% 7. VISUALIZAÇÃO DOS RESULTADOS (Grid 6x2 Fluxo Lógico)
figure('Name', 'Pipeline: Filtragem -> Amostragem -> Reconstrução', 'Position', [50, 50, 1600, 1050]);

% --- LINHA 1: SINAIS ORIGINAIS (A e B separados) ---
subplot(6, 2, 1);
plot(t, sinalA, 'b', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('1. Sinal A (1 kHz) - Desejado');
xlim([0 0.002]);

subplot(6, 2, 2);
plot(t, sinalB, 'r', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('2. Sinal B (5 kHz) - Indesejado');
xlim([0 0.002]);

% --- LINHA 2: SINAL COMPOSTO ---
subplot(6, 2, 3);
plot(t, sinalComposto, 'k', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('3. Sinal Composto x(t) = A + B');
xlim([0 0.002]);

subplot(6, 2, 4);
stem(f, abs(sinalCompostoFFT), 'k', 'filled'); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('4. Espectro do Sinal Composto');
xlim([-15000 15000]);

% --- LINHA 3: SINAL PRÉ-FILTRADO (Extração de A) ---
subplot(6, 2, 5);
plot(t, sinalA, 'b--', 'LineWidth', 1.5); hold on;
plot(t, sinalPreFiltrado, 'g', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('5. Sinal Filtrado vs Sinal A Original');
legend('Sinal A Original', 'Pré-Filtrado (Extraído)', 'Location', 'best');
xlim([0 0.002]);

subplot(6, 2, 6);
stem(f, abs(sinalPreFiltradoFFT), 'g', 'filled'); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('6. Espectro Filtrado (Apenas 1 kHz)');
xlim([-15000 15000]);

% --- LINHA 5: SINAL AMOSTRADO FLAT-TOP ---
subplot(6, 2, 7);
plot(t, sinalPreFiltrado, 'g--', 'LineWidth', 1); hold on;
stairs(t, sinalAmostrado, 'r', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('9. Sinal Amostrado Flat-Top');
legend('Pré-filtrado', 'Amostrado', 'Location', 'best');
xlim([0 0.002]);

subplot(6, 2, 8);
stem(f, imag(sinalAmostradoFFT), 'r', 'filled'); grid on; % Trocado abs por imag
xlabel('Frequência (Hz)'); ylabel('Parte Imaginária'); title('10. Espectro do Sinal Amostrado');
xlim([-55000 55000]);

% --- LINHA 6: RECONSTRUÇÃO FINAL ---
subplot(6, 2, 9);
plot(t, sinalA, 'b--', 'LineWidth', 1.5); hold on;
plot(t, sinalReconstruidoCompensado, 'm', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('11. Sinal Reconstruído Final');
legend('Sinal A Desejado', 'Reconstruído', 'Location', 'best');
xlim([0 0.002]);

subplot(6, 2, 10);
stem(f, abs(sinalReconstruidoFFT), 'm', 'filled'); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('12. Espectro Final Reconstruído');
xlim([-55000 55000]);
