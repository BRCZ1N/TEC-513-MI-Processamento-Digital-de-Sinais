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

%% 2. GERAÇÃO DOS SINAIS DE AMOSTRAGEM
kPasso = round(Ts/Ta);
d = (1/3) * Ts;             % Largura do pulso de retenção alterada para 1/3 de Ts
numPontosPulso = round(d/Ta);

% Trem de Pulsos Retangulares (Amostrador Prático)
tremPulsos = zeros(1, nPAnalog);
for i = 0 : floor(nPAnalog/kPasso) - 1
    idx = i*kPasso + 1;
    if idx <= nPAnalog
        idx_fim = min(idx + numPontosPulso - 1, nPAnalog);
        tremPulsos(idx:idx_fim) = 1;
    end
end

%% 3. AMOSTRAGEM
sinalFlatTop = sinalComposto .* tremPulsos;

%% 4. ANÁLISE ESPECTRAL (FFT)
f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog);

sinalCompostoFFT = fftshift(fft(sinalComposto)) / nPAnalog;
tremPulsosFFT = fftshift(fft(tremPulsos)) / nPAnalog;

% Corrigido: Usar a FFT direta do sinal em vez da convolução para manter o vetor do mesmo tamanho
sinalFlatTopFFT = fftshift(fft(sinalFlatTop)) / nPAnalog;

% Envoltória Teórica H(f) = (d/Ts) * sinc(f * d) do Efeito de Abertura
H_teorico = (d/Ts) * sinc(f * d);

%% 5. FILTRAGEM IDEAL NO DOMÍNIO DA FREQUÊNCIA
fCorte = 3e3;                         % Frequência de corte (3 kHz)

% 5.1. FFT do sinal amostrado (sem shift, para facilitar a aplicação do filtro)
X_k = fft(sinalFlatTop);
filter_ideal = zeros(1, nPAnalog);

% 5.2. Criar manualmente o filtro passa-baixa ideal
for i = 1:nPAnalog
    % Frequência correspondente ao índice atual da FFT
    f_atual = (i - 1) * (fAnalog / nPAnalog);

    % Permite as frequências abaixo do corte ou as frequências espelhadas (Nyquist a fs)
    if f_atual <= fCorte || f_atual >= (fAnalog - fCorte)
        filter_ideal(i) = 1;
    else
        filter_ideal(i) = 0;  % Bloqueia as frequências acima do corte
    end
end

% 5.3. Aplicar o filtro no domínio da frequência
X_k_filtered = X_k .* filter_ideal;

% 5.4. Regressar ao domínio do tempo com IFFT (retirando pequenos resíduos imaginários com 'real')
sinalFiltrado = real(ifft(X_k_filtered));

% 5.5. Compensação de Amplitude
% O duty cycle reduz a energia do sinal, então multiplicamos pelo inverso para recuperar a amplitude
ganhoReconstrucao = Ts / d;
sinalFiltradoCompensado = sinalFiltrado * ganhoReconstrucao;

%% 6. VISUALIZAÇÃO DOS RESULTADOS (Grid 5x2)
figure('Name', 'Análise Completa: Sinais e Amostragem', 'Position', [50, 50, 1600, 1000]);

% --- LINHA 1: SINAL A e SINAL B ---
subplot(5, 2, 1);
plot(t, sinalA, 'b', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('1. Sinal A (1 kHz) - Desejado');
xlim([0 0.002]);

subplot(5, 2, 2);
plot(t, sinalB, 'r', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('2. Sinal B (5 kHz) - Aliasing/Ruído');
xlim([0 0.002]);

% --- LINHA 2: SINAL COMPOSTO E ESPECTRO ---
subplot(5, 2, 3);
plot(t, sinalComposto, 'k', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('3. Sinal Composto x(t) = A + B');
xlim([0 0.002]);

subplot(5, 2, 4);
stem(f, abs(sinalCompostoFFT), 'k', 'filled'); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('4. Espectro Composto |X(f)|');
xlim([-15000 15000]);

% --- LINHA 3: TREM DE PULSOS ---
subplot(5, 2, 5);
stem(t, tremPulsos, 'k', 'Marker', 'none'); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('5. Trem de Pulsos p(t) (via stem)');
xlim([0 0.002]); ylim([-0.2 1.2]);

subplot(5, 2, 6);
plot(f, abs(tremPulsosFFT), 'k', 'LineWidth', 1); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('6. Espectro do Trem de Pulsos |P(f)|');
xlim([-55000 55000]);

% --- LINHA 4: SINAL AMOSTRADO FLAT-TOP ---
subplot(5, 2, 7);
plot(t, sinalComposto, 'k--', 'LineWidth', 1); hold on;
plot(t, sinalFlatTop, 'r', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('7. Sinal Amostrado (Multiplicação)');
legend('x(t) Composto', 'x_{amostrado}(t)', 'Location', 'best');
xlim([0 0.002]);

subplot(5, 2, 8);
stem(f, abs(sinalFlatTopFFT), 'r', 'filled'); hold on;
plot(f, abs(H_teorico * max(abs(sinalCompostoFFT))), 'k--', 'LineWidth', 1.5); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('8. Espectro do Sinal Amostrado e Sinc');
legend('Espectro', 'Envoltória Sinc', 'Location', 'best');
xlim([-55000 55000]);

% --- LINHA 5: RECONSTRUÇÃO E COMPARAÇÃO IDEAL ---
subplot(5, 2, 9);
plot(t, sinalA, 'b--', 'LineWidth', 1.5); hold on;
plot(t, sinalFiltradoCompensado, 'g', 'LineWidth', 1.5); grid on;
xlabel('Tempo (s)'); ylabel('Amplitude'); title('9. Sinal Reconstruído (Filtro Ideal Domínio da Frequência)');
legend('Sinal A (1 kHz)', 'Reconstruído', 'Location', 'best');
xlim([0 0.002]);

subplot(5, 2, 10);
stem(f, abs(fftshift(X_k_filtered) / nPAnalog), 'm', 'filled'); grid on;
xlabel('Frequência (Hz)'); ylabel('Magnitude'); title('10. Espectro do Sinal Filtrado');
xlim([-55000 55000]);
