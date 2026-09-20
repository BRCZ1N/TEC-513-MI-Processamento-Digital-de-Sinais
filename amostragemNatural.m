% =========================================================================
% SigmaDSP Inc. - Problema 01: Amostragem NATURAL e Reconstrução de Sinal
% Ferramenta: MATLAB / Octave
% Layout dos gráficos no mesmo padrão 5x2 da figura de referência:
%   Linha 1: Sinal A                       | Sinal B
%   Linha 2: Sinal composto                | Sinal Filtrado
%   Linha 3: Trem de pulsos                | Sinal amostrado
%   Linha 4: Espectro do Sinal Filtrado    | Espectro do Trem de Pulsos
%   Linha 5: Espectro Amostrado (Im{S(f)}) | Espectro Amostrado (|S(f)|)
% Figura 2: reconstrução (tempo e frequência).
% =========================================================================
clear; clc; close all;

cor = [0 0.447 0.741];      % cor única usada em TODOS os gráficos

%% 1. PARÂMETROS DA SIMULAÇÃO (TEMPO CONTÍNUO EMULADO)
% Fs_sim = 400 kHz e T_total = 0.1 s  ->  N = 40000 pontos, df = 10 Hz.
% fA, fB e fs (divisores de 100 kHz) caem exatamente em bins da FFT.
Fs_sim  = 400e3;
T_total = 0.1;
N       = round(T_total * Fs_sim);
n       = 0:N-1;
t       = n / Fs_sim;
f_eixo  = (-N/2 : N/2 - 1) * (Fs_sim / N);

%% 2. SINAIS ANALÓGICOS (senoides geradas com imag)
fA = 1000;  AA = 1.0;      % Sinal A: senoide útil
fB = 6000;  AB = 1.2;      % Sinal B: componente acima da banda permitida

% imag(exp(j*w*t)) = sen(w*t)
gA     = AA * imag(exp(1j*2*pi*fA*t));
gB     = AB * imag(exp(1j*2*pi*fB*t));
g_comp = gA + gB;                          % Sinal composto

%% 3. FILTRO DE LIMITAÇÃO EM BANDA (passa-baixa ideal, no domínio da frequência)
% true : B é removido ANTES de amostrar (como na figura de referência).
% false: o composto é amostrado direto e B pode causar aliasing.
usar_filtro_banda = true;
f_max_filtro      = 2000;                  % fA < f_max_filtro < fB

Gcomp_f = fftshift(fft(g_comp)) / N;
if usar_filtro_banda
    H_banda = double(abs(f_eixo) <= f_max_filtro);
    g_filt  = real(ifft(ifftshift(Gcomp_f .* H_banda))) * N;
else
    g_filt  = g_comp;
end
Gfilt_f = fftshift(fft(g_filt)) / N;

%% 4. FREQUÊNCIA DE AMOSTRAGEM (AJUSTÁVEL)
%   Sem aliasing: Fs_amostragem = 10000 Hz  (> 2*fA)
%   Com aliasing: Fs_amostragem = 1250  Hz  (< 2*fA)  -> fA aparece em 250 Hz
% Use divisores de 100 kHz (ex.: 10000, 5000, 4000, 2500, 2000, 1250, 1000).
Fs_amostragem = 10000;
Ts   = 1 / Fs_amostragem;
duty = 0.25;
d    = duty * Ts;

%% 5. TREM DE PULSOS p(t) E SINAL AMOSTRADO s(t)
% Pulso construído por contagem de amostras (duty cycle exato).
Ns = Fs_sim / Fs_amostragem;               % amostras por período de amostragem
Nd = round(duty * Ns);                     % amostras por pulso
if abs(Ns - round(Ns)) > 1e-9 || abs(duty*Ns - Nd) > 1e-9
    warning(['Fs_sim/Fs_amostragem nao gera pulso com numero inteiro de ' ...
             'amostras: o duty cycle real sera diferente de %.2f.'], duty);
end
Ns  = round(Ns);
p_t = double(mod(n, Ns) < Nd);
duty_real = mean(p_t);                     % c0 medido

s_t = g_filt .* p_t;                       % amostragem natural: s = g * p

%% 6. ESPECTROS
P_f = fftshift(fft(p_t)) / N;              % linhas em k*fs com peso c_k
S_f = fftshift(fft(s_t)) / N;

%% 7. RECONSTRUÇÃO (PASSA-BAIXA IDEAL, corte em fs/2)
fc = Fs_amostragem / 2;
H_rec = double(abs(f_eixo) <= fc);

G_rec_f = S_f .* H_rec / duty_real;        % 1/c0 = Ts/d
g_rec_t = real(ifft(ifftshift(G_rec_f))) * N;

erro_rms = sqrt(mean((g_rec_t - gA).^2));  % referência: sinal útil A
fprintf('fs = %g Hz | fc = %g Hz | duty real = %.4f | erro RMS vs A = %.4f\n', ...
        Fs_amostragem, fc, duty_real, erro_rms);

%% 8. FIGURA 1 - PADRÃO 5x2
it  = t <= 0.01;                           % janela de 10 ms nos gráficos do tempo
tt  = t(it);

lim_filt = 1.5 * fB;                       % eixo do "Espectro do Sinal Filtrado"
lim_trem = 2.1 * Fs_amostragem;            % eixo do espectro do trem
lim_samp = 3.5 * Fs_amostragem;            % eixo do espectro amostrado
m_filt = abs(f_eixo) <= lim_filt;
m_trem = abs(f_eixo) <= lim_trem;
m_samp = abs(f_eixo) <= lim_samp;

figure('Name', 'Amostragem Natural - Padrão 5x2', 'NumberTitle', 'off');

% --- Linha 1 ---
subplot(5,2,1);
plot(tt, gA(it), 'Color', cor); title('Sinal A');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);

subplot(5,2,2);
plot(tt, gB(it), 'Color', cor); title('Sinal B');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);

% --- Linha 2 ---
subplot(5,2,3);
plot(tt, g_comp(it), 'Color', cor); title('Sinal composto');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);

subplot(5,2,4);
plot(tt, g_filt(it), 'Color', cor); title('Sinal Filtrado');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);

% --- Linha 3 ---
subplot(5,2,5);
plot(tt, p_t(it), 'Color', cor); title(sprintf('Trem de pulsos p(t)  [f_s = %g Hz, d = %.3f ms]', ...
                                  Fs_amostragem, d*1e3));
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]); ylim([-0.1 1.1]);

subplot(5,2,6);
plot(tt, s_t(it), 'Color', cor); title('Sinal amostrado s(t)');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);

% --- Linha 4 ---
subplot(5,2,7);
stem(f_eixo(m_filt), abs(Gfilt_f(m_filt)), 'MarkerSize', 3, 'Color', cor);
title('Espectro do Sinal Filtrado');
xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on; xlim([-lim_filt lim_filt]);

subplot(5,2,8);
stem(f_eixo(m_trem), abs(P_f(m_trem)), 'MarkerSize', 3, 'Color', cor);
title('Espectro do Trem de Pulsos');
xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on; xlim([-lim_trem lim_trem]);

% --- Linha 5 ---
subplot(5,2,9);
stem(f_eixo(m_samp), imag(S_f(m_samp)), 'MarkerSize', 3, 'Color', cor);
title('Espectro do Sinal Amostrado - Im\{S(f)\}');
xlabel('Frequência (Hz)'); ylabel('Im\{S(f)\}'); grid on; xlim([-lim_samp lim_samp]);

subplot(5,2,10);
stem(f_eixo(m_samp), abs(S_f(m_samp)), 'MarkerSize', 3, 'Color', cor);
title('Espectro do Sinal Amostrado - |S(f)|');
xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on; xlim([-lim_samp lim_samp]);

%% 9. FIGURA 2 - RECONSTRUÇÃO
figure('Name', 'Amostragem Natural - Reconstrução', 'NumberTitle', 'off');

subplot(3,2,1);
plot(tt, g_rec_t(it), '-', 'Color', cor, 'LineWidth', 1.2); hold on;
im = find(it); im = im(1:40:end);          % marcadores esparsos p/ o sinal original
plot(t(im), gA(im), 'o', 'Color', cor, 'MarkerSize', 4);
title(sprintf('Reconstruído vs Sinal A  [erro RMS = %.3f]', erro_rms));
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);
legend('Reconstruído', 'Sinal A');

subplot(3,2,2);
plot(tt, g_rec_t(it) - gA(it), 'Color', cor);
title('Erro de reconstrução  g_{rec}(t) - A(t)');
xlabel('Tempo (s)'); ylabel('Amplitude'); grid on; xlim([0 0.01]);

subplot(3,2,3);
stem(f_eixo(m_samp), abs(S_f(m_samp)), 'MarkerSize', 3, 'Color', cor); hold on;
yl = [0 max(abs(S_f(m_samp)))*1.05];
line([-fc -fc], yl, 'Color', cor, 'LineStyle', '--', 'LineWidth', 1.5);
line([ fc  fc], yl, 'Color', cor, 'LineStyle', '--', 'LineWidth', 1.5);
title(sprintf('|S(f)| (linhas verdes: f_c = %g Hz)', fc));
xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on; xlim([-lim_samp lim_samp]);

subplot(3,2,4);
m_rec = abs(f_eixo) <= lim_filt;
stem(f_eixo(m_rec), abs(G_rec_f(m_rec)), 'MarkerSize', 3, 'Color', cor);
title('Espectro Reconstruído |G_{rec}(f)|');
xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on; xlim([-lim_filt lim_filt]);

subplot(3,2,5);
stem(f_eixo(m_rec), imag(G_rec_f(m_rec)), 'MarkerSize', 3, 'Color', cor);
title('Im\{G_{rec}(f)\}');
xlabel('Frequência (Hz)'); ylabel('Im\{G_{rec}(f)\}'); grid on; xlim([-lim_filt lim_filt]);

subplot(3,2,6);
stem(f_eixo(m_rec), imag(Gfilt_f(m_rec)), 'MarkerSize', 3, 'Color', cor);
title('Im\{G(f)\} do Sinal Filtrado (antes de amostrar)');
xlabel('Frequência (Hz)'); ylabel('Im\{G(f)\}'); grid on; xlim([-lim_filt lim_filt]);