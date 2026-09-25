clear;
clc;
close all;

pkg load signal;

f1 = 2e3;        % Frequencia do Sinal A (util)
f2 = 5e3;        % Frequencia do Sinal B (fora da banda permitida)
A1 = 1.0;        % Amplitude do Sinal A
A2 = 1.5;        % Amplitude do Sinal B

fs = 4.5e3;        % Frequencia de amostragem
Ts = 1/fs;       % Periodo de amostragem
fAnalog = 360e3; % Frequencia da grade computacional (emula tempo continuo)
Ta = 1/fAnalog;  % Passo da grade computacional
tempoTotal = 1;  % Tempo total de simulação
t = 0:Ta:tempoTotal-Ta; % Grade computacional
nPAnalog = length(t);   % Total de pontos da grade
fCorte = 3e3;    % Frequencia de corte do filtro passa baixas anti aliasing

duty = 1/3;      % Duty cycle do pulso
d = duty*Ts;     % Largura do pulso em segundos

sinalA = A1*sin(2*pi*f1*t); % Sinal a ser amostrado
sinalB = A2*sin(2*pi*f2*t); % Sinal fora da banda
sinalComposto = sinalA + sinalB;

% Emulando um filtro ideal anti-aliasing
filtroIdeal = zeros(1,nPAnalog);

for i = 1:nPAnalog
    fAtual = (i-1)*(fAnalog/nPAnalog);

    if fAtual <= fCorte || fAtual >= (fAnalog-fCorte)
        filtroIdeal(i) = 1;
    else
        filtroIdeal(i) = 0;
    end
end

sinalCompostoFFT = fft(sinalComposto);
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal;
sinalFiltrado = real(ifft(sinalFiltradoFFT));

kPasso = round(fAnalog / fs);
kPulso = round(duty * kPasso);

% Trem de pulsos retangulares para amostragem natural
tremPulsos = zeros(1,nPAnalog);
tremPulsos(mod(0:nPAnalog-1, kPasso) < kPulso) = 1;
dutyReal = mean(tremPulsos); % c0 medido do trem de pulsos

sinalAmostrado = sinalFiltrado .* tremPulsos; % Amostragem natural: x(t) * p(t)

% Figura 1: Sinais no tempo (4 subplots)
figure;

subplot(4,1,1);
plot(t,sinalA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.01]);
yticks(-1:0.2:1);

subplot(4,1,2);
plot(t,sinalB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.01]);

subplot(4,1,3);
plot(t,sinalComposto);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Composto');
xlim([0 0.01]);

subplot(4,1,4);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);
yticks(-1:0.2:1);

saveas(gcf, sprintf('Natural_%gHz_01_Sinais_No_Tempo.png', fs));

% Figura 2: Processo de amostragem no tempo (4 subplots: com a amostragem isolada E a sobreposta)
figure;

subplot(4,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

subplot(4,1,2);
plot(t,tremPulsos);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de Pulsos');
xlim([0 0.01]);
ylim([-0.1 1.1]);

subplot(4,1,3);
plot(t,sinalAmostrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Amostrado Natural ');
xlim([0 0.01]);

subplot(4,1,4);
plot(t, sinalFiltrado, '--');
hold on;
plot(t, sinalAmostrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Comparação: Sinal Filtrado vs Sinal Amostrado');
legend('Sinal Filtrado', 'Sinal Amostrado');
xlim([0 0.01]);
hold off;

saveas(gcf, sprintf('Natural_%gHz_02_Processo_De_Amostragem.png', fs));

f = (-nPAnalog/2:nPAnalog/2-1)*(fAnalog/nPAnalog);

sinalFiltradoFFT_plot = fftshift(fft(sinalFiltrado))/nPAnalog;
tremPulsosFFT_plot = fftshift(fft(tremPulsos))/nPAnalog;
sinalAmostradoFFT_plot = fftshift(fft(sinalAmostrado))/nPAnalog;

magFiltrado = abs(sinalFiltradoFFT_plot);
faseFiltrado = angle(sinalFiltradoFFT_plot) * (180/pi);
faseFiltrado(magFiltrado < 0.01 * max(magFiltrado)) = 0;

magTrem = abs(tremPulsosFFT_plot);

magAmostrado = abs(sinalAmostradoFFT_plot);
faseAmostrado = angle(sinalAmostradoFFT_plot) * (180/pi);
faseAmostrado(magAmostrado < 0.01 * max(magAmostrado)) = 0;

% Figura 3: Espectro do sinal - Parte Imaginária (2 subplots)
figure;

subplot(2,1,1);
stem(f,imag(sinalFiltradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(2,1,2);
stem(f,imag(sinalAmostradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

saveas(gcf, sprintf('Natural_%gHz_03_Espectros_Do_Sinal_Imag.png', fs));

% Figura 4: Espectros de Magnitude (3 subplots)
figure;

subplot(3,1,1);
stem(f,magFiltrado,'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(3,1,2);
stem(f,magTrem,'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Trem de Pulsos');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(3,1,3);
stem(f,magAmostrado,'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

saveas(gcf, sprintf('Natural_%gHz_04_Espectros_De_Magnitude.png', fs));

% Figura 5: Espectros de Fase (2 subplots)
figure;

subplot(2,1,1);
stem(f,faseFiltrado,'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Espectro de Fase - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);
ylim([-180 180]);
yticks(-180:90:180);

subplot(2,1,2);
stem(f,faseAmostrado,'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Espectro de Fase - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-180 180]);
yticks(-180:90:180);

saveas(gcf, sprintf('Natural_%gHz_05_Espectros_De_Fase.png', fs));

% Reconstrução por filtro passa-baixas ideal
fCorteReconstrucao = fs/2;
filtroReconstrucao = zeros(1,nPAnalog);

for i = 1:nPAnalog
    fAtual = (i-1)*(fAnalog/nPAnalog);
    if fAtual <= fCorteReconstrucao || fAtual >= (fAnalog-fCorteReconstrucao)
        filtroReconstrucao(i) = 1;
    else
        filtroReconstrucao(i) = 0;
    end
end

sinalReconstruidoFFT = fft(sinalAmostrado) .* filtroReconstrucao / dutyReal;
sinalReconstruido = real(ifft(sinalReconstruidoFFT));

filtroReconstrucaoPlot = fftshift(filtroReconstrucao);
sinalReconstruidoFFT_plot = fftshift(sinalReconstruidoFFT)/nPAnalog;

% Figura 6: Etapa de Reconstrução e Comparação (4 subplots)
figure;

subplot(4,1,1);
stem(f,imag(sinalAmostradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(4,1,2);
plot(f,filtroReconstrucaoPlot);
grid on;
xlabel('Frequência (Hz)');
ylabel('Ganho');
title('Filtro Ideal de Reconstrução');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-0.2 1.2]);

subplot(4,1,3);
stem(f,imag(sinalReconstruidoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Reconstruído');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(4,1,4);
plot(t, sinalFiltrado,'DisplayName', 'Sinal Filtrado');
hold on;
plot(t, sinalReconstruido, '--','DisplayName', 'Sinal Reconstruído');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Comparação: Sinal Filtrado vs Sinal Reconstruído');
xlim([0 0.01]);
legend('Location', 'northeast');
hold off;

saveas(gcf, sprintf('Natural_%gHz_06_Reconstrucao_E_Comparacao.png', fs));
