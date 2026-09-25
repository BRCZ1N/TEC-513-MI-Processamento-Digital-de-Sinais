clear;
clc;
close all;

pkg load signal;

f1 = 2e3;        % Frequencia do Sinal A (util)
f2 = 5e3;        % Frequencia do Sinal B (fora da banda permitida)
A1 = 1.0;        % Amplitude do Sinal A
A2 = 1.5;        % Amplitude do Sinal B

fs = 3e3;      % Frequencia de amostragem
Ts = 1/fs;       % Periodo de amostragem
fAnalog = 360e3; % Frequencia da grade computacional (emula tempo continuo)
Ta = 1/fAnalog;  % Passo da grade computacional
tempoTotal = 1;% Tempo total de simulação
t = 0:Ta:tempoTotal-Ta;% Grade computacional
nPAnalog = length(t);% Total de pontos da grade
fCorte = 3e3;    % Frequencia de corte do filtro passa baixas anti aliasing

duty = 1/3;      % Duty cycle do pulso
d = duty*Ts;     % Largura do pulso em segundos

sinalA = A1*sin(2*pi*f1*t);% Sinal a ser amostrado
sinalB = A2*sin(2*pi*f2*t);% Sinal fora da banda
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

sinalAmostrado = sinalFiltrado .* tremPulsos;% Amostragem natural: x(t) * p(t)

% Figura 1: Sinais no tempo (4 subplots)
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

saveas(gcf, sprintf('Natural_fs_%gHz_01_Sinais_No_Tempo.png', fs));

% Figura 2: Processo de amostragem no tempo (3 subplots - sem o 4° gráfico)
figure;

subplot(3,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.05]);

subplot(3,1,2);
plot(t,tremPulsos);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de Pulsos');
xlim([0 0.05]);
ylim([-0.1 1.1]);

subplot(3,1,3);
plot(t,sinalAmostrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amostragem (Sinal Amostrado)');
xlim([0 0.05]);

saveas(gcf, sprintf('Natural_fs_%gHz_02_Processo_De_Amostragem.png', fs));

f = (-nPAnalog/2:nPAnalog/2-1)*(fAnalog/nPAnalog);

sinalFiltradoFFT_plot = fftshift(fft(sinalFiltrado))/nPAnalog;
tremPulsosFFT_plot = fftshift(fft(tremPulsos))/nPAnalog;
sinalAmostradoFFT_plot = fftshift(fft(sinalAmostrado))/nPAnalog;

% Figura 3: Espectros de Magnitude (3 subplots)
figure;

subplot(3,1,1);
stem(f,abs(sinalFiltradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(3,1,2);
stem(f,abs(tremPulsosFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Trem de Pulsos');
xlim([-2.5*fs 2.5*fs]);
xticks(-2.5*fs : fs : 2.5*fs);

subplot(3,1,3);
stem(f,abs(sinalAmostradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Amostrado');
xlim([-2.5*fs 2.5*fs]);
xticks(-2.5*fs : fs : 2.5*fs);

saveas(gcf, sprintf('Natural_fs_%gHz_03_Espectros_De_Magnitude.png', fs));

% Figura 4: Espectros com foco na parte imaginária (2 subplots)
figure;

subplot(2,1,1);
stem(f,imag(sinalFiltradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(2,1,2);
stem(f,imag(sinalAmostradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Amostrado');
xlim([-2.5*fs 2.5*fs]);
xticks(-2.5*fs : fs : 2.5*fs);

saveas(gcf, sprintf('Natural_fs_%gHz_04_Espectros_Do_Sinal_Imag.png', fs));

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

% Compensação da atenuação do pulso dividindo pelo duty cycle real
sinalReconstruidoFFT = fft(sinalAmostrado) .* filtroReconstrucao / dutyReal;
sinalReconstruido = real(ifft(sinalReconstruidoFFT));

filtroReconstrucaoPlot = fftshift(filtroReconstrucao);
sinalReconstruidoFFT_plot = fftshift(sinalReconstruidoFFT)/nPAnalog;

% Etapa de Reconstrução
figure;

subplot(4,1,1);
stem(f,imag(sinalAmostradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Amostrado');
xlim([-2.5*fs 2.5*fs]);
xticks(-2.5*fs : fs : 2.5*fs);

subplot(4,1,2);
plot(f,filtroReconstrucaoPlot);
grid on;
xlabel('Frequência (Hz)');
ylabel('Ganho');
title('Filtro Ideal de Reconstrução');
xlim([-2.5*fs 2.5*fs]);
xticks(-2.5*fs : fs : 2.5*fs);
ylim([-0.2 1.2]);

subplot(4,1,3);
stem(f,imag(sinalReconstruidoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Reconstruído');
xlim([-2.5*fs 2.5*fs]);
xticks(-2.5*fs : fs : 2.5*fs);

subplot(4,1,4);
plot(t, sinalFiltrado, 'DisplayName', 'Sinal Filtrado');
hold on;
plot(t, sinalReconstruido, '--', 'DisplayName', 'Sinal Reconstruído');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Comparação: Sinal Filtrado vs Sinal Reconstruído');
xlim([0 0.05]);
legend('Location', 'northeast');
hold off;

saveas(gcf, sprintf('Natural_fs_%gHz_05_Reconstrucao_E_Comparacao.png', fs));
