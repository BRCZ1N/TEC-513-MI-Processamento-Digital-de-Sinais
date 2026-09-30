clear;
clc;
close all;

% pkg load signal;

f1 = 2e3;        % Frequencia da senoide a ser amostrada
f2 = 5e3;        % Frequencia da senoide do impulso na frequencia
A1 = 1.0;        % Amplitude da senoide a ser amostrada
A2 = 1.5;        % Amplitude da senoide do impulso na frequencia

fs = 4.5e3;      % Frequencia de amostragem
Ts = 1/fs;       % Periodo de amostragem
fAnalog = 360e3; % Frequencia da grade computacional
Ta = 1/fAnalog;  % Periodo da grade computacional
tempoTotal = 1;  % Tempo total de simulação
t = 0:Ta:tempoTotal-Ta; % Grade computacional
nPAnalog = length(t);   % Total de pontos da grade
fCorte = 3e3;    % Frequencia de corte do filtro passa baixas anti aliasing

sinalA = A1*sin(2*pi*f1*t); % Sinal a ser amostrado
sinalB = A2*sin(2*pi*f2*t); % Sinal que fornece o impulso em frequencia
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

sinalCompostoFFT = fft(sinalComposto, nPAnalog); % Jogando na frequencia para retirar o espectro da segunda senoide
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal; % Filtro no tempo pra retirar o espectro da segunda senoide
sinalFiltrado = real(ifft(sinalFiltradoFFT)); % Volto pro dominio do tempo pra ter o meu sinal filtrado e no tempo

% Figura 1: Sinais no Tempo
figure;
datacursormode on; % Habilita o cursor de dados (Data Tips) para interação

subplot(3,1,1);
plot(t,sinalA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.01]);

subplot(3,1,2);
plot(t,sinalB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.01]);

subplot(3,1,3);
plot(t,sinalComposto);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Composto');
xlim([0 0.01]);

saveas(gcf, sprintf('Ideal_%gHz_01_Sinais_No_Tempo_A.png', fs));

figure;
datacursormode on;

subplot(1,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

saveas(gcf, sprintf('Ideal_%gHz_01_Sinais_No_Tempo_B.png', fs));

kPasso = round(fAnalog / fs);
tremImpulsos = zeros(1,nPAnalog); % Grade computacional do trem de impulsos
tremImpulsos(1:kPasso:end) = 1;    % Faço a área dele tender a 1 discretamente

sinalAmostrado = sinalFiltrado .* tremImpulsos; % Amostrando o sinal x(t).s(t)

% Figura 2: Processo de amostragem no tempo
figure;
datacursormode on;

subplot(3,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

subplot(3,1,2);
stem(t,tremImpulsos,'.');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de Impulsos');
xlim([0 0.01]);

subplot(3,1,3);
stem(t,sinalAmostrado,'.');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Amostrado Ideal');
xlim([0 0.01]);

saveas(gcf, sprintf('Ideal_%gHz_02_Processo_De_Amostragem_A.png', fs));

figure;
datacursormode on;

subplot(1,1,1);
plot(t,sinalFiltrado,'--');
hold on;
stem(t,sinalAmostrado,'.');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado e Sinal Amostrado');
legend('Sinal Filtrado','Amostrado');
xlim([0 0.01]);
hold off;

saveas(gcf, sprintf('Ideal_%gHz_02_Processo_De_Amostragem_B.png', fs));

f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog); % Criando a grade computacional da frequencia

sinalFiltradoFFT_plot = fftshift(fft(sinalFiltrado)) / nPAnalog;
tremImpulsosFFT_plot = fftshift(fft(tremImpulsos)) / nPAnalog * fAnalog;
sinalAmostradoFFT_plot = fftshift(fft(sinalAmostrado)) / nPAnalog * fAnalog;

magFiltrado = abs(sinalFiltradoFFT_plot);
faseFiltrado = angle(sinalFiltradoFFT_plot) * (180/pi);
faseFiltrado(magFiltrado < 0.01 * max(magFiltrado)) = 0;

magTrem = abs(tremImpulsosFFT_plot);

magAmostrado = abs(sinalAmostradoFFT_plot);
faseAmostrado = angle(sinalAmostradoFFT_plot) * (180/pi);
faseAmostrado(magAmostrado < 0.01 * max(magAmostrado)) = 0;

% Figura 3: Espectro Imaginário
figure;
datacursormode on;

subplot(2,1,1);
stem(f, imag(sinalFiltradoFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Filtrado');
xlim([-f1*5 f1*5]);
xticks(-f1*5:f1:f1*5);

subplot(2,1,2);
stem(f, imag(sinalAmostradoFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

saveas(gcf, sprintf('Ideal_%gHz_03_Espectros_Do_Sinal_Imag.png', fs));

% Figura 4: Espectros de Magnitude
figure;
datacursormode on;

subplot(3,1,1);
stem(f, magFiltrado, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(3,1,2);
stem(f, magTrem, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Trem de Impulsos');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(3,1,3);
stem(f, magAmostrado, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

saveas(gcf, sprintf('Ideal_%gHz_04_Espectros_De_Magnitude.png', fs));

% Figura 5: Espectros de Fase
figure;
datacursormode on;

subplot(2,1,1);
stem(f, faseFiltrado, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Espectro de Fase - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);
ylim([-180 180]);
yticks(-180:90:180);

subplot(2,1,2);
stem(f, faseAmostrado, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Espectro de Fase - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-180 180]);
yticks(-180:90:180);

saveas(gcf, sprintf('Ideal_%gHz_05_Espectros_De_Fase.png', fs));

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

sinalReconstruidoFFT = fft(sinalAmostrado) .* filtroReconstrucao * kPasso;
sinalReconstruido = real(ifft(sinalReconstruidoFFT));

filtroReconstrucaoPlot = fftshift(filtroReconstrucao);
sinalReconstruidoFFT_plot = fftshift(sinalReconstruidoFFT) / nPAnalog;

% Figura 6: Reconstrução e Comparação
figure;
datacursormode on;

subplot(3,1,1);
stem(f,imag(sinalAmostradoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(3,1,2);
plot(f,filtroReconstrucaoPlot);
grid on;
xlabel('Frequência (Hz)');
ylabel('Ganho');
title('Filtro Ideal de Reconstrução');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-0.2 1.2]);

subplot(3,1,3);
stem(f,imag(sinalReconstruidoFFT_plot),'.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Reconstruído');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

saveas(gcf, sprintf('Ideal_%gHz_06_Reconstrucao_E_Comparacao_A.png', fs));

figure;
datacursormode on;

subplot(1,1,1);
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

saveas(gcf, sprintf('Ideal_%gHz_06_Reconstrucao_E_Comparacao_B.png', fs));