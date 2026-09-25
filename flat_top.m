clear;
clc;
close all;

pkg load signal;

f1 = 2e3;        % Frequencia da senoide a ser amostrada
f2 = 5e3;        % Frequencia da segunda senoide
A1 = 1.0;        % Amplitude da senoide a ser amostrada
A2 = 1.5;        % Amplitude da segunda senoide

fs = 4.5e3;        % Frequencia de amostragem
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

saveas(gcf, sprintf('FlatTop_%gHz_01_Sinais_No_Tempo.png', fs));

kPasso = round(Ts/Ta);

d = Ts/3; % Tempo de abertura do pulso
numPontosPulso = round(d/Ta);
dReal = numPontosPulso * Ta;

amostras = sinalFiltrado(1:kPasso:end); % Obtendo as amostras do sinal filtrado
pulsoFlatTop = [ones(1, numPontosPulso), zeros(1, kPasso - numPontosPulso)]; % Pulso retangular Flat-Top

sinalAmostradoFlatTop = kron(amostras, pulsoFlatTop); % Construindo o sinal Flat-Top
sinalAmostradoFlatTop = sinalAmostradoFlatTop(1:nPAnalog);

% Figura 2: Processo de amostragem no tempo (3 subplots - com os 2 gráficos de sobreposição no 3º subplot)
figure;

subplot(3,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

subplot(3,1,2);
stairs(t,sinalAmostradoFlatTop);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Amostrado Flat-Top');
xlim([0 0.01]);

subplot(3,1,3);
plot(t, sinalFiltrado);
hold on;
stairs(t, sinalAmostradoFlatTop);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Contínuo e Amostragem Flat-Top sobreposta');
legend('Sinal Contínuo', 'Flat-Top Amostrado');
xlim([0 0.01]);
hold off;

saveas(gcf, sprintf('FlatTop_%gHz_02_Processo_De_Amostragem.png', fs));

f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog); % Grade computacional da frequencia

sinalFiltradoFFT_plot = fftshift(fft(sinalFiltrado)) / nPAnalog;
sinalAmostradoFlatTopFFT_plot = fftshift(fft(sinalAmostradoFlatTop)) / nPAnalog;

magFiltrado = abs(sinalFiltradoFFT_plot);
faseFiltrado = angle(sinalFiltradoFFT_plot) * (180/pi);
faseFiltrado(magFiltrado < 0.01 * max(magFiltrado)) = 0;

magFlatTop = abs(sinalAmostradoFlatTopFFT_plot);
faseFlatTop = angle(sinalAmostradoFlatTopFFT_plot) * (180/pi);
faseFlatTop(magFlatTop < 0.01 * max(magFlatTop)) = 0;

% Figura 3: Espectro do sinal - Parte Imaginária (2 subplots)
figure;

subplot(2,1,1);
stem(f, imag(sinalFiltradoFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(2,1,2);
stem(f, imag(sinalAmostradoFlatTopFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Amostrado Flat-Top');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

saveas(gcf, sprintf('FlatTop_%gHz_03_Espectros_Do_Sinal_Imag.png', fs));

% Figura 4: Espectros de Magnitude (2 subplots)
figure;

subplot(2,1,1);
stem(f, magFiltrado, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude - Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);
xtickangle(45);

subplot(2,1,2);
stem(f, magFlatTop, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Flat-Top');
xlim([-15000 15000]);
xticks(-15000:3000:15000);

saveas(gcf, sprintf('FlatTop_%gHz_04_Espectros_De_Magnitude.png', fs));

% Figura 5: Espectros de Fase (2 subplots)
figure;

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
stem(f, faseFlatTop, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Espectro de Fase - Sinal Amostrado Flat-Top');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-180 180]);
yticks(-180:90:180);

saveas(gcf, sprintf('FlatTop_%gHz_05_Espectros_De_Fase.png', fs));

% Espectro Flat-Top e Envolvente Sinc adicional
fSinc = f;
fSinc(fSinc == 0) = 1e-10;
envolventeSinc = abs(sin(pi * fSinc * dReal) ./ (pi * fSinc * dReal));
amplitudeMaxima = max(magFlatTop);
envolventeSinc = envolventeSinc * amplitudeMaxima;

figure;

stem(f, magFlatTop, '.');
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

saveas(gcf, sprintf('FlatTop_%gHz_06_Envolvente_Sinc.png', fs));

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

% Figura 7: Reconstrução e Comparação (4 subplots)
figure;

subplot(4,1,1);
stem(f, imag(sinalAmostradoFlatTopFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Amostrado Flat-Top');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(4,1,2);
plot(f, filtroReconstrucaoPlot);
grid on;
xlabel('Frequência (Hz)');
ylabel('Ganho');
title('Filtro Ideal de Reconstrução');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-0.2 1.2]);

subplot(4,1,3);
stem(f, imag(sinalReconstruidoFFT_plot), '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do sinal - Sinal Reconstruído');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(4,1,4);
plot(t, sinalFiltrado, 'DisplayName', 'Sinal Filtrado');
hold on;
plot(t, sinalReconstruido, '--', 'DisplayName', 'Sinal Reconstruído');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Comparação: Sinal Filtrado vs Sinal Reconstruído');
xlim([0 0.01]);
legend('Location', 'northeast');
hold off;

saveas(gcf, sprintf('FlatTop_%gHz_07_Reconstrucao_E_Comparacao.png', fs));
