clear;
clc;
close all;

pkg load signal;

f1 = 2e3;
f2 = 5e3;
A1 = 1.0;
A2 = 1.5;

fs = 3e3;
Ts = 1/fs;
fAnalog = 360e3;
Ta = 1/fAnalog;
tempoTotal = 1;
t = 0:Ta:tempoTotal-Ta;
nPAnalog = length(t);

fCorte = 3e3;

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

sinalCompostoFFT = fft(sinalComposto, nPAnalog);
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal;
sinalFiltrado = real(ifft(sinalFiltradoFFT));

figure;

subplot(4,1,1);
plot(t,sinalA);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.05]);
yticks(-1:0.2:1);

subplot(4,1,2);
plot(t,sinalB);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.05]);
yticks(-1.5:0.3:1.5);

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
yticks(-1.5:0.3:1.5);

saveas(gcf, sprintf('Ideal_fs_%gHz_01_Sinais_No_Tempo.png', fs));

kPasso = round(fAnalog / fs);
tremImpulsos = zeros(1,nPAnalog);
tremImpulsos(1:kPasso:end) = 1;

sinalAmostrado = sinalFiltrado .* tremImpulsos;

figure;

subplot(3,1,1);
plot(t,sinalFiltrado);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.05]);

subplot(3,1,2);
stem(t,tremImpulsos,'.');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de Impulsos');
xlim([0 0.05]);

subplot(3,1,3);
stem(t,sinalAmostrado,'.');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amostragem');
xlim([0 0.05]);
hold off;

saveas(gcf, sprintf('Ideal_fs_%gHz_02_Processo_De_Amostragem.png', fs));

f = (-nPAnalog/2 : nPAnalog/2 - 1) * (fAnalog / nPAnalog);

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

figure;

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

saveas(gcf, sprintf('Ideal_fs_%gHz_03_Espectros_Do_Sinal_Imag.png', fs));

figure;

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

saveas(gcf, sprintf('Ideal_fs_%gHz_04_Espectros_De_Magnitude.png', fs));

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
stem(f, faseAmostrado, '.');
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Espectro de Fase - Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-180 180]);
yticks(-180:90:180);

saveas(gcf, sprintf('Ideal_fs_%gHz_05_Espectros_De_Fase.png', fs));

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
xlim([0 0.05]);
legend('Location', 'northeast');
hold off;

saveas(gcf, sprintf('Ideal_fs_%gHz_06_Reconstrucao_E_Comparacao.png', fs));
