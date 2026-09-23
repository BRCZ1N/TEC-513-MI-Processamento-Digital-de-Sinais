clear;
clc;
close all;

pkg load signal;

f1 = 2e3;%Frequencia da senoide a ser amostrada
f2 = 5e3;%Frequencia da senoide do impulso na frequencia
A1 = 1.0;%Amplitude da senoide a ser amostrada
A2 = 1.5;%Amplitude da senoide do impulso na frequencia

fs = 3e3; %Frequencia de amostragem
Ts = 1/fs; %Periodo de amostragem
fAnalog = 360e3; %Frequencia da grade computacional
Ta = 1/fAnalog;%Frequencia da grade computacional
tempoTotal = 1;%Tempo total de simulação
t = 0:Ta:tempoTotal-Ta;%Grade computacional
nPAnalog = length(t);%Total de pontos da grade
fCorte = 3e3; %Frequencia de corte do filtro passa baixas anti aliasing

sinalA = A1*sin(2*pi*f1*t);%Sinal a ser amostrado
sinalB = A2*sin(2*pi*f2*t);%Sinal que fornece o impulso em frequencia
sinalComposto = sinalA + sinalB;

%Emulando um filtro ideal criando uma grade no modelo da computacional
filtroIdeal = zeros(1,nPAnalog);

for i = 1:nPAnalog
    fAtual = (i-1)*(fAnalog/nPAnalog);

    if fAtual <= fCorte || fAtual >= (fAnalog-fCorte)
        filtroIdeal(i) = 1;
    else
        filtroIdeal(i) = 0;
    end
end

sinalCompostoFFT = fft(sinalComposto);%Jogando na frequencia para retirar o espectro da segunda senoide
sinalFiltradoFFT = sinalCompostoFFT .* filtroIdeal;%Filtro no tempo pra retirar o espectro da segunda senoide
sinalFiltrado = real(ifft(sinalFiltradoFFT));%Volto pro dominio do tempo pra ter o meu sinal filtrado e no tempo

figure;

subplot(4,1,1);
plot(t,sinalA,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.01]);

subplot(4,1,2);
plot(t,sinalB,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.01]);

subplot(4,1,3);
plot(t,sinalComposto,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Composto');
xlim([0 0.01]);

subplot(4,1,4);
plot(t,sinalFiltrado,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

kPasso = round(fAnalog / fs);

tremImpulsos = zeros(1,nPAnalog);%Grade computacional do trem de impulsos
tremImpulsos(1:kPasso:end) = 1;%Faço a área dele tender a 1 discretamente

sinalAmostrado = sinalFiltrado .* tremImpulsos;%Amostrando o sinal x(t).s(t)

figure;

subplot(3,1,1);
plot(t,sinalFiltrado,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

subplot(3,1,2);
stem(t,tremImpulsos,'filled');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de Impulsos');
xlim([0 0.01]);

subplot(3,1,3);
stem(t,sinalAmostrado,'filled');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amostragem: x(t) \cdot p(t)');
xlim([0 0.01]);

hold off;

f = (-nPAnalog/2:nPAnalog/2-1)*(fAnalog/nPAnalog);%Criando a grade computacional da frequencia

sinalFiltradoFFT_plot = fftshift(fft(sinalFiltrado))/nPAnalog;
tremImpulsosFFT_plot = (fftshift(fft(tremImpulsos))/nPAnalog) * fAnalog;
sinalAmostradoFFT_plot = (fftshift(fft(sinalAmostrado))/nPAnalog) * fAnalog;

figure;

subplot(3,1,1);
stem(f,abs(sinalFiltradoFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(3,1,2);
stem(f,abs(tremImpulsosFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Trem de Impulsos');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(3,1,3);
stem(f,abs(sinalAmostradoFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro de Magnitude do Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

figure;

subplot(2,1,1);
stem(f,imag(sinalFiltradoFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Filtrado');
xlim([-10000 10000]);
xticks(-10000:2000:10000);

subplot(2,1,2);
stem(f,imag(sinalAmostradoFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

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

sinalReconstruidoFFT = fft(sinalAmostrado) .* filtroReconstrucao * kPasso;
sinalReconstruido = real(ifft(sinalReconstruidoFFT));

filtroReconstrucaoPlot = fftshift(filtroReconstrucao);
sinalReconstruidoFFT_plot = fftshift(sinalReconstruidoFFT)/nPAnalog;

figure;

subplot(4,1,1);
stem(f,imag(sinalAmostradoFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Amostrado');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(4,1,2);
plot(f,filtroReconstrucaoPlot,'LineWidth',2,'Color','r');
grid on;
xlabel('Frequência (Hz)');
ylabel('Ganho');
title('Filtro Ideal de Reconstrução');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);
ylim([-0.2 1.2]);

subplot(4,1,3);
stem(f,imag(sinalReconstruidoFFT_plot),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro do Sinal Reconstruído');
xlim([-5*fs 5*fs]);
xticks(-5*fs : fs : 5*fs);

subplot(4,1,4);
plot(t, sinalFiltrado, 'LineWidth', 2, 'DisplayName', 'Sinal Filtrado');
hold on;
plot(t, sinalReconstruido, '--', 'LineWidth', 2, 'DisplayName', 'Sinal Reconstruído');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Comparação: Sinal Filtrado vs Sinal Reconstruído');
xlim([0 0.01]);
legend('Location', 'northeast');
hold off;
