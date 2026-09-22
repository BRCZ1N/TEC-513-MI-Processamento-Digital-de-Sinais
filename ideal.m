clear;
clc;
close all;

pkg load signal;

f1 = 2e3;%Frequencia da senoide a ser amostrada
f2 = 5e3;%Frequencia da senoide do impulso na frequencia
A1 = 1.0;%Amplitude da senoide a ser amostrada
A2 = 1.2;%Amplitude da senoide do impulso na frequencia

fs = 4.5e3; %Frequencia de amostragem
Ts = 1/fs; %Periodo de amostragem
fatorCompassagem = 100; %Fator de compassagem pra permitir escolher qualquer frequencia de amostragem e criar a grade computacional baseada nela
fAnalog = fatorCompassagem * fs; %Frequencia da grade computacional
Ta = 1/fAnalog;%Frequencia da grade computacional
tempoTotal = 1;%Tempo total de simulação
t = 0:Ta:tempoTotal-Ta;%Grade computacional
nPAnalog = length(t);%Total de pontos da grade
fCorte = 3e3; %Frequencia de corte do filtro passa baixas anti aliasing

sinalA = A1*sin(2*pi*f1*t);%Sinal a ser amostrado
sinalB = A2*sin(2*pi*f2*t);%Sinal que fornece o impulso em frequencia
sinalComposto = sinalA + sinalB;

figure('Position',[100 100 1200 750]);

subplot(3,1,1);
plot(t,sinalA,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal A');
xlim([0 0.01]);

subplot(3,1,2);
plot(t,sinalB,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal B');
xlim([0 0.01]);

subplot(3,1,3);
plot(t,sinalComposto,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Composto');
xlim([0 0.01]);

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

figure('Position',[100 100 1200 600]);

subplot(2,1,1);
plot(t,sinalComposto,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Composto Antes da Filtragem');
xlim([0 0.01]);

subplot(2,1,2);
plot(t,sinalFiltrado,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Após Filtragem');
xlim([0 0.01]);

kPasso = round(Ts/Ta);

tremImpulsos = zeros(1,nPAnalog);%Grade computacional do trem de impulsos
tremImpulsos(1:kPasso:end) = 1;%Faço a área dele tender a 1

sinalAmostrado = sinalFiltrado .* (tremImpulsos);%Amostrando o sinal x(t).s(t)

figure('Position',[100 100 1200 750]);

subplot(4,1,1);
plot(t,sinalFiltrado,'LineWidth',1.2);
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Sinal Filtrado');
xlim([0 0.01]);

subplot(4,1,2);
stem(t,tremImpulsos,'filled');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Trem de Impulsos');
xlim([0 0.01]);

subplot(4,1,3);
stem(t,sinalAmostrado,'filled');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amostragem: x(t) \cdot p(t)');
xlim([0 0.01]);

subplot(4,1,4);
plot(t,sinalFiltrado,'LineWidth',1.2);
hold on;
stem(t,sinalAmostrado,'filled');
grid on;
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amostragem Sobre o Sinal Original');
xlim([0 0.01]);
legend('Sinal Filtrado','Amostras');
hold off;

f = (-nPAnalog/2:nPAnalog/2-1)*(fAnalog/nPAnalog);%Criando a grade computacional da frequencia

sinalFiltradoFFT = fftshift(fft(sinalFiltrado))/nPAnalog;%Espectro de frequencia do sinal
tremImpulsosFFT = fftshift(fft(tremImpulsos))/nPAnalog;%Espectro de frequencia do trem
sinalAmostradoFFT = fftshift(fft(sinalAmostrado))/nPAnalog;%Espectro da amostragem

figure('Position',[100 100 1200 750]);

subplot(3,1,1);
plot(f,imag(sinalFiltradoFFT),'LineWidth',1.2);
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Sinal Filtrado');
xlim([-10000 10000]);

subplot(3,1,2);
stem(f,abs(tremImpulsosFFT),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Trem de Impulsos');
xlim([-25000 25000]);
xticks(-25000:5000:25000);

subplot(3,1,3);
stem(f,imag(sinalAmostradoFFT),'filled');
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Espectro do Sinal Amostrado');
xlim([-25000 25000]);
xticks(-25000:5000:25000);

