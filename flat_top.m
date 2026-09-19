clear;
clc;
close all;

pkg load signal;


% ============================================================
% 1. FREQUÊNCIAS E PARÂMETROS DO SINAL
% ============================================================

f1 = 1e3;                   % Senoide desejada = 1 kHz
f2 = 5e3;                   % Segunda senoide = 5 kHz

A1 = 1.0;                   % Amplitude da senoide 1
A2 = 1.2;                   % Amplitude da senoide 2


% ============================================================
% 2. SIMULAÇÃO DO SINAL ANALÓGICO
% ============================================================

fAnalog = 1e5;              % Frequência usada para simular o analógico
Ta = 1/fAnalog;             % Período da simulação analógica

tempoTotal = 1;             % Tempo total da simulação
t = 0:Ta:tempoTotal-Ta;

nPAnalog = length(t);


% ============================================================
% 3. FREQUÊNCIA DE AMOSTRAGEM
% ============================================================

fs = 10e3;                  % Frequência de amostragem
Ts = 1/fs;                  % Período de amostragem

fNyquist = fs/2;

fCorte = 2e3;


% ============================================================
% 4. SENOIDE A
% ============================================================

sinalA = A1*sin(2*pi*f1*t);


figure;

subplot(5,3,1);

plot(t, sinalA);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Sinal A - 1 kHz');

xlim([0 0.01]);


% ============================================================
% 5. SENOIDE B
% ============================================================

sinalB = A2*sin(2*pi*f2*t);


subplot(5,3,2);

plot(t, sinalB);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Sinal B - 5 kHz');

xlim([0 0.01]);


% ============================================================
% 6. SINAL COMPOSTO
% ============================================================

sinalComposto = sinalA + sinalB;


subplot(5,3,3);

plot(t, sinalComposto);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Sinal Composto');

xlim([0 0.01]);


% ============================================================
% 7. FILTRO IDEAL PASSA-BAIXA
% ============================================================

filtroIdeal = zeros(1, nPAnalog);

for i = 1:nPAnalog

    fAtual = (i - 1) * (fAnalog / nPAnalog);

    % Mantém:
    % frequências positivas abaixo do corte
    % e suas correspondentes negativas

    if fAtual <= fCorte || ...
       fAtual >= (fAnalog - fCorte)

        filtroIdeal(i) = 1;

    else

        filtroIdeal(i) = 0;

    end

end


% ============================================================
% 8. SINAL FILTRADO
% ============================================================

sinalCompostoFFT = fft(sinalComposto);

sinalFiltradoFFT = ...
    sinalCompostoFFT .* filtroIdeal;

sinalFiltrado = ...
    real(ifft(sinalFiltradoFFT));


subplot(5,3,4);

plot(t, sinalFiltrado);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Sinal Filtrado - 1 kHz');

xlim([0 0.01]);


% ============================================================
% 9. TREM DE IMPULSOS - AMOSTRAGEM IDEAL
% ============================================================

kPasso = round(Ts/Ta);

tremImpulsos = zeros(1, nPAnalog);

tremImpulsos(1:kPasso:end) = 1/Ta;


subplot(5,3,5);

stem(t, tremImpulsos*Ta);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Trem de Impulsos');

xlim([0 0.01]);


% ============================================================
% 10. AMOSTRAGEM IDEAL
% ============================================================

sinalAmostradoIdeal = ...
    sinalFiltrado .* (tremImpulsos*Ta);


subplot(5,3,6);

stem(t, sinalAmostradoIdeal);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Amostragem Ideal');

xlim([0 0.01]);


% ============================================================
% 11. ESPECTRO DO SINAL FILTRADO
% ============================================================

f = (-nPAnalog/2 : nPAnalog/2 - 1) ...
    * (fAnalog/nPAnalog);

sinalFiltradoFFT = ...
    fftshift(fft(sinalFiltrado)) / nPAnalog;


subplot(5,3,7);

stem(f, abs(sinalFiltradoFFT));

grid on;

xlabel('Frequência (Hz)');
ylabel('Magnitude');

title('Espectro do Sinal Filtrado');

xlim([-5000 5000]);

xticks(-5000:1000:5000);


% ============================================================
% 12. ESPECTRO DO TREM DE IMPULSOS
% ============================================================

tremImpulsosFFT = ...
    fftshift(fft(tremImpulsos)) / nPAnalog;


subplot(5,3,8);

stem(f, abs(tremImpulsosFFT));

grid on;

xlabel('Frequência (Hz)');
ylabel('Magnitude');

title('Espectro do Trem de Impulsos');

xlim([-20000 20000]);

xticks(-20000:5000:20000);


% ============================================================
% 13. ESPECTRO DA AMOSTRAGEM IDEAL
% ============================================================

sinalAmostradoIdealFFT = ...
    fftshift(fft(sinalAmostradoIdeal)) / nPAnalog;


subplot(5,3,9);

stem(f, abs(sinalAmostradoIdealFFT));

grid on;

xlabel('Frequência (Hz)');
ylabel('Magnitude');

title('Espectro da Amostragem Ideal');

xlim([-35000 35000]);

xticks(-30000:5000:30000);


% ============================================================
% ============================================================
%                  FLAT-TOP
% ============================================================
% ============================================================


% ============================================================
% 14. LARGURA DO TOPO PLANO
% ============================================================

% Vamos utilizar um pulso com duração de 1/3 de Ts.

d = Ts/3;

% Número de pontos da simulação analógica correspondentes
% à largura do pulso

numPontosPulso = round(d/Ta);

% Largura REAL obtida após a discretização

dReal = numPontosPulso * Ta;


fprintf('\n========================================\n');
fprintf('PARAMETROS DA AMOSTRAGEM FLAT-TOP\n');
fprintf('========================================\n');

fprintf('fs = %.2f Hz\n', fs);
fprintf('Ts = %.6e s\n', Ts);

fprintf('d ideal = %.6e s\n', d);

fprintf('d real  = %.6e s\n', dReal);

fprintf('Pontos por periodo Ts = %d\n', kPasso);

fprintf('Pontos do pulso = %d\n', numPontosPulso);

fprintf('Duty Cycle real = %.2f %%\n', ...
        100*dReal/Ts);


% ============================================================
% 15. OBTENDO AS AMOSTRAS
% ============================================================

% Os valores das amostras são exatamente aqueles obtidos
% nos instantes nTs.

amostras = sinalFiltrado(1:kPasso:end);


% ============================================================
% 16. CONSTRUÇÃO DO PULSO RETANGULAR
% ============================================================

% Durante d:
% valor = 1
%
% Após d até completar Ts:
% valor = 0

pulsoFlatTop = [ ...
    ones(1, numPontosPulso), ...
    zeros(1, kPasso - numPontosPulso) ...
];


% ============================================================
% 17. CONSTRUÇÃO DO SINAL FLAT-TOP
% ============================================================

% Cada amostra é mantida constante durante d.

sinalAmostradoFlatTop = ...
    kron(amostras, pulsoFlatTop);


% Ajuste para garantir que o vetor tenha
% exatamente o mesmo tamanho do sinal original.

sinalAmostradoFlatTop = ...
    sinalAmostradoFlatTop(1:nPAnalog);


% ============================================================
% 18. GRÁFICO DO FLAT-TOP
% ============================================================

subplot(5,3,10);

plot(t, sinalFiltrado, 'b--', 'LineWidth', 1);

hold on;

stairs(t, ...
       sinalAmostradoFlatTop, ...
       'r', ...
       'LineWidth', 1.5);

grid on;

xlabel('Tempo (s)');
ylabel('Amplitude');

title('Amostragem Flat-Top');

legend('Sinal Filtrado', ...
       'Flat-Top');

xlim([0 0.01]);


% ============================================================
% 19. MOSTRAR AS AMOSTRAS NUMERICAMENTE
% ============================================================

temposAmostra = t(1:kPasso:end);

% Mostrar somente as primeiras 10 amostras

numMostrar = min(10, length(amostras));

for i = 1:numMostrar

    if amostras(i) >= 0

        deslocamento = 0.15;

    else

        deslocamento = -0.15;

    end

    text( ...
        temposAmostra(i), ...
        amostras(i) + deslocamento, ...
        sprintf('%.2f', amostras(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 8, ...
        'FontWeight', 'bold');

end


% ============================================================
% 20. ESPECTRO DO SINAL FLAT-TOP
% ============================================================

sinalAmostradoFlatTopFFT = ...
    fftshift(fft(sinalAmostradoFlatTop)) / nPAnalog;


subplot(5,3,11);

stem(f, ...
     imag(sinalAmostradoFlatTopFFT));

grid on;

xlabel('Frequência (Hz)');
ylabel('Magnitude');

title('Espectro do Sinal Flat-Top');

xlim([-35000 35000]);

xticks(-30000:5000:30000);

% ============================================================
% 21. MAGNITUDE DO SINAL FLAT-TOP E RÉPLICAS (COM VALORES)
% ============================================================

subplot(5,3,12);

% 1. Plota a magnitude (abs) das réplicas
espectroMag = abs(sinalAmostradoFlatTopFFT);
stem(f, espectroMag, 'r');
hold on;

% 2. Cálculo da envolvente teórica (Função Sinc)
f_sinc = f;
f_sinc(f_sinc == 0) = 1e-10; % Evita divisão por zero
envolventeSinc = abs(sin(pi * f_sinc * dReal) ./ (pi * f_sinc * dReal));

% Ajusta a escala da envolvente para coincidir com a amplitude máxima
amplitudeMaxima = max(espectroMag);
envolventeSinc = envolventeSinc * amplitudeMaxima;

% Plota a envolvente
plot(f, envolventeSinc, 'k--', 'LineWidth', 1.5);

% 3. PLOTAR OS VALORES NUMÉRICOS EM CIMA DOS IMPULSOS
% Define um limiar (5% do máximo) para ignorar ruídos no chão do espectro
limiar = amplitudeMaxima * 0.05;

for i = 2:length(espectroMag)-1
    % Condição: O ponto atual tem de ser maior que o anterior, maior que o seguinte e maior que o limiar
    if espectroMag(i) > limiar && espectroMag(i) > espectroMag(i-1) && espectroMag(i) > espectroMag(i+1)

        % Escreve o valor do pico ligeiramente acima dele
        deslocamento = amplitudeMaxima * 0.08; % Espaço dinâmico para não colar na linha

        text(f(i), espectroMag(i) + deslocamento, ...
             sprintf('%.3f', espectroMag(i)), ... % %.3f para 3 casas decimais
             'HorizontalAlignment', 'center', ...
             'FontSize', 7, ...
             'FontWeight', 'bold', ...
             'Color', 'b'); % Cor azul para fácil leitura
    end
end

grid on;

xlabel('Frequência (Hz)');
ylabel('Magnitude');

title('Réplicas no Espectro com Valores Demarcados');

xlim([-35000 35000]);
xticks(-30000:5000:30000);

% Ajuste no limite do eixo Y para garantir que o texto não é cortado em cima
ylim([0 amplitudeMaxima * 1.2]);

legend('Réplicas Flat-Top', 'Envolvente Sinc', 'Location', 'best');
% ============================================================
% 23. INFORMAÇÕES DAS AMOSTRAS
% ============================================================

fprintf('\n========================================\n');
fprintf('PRIMEIRAS AMOSTRAS\n');
fprintf('========================================\n');

numMostrar = min(10, length(amostras));

for i = 1:numMostrar

    fprintf( ...
        'n = %d | t = %.6f ms | x[n] = %.4f\n', ...
        i-1, ...
        temposAmostra(i)*1000, ...
        amostras(i));

end


fprintf('\n========================================\n');
fprintf('FIM DA SIMULACAO\n');
fprintf('========================================\n');
