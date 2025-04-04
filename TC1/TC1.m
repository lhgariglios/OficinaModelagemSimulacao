% Luiz Henrique Gariglio dos Santos - 2022421137

clear
close all
clc
%% 1
disp('Questão 1:');
z1 = 2-1j*3;

r = 3;
theta = pi/4;
[a,b] = pol2cart(r,theta);

z2 = a+1j*b;

figure;
compass(z1);
hold on;
compass(z2);
title('Questão 1');
legend('z1','z2');

z_mult = z1*z2;
fprintf('z1*z2 = %.3f + %.3fi\n', real(z_mult), imag(z_mult));
z_div = z1/z2;
fprintf('z1/z2 = %.3f + %.3fi\n', real(z_div), imag(z_div));

%% 2
omega = 2 * pi * 3; % frequencia = 3;
alpha = log(0.5)/2; % decaimento de 50% a cada 2s
t = linspace(0, 5, 1000);
y = exp(alpha*t).*sin(omega*t);

figure;
plot(t,y);
title('Questão 2');
xlabel('Tempo[s]');
ylabel('Amplitude');

%% 3
t = linspace(0, 3, 1000);
x1 = real(2*exp((-1+1j*2*pi)*t));
x2 = imag(3-exp((1-1j*2*pi)*t));
x3 = 3-imag(exp((1-1j*2*pi)*t));
figure
plot(t,x1,'b');
hold on;
plot(t,x2,'r');
plot(t,x3,'g');
legend('x1','x2','x3');
title('Questão 3');

%% 4
t = linspace(0, 5, 1000);
y = cos(t).*sin(20*t);
figure;
plot(t,y);
title('Questão 4');
xlabel('Tempo[s]');
ylabel('Amplitude');

%% 5

N = 6*[1 1];
D = conv([1 0], conv([1 1.46], [1 0.13]));
[R1,P1,K1] = residue(N,D);


disp('Questão 5.a:');
for i = 1:length(R1)
    fprintf('%.4f/(s - (%.4f)) + ', R1(i), P1(i));
end
if ~isempty(K1)
    for i = 1:length(K1)
        fprintf('%.4f*s^%d + ', K1(i), length(K1)-i);
    end
end
fprintf('\b\b \n'); % Remove o último ' + '

N = [1 2 3];
D = [1 2 1];
[R2,P2,K2] = residue(N,D);

disp('Questão 5.b:');
for i = 1:length(R2)
    fprintf('%.4f/(s - (%.4f)) + ', R2(i), P2(i));
end
if ~isempty(K2)
    for i = 1:length(K2)
        fprintf('%.4f*s^%d + ', K2(i), length(K2)-i);
    end
end
fprintf('\b\b \n'); % Remove o último ' + '

N = 6*[1 34];
D = [1 10 34 0];
[R3,P3,K3] = residue(N,D);

disp('Questão 5.c:');
for i = 1:length(R3)
    fprintf('%.4f/(s - (%.4f)) + ', R3(i), P3(i));
end
if ~isempty(K3)
    for i = 1:length(K3)
        fprintf('%.4f*s^%d + ', K3(i), length(K3)-i);
    end
end
fprintf('\b\b \n'); % Remove o último ' + '

%% 6

t = linspace(0, 5, 1000);
x1 = cos(t).*sin(20*t);
x2 = cos(t);
x3 = sin(20*t);

figure;
plot(t,x1,'b');
hold on;
plot(t,x2,'r');
plot(t,x3,'g');
legend('x1','x2','x3');
title('Questão 6');
xlabel('Tempo[s]');
ylabel('Amplitude');

%% 7

t = linspace(0, 5, 1000);
w = 5;
f = -3*cos(w*t)+4*sin(w*t);
A = sqrt((-3)^2+4^2);
phi = atan2(4,-3);
f_eq = A*cos(w*t+phi);
figure;
plot(t,f,'k');
title('Questão 7');
xlabel('Tempo[s]');
ylabel('Amplitude');
figure;
plot(t,f_eq,'r');
title('Questão 7 - função equivalente');
xlabel('Tempo[s]');
ylabel('Amplitude');

%% 8
disp('Questão 8');
A = [1 1 6; 5 -2 1; -8 2 -3];
B = [2 9; -5 -1; 9 2];

disp('A');
if (size(A,1) == size(A, 2))
    disp('É uma matriz quadrada');
else
    disp('Não é uma matriz quadrada');
end

disp('Elementos com valor 2:');
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if(A(i,j) == 2)
            a = [i,j];
            disp(a);
        end
    end
end

disp('Elementos com valores negativos:');
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if(A(i,j) < 0)
            a = [i,j];
            disp(a);
        end
    end
end

disp('B');
if (size(B,1) == size(B, 2))
    disp('É uma matriz quadrada');
else
    disp('Não é uma matriz quadrada');
end

disp('Elementos com valor 2:');
for i = 1:size(B,1)
    for j = 1:size(B,2)
        if(B(i,j) == 2)
            a = [i,j];
            disp(a);
        end
    end
end

disp('Elementos com valores negativos:');
for i = 1:size(B,1)
    for j = 1:size(B,2)
        if(B(i,j) < 0)
            a = [i,j];
            disp(a);
        end
    end
end

%% 9

syms a b c d;
M = [a b; c d];
d = det(M);
disp('Questão 9');
fprintf('Determinante: %s\n', d);
I = inv(M);
disp('Inversa:');
disp(I);
t = trace(M);
fprintf('Traço: %s\n', t);

%% 10

disp('Questão 10');
syms s;
f = s^4 + 5*s^3 - s^2 + 3*s + 2;
df = diff(f,s);
df2 = diff(f,s,2);
disp('Primeira derivada:');
pretty(df);
disp('Segunda derivada:');
pretty(df2);

%% 11

syms a b c d x y;
M = [a*x b*x^2; c*x^3 d*y];
disp('Questão 11');
df = diff(M,x);
pretty(df);

hold off;

%% 12

syms x;
p = (x^2-1)*(x-2)*(x-3);
ex = expand(p); % Expande expressões simbólicas. Aplica propriedades distributivas e produtos notáveis
fac = factor(p); % Fatora expressões simbólicas
disp('Questão 12');
pretty(ex);
pretty(fac);