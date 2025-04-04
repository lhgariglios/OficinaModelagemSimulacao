% Luiz Henrique Gariglio dos Santos - 2022421137

clear
close all
clc
%% 1

t = linspace(0, 10, 1000);

x = -10*exp(-t) + cos(5*t);  % Resposta total
x_t = -10*exp(-t);           % Resposta transitória
x_p = cos(5*t);              % Resposta permanente

figure;
hold on;
plot(t, x, 'b');
plot(t, x_t, 'r');
plot(t, x_p, 'g');
xlabel('Tempo');
ylabel('Amplitude');
legend('Resposta Total (x(t))','Resposta Transitória (-10e^{-t})','Resposta Permanente (cos(5t))');
grid on;
hold off;

%% 2

syms m c k s
G = 1/(m*s^2 + c*s + k);
pretty(G)

m = 1;
c = 2;
k = 4;

num = 1;                     
den = [m c k];               
G = tf(num, den);           

figure;
step(G);
xlabel('Tempo');
ylabel('Deslocamento');
grid on;

%% 3

syms m c k s
G = (c*s + k) / (m*s^2 + c*s + k);
pretty(G)

m = 1;
c = 2;
k = 4;
num = [c k];        
den = [m c k];     
G = tf(num, den);

t = 0:0.01:20;
y = 2 * sin(t);
x = lsim(G, y, t);

figure;
plot(t, x, 'b');
hold on;
plot(t, y, 'r');
legend('resposta','excitação');
xlabel('Tempo (s)');
ylabel('Deslocamento');
hold off;

%% 4
t = linspace(0, 1, 100);

f = zeros(size(t));
f(t <= 0.5) = 2 * t(t <= 0.5);
f(t > 0.5) = 0;

figure;
plot(t, f);
xlabel('Tempo (s)');
ylabel('f(t)');

%% 5
t = 0:0.1:3;

f = zeros(size(t));
f(t <= 1) = 2 * t(t <= 1);
f(t > 1 & t <= 2) = 4 - 2 * t(t > 1 & t <= 2);

figure;
plot(t, f);
xlabel('Tempo (s)');
ylabel('f(t)');

%% 6

syms t s
f = -4*t^2+4*t;

F = laplace(f, t, s);
pretty(F);

t_var = linspace(0, 1, 100);
f_var = -4*t_var.^2+4*t_var;

figure;
plot(t_var, f_var);
xlabel('tempo [s]');
ylabel('f(t)');

%% 7

t = linspace(-1,5,500);
f = zeros(500);

for i=1:500
    if t(i) > 0 && t(i) <= 1
        f(i) = t(i)* 10;
    elseif t(i) > 1 && t(i) <= 2
        f(i) = 10;
    elseif t(i) > 2 && t(i) < 3
        f(i) = 10 - (10*(t(i)-2));
    else
        f(i) = 0;
    end
end

plot(t,f);