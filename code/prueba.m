%% Declaración de variables
clear
load('../data/matlab.mat')

rng(48151623);

% Tamaño de la simulación
tam = zeros(4,1);
tam(1) = length(industrial);
tam(2) = length(construccion);
tam(3) = length(comercio);
tam(4) = length(servicios);

z = norminv(0.975);
n_sim = zeros(4,1);
n_sim(1) = ceil((z*std(industrial)/10)^2);
n_sim(2) = ceil((z*std(construccion)/10)^2);
n_sim(3) = ceil((z*std(comercio)/10)^2);
n_sim(4) = ceil((z*std(servicios)/10)^2);

m = 4;
p = [0.007,0.009,0.0065,0.0060]';
rho_dif = 0.005;
r0 = .5;
rho = [.0009,.0004,.0005,.0007]';

% calibrar u
for i=1:m
    u(i) = norminv(p(i));
end

% calibrar ro
for j = 1:m
    fun = @(r) (mvncdf([u(j); u(j)],[0; 0],[1 r; r 1]) - (p(j)^2 + rho(j)*p(j)*(1-p(j))));
    rho_tilde(j) = fzero(fun,r0);
end
%% Individual 
L_ind = zeros(n_sim(1),1);
L_cons = zeros(n_sim(2),1);
L_com = zeros(n_sim(3),1);
L_serv = zeros(n_sim(4),1);

tic;
for i = 1:max(n_sim)
    [ I_ind, I_cons, I_com, I_serv ] = individual(u,rho_tilde, tam, rho_dif);
    L_ind(i) = I_ind'*industrial;
    L_cons(i) = I_cons'*construccion;
    L_com(i) = I_com'*comercio;
    L_serv(i) = I_serv'*servicios;
end
K = [sum(I_ind);sum(I_cons);sum(I_com);sum(I_serv)];
L_ind = L_ind(1:n_sim(1));
L_cons = L_cons(1:n_sim(2));
L_com = L_com(1:n_sim(3));
L_serv = L_serv(1:n_sim(4));
L = sum(L_ind) + sum(L_cons) + sum(L_com) + sum(L_serv)

% ------ PLOTS ------

h = figure(1)
subplot(2,2,1)
hist(L_ind, 20)
xlim([0,max(L_ind)])
title('Industrial')

subplot(2,2,2)
hist(L_cons, 20)
xlim([0,max(L_cons)])
title('Construcción')

subplot(2,2,3)
hist(L_com, 20)
xlim([0,max(L_com)])
title('Comercio')

subplot(2,2,4)
hist(L_serv, 20)
xlim([0,max(L_serv)])
title('Servicios')

saveas(h,'../output/histogramas_individual.jpg')

h = figure(2)
subplot(2,2,1)
boxplot(L_ind)
title('Industrial')

subplot(2,2,2)
boxplot(L_cons)
title('Construcción')

subplot(2,2,3)
boxplot(L_com)
title('Comercio')

subplot(2,2,4)
boxplot(L_serv)
title('Servicios')

saveas(h,'../output/boxplots_individual.jpg')
close all;

% ------ ESTADÍSTICOS ------
mu(1) = mean(L_ind); mu(2) = mean(L_cons); mu(3) = mean(L_com); mu(4) = mean(L_serv);
vari(1) = var(L_ind); vari(2) = var(L_cons); vari(3) = var(L_com); vari(4) = var(L_serv);
k(1) = kurtosis(L_ind); k(2) = kurtosis(L_cons); k(3) = kurtosis(L_com); k(4) = kurtosis(L_serv);
sesgo(1) = skewness(L_ind); sesgo(2) = skewness(L_cons); sesgo(3) = skewness(L_com); sesgo(4) = skewness(L_serv);
quart_ind = quantile(L_ind,[0.25, 0.5, 0.75, 1]);
quart_cons = quantile(L_cons,[0.25, 0.5, 0.75, 1]);
quart_com = quantile(L_com,[0.25, 0.5, 0.75, 1]);
quart_serv = quantile(L_serv,[0.25, 0.5, 0.75, 1]);
quart = [quart_ind; quart_cons; quart_com; quart_serv];

mu_gen = mean([L_ind;L_cons;L_com;L_serv])
var_gen = var([L_ind;L_cons;L_com;L_serv])
k_gen = kurtosis([L_ind;L_cons;L_com;L_serv])
sesgo = skewness([L_ind;L_cons;L_com;L_serv])

t = toc

%% Colectivo

Z0 = randn;
Z = randn(m,1);
p2 = zeros(m,1);
for i = 1:m
    z = (u(i) - sqrt(rho_dif)*Z0 - sqrt(rho_tilde(i)-rho_dif)*Z(i))/(sqrt(1-rho_tilde(i)));
    p2(i) = normcdf(z);
end
p2

L_ind = zeros(n_sim(1),1);
L_cons = zeros(n_sim(2),1);
L_com = zeros(n_sim(3),1);
L_serv = zeros(n_sim(4),1);
tic;
for i=1:n_sim(1)
    K = binornd(tam(1),p2(1));
    f_ind = randsample(industrial,K);
    L_ind(i) = sum(f_ind);
end

for i=1:n_sim(2)
    K = binornd(tam(2),p2(2));
    f_cons = randsample(construccion,K);
    L_cons(i) = sum(f_cons);
end

for i=1:n_sim(3)
    K = binornd(tam(3),p2(3));
    f_com = randsample(comercio,K);
    L_com(i) = sum(f_com);
end

for i=1:n_sim(4)
    K = binornd(tam(4),p2(4));
    f_serv = randsample(servicios,K);
    L_serv(i) = sum(f_serv);
end
L = sum(L_ind) + sum(L_cons) + sum(L_com) + sum(L_serv)


%------ PLOTS ------

h = figure(1)
subplot(2,2,1)
hist(L_ind, 20)
xlim([0,max(L_ind)])
title('Industrial')

subplot(2,2,2)
hist(L_cons, 20)
xlim([0,max(L_cons)])
title('Construcción')

subplot(2,2,3)
hist(L_com, 20)
xlim([0,max(L_com)])
title('Comercio')

subplot(2,2,4)
hist(L_serv, 20)
xlim([0,max(L_serv)])
title('Servicios')

saveas(h,'../output/histogramas_colectivo.jpg')

h = figure(2)
subplot(2,2,1)
boxplot(L_ind)
title('Industrial')

subplot(2,2,2)
boxplot(L_cons)
title('Construcción')

subplot(2,2,3)
boxplot(L_com)
title('Comercio')

subplot(2,2,4)
boxplot(L_serv)
title('Servicios')

saveas(h,'../output/boxplots_colectivo.jpg')
close all;

% ------ ESTADÍSTICOS ------
mu(1) = mean(L_ind); mu(2) = mean(L_cons); mu(3) = mean(L_com); mu(4) = mean(L_serv);
vari(1) = var(L_ind); vari(2) = var(L_cons); vari(3) = var(L_com); vari(4) = var(L_serv);
k(1) = kurtosis(L_ind); k(2) = kurtosis(L_cons); k(3) = kurtosis(L_com); k(4) = kurtosis(L_serv);
sesgo(1) = skewness(L_ind); sesgo(2) = skewness(L_cons); sesgo(3) = skewness(L_com); sesgo(4) = skewness(L_serv);
quart_ind = quantile(L_ind,[0.25, 0.5, 0.75, 1]);
quart_cons = quantile(L_cons,[0.25, 0.5, 0.75, 1]);
quart_com = quantile(L_com,[0.25, 0.5, 0.75, 1]);
quart_serv = quantile(L_serv,[0.25, 0.5, 0.75, 1]);
quart = [quart_ind; quart_cons; quart_com; quart_serv];
t = toc

mu_gen = mean([L_ind;L_cons;L_com;L_serv])
var_gen = var([L_ind;L_cons;L_com;L_serv])
k_gen = kurtosis([L_ind;L_cons;L_com;L_serv])
sesgo = skewness([L_ind;L_cons;L_com;L_serv])

%% Variables antitéticas
L_ind = zeros(n_sim(1),1);
L_cons = zeros(n_sim(2),1);
L_com = zeros(n_sim(3),1);
L_serv = zeros(n_sim(4),1);

tic;
for i = 1:max(n_sim)
    [ I_ind, I_cons, I_com, I_serv ] = individual_antitetico(u,rho_tilde, tam, rho_dif);
    L_ind(i) = I_ind'*industrial;
    L_cons(i) = I_cons'*construccion;
    L_com(i) = I_com'*comercio;
    L_serv(i) = I_serv'*servicios;
end
K = [sum(I_ind);sum(I_cons);sum(I_com);sum(I_serv)];
L_ind = L_ind(1:n_sim(1));
L_cons = L_cons(1:n_sim(2));
L_com = L_com(1:n_sim(3));
L_serv = L_serv(1:n_sim(4));
L = sum(L_ind) + sum(L_cons) + sum(L_com) + sum(L_serv)

% ------ PLOTS ------


h = figure(1)
subplot(2,2,1)
hist(L_ind, 20)
xlim([0,max(L_ind)])
title('Industrial')

subplot(2,2,2)
hist(L_cons, 20)
xlim([0,max(L_cons)])
title('Construcción')

subplot(2,2,3)
hist(L_com, 20)
xlim([0,max(L_com)])
title('Comercio')

subplot(2,2,4)
hist(L_serv, 20)
xlim([0,max(L_serv)])
title('Servicios')

saveas(h,'../output/histogramas_individual_antitetico.jpg')

h = figure(2)
subplot(2,2,1)
boxplot(L_ind)
title('Industrial')

subplot(2,2,2)
boxplot(L_cons)
title('Construcción')

subplot(2,2,3)
boxplot(L_com)
title('Comercio')

subplot(2,2,4)
boxplot(L_serv)
title('Servicios')

saveas(h,'../output/boxplots_individual_antitetico.jpg')
close all;



