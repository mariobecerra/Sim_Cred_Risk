function [ I_ind, I_cons, I_com, I_serv ] = individual(u, rho_tilde, tam, rho_dif )
% Este programa simula las variables binarias de impago en un modelos de
% pérdida agregada en riesgo de crédito con varios grupos homogéneos de
% riesgo
%
% In: 
%       u, vector que contiene las u's calibradas
%       rho_tilde, vector que contiene las rho_tildes's calibradas
%       tam, vector que contiene el número total de créditos de cada grupo
%       rho_dif, correlación extragrupo
% Out: 
%       I_ind, vector de simulaciones para el grupo industrial
%       I_cons, vector de simulaciones para el grupo construcción
%       I_com, vector de simulaciones para el grupo comercio
%       I_serv, vector de simulaciones para el grupo servicios
% 
% -------------------------------------------------------------------------

%calcular las xi
m = 4;
Z0=randn;
Z=randn(1,m);
for k=1:m
    x=zeros(tam(k),1);
    for l=1:tam(k)
        eps=randn;
        x(l)=sqrt(rho_dif)*Z0+sqrt(rho_tilde(k)-rho_dif)*Z(k)+sqrt(1- rho_tilde(k))*eps;
    end
    % ---------------------------------------------------------------------
    switch k
        case 1
            I_ind = x < u(k);
        case 2
            I_cons = x < u(k);
        case 3
            I_com = x < u(k);
        case 4
            I_serv = x < u(k);
    end
end












