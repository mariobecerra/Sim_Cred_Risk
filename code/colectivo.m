function [ I_ind, I_cons, I_com, I_serv ] = colectivo(m, u, rho_tilde, tam, rho_dif )

%calcular las xi
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












