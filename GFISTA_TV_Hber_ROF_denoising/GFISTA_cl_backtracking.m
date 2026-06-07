function [ u_k2, L1_2, e1_2 ] = GFISTA_cl_backtracking( u0, Fu_star, iter2, mu, mu_f, mu_g, lambda, ep, L_0, rho )
disp(['GFISTA_cl_backtracking L_0=',num2str(L_0),' ...']);
L(1)=L_0;
tau(1)=1/L(1);
%-------------------------------------------------------
q(1)=(mu*tau(1))/(1+tau(1)*mu_g);
x0=D_gradient( u0 );
x1=x0;
%-------------------------------------------------------
t(1)=1;
%------------------------------------------------------
y2=x0;
y_tilde=y2-tau(1)*grad_f(y2,u0);
x2=prox_tau_g( y_tilde,lambda,tau(1),mu_g );
F_AD=F_dual( x2,u0,lambda,ep );
F_u2(1)=F_dual( x2,u0,lambda,ep );
%-------------------------------------------------------
%%
k=1;
tau0(k+1)=tau(k);
i=0;
while (i>=0)
tau(k+1)=(rho^i)*tau0(k+1);
q(k+1)=(mu*tau(k+1))/(1+tau(k+1)*mu_g);
t(k+1)=(1-q(k)*(t(k)^2)+sqrt(((1-q(k)*(t(k)^2))^2)+4*(q(k)/q(k+1))*(t(k)^2)))/2;
beta(k+1)=((t(k)-1)/t(k+1))*((1+tau(k+1)*mu_g-t(k+1)*tau(k+1)*mu)/(1-tau(k+1)*mu_f));
y3=x2+(beta(k+1)*(x2-x1));
y_tilde=y3-tau(k+1)*grad_f(y3,u0);
p=prox_tau_g( y_tilde,lambda,tau(k+1),mu_g );
if ((D_f(p,y3,u0)>(norm_2_p(p-y3)/(2*tau(k+1)))))
    i=i+1;
else
    break
end
end
%----------------------------------------------------------
F_ad=F_dual( p,u0,lambda,ep );
if F_ad<=F_AD
   x3=p;
   F_AD=F_ad;
else
   x3=x2;
   x2=x1;
   t(k+1)=1;
end
%----------------------------------------------------------
F_u2(k+1)=F_dual( x3,u0,lambda,ep );
x1=x2;
x2=x3;
y2=y3;
T=p;
%%
for k=2:iter2
    tau0(k+1)=tau(k);
    i=0;
    while (i>=0)
    tau(k+1)=(rho^i)*tau0(k+1);
    q(k+1)=(mu*tau(k+1))/(1+tau(k+1)*mu_g);
    t(k+1)=(1-q(k)*(t(k)^2)+sqrt(((1-q(k)*(t(k)^2))^2)+4*(q(k)/q(k+1))*(t(k)^2)))/2;
    beta(k+1)=((t(k)-1)/t(k+1))*((1+tau(k+1)*mu_g-t(k+1)*tau(k+1)*mu)/(1-tau(k+1)*mu_f));
    y3=x2+(beta(k+1)*((x2-x1)+((t(k)/(t(k)-1))*(T-x2))));
    y_tilde=y3-tau(k+1)*grad_f(y3,u0);
    p=prox_tau_g( y_tilde,lambda,tau(k+1),mu_g );
    if ((D_f(p,y3,u0)>(norm_2_p(p-y3)/(2*tau(k+1)))))
        i=i+1;
    else
        break
    end
    end
    %------------------------------------------------------
    F_ad=F_dual( p,u0,lambda,ep );
    if F_ad<=F_AD
        x3=p;
        F_AD=F_ad;
    else
        x3=x2;
        x2=x1;
        t(k+1)=1;
    end
    %------------------------------------------------------
    F_u2(k+1)=F_dual( x3,u0,lambda,ep );
    x1=x2;
    x2=x3;
    y2=y3;
    T=p;
end
%%
u_k2=u0-D_adjoint( x3 ); 
%------
for k=1:(iter2+1)
    e1_2(k)=F_u2(k)-Fu_star;
end
%------
for k=1:iter2+1;
    L1_2(k)=1/tau(k);
end
end

