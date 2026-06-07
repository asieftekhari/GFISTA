function [ u_k1, L1_1, e1_1 ] = GFISTA_modified( u0, Fu_star, iter2, L_f, mu, mu_f, mu_g, lambda, ep )
disp(['GFISTA_modified L_f=', num2str(L_f), ' ...']);
tau=1/L_f;
x0=D_gradient( u0 );
x1=x0;
q=((tau*mu)/(1+tau*mu_g));
t(1)=1;
%%
k=1;
t(k+1)=(1-q*(t(k)^2)+sqrt(((1-q*(t(k)^2))^2)+4*(t(k)^2)))/2;
beta(k)=((t(k)-1)/t(k+1))*((1+tau*mu_g-t(k+1)*tau*mu)/(1-tau*mu_f));
%-----------------------------------------------------
y=x1+(beta(k)*(x1-x0));
y_tilde=y-tau*grad_f(y,u0);
%-----------------------------------------------------
p=prox_tau_g( y_tilde,lambda,tau,mu_g );
x2=p;
F_AD=F_dual( p,u0,lambda,ep );
%-----------------------------------------------------
x0=x1;
x1=x2;
F_u1(k)=F_dual( x2,u0,lambda,ep );
%%
for k=2:iter2+1
    t(k+1)=(1-q*(t(k)^2)+sqrt(((1-q*(t(k)^2))^2)+4*(t(k)^2)))/2;
    beta(k)=((t(k)-1)/t(k+1))*((1+tau*mu_g-t(k+1)*tau*mu)/(1-tau*mu_f));
    y=x1+beta(k)*((x1-x0)+(t(k)/(t(k)-1))*(p-x1));
    y_tilde=y-tau*grad_f(y,u0);
    %-----------------------------------------------------------------
    p=prox_tau_g( y_tilde,lambda,tau,mu_g );
    %-----------------------------------------------------------------
    F_ad=F_dual( p,u0,lambda,ep );
    if ( F_ad<=F_AD )
        x2=p;
        F_AD=F_ad;
    else
        x2=x1;
    end
    %-----------------------------------------------------------------
    F_u1(k)=F_dual( x2,u0,lambda,ep );
    x0=x1;
    x1=x2;
end
%%
a=Fu_star;
u_k1=u0-D_adjoint( x2 );
for i=1:(iter2+1);
    e1_1(i)=F_u1(i)-Fu_star;
end
%%
for k=0:iter2
    L1_1(k+1)=1/tau;
end
end

