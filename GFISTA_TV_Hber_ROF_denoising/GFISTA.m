function [ u_star, F_u, improvement ] = GFISTA( u, u0, iter1, L_f, mu, mu_f, mu_g, lambda, ep )
disp(['GFISTA L_f=', num2str(L_f), ' ...']);
tau=1/L_f;
x0=D_gradient( u0 );
x1=x0;
q=((tau*mu)/(1+tau*mu_g));
t(1)=1;
%%
for k=1:iter1+1
    t(k+1)=(1-q*(t(k)^2)+sqrt(((1-q*(t(k)^2))^2)+4*(t(k)^2)))/2;
    beta(k)=((t(k)-1)/t(k+1))*((1+tau*mu_g-t(k+1)*tau*mu)/(1-tau*mu_f));
    y=x1+(beta(k))*(x1-x0);
    y_tilde=y-tau*grad_f(y,u0);
    %-----------------------------------------------------------------
    x2=prox_tau_g( y_tilde,lambda,tau,mu_g );
    %-----------------------------------
    u2=u0-D_adjoint( x2 );
    improvement(k)=(1 - var(u2(:)-u(:)) / var(u0(:)-u(:))) * 100;
    %-----------------------------------
    F_u(k)=F_dual( x2,u0,lambda,ep );
    x0=x1;
    x1=x2;
    %-----------------------------------
    if mod(k,100)==0
        disp(['iter = ', num2str(k)]);
    end
end
%%
u_star=u0-D_adjoint( x2 );
end

