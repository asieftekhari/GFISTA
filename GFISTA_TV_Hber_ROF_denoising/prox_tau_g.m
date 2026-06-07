function [ p ] = prox_tau_g( x_tilde,lambda,tau,mu_g )
[m,n,k]=size(x_tilde);
d=zeros(m,n);
p=zeros(m,n,k);
for i=1:m
    for j=1:n
        d(i,j)=sqrt( (x_tilde(i,j,1)^2) + (x_tilde(i,j,2)^2) );
        a=(1+(tau*mu_g))^(-1);
        b=( lambda*(1+(tau*mu_g)) )^(-1);
        c=max(1,b*d(i,j));
        p(i,j,1)= (a*x_tilde(i,j,1))/c;
        p(i,j,2)= (a*x_tilde(i,j,2))/c;
    end
end
end

