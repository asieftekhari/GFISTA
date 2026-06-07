function [ F_ad ] = F_dual( p,u0,lambda,ep )
[m,n,r]=size(p);
A=norm_2_p( p );
B=zeros(m,n);
%--------------------------------------
for i=1:m
    for j=1:n
        B(i,j)=sqrt((p(i,j,1)^2)+(p(i,j,2)^2));
    end
end
%--------------------------------------
if (B<=(lambda+1e-16))
    delta=0;
else
    delta=+Inf;
end
%--------------------------------------
f_p=(0.5)*( (N_2( D_adjoint(p)-u0 ))^2 );
g_p=((ep/(2*lambda))*A)+delta;
F_ad=f_p+g_p;
end

