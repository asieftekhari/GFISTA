function [ D_f ] = D_f( x,y,u0 )
g=grad_f(y,u0);
d=x-y;
a=(0.5)*( (N_2( D_adjoint(x)-u0 ))^2 );
b=(0.5)*( (N_2( D_adjoint(y)-u0 ))^2 );
D_f = a - b - prodd(g,d);
end