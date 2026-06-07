function [ g ] = grad_f( p,u0 )
a=D_adjoint( p );
g=D_gradient( a - u0 );
end

