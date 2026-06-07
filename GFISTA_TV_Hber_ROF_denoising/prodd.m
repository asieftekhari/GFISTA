function [ p ] = prodd( x,y )
A=x(:,:,1).*y(:,:,1)+x(:,:,2).*y(:,:,2);
p=sum(sum(A));
end

