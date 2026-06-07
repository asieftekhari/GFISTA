function [ N ] = N_2( u )
[m,n]=size(u);
M=0;
for i=1:m
    for j=1:n
        M=M+u(i,j)^2;
    end
end
N=sqrt(M);
end

