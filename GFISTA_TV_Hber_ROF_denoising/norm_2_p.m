function [ N ] = norm_2_p( p )
[m,n,r]=size(p);
N=0;
for i=1:m
    for j=1:n
        for k=1:r
            N=N+(p(i,j,k)^2);
        end
    end
end
end

