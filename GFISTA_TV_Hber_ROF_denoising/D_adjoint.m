function [ D_star ] = D_adjoint( p )
[m,n,k]=size(p);
div=zeros(m,n,k);
div_p=zeros(m,n);
for i=1:m
    for j=1:n
        if (i>1 && i<m)
            div(i,j,1)=p(i,j,1)-p(i-1,j,1);
        elseif (i==1)
            div(i,j,1)=p(i,j,1);
        elseif (i==m)
            div(i,j,1)=-p(i-1,j,1);
        end
        %-----------------------------------
        if (j>1 && j<n)
            div(i,j,2)=p(i,j,2)-p(i,j-1,2);
        elseif (j==1)
            div(i,j,2)=p(i,j,2);
        elseif (j==n)
            div(i,j,2)=-p(i,j-1,2);
        end   
        div_p(i,j)=div(i,j,1)+div(i,j,2);
    end
end
D_star=-div_p;
end

