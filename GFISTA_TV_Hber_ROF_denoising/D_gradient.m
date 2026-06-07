function [ D ] = D_gradient( u )
[m,n]=size(u);
D=zeros(m,n,2);
%%
for i=1:m
    for j=1:n
       if (i<m)
          D(i,j,1)=u(i+1,j)-u(i,j);
       elseif (i==m)
           D(i,j,1)=0;
       end
       %----------------------
       if (j<n)
          D(i,j,2)=u(i,j+1)-u(i,j);
       elseif (j==n)
           D(i,j,2)=0;
       end
    end
end
end

