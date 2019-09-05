function [y] = NEO(x,L,offset)  
% Function to return the NEO of an input signal x 

x=x(:)'+offset;
xp = [x(L+1:end) x(end-L+1:end)]; 
xn = [x(1:L) x(1:end-L)];   
y = x.^2 - xp.*xn;   
end 
