function [y] = ASO(x,L,offset) 
% Function to return the proposed of an input signal x   
x=x(:)'+offset;
xp = [x(1:L) x(1:end-L)];
for i = 1:L
    if abs(x(i))>=64
        x(i)=sign(x(i))*pow2(nextpow2(x(i)));
    else
        x(i)=sign(x(i))*pow2(nextpow2(x(i))-1);
    end
end
y = ((x - xp).*x);
end 
