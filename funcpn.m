function [ Pk ] = funcpn( n,f )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gnk(1)=0;

for k=1:n
  gnk(k+1)=gnk(k)+(1/n)*(1/2)*(f(k+1)+f(k));  
end

Ig=0;

for k=1:n
Ig=Ig+(1/n)*(1/2)*(gnk(k+1)+gnk(k)); 
end

for k=1:n+1
 Pk(k)= gnk(k)-Ig;
end

end