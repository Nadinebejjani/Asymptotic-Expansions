function [ moy ] = average(n,g )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

% NCouches = Plaque.NCouches ;
%         NElmts = Plaque.NElmts ;
%         e = Plaque.e ; % Epaisseurs
% h = sum(e) ;
% 
%         Nodes = [] ;
%         Nodes(1) = -1/2 ;
% for k = 1:n
%                 Nodes(k+1)=-0.5+k/n;
% end

  
% moy=0;
% for k = 1:n
%      moy=moy+0.5*(1/(n))*(g(k)+g(k+1))  ;
% end   

moy=0;
for k = 1:n
     moy=moy+(g(k)+g(k+1))*0.5*(1/(n));
end  


end

