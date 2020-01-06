function [ res ] = calculus(Plaque,PHI,n,v )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

NCouches = Plaque.NCouches ;
NElmts = Plaque.NElmts ;
THETA = Plaque.THETA ;
e = Plaque.e ; % Epaisseurs
EL = Plaque.EL ;
ET = Plaque.ET ;
EN = Plaque.EN ;
nuLN = Plaque.nuLN ;
nuTN = Plaque.nuTN ;
nuLT = Plaque.nuLT ;
GLT = Plaque.GLT ;
GNL = Plaque.GNL ;
GTN = Plaque.GTN ;
rho = Plaque.rho ;
        
        
% Vecteur angle de propagation de l'onde
NPHI = length(PHI) ; 
h = sum(e) ;
p=1;
       
Nodes = [] ; % points de -1/2 à 1/2
        Nodes(1) = -h/2 ;
        for k = 1:NCouches
                Nodes = [Nodes ; Nodes(end)+e(k)] ;  %%%position des points
        end
        
L = Nodes(2:end)-Nodes(1:end-1) ; 

%moyenne de rho (-1/2,1/2)
rhobar=0;
for k=1:NCouches
rhobar=rhobar+ rho(k)*L(k);
end


%moyenne de G_LN
Gbar=0;

for k=1:NCouches
   Gbar=Gbar+Plaque.GNL(k)*L(k); 
end


res.Nodes=Nodes;
res.rhobar=rhobar;
res.Gbar=Gbar;
