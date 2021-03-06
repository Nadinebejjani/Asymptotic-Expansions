function [ res ] = fct_cisaillement(Plaque,PHI,n,v )
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
       
       % Calcul de C
        O = zeros(1,1,NCouches) ;
        nuTL = nuLT.*real(ET)./real(EL) ;
        nuNL = nuLN.*real(EN)./real(EL) ;
        nuNT = nuTN.*real(EN)./real(ET) ;  
        
         d = 1-nuTN.*nuNT-nuNL.*nuLN-nuLT.*nuTL-nuTL.*nuLN.*nuNT-nuLT.*nuTN.*nuNL ;
        C11 = EL.*(1-nuTN.*nuNT)./d ;
        C22 = ET.*(1-nuNL.*nuLN)./d ; 
        C33 = EN.*(1-nuLT.*nuTL)./d ; 
        C12 = real(ET).*(nuLT+nuNT.*nuLN)./d ; 
        C13 = real(EN).*(nuLN+nuTN.*nuLT)./d ; 
        C23 = real(EN).*(nuTN+nuLN.*nuTL)./d ; 
        C21 = C12 ;%EL.*(nuTL+nuNL.*nuTN)./d ; 
        C31 = C13 ;%EL.*(nuNL+nuTL.*nuNT)./d ; 
        C32 = C23 ;%ET.*(nuNT+nuLT.*nuNL)./d ;
            C11 = reshape(C11,[1 1 NCouches]) ;
            C22 = reshape(C22,[1 1 NCouches]) ;
            C33 = reshape(C33,[1 1 NCouches]) ;
            C12 = reshape(C12,[1 1 NCouches]) ;
            C13 = reshape(C13,[1 1 NCouches]) ;
            C23 = reshape(C23,[1 1 NCouches]) ;
            C21 = reshape(C21,[1 1 NCouches]) ;
            C31 = reshape(C31,[1 1 NCouches]) ;
            C32 = reshape(C32,[1 1 NCouches]) ;
            C44 = reshape(GTN,[1 1 NCouches]) ;
            C55 = reshape(GNL,[1 1 NCouches]) ;
            C66 = reshape(GLT,[1 1 NCouches]) ;
        C = [C11 C12 C13 O O O ; ... Raideur
                C21 C22 C23 O O O ; ...
                C31 C32 C33 O O O ; ...
                O O O C44 O O ; ...
                O O O O C55 O ; ...
                O O O O O C66 ];
    
            % Calcul de S
            for k=1:NCouches
                S(:,:,k,1)=inv(C(:,:,k,1)) ;
            end
                
            
            % Matrice de rotation des contraintes
   Os = @(c,s)[c^2, s^2, 0, 0, 0, 2*c*s ;...
                    s^2, c^2, 0, 0, 0, -2*c*s ;...
                    0, 0, 1, 0, 0, 0 ;...
                    0, 0, 0, c, -s, 0 ;...
                    0, 0, 0, s, c, 0 ;...
                    -c*s, c*s, 0, 0, 0, c^2-s^2 ] ;
         

               % Matrice de rotation des deformations 
                Oe = @(c,s)[c^2, s^2, 0, 0, 0, c*s ;...
                    s^2, c^2, 0, 0, 0, -c*s ;...
                    0, 0, 1, 0, 0, 0 ;...
                    0, 0, 0, c, -s, 0 ;...
                    0, 0, 0, s, c, 0 ;...
                    -2*c*s, 2*c*s, 0, 0, 0, c^2-s^2 ] ;
                
          
        
          % Rotation de C et S
                CM = zeros(6,6,NCouches,NPHI) ;
                SM=zeros(6,6,NCouches,NPHI) ;
                for p = 1:length(PHI)
        for k = 1:NCouches
            CM(:,:,k,p) = Os(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p)))*C(:,:,k)*Os(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p))).' ;
            SM(:,:,k,p) = Oe(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p)))*S(:,:,k)*Oe(cos(THETA(k)-PHI(p)),sin(THETA(k)-PHI(p))).' ;
        end
                end
    
        
   Csisgma = zeros(3,3,NCouches,NPHI) ;   %%%%%tenseur d elasticite plane
   for p = 1:length(PHI)
        for k = 1:NCouches
   Csigma(1,1,k,p)=CM(1,1,k,p)-(CM(1,3,k,p)*CM(1,3,k,p))/CM(3,3,k,p);
   Csigma(1,2,k,p)=CM(1,2,k,p)-(CM(1,3,k,p)*CM(2,3,k,p))/CM(3,3,k,p);
   Csigma(1,3,k,p)=CM(1,6,k,p)-(CM(1,3,k,p)*CM(3,6,k,p))/CM(3,3,k,p);
   Csigma(2,1,k,p)=Csigma(1,2,k,p);
   Csigma(2,2,k,p)=CM(2,2,k,p)-(CM(2,3,k,p)*CM(2,3,k,p))/CM(3,3,k,p);
   Csigma(2,3,k,p)=CM(2,6,k,p)-(CM(2,3,k,p)*CM(6,3,k,p))/CM(3,3,k,p);
   Csigma(3,1,k,p)=Csigma(1,3,k,p);
   Csigma(3,2,k,p)=Csigma(2,3,k,p);
   Csigma(3,3,k,p)=CM(6,6,k,p)-(CM(3,6,k,p)*CM(3,6,k,p))/CM(3,3,k,p);
        end
   end 
   
     NPHI = length(PHI) ; 
        h = sum(e) ;
          Codes = [] ;
        Codes(1) = -1/2 ;
 for k = 1:NCouches      
                Codes = [Codes ; Codes(end)+e(k)/h] ;
 end
    
 %Nodes 
        Nodes = zeros(1,n+1) ;
        Nodes(1,1) = -1/2 ;
for k = 1:n
                Nodes(1,k+1)=-0.5+k/n;
end

%tableau des rho
RHO=zeros(1,n+1);
i=1;

    while i<= NCouches
        for k=1:n+1
if (Nodes(k) >= Codes(i)) && (Nodes(k) <= Codes(i+1))
    RHO(1,k)=rho(i);
end
        end,

i=i+1;
    end
RHO(1,n+1)=rho(NCouches);

% tableaux des C et S
 % tableaux des C et S
  p=1;
c11=zeros(1,n+1);  
c16=zeros(1,n+1);   
c13=zeros(1,n+1);
c33=zeros(1,n+1);   
c23=zeros(1,n+1);   
c66=zeros(1,n+1);
s44=zeros(1,n+1);
s45=zeros(1,n+1);
s55=zeros(1,n+1);
cs11=zeros(1,n+1);
cs12=zeros(1,n+1);
cs22=zeros(1,n+1);
i=1;

    while i<= NCouches
        for k=1:n+1
if (Nodes(k) >= Codes(i)) && (Nodes(k) <= Codes(i+1))
    c11(k)=CM(1,1,i,p);
    c13(k)=CM(1,3,i,p);
    c33(k)=CM(3,3,i,p);
    c16(k)=CM(1,6,i,p);
    c23(k)=CM(3,6,i,p);
    c66(k)=CM(6,6,i,p);    
    s44(k)=SM(4,4,i,p);
    s45(k)=SM(4,5,i,p);
    s55(k)=SM(5,5,i,p);
    cs11(k)=Csigma(1,1,i,p);
    cs12(k)=Csigma(1,3,i,p);
    cs22(k)=Csigma(3,3,i,p);
end
        end
i=i+1;
    end
    c11(n+1)=CM(1,1,NCouches,p);
    c13(n+1)=CM(1,3,NCouches,p);
    c33(n+1)=CM(3,3,NCouches,p);
    c16(n+1)=CM(1,6,NCouches,p);
    c23(n+1)=CM(3,6,NCouches,p);
    c66(n+1)=CM(6,6,NCouches,p);    
    s44(n+1)=SM(4,4,NCouches,p);
    s45(n+1)=SM(4,5,NCouches,p);
    s55(n+1)=SM(5,5,NCouches,p);
    cs11(n+1)=Csigma(1,1,NCouches,p);
    cs12(n+1)=Csigma(1,3,NCouches,p);
    cs22(n+1)=Csigma(3,3,NCouches,p);
    
 %%%%%%%%%%%%%%%%%%%
    %Matrice A %
%%%%%%%%%%%%%%%%%%
A=zeros(2,2); 
A(1,1)=average(n,cs11);
A(1,2)=average(n,cs12);
A(2,1)=A(1,2);
A(2,2)=average(n,cs22);
%%%%%%%%%%%%%%%%%%%
    %D11%
%%%%%%%%%%%%%%%%%%
D11=average(n,cs11.*(Nodes.^2));
  

%%%%%%%%%%%%%%%%%%
%Moyenne de Rho%%
%%%%%%%%%%%%%%%%%%

intrho=average(n,RHO);

  %%%%%%%%%%%%%%%%Valeurs propres et vecteurs propres de A
X1=zeros(2,1);
X2=zeros(2,1);
[Vecttemp,Valtemp]=eig(A);
%%%%%%%%%A1%%%%%%%%%%%
%Vecteurs propres de A%
%%%%%%%%%A2%%%%%%%%%%%

X1(1,1)=Vecttemp(1,1); 
X1(2,1)=Vecttemp(2,1);


X2(1,1)=Vecttemp(1,2);
X2(2,1)=Vecttemp(2,2);
        
%%%%%%%%%A1%%%%%%%%%%%
%Valeurs propres de A%
%%%%%%%%%A2%%%%%%%%%%%
A1=Valtemp(1,1);
A2=Valtemp(2,2);
        
%%%%%%%%%%%%%%%%%%
%alpha0 et alpha0bar%%
%%%%%%%%%%%%%%%%%%
alpha0= A1/intrho;
alpha0bar= A2/intrho;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ordre p=0
sigma33ord0=0;
 
s=1;
v1ord0=s*X1(1,1);
v2ord0=s*X1(2,1);
V0=[v1ord0;v2ord0];

eps11ord0=(-1i)*v1ord0;
eps12ord0=(-1i)*v2ord0;

sigma11ord0=cs11.*eps11ord0+cs12.*eps12ord0;
sigma12ord0=cs12.*eps11ord0+cs22.*eps12ord0;

eps33ord0=-(c13./c33)*eps11ord0-(c23./c33)*eps12ord0;
v3barord1=funcpn(n,eps33ord0);
V3ord1=-average(n,RHO.*v3barord1)/intrho;
v3ord1=v3barord1+V3ord1;

sigma13ord1=funcqn(n,(1i)*sigma11ord0-RHO*alpha0*v1ord0);
sigma23ord1=funcqn(n,(1i)*sigma12ord0-RHO*alpha0*v2ord0);

epsilon13ord1=(s55).*sigma13ord1+(s45).*sigma23ord1;
v1barord2=funcpn(n,epsilon13ord1)+funcpn(n,(1i)*v3barord1);


epsilon23ord1=(s45).*sigma13ord1+(s44).*sigma23ord1;
v2barord2=funcpn(n,epsilon23ord1);

eps11barord2=(-1i)*v1barord2;
eps12barord2=(-1i)*v2barord2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ordre p=1 
sigma33ord2=funcqn(n,(1i)*sigma13ord1-alpha0*RHO.*v3ord1);
 
sigma11barord2=cs11.*eps11barord2+cs12.*eps12barord2+(c13./c33).*sigma33ord2;
sigma12barord2=cs12.*eps11barord2+cs22.*eps12barord2+(c23./c33).*sigma33ord2;

N11barord2=average(n,sigma11barord2);
N12barord2=average(n,sigma12barord2);

U0ord2=[-(1i)*N11barord2+alpha0*average(n,RHO.*v1barord2);-(1i)*N12barord2+alpha0*average(n,RHO.*v2barord2)];
 
alpha2=-prodscal(U0ord2,X1)/(intrho*prodscal(V0,X1));%-prodscal(U0ord2,X1)/(intrho*prodscal(V0,X1));

Vord2=(1/(A2-A1))*(U0ord2-prodscal(U0ord2,X1)*X1);
Vord2test=(1/(A2-A1))*(alpha2*intrho*X1+U0ord2);

 
v1ord2=Vord2(1)+v1barord2+(1i)*Nodes*V3ord1;
v2ord2=Vord2(2)+v2barord2;
 
eps11ord2=(-1i)*v1ord2;
eps12ord2=(-1i)*v2ord2;
 
sigma11ord2=cs11.*eps11ord2+cs12.*eps12ord2+(c13./c33).*sigma33ord2;
sigma12ord2=cs12.*eps11ord2+cs22.*eps12ord2+(c23./c33).*sigma33ord2;

eps33ord2=(1./c33).*sigma33ord2-(c13./c33).*eps11ord2-(c23./c33).*eps12ord2; 
v3barord3=funcpn(n,eps33ord2);

M11ord2=average(n,Nodes.*sigma11ord2);

V3ord3=(M11ord2+(1i)*alpha0*average(n,RHO.*Nodes.*v1ord2)-alpha2*average(n,RHO.*v3ord1)-alpha0*average(n,RHO.*v3barord3))/(alpha0*intrho);
v3ord3=v3barord3+V3ord3;

sigma13ord3=funcqn(n,(1i)*sigma11ord2-RHO.*(alpha2*v1ord0+alpha0*v1ord2));
sigma23ord3=funcqn(n,(1i)*sigma12ord2-RHO.*(alpha2*v2ord0+alpha0*v2ord2));

epsilon13ord3=(s55).*sigma13ord3+(s45).*sigma23ord3;
v1barord4=funcpn(n,epsilon13ord3)+funcpn(n,(1i)*v3barord3);

epsilon23ord3=(s45).*sigma13ord3+(s44).*sigma23ord3;
v2barord4=funcpn(n,epsilon23ord3);

eps11barord4=(-1i)*v1barord4;
eps12barord4=(-1i)*v2barord4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calcul de alpha4
sigma33ord4=funcqn(n,(1i)*sigma13ord3)-alpha2*funcqn(n,RHO.*v3ord1)-alpha0*funcqn(n,RHO.*v3ord3);

    
 sigma11bar4=cs11.*eps11barord4+cs12.*eps12barord4+(c13./c33).*sigma33ord4;
 sigma12bar4=cs12.*eps11barord4+cs22.*eps12barord4+(c23./c33).*sigma33ord4;
    
N11bar4=average(n,sigma11bar4);
N12bar4=average(n,sigma12bar4);   

p14=alpha2*average(n,RHO.*v1ord2)+alpha0*average(n,RHO.*v1barord4);
p24=alpha2*average(n,RHO.*v2ord2)+alpha0*average(n,RHO.*v2barord4);

Ualpha4=[-(1i)*N11bar4+p14;-(1i)*N12bar4+p24];
    
alpha4=-prodscal(Ualpha4,X1)/(intrho*prodscal(X1,X1));  

%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ordre qlq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tableaux initiaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
sigma33=zeros(1,n+1,v);
sigma33(:,:,1)=0;
sigma33(:,:,2)=sigma33ord2;

alpha=zeros(1,v);
alpha(1)=alpha0;
alpha(2)=alpha2;

v3=zeros(1,n+1,v);
v3(:,:,1)=v3ord1;
v3(:,:,2)=v3ord3;

sigma13=zeros(1,n+1,v);
sigma13(:,:,1)=sigma13ord1;
sigma13(:,:,2)=sigma13ord3;

sigma23=zeros(1,n+1,v);
sigma23(:,:,1)=sigma23ord1;
sigma23(:,:,2)=sigma23ord3;

M11=zeros(1,v);
N11bar=zeros(1,v);
N12bar=zeros(1,v);

eps11bar=zeros(1,n+1,v);
eps11bar(:,:,1)=0;
eps11bar(:,:,2)=eps11barord2;
eps11bar(:,:,3)=eps11barord4;

eps12bar=zeros(1,n+1,v);
eps12bar(:,:,1)=0;
eps12bar(:,:,2)=eps12barord2;
eps12bar(:,:,3)=eps12barord4;

sigma11bar=zeros(1,n+1,v);
sigma11bar(:,:,1)=0;
sigma11bar(:,:,2)=sigma11barord2;

sigma12bar=zeros(1,n+1,v);
sigma12bar(:,:,1)=0;
sigma12bar(:,:,2)=sigma12barord2;


v1=zeros(1,n+1,v);
v2=zeros(1,n+1,v);
v1(:,:,1)=v1ord0;
v1(:,:,2)=v1ord2; 
v2(:,:,1)=v2ord0;
v2(:,:,2)=v2ord2;

v1bar=zeros(1,n+1,v);
v2bar=zeros(1,n+1,v);
v1bar(:,:,1)=0;
v1bar(:,:,2)=v1barord2;
v1bar(:,:,3)=v1barord4;


v2bar(:,:,1)=0;
v2bar(:,:,2)=v2barord2;
v2bar(:,:,3)=v2barord4;

p1=zeros(1,v);
p2=zeros(1,v);

eps11=zeros(1,n+1,v);
eps12=zeros(1,n+1,v);
% 
eps11(:,:,1)=eps11ord0;
eps11(:,:,2)=eps11ord2;

eps12(:,:,1)=eps12ord0;
eps12(:,:,2)=eps12ord2;

eps33=zeros(1,n+1,v);
eps33(:,:,1)=eps33ord0;
eps33(:,:,2)=eps33ord2;

sigma11=zeros(1,n+1,v);
sigma11(:,:,1)=sigma11ord0;
sigma11(:,:,2)=sigma11ord2;

sigma12=zeros(1,n+1,v);
sigma12(:,:,1)=sigma12ord0;
sigma12(:,:,2)=sigma12ord2;

V3=zeros(1,v);
V3(1)=0;
V3(2)=V3ord1;
V3(3)=V3ord3;

v3bar=zeros(1,n+1,v);
v3bar(:,:,1)=0;
v3bar(:,:,2)=v3barord1;
v3bar(:,:,3)=v3barord3;

epsilon23=zeros(1,n+1,v);
epsilon23(:,:,1)=epsilon23ord1;
epsilon23(:,:,2)=epsilon23ord3;

epsilon13=zeros(1,n+1,v);
epsilon13(:,:,1)=epsilon13ord1;
epsilon13(:,:,2)=epsilon13ord3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcul
for k=3:v
    
 sigma33(:,:,k)=funcqn(n,(1i)*sigma13(:,:,k-1)-alpha(k-1)*RHO.*v3(:,:,1)-alpha(1)*RHO.*v3(:,:,k-1));
    for u=2:k-2
sigma33(:,:,k)=sigma33(:,:,k)-(alpha(u)*funcqn(n,RHO.*v3(:,:,k-u)));
    end
    
sigma11bar(:,:,k)=cs11.*eps11bar(:,:,k)+cs12.*eps12bar(:,:,k)+(c13./c33).*sigma33(:,:,k);
sigma12bar(:,:,k)=cs12.*eps11bar(:,:,k)+cs22.*eps12bar(:,:,k)+(c23./c33).*sigma33(:,:,k);  
    
    
N11bar(k)=average(n,sigma11bar(:,:,k));
N12bar(k)=average(n,sigma12bar(:,:,k));   

p1(k)=alpha(k-1)*average(n,RHO.*v1(:,:,2))+alpha(1)*average(n,RHO.*v1bar(:,:,k));
for u=2:k-2
p1(k)=p1(k)+alpha(u)*average(n,RHO.*v1(:,:,k-u+1));
end

p2(k)=alpha(k-1)*average(n,RHO.*v2(:,:,2))+alpha(1)*average(n,RHO.*v2bar(:,:,k));
for u=2:k-2
p2(k)=p2(k)+alpha(u)*average(n,RHO.*v2(:,:,k-u+1));
end

Ualpha2p=[-(1i)*N11bar(k)+p1(k);-(1i)*N12bar(k)+p2(k)];
    
alpha(k)=-prodscal(Ualpha2p,X1)/(intrho*prodscal(X1,X1));  

Vord2k=(1/(A2-A1))*(Ualpha2p+alpha(k)*intrho*X1);

v1(:,:,k)=Vord2k(1)+v1bar(:,:,k)+(1i)*Nodes*V3(k);
v2(:,:,k)=Vord2k(2)+v2bar(:,:,k);
 

eps11(:,:,k)=(-1i)*v1(:,:,k);
eps12(:,:,k)=(-1i)*v2(:,:,k);

sigma11(:,:,k)=cs11.*eps11(:,:,k)+cs12.*eps12(:,:,k)+(c13./c33).*sigma33(:,:,k);
sigma12(:,:,k)=cs12.*eps11(:,:,k)+cs22.*eps12(:,:,k)+(c23./c33).*sigma33(:,:,k);

eps33(:,:,k)=(1./c33).*sigma33(:,:,k)-(c13./c33).*eps11(:,:,k)-(c23./c33).*eps12(:,:,k);
v3bar(:,:,k+1)=funcpn(n,eps33(:,:,k));

M11(k)=average(n,Nodes.*sigma11(:,:,k));

r1(k)=(1i)*alpha(k-1)*average(n,RHO.*Nodes.*v1(:,:,2))+(1i)*alpha(1)*average(n,RHO.*Nodes.*v1(:,:,k))-alpha(k)*average(n,RHO.*v3(:,:,1))-alpha(k-1)*average(n,RHO.*v3(:,:,2))-alpha(1)*average(n,RHO.*v3bar(:,:,k+1));
for u=2:k-2
r1(k)=r1(k)+(1i)*alpha(u)*average(n,RHO.*Nodes.*v1(:,:,k-u+1))-alpha(u)*average(n,RHO.*v3(:,:,k-u+1));
end


V3(k+1)=(M11(k)+r1(k))/(alpha0*intrho);

v3(:,:,k)=v3bar(:,:,k+1)+V3(k+1);

sigma13(:,:,k)=funcqn(n,(1i)*sigma11(:,:,k)-RHO.*(alpha(k)*v1(:,:,1)+alpha(k-1)*v1(:,:,2)+alpha(1)*v1(:,:,k)));
    for u=2:k-2
sigma13(:,:,k)=sigma13(:,:,k)-funcqn(n,RHO.*(alpha(u)*v1(:,:,k-u+1)));
    end
    
    sigma23(:,:,k)=funcqn(n,(1i)*sigma12(:,:,k)-RHO.*(alpha(k)*v2(:,:,1)+alpha(k-1)*v2(:,:,2)+alpha(1)*v2(:,:,k)));
    for u=2:k-2
sigma23(:,:,k)=sigma23(:,:,k)-funcqn(n,RHO.*(alpha(u)*v2(:,:,k-u+1)));
    end
    
    
epsilon13(:,:,k)=(s55).*sigma13(:,:,k)+(s45).*sigma23(:,:,k);
epsilon23(:,:,k)=(s45).*sigma13(:,:,k)+(s44).*sigma23(:,:,k);   
    
    
v1bar(:,:,k+1)=funcpn(n,epsilon13(:,:,k)+(1i)*v3(:,:,k));
v2bar(:,:,k+1)=funcpn(n,epsilon23(:,:,k)); 
 
eps11bar(:,:,k+1)=(-1i)*v1bar(:,:,k+1);
eps12bar(:,:,k+1)=(-1i)*v2bar(:,:,k+1);    
       
end

res.Nodes=Nodes;
res.Codes=Codes;
res.RHO=RHO;
res.c13=c13;
res.CM=CM;
res.A=A;
res.D11=D11;
res.intrho=intrho;
res.X1=X1;
res.X2=X2;
res.A1=A1;
res.A2=A2;
res.alpha0bar=alpha0bar;
res.alpha0=alpha0;
res.alpha2=alpha2;
res.alpha4=alpha4;
res.alpha=alpha;
end

