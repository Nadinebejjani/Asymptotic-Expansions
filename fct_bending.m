function [ res ] = fct_bending(Plaque,PHI,n,v )
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
intrhoz2=average(n,RHO.*Nodes.^2);

%%%%%%%%%%%%%%%%Valeurs propres et vecteurs propres de A
a2=real(D11)/real(intrho);

%Caldul de a4%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%p=0%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
v3moins1=1;
V1ord0=0;
V2ord0=0;
v1ord0=V1ord0+(1i)*Nodes*v3moins1;
v2ord0=V2ord0;
eps33ord0=-(c13./c33).*Nodes;
v3barord1=funcpn(n,eps33ord0);

eps11ord0=(-1i)*v1ord0;
eps12ord0=(-1i)*v2ord0;

sigma11ord0=cs11.*eps11ord0;
sigma12ord0=cs12.*eps11ord0;

%sigma13ord1=funcqn(n,(1i)*cs11.*Nodes);
sigma13ord1=funcqn(n,(1i)*sigma11ord0);
sigma23ord1=funcqn(n,(1i)*sigma12ord0);
%sigma23ord1=funcqn(n,(1i)*cs12.*Nodes);

epsilon13ord1=(s55).*sigma13ord1+(s45).*sigma23ord1;
v1barord2=funcpn(n,epsilon13ord1+(1i)*v3barord1);
%v1barord2=funcpn(n,epsilon13ord1+(1i)*v3ord1);

epsilon23ord1=(s45).*sigma13ord1+(s44).*sigma23ord1;
v2barord2=funcpn(n,epsilon23ord1);

eps11barord2=(-1i)*v1barord2;
eps12barord2=(-1i)*v2barord2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%p=1%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
sigma33ord2=funcqn(n,(1i)*sigma13ord1-a2*RHO);

sigma11barord2=(cs11.*((-1i)*v1barord2)+cs12.*((-1i)*v2barord2)+(c13./c33).*sigma33ord2);
sigma12barord2=(cs12.*((-1i)*v1barord2)+cs22.*((-1i)*v2barord2)+(c23./c33).*sigma33ord2);

M11barord2=average(n,Nodes.*sigma11barord2);

a4=(M11barord2-a2*average(n,RHO.*Nodes.^2)-a2*average(n,RHO.*v3barord1))/average(n,RHO);
num=M11barord2-a2*average(n,RHO.*Nodes.^2)-a2*average(n,RHO.*v3barord1);
w1=-a2*average(n,RHO.*Nodes.^2);
w2=-a2*average(n,RHO.*v3barord1);


V3ord1=0;% 1
v3ord1=v3barord1+V3ord1;

N11barord2=average(n,sigma11barord2);
N12barord2=average(n,sigma12barord2);

b=[(-1i)*N11barord2;(-1i)*N12barord2];
valphaord2=inv(A)*b;
V1ord2=valphaord2(1);
V2ord2=valphaord2(2);

v1ord2=v1barord2+V1ord2+(1i)*V3ord1*Nodes;
v2ord2=v2barord2+V2ord2;

eps11ord2=(-1i)*v1ord2;
eps12ord2=(-1i)*v2ord2;


sigma11ord2=cs11.*eps11ord2+cs12.*eps12ord2+(c13./c33).*sigma33ord2;
sigma12ord2=cs12.*eps11ord2+cs22.*eps12ord2+(c23./c33).*sigma33ord2;

eps33ord2=(sigma33ord2./c33)-(c13./c33).*eps11ord2-(c23./c33).*eps12ord2;
v3barord3=funcpn(n,eps33ord2);

sigma13ord3=funcqn(n,(1i)*sigma11ord2-RHO.*(a2*v1ord0));
sigma23ord3=funcqn(n,(1i)*sigma12ord2-RHO.*(a2*v2ord0));

epsilon13ord3=(s55).*sigma13ord3+(s45).*sigma23ord3;
v1barord4=funcpn(n,epsilon13ord3+(1i)*v3barord3);

epsilon23ord3=(s45).*sigma13ord3+(s44).*sigma23ord3;
v2barord4=funcpn(n,epsilon23ord3);

eps11barord4=(-1i)*v1barord4;
eps12barord4=(-1i)*v2barord4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%p=2% CALCUL a6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
sigma33ord4=funcqn(n,(1i)*sigma13ord3-a4*RHO*v3moins1-RHO.*(a2*v3ord1));
   
sigma11barord4=cs11.*eps11barord4+cs12.*eps12barord4+(c13./c33).*sigma33ord4;
sigma12barord4=cs12.*eps11barord4+cs22.*eps12barord4+(c23./c33).*sigma33ord4;
    
M11barord4=average(n,Nodes.*sigma11barord4);

r1=(1i)*a4*average(n,RHO.*Nodes.*v1ord0)+(1i)*a2*average(n,RHO.*Nodes.*v1ord2)-a4*average(n,RHO.*v3ord1)-a2*average(n,RHO.*v3barord3);

a6=(M11barord4+r1)/(intrho*v3moins1);
V3ord3=0;% 1
v3ord3=v3barord3+V3ord3;

N11barord4=average(n,sigma11barord4);
N12barord4=average(n,sigma12barord4);

b=[(-1i)*N11barord4+a4*average(n,RHO.*v1ord0)+a2*average(n,RHO.*v1ord2);(-1i)*N12barord4+a4*average(n,RHO.*v2ord0)+a2*average(n,RHO.*v2ord2)];
valphaord4=inv(A)*b;
V1ord4=valphaord4(1);
V2ord4=valphaord4(2);

v1ord4=v1barord4+V1ord4+(1i)*V3ord3*Nodes;
v2ord4=v2barord4+V2ord4;

eps11ord4=(-1i)*v1ord4;
eps12ord4=(-1i)*v2ord4;

sigma11ord4=cs11.*eps11ord4+cs12.*eps12ord4+(c13./c33).*sigma33ord4;
sigma12ord4=cs12.*eps11ord4+cs22.*eps12ord4+(c23./c33).*sigma33ord4;

eps33ord4=(sigma33ord4./c33)-(c13./c33).*eps11ord4-(c23./c33).*eps12ord4;
v3barord5=funcpn(n,eps33ord4);

sigma13ord5=funcqn(n,(1i)*sigma11ord4-RHO.*(a4*v1ord0+a2*v1ord2));
sigma23ord5=funcqn(n,(1i)*sigma12ord4-RHO.*(a4*v2ord0+a2*v2ord2));

epsilon13ord5=(s55).*sigma13ord5+(s45).*sigma23ord5;
v1barord6=funcpn(n,epsilon13ord5+(1i)*v3barord5);

epsilon23ord5=(s45).*sigma13ord5+(s44).*sigma23ord5;
v2barord6=funcpn(n,epsilon23ord5);

eps11barord6=(-1i)*v1barord6;
eps12barord6=(-1i)*v2barord6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%p=3% CALCUL a8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
sigma33ord6=funcqn(n,(1i)*sigma13ord5-a6*RHO*v3moins1-RHO.*(a4*v3ord1)-RHO.*(a2*v3ord3));
   
sigma11barord6=cs11.*eps11barord6+cs12.*eps12barord6+(c13./c33).*sigma33ord6;
sigma12barord6=cs12.*eps11barord6+cs22.*eps12barord6+(c23./c33).*sigma33ord6;
    
M11barord6=average(n,Nodes.*sigma11barord6);

r1=(1i)*a6*average(n,RHO.*Nodes.*v1ord0)+(1i)*a4*average(n,RHO.*Nodes.*v1ord2)+(1i)*a2*average(n,RHO.*Nodes.*v1ord4)-a6*average(n,RHO.*v3ord1)-a4*average(n,RHO.*v3ord3)-a6*average(n,RHO.*v3barord5);

a8=(M11barord6+r1)/(intrho*v3moins1);
V3ord5=0;% 1
v3ord5=v3barord5+V3ord5;

N11barord6=average(n,sigma11barord6);
N12barord6=average(n,sigma12barord6);

b=[(-1i)*N11barord6+a6*average(n,RHO.*v1ord0)+a4*average(n,RHO.*v1ord2)+a2*average(n,RHO.*v1ord4);(-1i)*N12barord4+a6*average(n,RHO.*v2ord0)+a4*average(n,RHO.*v2ord2)+a2*average(n,RHO.*v2ord4)];
valphaord6=inv(A)*b;
V1ord6=valphaord6(1);
V2ord6=valphaord6(2);

v1ord6=v1barord6+V1ord6+(1i)*V3ord5*Nodes;
v2ord6=v2barord6+V2ord6;

eps11ord6=(-1i)*v1ord6;
eps12ord6=(-1i)*v2ord6;

sigma11ord6=cs11.*eps11ord6+cs12.*eps12ord6+(c13./c33).*sigma33ord6;
sigma12ord6=cs12.*eps11ord6+cs22.*eps12ord6+(c23./c33).*sigma33ord6;

eps33ord6=(sigma33ord6./c33)-(c13./c33).*eps11ord6-(c23./c33).*eps12ord6;
v3barord7=funcpn(n,eps33ord6);

sigma13ord7=funcqn(n,(1i)*sigma11ord6-RHO.*(a6*v1ord0+a4*v1ord2+a2*v1ord4));
sigma23ord7=funcqn(n,(1i)*sigma12ord6-RHO.*(a6*v2ord0+a4*v2ord2+a2*v2ord4));

epsilon13ord7=(s55).*sigma13ord7+(s45).*sigma23ord7;
v1barord8=funcpn(n,epsilon13ord7+(1i)*v3barord7);

epsilon23ord7=(s45).*sigma13ord7+(s44).*sigma23ord7;
v2barord8=funcpn(n,epsilon23ord7);

eps11barord8=(-1i)*v1barord8;
eps12barord8=(-1i)*v2barord8;


%%%%%%%%%%
%Algo%
%%%%%%%%%%
sigma33=zeros(1,n+1,v);
sigma33(:,:,1)=0;
sigma33(:,:,2)=sigma33ord2;

sigma13=zeros(1,n+1,v);
sigma13(:,:,1)=sigma13ord1;
sigma13(:,:,2)=sigma13ord3;

v3=zeros(1,n+1,v);
v3(:,:,1)=v3moins1;
v3(:,:,2)=v3ord1;
v3(:,:,3)=v3ord3;

alpha=zeros(1,v);
alpha(1)=a2;
alpha(2)=a4;

sigma11bar=zeros(1,n+1,v);
sigma11bar(:,:,1)=0;
sigma11bar(:,:,2)=sigma11barord2;

sigma11=zeros(1,n+1,v);
sigma11(:,:,1)=sigma11ord0;
sigma11(:,:,2)=sigma11ord2;

sigma12bar=zeros(1,n+1,v);
sigma12bar(:,:,1)=0;
sigma12bar(:,:,2)=sigma12barord2;

sigma12=zeros(1,n+1,v);
sigma12(:,:,1)=sigma12ord0;
sigma12(:,:,2)=sigma12ord2;

eps11=zeros(1,n+1,v);
eps12=zeros(1,n+1,v);
% 
eps11(:,:,1)=eps11ord0;
eps11(:,:,2)=eps11ord2;

eps12(:,:,1)=eps12ord0;
eps12(:,:,2)=eps12ord2;


eps11bar=zeros(1,n+1,v);
% eps11bar(:,:,1)=0;
% eps11bar(:,:,2)=eps11barord2;
% eps11bar(:,:,3)=eps11barord4;
eps11bar(:,:,1)=0;
eps11bar(:,:,2)=eps11barord2;
eps11bar(:,:,3)=eps11barord4;

eps12bar=zeros(1,n+1,v);
% eps12bar(:,:,1)=0;
% eps12bar(:,:,2)=eps12barord2;
% eps12bar(:,:,3)=eps12barord4;
eps12bar(:,:,1)=0;
eps12bar(:,:,2)=eps12barord2;
eps12bar(:,:,3)=eps12barord4;


M11bar=zeros(1,v);
N11bar=zeros(1,v);
N12bar=zeros(1,v);
p1=zeros(1,v);
p2=zeros(1,v);
r1=zeros(1,v);

v1=zeros(1,n+1,v);
v2=zeros(1,n+1,v);
v1(:,:,1)=v1ord0;
v1(:,:,2)=v1ord2; 
v2(:,:,1)=v2ord0;
v2(:,:,2)=v2ord2;

V1=zeros(1,v);
V2=zeros(1,v);

V1(1)=0;
V1(2)=V1ord2;
V2(1)=0;
V2(2)=V2ord2;


v1bar=zeros(1,n+1,v);
v2bar=zeros(1,n+1,v);
v1bar(:,:,1)=0;
v1bar(:,:,2)=v1barord2;
v1bar(:,:,3)=v1barord4;
% 
v2bar(:,:,1)=0;
v2bar(:,:,2)=v2barord2;
v2bar(:,:,3)=v2barord4;
%v1bar(:,:,1)=v1barord2;
%v1bar(:,:,2)=v1barord4;

% 
% v2bar(:,:,1)=v2barord2;
% v2bar(:,:,2)=v2barord4;
v3bar=zeros(1,n+1,v);
% v3bar(:,:,1)=0;
% v3bar(:,:,2)=v3barord1;
% v3bar(:,:,3)=v3ord3;
v3bar(:,:,1)=0;
v3bar(:,:,2)=v3barord1;
v3bar(:,:,3)=v3barord3;


eps33=zeros(1,n+1,v);
eps33(:,:,1)=eps33ord0;
eps33(:,:,2)=eps33ord2;

sigma13=zeros(1,n+1,v);
% sigma13(:,:,1)=sigma13ord1;
% sigma13(:,:,2)=sigma13ord3;
sigma13(:,:,2)=sigma13ord1;
sigma13(:,:,3)=sigma13ord3;


sigma23=zeros(1,n+1,v);
% sigma23(:,:,1)=sigma23ord1;
% sigma23(:,:,2)=sigma23ord3;
sigma23(:,:,2)=sigma23ord1;
sigma23(:,:,3)=sigma23ord3;

epsilon23=zeros(1,n+1,v);
% epsilon23(:,:,1)=epsilon23ord1;
% epsilon23(:,:,2)=epsilon23ord3;
epsilon23(:,:,2)=epsilon23ord1;
epsilon23(:,:,3)=epsilon23ord3;

epsilon13=zeros(1,n+1,v);
% epsilon13(:,:,1)=epsilon13ord1;
% epsilon13(:,:,2)=epsilon13ord3;
epsilon13(:,:,2)=epsilon13ord1;
epsilon13(:,:,3)=epsilon13ord3;
%V3=ones(1,n+1);
V3=zeros(1,n+1);

for k=3:v
 % sigma33(:,:,k)=funcqn(n,(1i)*sigma13(:,:,k-1)-alpha(k-1)*RHO.*(v3(:,:,1)));
 sigma33(:,:,k)=funcqn(n,(1i)*sigma13(:,:,k)-alpha(k-1)*RHO.*(v3(:,:,1)));
    for u=1:k-2
sigma33(:,:,k)=sigma33(:,:,k)-funcqn(n,alpha(u)*RHO.*v3(:,:,k-u));
    end  
    
sigma11bar(:,:,k)=cs11.*eps11bar(:,:,k)+cs12.*eps12bar(:,:,k)+(c13./c33).*sigma33(:,:,k);
sigma12bar(:,:,k)=cs12.*eps11bar(:,:,k)+cs22.*eps12bar(:,:,k)+(c23./c33).*sigma33(:,:,k); 
  
M11bar(k)=average(n,Nodes.*sigma11bar(:,:,k));

% r1(k)=(1i)*alpha(k-1)*average(n,RHO.*Nodes.*v1(:,:,1))-alpha(k-1)*average(n,RHO.*v3bar(:,:,2));
% for u=1:k-2
% r1(k)=r1(k)+(1i)*alpha(u)*average(n,RHO.*Nodes.*v1(:,:,k-u))-alpha(u)*average(n,RHO.*v3bar(:,:,k-u+1));
% end
%    
r1(k)=(1i)*alpha(k-1)*average(n,RHO.*Nodes.*v1(:,:,1))+(1i)*alpha(1)*average(n,RHO.*Nodes.*v1(:,:,k-1))-alpha(k-1)*average(n,RHO.*v3(:,:,2))-alpha(1)*average(n,RHO.*v3bar(:,:,k));
for u=2:k-2
r1(k)=r1(k)+(1i)*alpha(u)*average(n,RHO.*Nodes.*v1(:,:,k-u))-alpha(u)*average(n,RHO.*v3(:,:,k-u+1));
end


alpha(k)=(M11bar(k)+r1(k))/(intrho*v3moins1); 

N11bar(k)=average(n,sigma11bar(:,:,k));
N12bar(k)=average(n,sigma12bar(:,:,k));

p1=alpha(k-1)*average(n,RHO.*v1(:,:,1));
for u=1:k-2%2:k-2
p1=p1+alpha(u)*average(n,RHO.*v1(:,:,k-u));
end

p2=alpha(k-1)*average(n,RHO.*v2(:,:,1));
for u=1:k-2
p2=p2+alpha(u)*average(n,RHO.*v2(:,:,k-u));
end

b=[(-1i)*N11bar(k)+p1;(-1i)*N12bar(k)+p2];
valpha=inv(A)*b;
V1(k)=valpha(1);
V2(k)=valpha(2);

v1(:,:,k)= V1(k)+v1bar(:,:,k)+(1i)*Nodes*V3(k);   
v2(:,:,k)= V2(k)+v2bar(:,:,k);

eps11(:,:,k)=(-1i)*v1(:,:,k);
eps12(:,:,k)=(-1i)*v2(:,:,k);

sigma11(:,:,k)=cs11.*eps11(:,:,k)+cs12.*eps12(:,:,k)+(c13./c33).*sigma33(:,:,k);
sigma12(:,:,k)=cs12.*eps11(:,:,k)+cs22.*eps12(:,:,k)+(c23./c33).*sigma33(:,:,k);

eps33(:,:,k)=(1./c33).*sigma33(:,:,k)-(c13./c33).*eps11(:,:,k)-(c23./c33).*eps12(:,:,k);
v3bar(:,:,k+1)=funcpn(n,eps33(:,:,k));
v3(:,:,k+1)=v3bar(:,:,k+1)+V3(k+1);

sigma13(:,:,k+1)=funcqn(n,(1i)*sigma11(:,:,k)-RHO.*(alpha(k-1)*v1(:,:,1)));
    for u=1:k-2
sigma13(:,:,k+1)=sigma13(:,:,k+1)-funcqn(n,RHO.*(alpha(u)*v1(:,:,k-u)));
    end
    
    sigma23(:,:,k+1)=funcqn(n,(1i)*sigma12(:,:,k)-RHO.*(alpha(k-1)*v2(:,:,1)));
    for u=1:k-2
sigma23(:,:,k+1)=sigma23(:,:,k+1)-funcqn(n,RHO.*(alpha(u)*v2(:,:,k-u)));
    end

epsilon13(:,:,k+1)=(s55).*sigma13(:,:,k+1)+(s45).*sigma23(:,:,k+1);
epsilon23(:,:,k+1)=(s45).*sigma13(:,:,k+1)+(s44).*sigma23(:,:,k+1);

v1bar(:,:,k+1)=funcpn(n,epsilon13(:,:,k+1))+funcpn(n,(1i)*v3bar(:,:,k+1));
v2bar(:,:,k+1)=funcpn(n,epsilon23(:,:,k+1)); 
 
eps11bar(:,:,k+1)=(-1i)*v1bar(:,:,k+1);
eps12bar(:,:,k+1)=(-1i)*v2bar(:,:,k+1);
end



res.Nodes=Nodes;
res.a2=a2;
res.a4=a4;
res.a6=a6;
res.alpha=alpha;
res.RHO=RHO;
res.A=A;
res.num=num;
res.w1=w1;
res.w2=w2;
res.M11barord2=M11barord2;
res.v3barord1=v3barord1;
res.c13=c13;
res.c33=c33;
res.sigma33ord2=sigma33ord2;
res.cs11=cs11;
res.v1barord2=v1barord2;
res.s55=s55;
res.d=d;
res.CM=CM;
res.a6=a6;
res.a8=a8;
res.intrho=intrho;
res.intrhoz2=intrhoz2;
res.D11=D11;
end