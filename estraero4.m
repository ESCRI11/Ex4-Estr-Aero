clear all

% Geometria

H1 = 0.45;
H2 = 0.9;
L  = 1.2;
Aab= 12e-6;

% Parámetros

P  = -9500;
Ibc= 1.5e18;
Iac= 3e18;
Icd= 3e18;
Iab= Aab/4;
E  = 60e9;
I  = [Iab Ibc Icd Iac];

% Grados de libertad

nd = 3;
nel = 4;
n = nel; 
nnod = 2;
nig = 3;
ngl = n * nig;

% Discretización

x = zeros(nel,nnod);
x = [0     L  0  0;
     H1+H2 H2 H2 0];

T=[1 2 3 1;
   2 3 4 3];
 
T1=[1 4 7  1;
    2 5 8  2;
    3 6 9  3;
    4 7 10 7;
    5 8 11 8;
    6 9 12 9];

% Computation of element stiffness matrices

for e = 1:nel
    x1e = x(1,T(1,e));
    x2e = x(1,T(2,e));
    y1e = x(2,T(1,e));
    y2e = x(2,T(2,e));
    le = sqrt((x2e-x1e)^2+(y2e-y1e)^2);
    ke = ((E*I(e))/(le^3))*[1 0    0      1 0     0;
                        0 12   6*le   0 -12   6*le;
                        0 6*le 4*le^2 0 -6*le 2*le^2;
                        1 0    0      1 0     0; 
                        0 -12 -6*le   0 12    -6*le;
                        0 6*le 2*le^2 0 -6*le 4*le^2];
                    
    for r=1:nnod*nig
        for s=1:nnod*nig
            Kel(r,s,e)=ke(r,s);
        end
    end       
end

% Global matrix

KG = zeros(ngl,ngl);
for e = 1:nel
    for i = 1:nnod*nig
        I = T1(i,e);
        for j = 1:nnod*nig
            J = T1(j,e);
            KG(I,J) = KG(I,J) + Kel(i,j,e);
        end
    end
end

% Global system of equations

Vl = zeros(3*n,1);
for i = 1:3*n
    Vl(i) = i;
end
Vr = [10 11 12];
Vl(Vr)=[];
Kll = KG(Vl,Vl);
Klr = KG(Vl,Vr);
Krl = KG(Vr,Vl);
Krr = KG(Vr,Vr);

Fext = zeros(3*n,1);

Fext(5) = P; 

