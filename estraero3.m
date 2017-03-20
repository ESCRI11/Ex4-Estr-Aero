clear all
% Material data

EI=3E4;
P1=-40;
P2=10;
q=-15;
m=10;
sm=-50;
d1=7.5;
d2=7.5;
d3=5;

% Geometry

nd = 1;
nel = 40;
n = nel+1; 
nnod = 2;
nig = 2;
ngl = n * nig;
dnod=(d1+d2+d3)/nel;
le=dnod;

x = zeros(n,1);
y = zeros(n,1);
for i=1:n
x(i) = (i-1)*dnod;
end
for i=1:nel
T1(1,i) = i;
T1(2,i) = i+1;
end 
T2 = zeros(nnod*nig,nel);
for i = 1:nel
    for j = 1:nnod
    T2((2*j)-1,i) = (T1(j,i)*2)-1;
    T2((2*j),i) = (T1(j,i)*2);
    end
end

% Computation of element stiffness matrices

for e = 1:nel
    ke = ((EI)/(le^3))*[12 6*le -12 6*le;
                          6*le 4*le^2 -6*le 2*le^2;
                          -12 -6*le    12 -6*le;
                          6*le 2*le^2 -6*le 4*le^2];
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
        I = T2(i,e);
        for j = 1:nnod*nig
            J = T2(j,e);
            KG(I,J) = KG(I,J) + Kel(i,j,e);
        end
    end
end

% Global system of equations

for i=1:n*nig
Vl(i) = i;
end
Vr = [1 2];
Vl(Vr)=[];
Kll = KG(Vl,Vl);
Klr = KG(Vl,Vr);
Krl = KG(Vr,Vl);
Krr = KG(Vr,Vr);
Fext=zeros(length(KG),1);
for i=1:((d1+d2)/dnod)+1
    Fext(i*2-1)=Fext(i*2-1)+(1/2)*q*dnod;
end
for i=2:((d1+d2)/dnod)
    Fext(i*2-1)=Fext(i*2-1)+(1/2)*q*dnod;
end
for i=1:((d1+d2)/dnod)
    Fext(i*2)=Fext(i*2)+(1/12)*q*dnod^2;
end
for i=2:((d1+d2)/dnod)+1
    Fext(i*2)=Fext(i*2)-(1/12)*q*dnod^2;
end

Fext(((d1/dnod)+1)*2-1)=Fext(((d1/dnod)+1)*2-1)+P1;
Fext((((d1+d2)/dnod)*2)+1)=Fext((((d1+d2)/dnod)*2)+1)+P2;
Fext(81)=Fext(81)+sm;
Fext(82)=Fext(82)+m;
Fextsol = Fext(Vl);
ul = Kll\(Fextsol);
R = Krl*ul;
u = zeros(1,ngl);
u(Vl) = ul;

% Shear force

Fshear = zeros((length(Fext)/2)+1,1);
a = 1;
for i = 1:2:(length(Fext)-1)
    Fshear(a,1) = Fext(i);
    a = a + 1;
end

shforce = zeros(n,1);
for i = 1:n
    for j = i:n
        shforce(i) = shforce(i) + Fshear(j+1);
    end
end

% Moment diagram

Mdiag = zeros((length(Fext)/2)+1,1);
a = 1;
for i = 2:2:(length(Fext))
    Mdiag(a,1) = Fext(i);
    a = a + 1;
end

moment = zeros(n,1);
for i = 1:n
    for j = i:n
        moment(i) = moment(i) + ((j-i)*dnod)*Fshear(j+1) + Mdiag(j+1);
    end
end

% Plots

hold all
plot(x(:,1),y);
desplP = zeros(length(u)/2,1);
a = 1;
for i = 1:2:(length(u)-1)
    desplP(a,1) = u(i);
    a = a + 1;
end

scatter(x,desplP);
scatter(x,shforce);
scatter(x,moment);
xlabel('Longitud [m]')
ylabel('Desplazamiento [m]')
ylabel('Momento [Nm]')
ylabel('Esfuerzo cortante [N]')