clear; 
clc;
L=1.5;
A=1.5e-5;
H=0.75;
nnod=2;
nn=4;
ndim=3;
nel=3;
ngl=nn*ndim;
E=70e9; 
I=2.5e-6; 
P=12000; 
X=[0 0 L 2*L;
   0 H 0 0];
T=[1 2 3; 
   3 3 4]; 
for i=1:nel
    x1e=X(1,T(1,i)); x2e=X(1,T(2,i));
    y1e=X(2,T(1,i)); y2e=X(2,T(2,i));
    l(:,i)=[x2e-x1e; y2e-y1e];
    Le(i)=norm(l(:,i));
    Lunit(:,i)=l(:,i)/Le(i);
    Re(:,:,i)=[Lunit(1,i) Lunit(2,i) 0 0 0 0 ;
        -Lunit(2,i) Lunit(1,i) 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 Lunit(1,i) Lunit(2,i) 0;
        0 0 0 -Lunit(2,i) Lunit(1,i) 0;
        0 0 0 0 0 1];
    Kax=A*E/Le(i)*[1 0 0 -1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        -1 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    Kmom=(E*I/Le(i)^3).*[0 0 0 0 0 0;
        0 12 6*Le(i) 0 -12 6*Le(i);
        0 6*Le(i) 4*Le(i)^2 0 -6*Le(i) 2*Le(i)^2;
        0 0 0 0 0 0 ;
        0 -12 -6*Le(i) 0 12 -6*Le(i);
        0 6*Le(i) 2*Le(i)^2 0 -6*Le(i) 4*Le(i)^2];
    Ke_=Kax+Kmom; 
    Ke=Re(:,:,i)'*Ke_*Re(:,:,i);
    for r=1:nnod*ndim
        for s=1:nnod*ndim
Kel(r,s,i)=Ke(r,s);
        end
    end
end

KG=zeros(ngl,ngl); 
for el=1:nel
    for a=1:nnod
        for i=1:ndim
            r=ndim*(a-1)+i; 
            A=T(a,el); 
            p=ndim*(A-1)+i; 
            for b=1:nnod
                for j=1:ndim 
                    s=ndim*(b-1)+j;
                    B=T(b,el);
                    q=ndim*(B-1)+j;
                    KG(p,q)=KG(p,q)+Kel(r,s,el);
                end
            end
        end
    end
end

Fext=zeros(ngl,1); 
Fext(11)=-P;

vL=[6 7 8 9 10 11 12]; 
vR=[1 2 3 4 5]; 
FextL=Fext(vL); 
FextR=Fext(vR);

KLL=KG(vL,vL); 
KLR=KG(vL,vR); 
KRL=KG(vR,vL); 
KRR=KG(vR,vR);
uR=[0 0 0 0 0]'; 
uL=(KLL)\(FextL-KLR*uR); 
u=vertcat(uR,uL);
Rx=KG*u-Fext;
Esfglob=zeros(6,nel);
Esfglob(:,1)=Kel(:,:,1)*vertcat(u(1:3,1),u(7:9,1)); 
Esfglob(:,2)=Kel(:,:,2)*vertcat(u(4:6,1),u(7:9,1)); 
Esfglob(:,3)=Kel(:,:,3)*u(7:12,1);
Esfloc=zeros(6,nel);


for i=1:nel 
    Esfloc(:,i)=Re(:,:,i)*Esfglob(:,i);
end

Axil=zeros(4,100);
Ty=zeros(4,100);
MZ=zeros(4,100);
Long=zeros(4,100);


for i=1:nel 
    Axil(i,:)=linspace(-Esfloc(1,i),Esfloc(4,i)); 
    Ty(i,:)=linspace(-Esfloc(2,i),Esfloc(5,i)); 
    MZ(i,:)=-linspace(Esfloc(3,i),Esfloc(6,i)); 
    Long(i,:)=linspace(0,Le(i));
end


for i=1:nel
    figure(i);
    
    subplot(2,2,1); area(Long(i,:),Axil(i,:));ylim auto; xlim([0 Le(i)]); 
    xlabel('Longitud');
    ylabel ('Axil [N]');
    str1=sprintf('N (elemento %d)',i);
    title(str1); 
    
    subplot(2,2,2); area(Long(i,:),Ty(i,:));ylim auto; xlim([0 Le(i)]); 
    xlabel('Longitud');
    str2=sprintf('Ty (elemento %d)',i);
    title(str2);
    ylabel ('Cortante [N]'); 
    
    if i==2
        subplot(2,2,3); area(Long(i,:),-MZ(i,:));ylim auto; xlim([0 Le(i)]);
        xlabel('Longitud');
        ylabel ('Momento flector [Nm]'); 
        str3=sprintf('Mz (elemento %d)',i); 
        title(str3);
    else
        subplot(2,2,3); area(Long(i,:),MZ(i,:));ylim auto; xlim([0 Le(i)]);
        xlabel('Longitud');
        ylabel ('Momento flector [Nm]'); 
        str4=sprintf('Mz (elemento %d)',i); 
        title(str4);
    end
end











