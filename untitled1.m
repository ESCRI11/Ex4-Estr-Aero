clear;	
clc;

H2= 0.45;   
H1= 0.9;    
L=1.2;      
P=9500;     
A1=12e-6;   
E=60e9;     
nnod=2;
ndim=3;
nel = 4;
N=1;
M=1;
Z=1;
e=N+M+Z+1;
ngl=(N+M+Z+1)*ndim;

x=zeros(1,e);
y=zeros(1,e);
 
for i=2:N+1
    y(i)=(i-1)*H1/(N);
    x(i)=0;
end
 
for i=2:M+1
     y(N+i)=H1+(i-1)*H2/M;
    x(N+i)=0;
end
 
y(N+M+2)=H1;
x(N+M+2)=L;
 
for i=2:Z
    y(N+M+1+i)=H1;
    x(N+M+1+i)=L-(L/Z)*(i-1);
end

T=zeros(2,e);
 
for j=1:e-1
    T(1,j)=j;
    T(2,j)=j+1;
end
    T(1,e)=e;
    T(2,e)=N+1;

    Le=zeros(1,e);
Kel=zeros(ndim*nnod,ndim*nnod,e);
 
for i=1:e-1
    Le(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
end
 
Le(e)=sqrt((x(N+1)-x(e))^2+(y(N+1)-y(e))^2);
Re=zeros(6,6,e);
 
 
for i=1:e
    x1e=x(1,T(1,i)); x2e=x(1,T(2,i));
    y1e=y(1,T(1,i)); y2e=y(1,T(2,i));
    Re(:,:,i)=[-(x1e-x2e)/Le(i) -(y1e-y2e)/Le(i) 0 0 0 0 ;-(y2e-y1e)/Le(i) -(x1e-x2e)/Le(i) 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -(x1e-x2e)/Le(i) -(y1e-y2e)/Le(i) 0;0 0 0 -(y2e-y1e)/Le(i) -(x1e-x2e)/Le(i) 0; 0 0 0 0 0 1];
 
    if i<=N+M
        
        I=3e6/(1e12);
        A=sqrt(12*I);
        Kax=A*E/Le(i)*[1 0 0 -1 0 0;
                       0 0 0  0 0 0;
                       0 0 0  0 0 0;
                      -1 0 0  1 0 0; 
                       0 0 0  0 0 0;
                       0 0 0  0 0 0];
        
       Kmom=(E*I/Le(i)^3).*[0           0           0            0           0           0;
                            0           12          6*Le(i)      0           -12         6*Le(i);
                            0           6*Le(i)     4*Le(i)^2    0           -6*Le(i)    2*Le(i)^2;
                            0           0           0            0           0           0
                            0           -12        -6*Le(i)      0           12          -6*Le(i);
                            0           6*Le(i)     2*Le(i)^2    0           -6*Le(i)    4*Le(i)^2];

                        elseif i>N+M+1
        I=1.5e6/(1e12);
        A=sqrt(12*I);
        Kax=A*E/Le(i)*[1 0 0 -1 0 0;
                       0 0 0  0 0 0;
                       0 0 0  0 0 0;
                      -1 0 0  1 0 0; 
                       0 0 0  0 0 0;
                       0 0 0  0 0 0];
        
       Kmom=(E*I/Le(i)^3).*[0           0           0            0           0           0;
                            0           12          6*Le(i)      0           -12         6*Le(i);
                            0           6*Le(i)     4*Le(i)^2    0           -6*Le(i)    2*Le(i)^2;
                            0           0           0            0           0           0
                            0           -12        -6*Le(i)      0           12          -6*Le(i);
                            0           6*Le(i)     2*Le(i)^2    0           -6*Le(i)    4*Le(i)^2];
 
    else  
        A=A1;
         Kax=A*E/Le(i)*[1 0 0 -1 0 0;
                   0 0 0 0 0 0;
                   0 0 0 0 0 0 ;
                   -1 0 0 1 0 0; 
                   0 0 0 0 0 0;
                   0 0 0 0 0 0];
         I=0;
         
        Kmom=(E*I/Le(i)^3).*[0           0           0            0           0           0;
                             0           12          6*Le(i)      0           -12         6*Le(i);
                             0           6*Le(i)     4*Le(i)^2    0           -6*Le(i)    2*Le(i)^2;
                             0           0           0            0           0           0
                             0           -12        -6*Le(i)      0           12          -6*Le(i);
                             0           6*Le(i)     2*Le(i)^2    0           -6*Le(i)    4*Le(i)^2];
    
    end
    Ke_=Kax+Kmom;
    Ke=Re(:,:,i)'*Ke_*Re(:,:,i);
    
    for r=1:nnod*ndim
        for s=1:nnod*ndim
            
            Kel(r,s,i)=Ke(r,s);
        end
    end
    
end

vL=[4 5 6 7 8 9 10 11 12];
vR=[1 2 3]; 

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

 
KLL=KG(vL,vL);
KLR=KG(vL,vR);
KRL=KG(vR,vL);
KRR=KG(vR,vR);

FextL=[0 0 0 -2340.8 877.8 0 2340.8 -877.8-P 0]'; 
uR=[0 0 0]'; 
uL=inv(KLL)*(FextL-KLR*uR); 
 
u=zeros(ngl,1);
Fext=zeros(ngl,1);
Fext(11)=-P;
u(vL)=uL; 
rx=KG*u-Fext;
 
Esfglob=zeros(6,e);
 
Esfglob(:,1)=Kel(:,:,1)*u(1:6,1);
Esfglob(:,2)=Kel(:,:,2)*u(4:9,1);
Esfglob(:,3)=Kel(:,:,3)*u(7:12,1);
Esfglob(:,4)=Kel(:,:,4)*vertcat(u(10:12,1),u(4:6,1));
 
Esfloc=zeros(6,e);
 
for i=1:e
    Esfloc(:,i)=Re(:,:,i)*Esfglob(:,i);
end
 
 
Axil=zeros(4,100);
Ty=zeros(4,100);
MZ=zeros(4,100);
Long=zeros(4,100);
 
for i=1:e
Axil(i,:)=linspace(-Esfloc(1,i),Esfloc(4,i));
Ty(i,:)=linspace(-Esfloc(2,i),Esfloc(5,i));
if i==1
MZ(i,:)=linspace(-Esfloc(3,i),-Esfloc(3,i));
else
MZ(i,:)=linspace(Esfloc(3,i),Esfloc(6,i));
end
Long(i,:)=linspace(0,Le(i));
end

for i=1:e
    
    figure(i);
    
    subplot(2,2,1); area(Long(i,:),Axil(i,:));ylim auto; xlim([0 Le(i)]);
    xlabel('Longitud');
    ylabel ('Axial [N]');
    str1=sprintf('N  (elemento %d)',i);
    title(str1);
    
    
    subplot(2,2,2); 
    area(Long(i,:),Ty(i,:));
    ylim auto; xlim([0 Le(i)]);
    xlabel('Longitud');
    ylabel ('Cortante [N]');
    str2=sprintf('Ty (elemento %d)',i);
    title(str2);
    
    
    subplot(2,2,3); area(Long(i,:),MZ(i,:));ylim auto; xlim([0 Le(i)]);
    xlabel('Longitud');
    ylabel ('Momento flector [Nm]');
    str3=sprintf('Mz (elemento %d)',i);
    title(str3);
end

