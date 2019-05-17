clear all
clc
% A1=[    0.0100    0.0029
%        -0.0600   -0.0029];
% B1=[   -0.0198         0
%         0.1190    0.0100];
% A2 =[0.0160    0.0010
%     -0.0660   -0.0010];
% B2=[-0.0148         0
%     0.0612    0.0160];
% A3=[0.0050    0.0027
%    -0.0550   -0.0027];
% B3=[-0.0135         0
%     0.1478    0.0050];
A1=[   0.006497175141243   0.002987048740783
      -0.056497175141243  -0.002987048740783];
B1=[  -0.016100000000000                   0
       0.140000000000000   0.006480000000000];
A2=[  -0.018000000000000                   0
       0.080000000000000   0.014550000000000];
B2=[  -0.019800000000000                   0
       0.118900000000000   0.010000000000000];
A3=[   0.011250805931657   0.002682863866452
      -0.061250805931657  -0.002682863866452];
B3=[  -0.020200000000000                   0
       0.110000000000000   0.011250000000000];
A4=[   0.012500000000000   0.002343750000000
      -0.062500000000000  -0.002343750000000];
B4=[  -0.020000000000000                   0
       0.100000000000000   0.012500000000000];
A5=[   0.014553990610329   0.001599385042650
      -0.064553990610329  -0.001599385042650];
B5=[  -0.018000000000000                   0
       0.080000000000000   0.014550000000000];
C=[0 1];
D=[0 0];
Ts=0.05;
[n1,d1]=ss2tf(A1,B1,C,D,1);
[n2,d2]=ss2tf(A1,B1,C,D,2);
g11 = tf(n1,d1);
gz11 = c2d(g11,Ts,'impulse');
g21 = tf(n2,d2);
gz21 = c2d(g21,Ts,'impulse');

[n12,d12]=ss2tf(A2,B2,C,D,1);
[n22,d22]=ss2tf(A2,B2,C,D,2);
g12 = tf(n12,d12);
gz12 = c2d(g12,Ts,'impulse');
g22 = tf(n22,d22);
gz22 = c2d(g22,Ts,'impulse');

[n13,d13]=ss2tf(A3,B3,C,D,1);
[n23,d23]=ss2tf(A3,B3,C,D,2);
g13 = tf(n13,d13);
gz13 = c2d(g13,Ts,'impulse');
g23 = tf(n23,d23);
gz23 = c2d(g23,Ts,'impulse');

[n14,d14]=ss2tf(A4,B4,C,D,1);
[n24,d24]=ss2tf(A4,B4,C,D,2);
g14 = tf(n14,d14);
gz14 = c2d(g14,Ts,'impulse');
g24 = tf(n24,d24);
gz24 = c2d(g24,Ts,'impulse');

[n15,d15]=ss2tf(A5,B5,C,D,1);
[n25,d25]=ss2tf(A5,B5,C,D,2);
g15 = tf(n15,d15);
gz15 = c2d(g15,Ts,'impulse');
g25 = tf(n25,d25);
gz25 = c2d(g25,Ts,'impulse');

[num11,den11] = tfdata(gz11,'v');
[num21,den21] = tfdata(gz21,'v');
[num12,den12] = tfdata(gz12,'v');
[num22,den22] = tfdata(gz22,'v');
[num13,den13] = tfdata(gz13,'v');
[num23,den23] = tfdata(gz23,'v');
[num14,den14] = tfdata(gz14,'v');
[num24,den24] = tfdata(gz24,'v');
[num15,den15] = tfdata(gz15,'v');
[num25,den25] = tfdata(gz25,'v');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for model=1:5
P=20;
M=8;
na=2;
nb1=1;
nb2=1;
nb=nb1;
d=0;
N1=d+1;
N2=d+P;
a_tilt = [1 -2 1];
f=zeros(P+d,na+1);
f(1,1:2)=-1*a_tilt(2:3);
for j=1:P+d-1
   for i=1:na
    f(j+1,i)=f(j,i+1)-f(j,1)*a_tilt(i+1);
   end
end
F=f(N1:N2,1:na);
% %.......................................
E=zeros(P,P);
E(:,1)=1;
for j=1:P-1
    E(j+1:P,j+1)=f(j,1);
end
bb1=[0.007 0.005945 0.0055 0.005 0.004];
Eaux11=zeros(P,P+1);
M1=zeros(P,P+1);
Mi1=zeros(P,P);
Mprim1=zeros(P,1);
Eaux1=zeros(P,P);
Eaux1=E*bb1(model);
Eaux11(:,1:P)=Eaux1;%%b11
M1=Eaux11;
Mi1=toeplitz(M1(:,1),[M1(1,1) zeros(P-1,1)']);
if M~=P
    if M<P
    Mi1(:,M)=Mi1(:,M:P)*ones(P-M+1,1);
    Mi1=Mi1(:,1:M);
    else
        M=P;
    end
end
    
for i=1:P
 Mprim1(i)=M1(i,i+1);   
end

bb2=[0.000324 0.0005 0.0005625 0.000625 0.000727];
Eaux11=zeros(P,P+1);
M2=zeros(P,P+1);
Mi2=zeros(P,P);
Mprim2=zeros(P,1);
Eaux1=zeros(P,P);
Eaux1=E*bb2(model);
Eaux11(:,1:P)=Eaux1;%%b21
M2=Eaux11;
Mi2=toeplitz(M2(:,1),[M2(1,1) zeros(P-1,1)']);
if M~=P
    if M<P
    Mi2(:,M)=Mi2(:,M:P)*ones(P-M+1,1);
    Mi2=Mi2(:,1:M);
    else
        M=P;
    end
end
for i=1:P
 Mprim2(i)=M2(i,i+1); 
end
Mi=[Mi1 Mi2];
Mprim=[Mprim1 Mprim2];
if model==1
    Mi_1=Mi;
    Mprim_1=Mprim;
    
    Q=1*eye(P,P);
    gamma1=1;gamma2=1;
    R1=gamma1*(0.1)^2*eye(M,M);
    R2=gamma2*(0.1)^2*eye(M,M);  
    R=zeros(2*M,2*M);
    R(1:M,1:M)=R1;
    R(M+1:2*M,M+1:2*M)=R2;
    a=Mi_1'*Q*Mi_1+R;
    b=Mi_1'*Q;
    k_gpc1 = a\b; %%%%%%%constant
elseif model==2
    Mi_2=Mi;
    Mprim_2=Mprim;
    
    Q=1*eye(P,P);
    gamma1=0.1;gamma2=0.1;
    R1=gamma1*(1)^2*eye(M,M);
    R2=gamma2*(1)^2*eye(M,M);  
    R=zeros(2*M,2*M);
    R(1:M,1:M)=R1;
    R(M+1:2*M,M+1:2*M)=R2;
    a=Mi_2'*Q*Mi_2+R;
    b=Mi_2'*Q;
    k_gpc2 = a\b; %%%%%%%constant
elseif model==3
    Mi_3=Mi;
    Mprim_3=Mprim;
    
    Q=1*eye(P,P);
    gamma1=0.1;gamma2=0.1;
    R1=gamma1*(1)^2*eye(M,M);
    R2=gamma2*(1)^2*eye(M,M);  
    R=zeros(2*M,2*M);
    R(1:M,1:M)=R1;
    R(M+1:2*M,M+1:2*M)=R2;
    a=Mi_3'*Q*Mi_3+R;
    b=Mi_3'*Q;
    k_gpc3=a\b; %%%%%%%constant
elseif model==4
    Mi_4=Mi;
    Mprim_4=Mprim;
    
    Q=1*eye(P,P);
    gamma1=0.1;gamma2=0.1;
    R1=gamma1*(1)^2*eye(M,M);
    R2=gamma2*(1)^2*eye(M,M);  
    R=zeros(2*M,2*M);
    R(1:M,1:M)=R1;
    R(M+1:2*M,M+1:2*M)=R2;
    a=Mi_4'*Q*Mi_4+R;
    b=Mi_4'*Q;
    k_gpc4=a\b; %%%%%%%constant
else
    Mi_5=Mi;
    Mprim_5=Mprim;
    
    Q=1*eye(P,P);
    gamma1=0.1;gamma2=0.1;
    R1=gamma1*(1)^2*eye(M,M);
    R2=gamma2*(1)^2*eye(M,M);  
    R=zeros(2*M,2*M);
    R(1:M,1:M)=R1;
    R(M+1:2*M,M+1:2*M)=R2;
    a=Mi_5'*Q*Mi_5+R;
    b=Mi_5'*Q;
    k_gpc5=a\b; %%%%%%%constant
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaU1=zeros(2*M,1);
deltaU2=zeros(2*M,1);
deltaU3=zeros(2*M,1);
deltaU4=zeros(2*M,1);
deltaU5=zeros(2*M,1);
deltaUp=zeros(2*M,1);
deltaUp1=zeros(2*M,1);
deltaUp2=zeros(2*M,1);
deltaUp3=zeros(2*M,1);
deltaUp4=zeros(2*M,1);
deltaUp5=zeros(2*M,1);
U=zeros(2*M,1);
deltaUprim1=zeros(2,1);
deltaUprim2=zeros(2,1);
deltaUprim3=zeros(2,1);
deltaUprim4=zeros(2,1);
deltaUprim5=zeros(2,1);
Yp1(1)=0;
Yp2(1)=0;
Yp3(1)=0;
Yp4(1)=0;
Yp5(1)=0;
Yp(1)=0;
Ypp(1)=0;
Ym=zeros(P,1,5);
Yd1=zeros(P,1);
Yd2=zeros(P,1);
Yd3=zeros(P,1);
Yd4=zeros(P,1);
Yd5=zeros(P,1);
Yprim1=zeros(2,1);
Yprim2=zeros(2,1);
Yprim3=zeros(2,1);
Yprim4=zeros(2,1);
Yprim5=zeros(2,1);
Ypast=zeros(P,5);
d=0;
d1=0;
d2=0;
d3=0;
d4=0;
d5=0;
t1=1;

t=0:0.05:1;
yy=[0.0161    0.0198    0.0202    0.0200    0.0180
    0.1600    0.1811    0.1900    0.2000    0.2200];
IC=[0.019821 0.18107];
% IC=[0.0173    0.2099];
IC1=[0.0161 0.16];
IC2=[0.0198 0.1811];
IC3=[0.0202 0.19];
IC4=[0.02 0.2];
IC5=[0.018 0.22];

alpha=0.5;
ki=0.1;
ks=1;
miomax=0.5;
select=12;
Ym_avr(1)=0;
in_init=[0.00648 0.01 0.01125 0.0125 0.01455;0.3 0.3 0.3 0.3 0.3];
e=zeros(5,1);
e(1)=0.5;
e(2)=0.5;
p=zeros(5,1);
p(1)=0.5;
p(2)=0.5;
ee1=0;
ee2=0;
ee3=0;
ee4=0;
ee5=0;
tsim=1:t1:550;
flag1=0;
flag2=0;
flag3=0;
% [r,t3]=gensig('square',8,20,0.05);
% r=0.03*r+0.15;
r=[0.181*ones(1,150) 0.19*ones(1,150) 0.20*ones(1,150) 0.22*ones(1,150)];
for j=1:length(tsim)
    rd=r(j:j+P-1);
        if j==1
           Yd1(1)=0; 
           Yd2(1)=0;
           Yd3(1)=0; 
           Yd4(1)=0;
           Yd5(1)=0; 
        end
    for i=1:P-1 
        Yd1(i+1)=alpha*Yd1(i)+(1-alpha)*(rd(i)-IC1(2));
    end
    for i=1:P-1 
        Yd2(i+1)=alpha*Yd2(i)+(1-alpha)*(rd(i)-IC2(2));
    end
    for i=1:P-1 
        Yd3(i+1)=alpha*Yd3(i)+(1-alpha)*(rd(i)-IC3(2));
    end
    for i=1:P-1 
        Yd4(i+1)=alpha*Yd4(i)+(1-alpha)*(rd(i)-IC4(2));
    end
    for i=1:P-1 
        Yd5(i+1)=alpha*Yd5(i)+(1-alpha)*(rd(i)-IC5(2));
    end
    deltaU1(2:M)=deltaU1(1:M-1);
    deltaU1(M+2:2*M)=deltaU1(M+1:2*M-1);
    deltaU1(1)=deltaUp1(1);
    deltaU1(M+1)=deltaUp1(M+1);
    deltaUprim1(1)=deltaUp1(1);
    deltaUprim1(2)=deltaUp1(M+1);
    
    Yprim1(2)=Yprim1(1);
    Yprim1(1) = Ym(1,j,1);
    
    deltaU2(2:M)=deltaU2(1:M-1);
    deltaU2(M+2:2*M)=deltaU2(M+1:2*M-1);
    deltaU2(1)=deltaUp2(1);
    deltaU2(M+1)=deltaUp2(M+1);
    deltaUprim2(1)=deltaUp2(1);
    deltaUprim2(2)=deltaUp2(M+1);
    
    Yprim2(2)=Yprim2(1);
    Yprim2(1) =Ym(1,j,2);
    
    deltaU3(2:M)=deltaU3(1:M-1);
    deltaU3(M+2:2*M)=deltaU3(M+1:2*M-1);
    deltaU3(1)=deltaUp3(1);
    deltaU3(M+1)=deltaUp3(M+1);
    deltaUprim3(1)=deltaUp3(1);
    deltaUprim3(2)=deltaUp3(M+1);
    
    Yprim3(2)=Yprim3(1);
    Yprim3(1) = Ym(1,j,3);
    
    deltaU4(2:M)=deltaU4(1:M-1);
    deltaU4(M+2:2*M)=deltaU4(M+1:2*M-1);
    deltaU4(1)=deltaUp4(1);
    deltaU4(M+1)=deltaUp4(M+1);
    deltaUprim4(1)=deltaUp4(1);
    deltaUprim4(2)=deltaUp4(M+1);
    
    Yprim4(2)=Yprim4(1);
    Yprim4(1) = Ym(1,j,4);
    
    deltaU5(2:M)=deltaU5(1:M-1);
    deltaU5(M+2:2*M)=deltaU5(M+1:2*M-1);
    deltaU5(1)=deltaUp5(1);
    deltaU5(M+1)=deltaUp5(M+1);
    deltaUprim5(1)=deltaUp5(1);
    deltaUprim5(2)=deltaUp5(M+1);
    
    Yprim5(2)=Yprim5(1);
    Yprim5(1) = Ym(1,j,5);
    
    Ypast(:,1)=F*Yprim1+Mprim_1*deltaUprim1;
    Ypast(:,2)=F*Yprim2+Mprim_2*deltaUprim2;
    Ypast(:,3)=F*Yprim3+Mprim_3*deltaUprim3;
    Ypast(:,4)=F*Yprim4+Mprim_4*deltaUprim4;
    Ypast(:,5)=F*Yprim5+Mprim_5*deltaUprim5;
    
    E_1=Yd1-Ypast(:,1)-d1*ones(P,1);
    E_2=Yd2-Ypast(:,2)-d2*ones(P,1);
    E_3=Yd3-Ypast(:,3)-d3*ones(P,1);
    E_4=Yd4-Ypast(:,4)-d4*ones(P,1);
    E_5=Yd5-Ypast(:,5)-d5*ones(P,1);
    
    deltaUp1=k_gpc1*E_1;
    deltaUp2=k_gpc2*E_2;
    deltaUp3=k_gpc3*E_3;
    deltaUp4=k_gpc4*E_4;
    deltaUp5=k_gpc5*E_5;
    
         deltaUp=p(1)*deltaUp1+p(2)*deltaUp2+p(3)*deltaUp3+p(4)*deltaUp4+p(5)*deltaUp5;
         U=deltaUp+U;
         in(j,:)=[U(1)+p(1)*in_init(1,1)+p(2)*in_init(1,2)+p(3)*in_init(1,3)+p(4)*in_init(1,4)+p(5)*in_init(1,5) U(M+1)+in_init(2,1)];
    
%     if select==23
%        deltaUp=e3*deltaUp3+e2*deltaUp2;
%        U=deltaUp+U;
%        in(j,:)=[U(1)+e3*in_init(1,3)+e2*in_init(1,2) U(M+1)+in_init(2,1)];
%     end
%     if select==13
%         deltaUp=e1*deltaUp1+e3*deltaUp3;
%         U=deltaUp+U;
%         in(j,:)=[U(1)+e3*in_init(1,3)+e1*in_init(1,1) U(M+1)+in_init(2,1)];
%     end
%     if select==12
%         deltaUp=e1*deltaUp1+e2*deltaUp2;
%         U=deltaUp+U;
%         in(j,:)=[U(1)+e1*in_init(1,1)+e2*in_init(1,2) U(M+1)+in_init(2,1)];
%     end
    [T,X]=ode45(@(t,x) bioreact(t,x,in(j,:),miomax,ki,ks),[0 0.05],IC);
    Yp1(j+1) = X(end,2)-IC1(2);
    Yp2(j+1) = X(end,2)-IC2(2);
    Yp3(j+1) = X(end,2)-IC3(2);
    Yp4(j+1) = X(end,2)-IC4(2);
    Yp5(j+1) = X(end,2)-IC5(2);
    Ym(:,j+1,1) = Mi_1*deltaU1+Ypast(:,1);
    Ym(:,j+1,2) = Mi_2*deltaU2+Ypast(:,2);
    Ym(:,j+1,3) = Mi_3*deltaU3+Ypast(:,3);
    Ym(:,j+1,4) = Mi_4*deltaU4+Ypast(:,4);
    Ym(:,j+1,5) = Mi_5*deltaU5+Ypast(:,5);
    IC = X(end,:);
    d1 = Yp1(j+1)-Ym(1,j+1,1);
    d2 = Yp2(j+1)-Ym(1,j+1,2);
    d3 = Yp3(j+1)-Ym(1,j+1,3);
    d4 = Yp4(j+1)-Ym(1,j+1,4);
    d5 = Yp5(j+1)-Ym(1,j+1,5);
    if d1==0
     d1=d1+0.01;
     flag1=flag1+1;
    end
    if d2==0
     d2=d2+0.01;  
     flag2=flag2+1;
    end
    if d3==0
     d3=d3+0.01;
     flag3=flag3+1;
    end
    if d4==0
     d4=d+0.01;
     flag3=flag3+1;
    end
    if d5==0
     d5=d5+0.01;
     flag3=flag3+1;
    end
    sum= abs(1/d1)+abs(1/d2)+abs(1/d3)+abs(1/d4)+abs(1/d5);
    e(1) = abs(1/d1)/sum;
    e(2) = abs(1/d2)/sum;
    e(3) = abs(1/d3)/sum;
    e(4) = abs(1/d4)/sum;
    e(5) = abs(1/d5)/sum;
    Yd1=zeros(P,1); 
    Yd2=zeros(P,1); 
    Yd3=zeros(P,1); 
    Yd4=zeros(P,1); 
    Yd5=zeros(P,1); 
    [e_sort,e_ind]=sort([e(1) e(2) e(3) e(4) e(5)]);
        v2=e_ind(2);
        v1=e_ind(1);
        sum=e(1)+e(2)+e(3)+e(4)+e(5)-e(v1)-e(v2);
        p(1)=e(1)/sum;
        p(2)=e(2)/sum;
        p(3)=e(3)/sum;
        p(4)=e(4)/sum;
        p(5)=e(5)/sum;
        p(v1)=0;
        p(v2)=0;
       sw_max(j+1)=e_ind(5);
       sw_off1(j+1)=e_ind(1);
       sw_off2(j+1)=e_ind(2);
%     if e_sort(1)==e(1)
%         v=e_ind(2);
%         sum=e(2)+e(3)+e(4)+e(5)-e(v);
%         p(1)=0;
%         p(2)=e(2)/sum;
%         p(3)=e(3)/sum;
%         p(4)=e(4)/sum;
%         p(5)=e(5)/sum;
%         p(v)=0; 
%        sw_max(j+1)=e_ind(5);
%        sw_off(j+1)=e_ind(1);
%     else if e_sort(1)==e(2)
%             v=e_ind(2);
%             sum=e(1)+e(3)+e(4)+e(5)-e(v);
%             p(2)=0;
%             p(1)=e(1)/sum;
%             p(3)=e(3)/sum;
%             p(4)=e(4)/sum;
%             p(5)=e(5)/sum;
%             p(v)=0; 
%             sw_max(j+1)=e_ind(5);
%             sw_off(j+1)=e_ind(1);
%         else if e_sort(1)==e(3)
%                 v=e_ind(2);
%                 sum=e(1)+e(2)+e(4)+e(5)-e(v);
%                 p(3)=0;
%                 p(1)=e(1)/sum;
%                 p(2)=e(2)/sum;
%                 p(4)=e(4)/sum;
%                 p(5)=e(5)/sum;
%                 p(v)=0; 
%                 sw_max(j+1)=e_ind(5);
%                 sw_off(j+1)=e_ind(1);
%             else if e_sort(1)==e(4)
%                     v=e_ind(2);
%                     sum=e(1)+e(2)+e(3)+e(5)-e(v);
%                     p(4)=0;
%                     p(1)=e(1)/sum;
%                     p(2)=e(2)/sum;
%                     p(3)=e(3)/sum;
%                     p(5)=e(5)/sum;
%                     p(v)=0; 
%                     sw_max(j+1)=e_ind(5);
%                     sw_off(j+1)=e_ind(1);
%                 else if e_sort(1)==e(5)
%                         v=e_ind(2);
%                         sum=e(1)+e(2)+e(4)+e(3)-e(v);
%                         p(5)=0;
%                         p(1)=e(1)/sum;
%                         p(2)=e(2)/sum;
%                         p(4)=e(4)/sum;
%                         p(3)=e(3)/sum; 
%                         p(v)=0;
%                         sw_max(j+1)=e_ind(5);
%                         sw_off(j+1)=e_ind(1);
%                     end
%                 end
%             end
%         end
%     end
       Yp(j+1)=p(1)*(Yp1(j+1)+IC1(2))+p(2)*(Yp2(j+1)+IC2(2))+p(3)*(Yp3(j+1)+IC3(2))+p(4)*(Yp4(j+1)+IC4(2))+p(5)*(Yp5(j+1)+IC5(2));
       Ypp(j+1)=p(1)*(Yp1(j+1))+p(2)*(Yp2(j+1))+p(3)*(Yp3(j+1))+p(4)*(Yp4(j+1))+p(5)*(Yp5(j+1));
       Yd1(1)=p(1)*Yp1(j+1)+p(2)*Yp2(j+1)+p(3)*Yp3(j+1)+p(4)*Yp4(j+1)+p(5)*Yp5(j+1);
       Yd2(1)=p(1)*Yp1(j+1)+p(2)*Yp2(j+1)+p(3)*Yp3(j+1)+p(4)*Yp4(j+1)+p(5)*Yp5(j+1);
       Yd3(1)=p(1)*Yp1(j+1)+p(2)*Yp2(j+1)+p(3)*Yp3(j+1)+p(4)*Yp4(j+1)+p(5)*Yp5(j+1);
       Yd4(1)=p(1)*Yp1(j+1)+p(2)*Yp2(j+1)+p(3)*Yp3(j+1)+p(4)*Yp4(j+1)+p(5)*Yp5(j+1);
       Yd5(1)=p(1)*Yp1(j+1)+p(2)*Yp2(j+1)+p(3)*Yp3(j+1)+p(4)*Yp4(j+1)+p(5)*Yp5(j+1);
    ee1=[ee1;e(1)];
    ee2=[ee2;e(2)];
    ee3=[ee3;e(3)];
    ee4=[ee4;e(4)];
    ee5=[ee5;e(5)];
    ee=[ee1 ee2 ee3 ee4 ee5];
        
end
plot(Yp,'r--')
hold on
plot(r,'g')
plot(Ym_avr)
plot(0.1*sw_max)
plot(0.1*sw_off1,'r')
plot(0.1*sw_off2,'r--')
% plot(Ypp,'r--')
% hold on
% plot(Ym(1,:,1),'b--')
% plot(Ym(1,:,2),'c--')
% plot(Ym(1,:,3),'m--')
