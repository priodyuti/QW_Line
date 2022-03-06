%%Code Ref: Mahesh N Jayakody and Asiri Nanayakkara, Full state revivals in higher dimensional quantum walks, 
%Phys. Scr. 94 (2019) 045101.
% E-mail: asiri@ifs.ac.lk


clear;
clc;

NSTEPS = 5;
N = 2*NSTEPS+3; %size of the line
mid = NSTEPS+2; %middle point 

w=sqrt(-1);
%S=zeros(1,N);

%Coin operator
%al=45;
%alpha=(pi*al)/180;

%C(1,1)=cos(alpha);            C(1,2)=sin(alpha);
%C(2,1)=sin(alpha);            C(2,2)=-cos(alpha);

a = 45;
b = 0;
c = 0;

alpha = (pi*a)/180;
phi1 = (pi*b)/180;
phi2 = (pi*c)/180;
     
C=zeros(2,2);
C(1,1) = cos(alpha);  
C(1,2) = exp(w*phi1)*sin(alpha);
C(2,1) = exp(w*phi2)*sin(alpha);
C(2,2) = -exp(w*(phi1+phi2))*cos(alpha);

%Amplitude Matrix
At = zeros(2,N);
A = zeros(2,N);



%Initial Coinditions of Amplitude Matrix

At(1,mid) = sqrt(0.5);
At(2,mid) = -w*sqrt(0.5);
% % At(1,mid,1) = 1;
% At(1,mid,2) = 1;

%NSTEPS=1;
%Amplitude Calculation Loop
for k= 1:NSTEPS
  for m = 2:2*NSTEPS+2
    A(1,m)=C(1,1)*At(1,m-1)+C(1,2)*At(2,m-1);
    A(2,m)=C(2,1)*At(1,m+1)+C(2,2)*At(2,m+1);
  end
  At = A;
end

%Probability Density Function
P=zeros(1,N);

%Probability Calculation Loop
for m = 1:N         
P(m) = At(1,m)*conj(At(1,m))+ At(2,m)*conj(At(2,m)); 
x(m)= -((N+1)/2)+ m;
end

QRWtotal_probability = sum(sum(P));
plot(x,P,'ro-')
xlabel('x_i','fontweight','bold'); 
ylabel('P(t) ','fontweight','bold');
%title('Numerical perturbation tracking')
set(gca,'FontSize',30)


P = sort(P,'descend');
for m = 1:2         
%fprintf('%0.20f\n', P(m));  
end
%P(1)-P(2)

s = std(P);
%fprintf('std = %0.10f\n', s);