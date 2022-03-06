


clear;
clc;

NSTEPS = 5;
N = 2*NSTEPS+3; %size of the line
mid =NSTEPS+2; %middle point 

w=sqrt(-1);

a = 315;
b = 0;
c = 0;

output_path = "C:\Users\user\Documents\MATLAB\prio\quantum_codes\images\data\theta_vary_"+int2str(a)+"_"+int2str(b)+"_"+int2str(c)+".txt";
fd = fopen(output_path,'wt');

alpha = (pi*a)/180;
phi1 = (pi*b)/180;
phi2 = (pi*c)/180;
     
C=zeros(2,2);
C(1,1) = cos(alpha);  
C(1,2) = exp(w*phi1)*sin(alpha);
C(2,1) = exp(w*phi2)*sin(alpha);
C(2,2) = -exp(w*(phi1+phi2))*cos(alpha);



%Amplitude Matrix
At=zeros(1,N,2);
A=zeros(1,N,2);

%Initial Coinditions of Amplitude Matrix

At(1,mid,1) = sqrt(0.5);
At(1,mid,2) = -w*sqrt(0.5);
% % At(1,mid,1) = 1;
% At(1,mid,2) = 1;

%NSTEPS=1;
%Amplitude Calculation Loop
for k= 1:NSTEPS
  for m = 2:2*NSTEPS+2
    A(1,m,1)=C(1,1)*At(1,m-1,1)+C(1,2)*At(1,m-1,2);
    A(1,m,2)=C(2,1)*At(1,m+1,1)+C(2,2)*At(1,m+1,2);
  end
  At=A;
end

%Probability Density Function
P=zeros(1,N);

%Probability Calculation Loop
for m = 1:N         
P(m) = At(1,m,1)*conj(At(1,m,1))+ At(1,m,2)*conj(At(1,m,2)); 
x(m)= -((N+1)/2)+ m;
end


QRWtotal_probability =sum(sum(P));
plot(x,P,'r-')

n = size(P,2);

for i=1:n
   fprintf(fd,'%d %f\n',x(i),P(i)); 
end    
fclose(fd);

P = sort(P,'descend');
for m = 1:2         
fprintf('%0.20f\n', P(m));  
end
%P(1)-P(2)

%s = std(P);
%fprintf('std = %0.10f\n', s);