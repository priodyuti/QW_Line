clear;
clc;

NSTEPS = 100;
N = 2*NSTEPS+3; %size of the line
mid = NSTEPS+2; %middle point 

%tpath = '/home/priodyuti/Dropbox/Priodyuti/Bar-Ilan//others/Manuscript/quantum_walk/codes/plots/';
tpath = 'C:\Users\user\Documents\MATLAB\prio\quantum_codes\plots\';

out = 'C:\Users\user\Documents\MATLAB\prio\quantum_codes\output_init2.txt';
fd = fopen(out,'wt');

lst1 = 0:15:360;
lst2 = 0:15:180;
lst3 = 0:15:180;
s1 = size(lst1,2);                                                            
s2 = size(lst2,2);
s3 = size(lst3,2);

w = sqrt(-1);
for i1=1:s1
  color = rand(1,3);  
  for j1=1:s2
    for k1=1:s3
     a = lst1(i1); 
     b = lst2(j1); 
     c = lst3(k1);
     
     img = int2str(a)+"_"+int2str(b)+"_"+int2str(c)+".png";
     imgpath = strcat(tpath,img);
     
     %S = zeros(1,N);
     
     %%Coin operator
     %al = 45;
     alpha = (pi*a)/180;
     phi1 = (pi*b)/180;
     phi2 = (pi*c)/180;
     
     C=zeros(2,2);
     C(1,1) = cos(alpha);  
     C(1,2) = exp(w*phi1)*sin(alpha);
     C(2,1) = exp(w*phi2)*sin(alpha);
     C(2,2) = -exp(w*(phi1+phi2))*cos(alpha);

     %Amplitude Matrix
     At = zeros(1,N,2);
     A = zeros(1,N,2);

     %%Initial Coinditions of Amplitude Matrix
     %At(1,mid,1) = sqrt(0.5);
     %At(1,mid,2) = -w*sqrt(0.5);
     At(1,mid,1) = 0;
     At(1,mid,2) = 1;


     %%Amplitude Calculation Loop
     for k= 1:NSTEPS
       for m = 2:2*NSTEPS+2
         A(1,m,1) = C(1,1)*At(1,m-1,1) + C(1,2)*At(1,m-1,2);
         A(1,m,2) = C(2,1)*At(1,m+1,1) + C(2,2)*At(1,m+1,2);
       end 
       At=A;
     end

     %%Probability Density Function
     P = zeros(1,N);

     %%Probability Calculation Loop
     for m = 1:N         
      P(m) = At(1,m,1)*conj(At(1,m,1))+ At(1,m,2)*conj(At(1,m,2)); 
      x(m)= -((N+1)/2)+ m;
     end

     QRWtotal_probability =sum(sum(P));
     %plot(x,P)
     iterate_draw_plot_QW(x,P,figure(1),imgpath,a,b,c,color)
     st = std(P);
     P = sort(P,'descend');
     fprintf(fd,'%d %d %d %f %f\n',a,b,c,st,P(1)-P(2));
     %fprintf(fd,'%d %d %f %f\n',a,b+c,st,P(1)-P(2));      
     clear C At A P x
     
     
     
%      cls1=1; cls2=2; cls3=3; cls4=4; cls5=5; cls6=6; cls7=7; cls8=8;
%      
%      if a==0 || a==180
%        fprintf(fd,'%d %d %d %d\n',a,b,c,cls1);  
%      elseif  a==90 || a==270
%        fprintf(fd,'%d %d %d %d\n',a,b,c,cls2);
%      elseif  a==30 || a==150 || a==210 || a==330
%          if (a==30 && b==0) || (a==150 && b==0) || (a==210 && b==0) || (a==330 && b==0)
%             fprintf(fd,'%d %d %d %d\n',a,b,c,cls3);
%          elseif (a==30 && b>0) || (a==210 && b>0)
%             fprintf(fd,'%d %d %d %d\n',a,b,c,cls4);
%          elseif (a==150 && b>0) || (a==330 && b>0)
%             fprintf(fd,'%d %d %d %d\n',a,b,c,cls5);  
%          end    
%      elseif  a==60 || a==120 || a==240 || a==300
%          if (a==60 && b==0) || (a==120 && b==0) || (a==240 && b==0) || (a==300 && b==0)
%             fprintf(fd,'%d %d %d %d\n',a,b,c,cls6);
%          elseif (a==60 && b>0) || (a==240 && b>0)
%             fprintf(fd,'%d %d %d %d\n',a,b,c,cls7);
%          elseif (a==120 && b>0) || (a==300 && b>0)
%             fprintf(fd,'%d %d %d %d\n',a,b,c,cls8);  
%          end    
%        
%      end        
   end
  end
end
clear 