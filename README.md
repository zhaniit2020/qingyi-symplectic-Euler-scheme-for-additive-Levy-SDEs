# qingyi-symplectic-Euler-scheme-for-additive-Levy-SDEs

N=300;  dt=0.08;  T=20.0;        
sigma=0.2;
sigma1=5.0;
lambda=5.0;

tau=zeros(1,N);
R=zeros(1,N);
tau=exprnd(1/lambda,1,N);    
R=normrnd(0,sigma,[1,N]);    
lambda2=2.0;
P=poissrnd(lambda2*dt,[1,N]);
beta=0.20;

L=zeros(1,N);
L(1,1)=R(1,1);
for i=2:N
    L(1,i)=L(1,i-1)+R(1,i);
end
   t=zeros(1,N);
   t(1,1)=0;                    
    for i=2:N
        t(1,i)=t(1,i-1)+dt;
    end



x00=zeros(1,N);
y00=zeros(1,N);
x00(1,1)=0;
y00(1,1)=1.0;



 for j=2:N    
    for s=2:j   
     x00(1,j)=x00(1,1)*cos(j*dt)+y00(1,1)*sin(j*dt)+beta*trapz(sin((j-s)*dt)*R(1,s),0,j);    
     y00(1,j)=-x00(1,1)*sin(j*dt)+y00(1,1)*cos(j*dt)+beta*trapz(cos((j-s)*dt)*R(1,s),0,j); 
    end
 end

 t0=zeros(1,N);           
    t0(1,1)=0;
    for s=2:N
        t0(1,s)=t0(1,s-1)+tau(1,s-1);
    end
    
     x=cell(1,N);
    y=cell(1,N);           
    for i=1:N
        x{i}=zeros(1,length(0:dt:tau(1,i)));
        y{i}=zeros(1,length(0:dt:tau(1,i)));
    end
 
 x{1}(1,1)=1.0;
 y{1}(1,1)=0; 
for i=1:N
          if  t0(1,i)<=T;               
                    
                     
                      
                  for    j=2:length(0:dt:tau(1,i))
                                 x{i}(1,j)=x{i}(1,j-1)-y{i}(1,j-1)*dt; 
                                 y{i}(1,j)=y{i}(1,j-1)+x{i}(1,j-1)*dt; 
                  end
                             x{i+1}(1,1)=x{i}(1,length(0:dt:tau(1,i)))+beta*normrnd(0,sigma); 
                             y{i+1}(1,1)=y{i}(1,length(0:dt:tau(1,i)));                   
                     
                  end  
                  
                  sum=t0(1,i);             
                   
          if    t0(1,i)>T;       
                  break
         end
end
xe=cell2mat(x);
ye=cell2mat(y);     
    
t1=0:dt:sum;       
number=length(t1); 
x2=xe(1:number);
y2=ye(1:number);   

xx=x00(1:number);
yy=y00(1:number);   


 plot(x2,y2,'--r','LineWidth',0.8);
 xlabel('P','FontSize',16)
 ylabel('Q','FontSize',16)
 hold on

 plot(xx,yy,'-b','LineWidth',0.8);
 xlabel('P','FontSize',16)
 ylabel('Q','FontSize',16)
 legend('numerical solution by EEM', 'the exact solution'); 
 hold on 
 figure
 
plot3(t1,x2,yy,'--r','LineWidth',0.8);
xlabel('Time','FontSize',16)
ylabel('P','FontSize',16) 
zlabel ('Q','FontSize',16)
hold on
  
plot3(t1,xx,yy,'-b','LineWidth',0.8);
xlabel('Time','FontSize',16)
ylabel('P','FontSize',16) 
zlabel ('Q','FontSize',16)
legend('numerical solution by EEM', 'the exact solution'); 
hold on  
figure

plot(t1,x2,'--r');                  
xlabel('Time','FontSize',16)
ylabel('P','FontSize',16)
hold on

plot(t1,yy,'-b');                  
xlabel('Time','FontSize',16)
ylabel('P','FontSize',16)
legend('numerical solution by EEM', 'the exact solution');  
hold off
