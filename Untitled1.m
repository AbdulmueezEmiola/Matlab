clear; close all;

n = 2; M= 3;
K = 200;
m = [-2 3;10 1;4 -1]';
pw =[0.3; 0.6;0.1];

C = [5 -1;-1 4];
C_ = C^-1;
N = K*M;

for k = 1 : M - 1
    NN(k) = uint16(N * pw(k));
end
NN(M) = N - sum(NN);

label = {'bo', 'r+', 'k*', 'gx'};

IMS=[];
figure; hold on;title('X- Points');
for i=1:M
    ims = mvnrnd(m(:,i),C,NN(i))';
    plot(ims(1,:),ims(2,:),'+');
    IMS = [IMS,ims];
end;
U  =[];
for i = 1:M
    comp_1 = IMS'*C_*m(:,i)-0.5*m(:,i)'*C_*m(:,i) +log(pw(i));    
    U = [U,comp_1];
end
[ui,iai]=max(U,[],2);

figure; hold on;
for i = 1:N
    plot(IMS(1, i), IMS(2, i),label{iai(i)});
end


PIJ = zeros(M);
l0_ = zeros(M);
for i = 1 : M
    for j=i+1:M
        l0_(i,j)=log(pw(j)/pw(i)); 
        h=0.5*(m(:,i)-m(:,j))'*C_*(m(:,i)-m(:,j)); sD=sqrt(2*h);
        PIJ(i,j)=normcdf(l0_(i,j),h,sD); 
        PIJ(j,i)=1-normcdf(l0_(i,j),-h,sD);
    end
    PIJ(i,i)=1-sum(PIJ(i,:));
end;

Pc_=zeros(M);
for i=1:M
    X_test = mvnrnd(m(:,i),C,K)';    
    U  =[];
    for j = 1:M
        comp_1 = X_test'*C_*m(:,j)-0.5*m(:,j)'*C_*m(:,j) +log(pw(j));
        U = [U,comp_1];
    end
    [ui,iai]=max(U,[],2);
    Pc_(i,1)=Pc_(i,1)+ sum(iai==1);
    Pc_(i,2)=Pc_(i,2)+ sum(iai==2);
    Pc_(i,3)=Pc_(i,3)+ sum(iai==3);
end
Pc_ = Pc_/K;
disp("Theoretical Error rates");
disp(PIJ);
disp("Experimental Error rates");
disp(Pc_);


D = C(1,1);
xmin1=-4*sqrt(D)+min(m(1,:)); 
xmax1=4*sqrt(D)+max(m(1,:));
xmin2=-4*sqrt(D)+min(m(2,:));
xmax2=4*sqrt(D)+max(m(2,:));
x1=xmin1:0.05:xmax1; x2=xmin2:0.05:xmax2; 
axis([xmin1,xmax1,xmin2,xmax2]);
figure; hold on; grid on;
[X1,X2]=meshgrid(x1,x2);
x12=[X1(:),X2(:)]; 
for i = 1:M
    f2=mvnpdf(x12,m(:,i)',C);
    f3=reshape(f2,length(x2),length(x1));
    [Ch,h]=contour(x1,x2,f3,[0.01,0.5*max(f3(:))],'Color','b','LineWidth',0.75);    clabel(Ch,h);
    for j=i+1:M,%èçîáðàæåíèå ðàçäåëÿþùèõ ãðàíèö
              wij=C_*(m(:,i)-m(:,j)); 
              wij0=-0.5*(m(:,i)+m(:,j))'*C_*(m(:,i)-m(:,j));
              f4=wij'*x12'+wij0; 
              f5=reshape(f4,length(x2),length(x1));
              [Ch_,h_] = contour(x1,x2,f5,[l0_(i,j)+0.0001]); 
    end
end





