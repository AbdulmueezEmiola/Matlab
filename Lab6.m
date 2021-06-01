clear all; close all;
K = 200;
n=2;M=3;
dm = 2.0;
C = [5 -1;-1 4];
C_ = C^-1;
pw =[0.3; 0.6;0.1];
np=sum(pw); pw=pw/np;
m = [-2 3;10 1;4 -1]';

N = K*M;
for k = 1 : M - 1
    NN(k) = uint16(N * pw(k));
end
NN(M) = N - sum(NN);

for i=1:M
    ims = mvnrnd(m(:,i),C,NN(i))';
    IMS{i} = ims;
end;

%Определение вероятностей ошибок методом скользящего контроля
r = 0.25;
Pc1 = zeros(M);
p1_ = zeros(M,1);
for i = 1:M
    N = double(NN(i));
    XNi = IMS{i};
    XNi_ = zeros(n,N-1);
    indi=[1:i-1,i+1:M];
    for j=1:N
        x=XNi(:,j); indj=[1:j-1,j+1:N];
        XNi_(:,1:j-1)=XNi(:,1:j-1); XNi_(:,j:end)=XNi(:,j+1:end);
        kn=2*round(N^r)+1;
        p1_(i)=vknn(x,IMS{i},kn);
        for t=1:M-1
            ij=indi(t);
            p1_(ij)=vknn(x,IMS{ij},kn);
        end
        [ui1,iai1]=max(p1_); Pc1(i,iai1)=Pc1(i,iai1)+1;
    end
    Pc1(i,:)=Pc1(i,:)/N;
end

%Тестирование алгоритма методом статистических испытаний
Pcv=zeros(M); 
p=zeros(M,1); 
 
x=ones(n,1); 
u=zeros(M,1);
Pc_=zeros(M);
for k=1:K
    for i=1:M
        [x,px]=randncor(n,1,C);
        x=x+m(:,i);
        for j=1:M
            u(j)=-0.5*(x-m(:,j))'*C_*(x-m(:,j))-0.5*log(det(C))+log(pw(j)); 
            kn=2*round(N^r)+1;
            p(j)=vknn(x,IMS{j},kn);
        end
        [ui,iai]=max(u);
        Pc_(i,iai)=Pc_(i,iai)+1;%фиксация результата распознавания
        [ui,iai]=max(p); % + %определение максимума
        Pcv(i,iai)=Pcv(i,iai)+1;% + %фиксация результата распознавания
    end
end
Pc_=Pc_/K;
Pcv=Pcv/K;

%Расчет матриц вероятностей ошибок распознавания
PTheory = zeros(M);
l0_ = zeros(M);
for i = 1 : M
    for j=i+1:M
        l0_(i,j)=log(pw(j)/pw(i)); 
        h=0.5*(m(:,i)-m(:,j))'*C_*(m(:,i)-m(:,j)); sD=sqrt(2*h);
        PTheory(i,j)=normcdf(l0_(i,j),h,sD); 
        PTheory(j,i)=1-normcdf(l0_(i,j),-h,sD);
    end
    PTheory(i,i)=1-sum(PTheory(i,:));
end

PExp_=zeros(M);
for i=1:M
    U  =[];
    for j = 1:M
        comp_1 = IMS{i}'*C_*m(:,j)-0.5*m(:,j)'*C_*m(:,j) +log(pw(j));
        U = [U,comp_1];
    end
    [ui,iai]=max(U,[],2);
    PExp_(i,1)=PExp_(i,1)+ sum(iai==1);
    PExp_(i,2)=PExp_(i,2)+ sum(iai==2);
    PExp_(i,3)=PExp_(i,3)+ sum(iai==3);
end

disp('Экспериментальная матрица ошибок (гауссовский классификатор)');
disp(Pc_)
disp('Экспериментальная матрица ошибок (k ближайших соеседей)');
disp(Pcv)
disp('Матрица ошибок по методу скользящего контроля');
disp(Pc1)
disp('Теоретическая матрица вероятностей ошибок');
disp(PTheory)
disp('Матрица ошибок на основе границы Чернова');
disp(PExp_)