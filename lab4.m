clear all; close all;

n=45;
M=3;
s=zeros(n,M); 

letterA =  [0 0 0 0 1 0 0 0 0 ...
            0 0 0 1 0 1 0 0 0 ...
            0 0 1 0 0 0 1 0 0 ...
            0 1 0 1 0 0 0 1 0 ...
            1 0 0 0 0 0 0 0 1 ]';
letterO =  [0 0 0 1 1 1 0 0 0 ...
            0 0 1 0 0 0 1 0 0 ...
            0 0 1 0 0 0 1 0 0 ...
            0 0 1 0 0 0 1 0 0 ...
            0 0 0 1 1 1 0 0 0 ]';
letterE =  [1 1 1 1 1 1 1 1 1 ...
            1 0 0 0 0 0 0 0 0 ...
            1 1 1 1 1 1 1 1 1 ...
            1 0 0 0 0 0 0 0 0 ...
            1 1 1 1 1 1 1 1 1 ]';        
s(:,1)=letterA; 
s(:,2)=letterO;
s(:,3)=letterE;

pw = [0.3,0.5,0.2];

K = 1000;
Pc_ = zeros(M);
Pt = zeros(M);

pI=0.6;  
pI_ = 1-pI;
s_=1-s;

G1=zeros(1,n);  
G2=zeros(1,n); 

for i = 1:M-1
    for j = i+1:M
        ns=sum(abs(s(:,i)-s(:,j)));
        l0_ = log(pw(j)/pw(i));
        L0  = log(pw(j)/pw(i)) / (2*log(pI_)-2*log(pI)) + ns/2;
        L0r = floor(L0);                
        Pt(i,j) = 1 - binocdf(L0r,ns,1-pI);
        Pt(j,i) = binocdf(L0r,ns,pI);
    end
end

Pt = Pt + diag(ones(3, 1) - sum(Pt, 2));

for k=1:K
    for i=1:M
        x = s(:,i);
        r = rand(n,1); 
        ir = find(r < pI);
        x(ir) = 1-x(ir);
        x_ = 1-x;
        
        iais= [];
        for ii = 1:M-1
            for jj  = ii+1:M
                l0_ = log(pw(jj)/pw(ii));
                for kk=1:n
                   G1(1,kk)=log((s(kk,ii)*pI_+s_(kk,ii)*pI)/(s(kk,jj)*pI_+s_(kk,jj)*pI));
                   G2(1,kk)=log((s(kk,ii)*pI+s_(kk,ii)*pI_)/(s(kk,jj)*pI+s_(kk,jj)*pI_));
                end
                u=G1*x+G2*x_-l0_; 
                
                if u>0
                    iai=ii;
                else
                    iai=jj;
                end
                
                iais = [iais, iai];
            end
        end
        id = mode(iais);
        Pc_(i,id)=Pc_(i,id) + 1;
    end
end
Pc_ = Pc_ / K;
disp('Theoretical Result')
disp(Pt)
disp('Experimental Result')
disp(Pc_)



