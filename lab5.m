clear all; close all;



RR = 0.1 : 0.1 : 0.9;
err = RR * 0; 

for t = 1 : numel(RR)
    n=1;
    N =1000;
    r = RR(t);
    h_N=N^(-r/n);
    x0=-3:0.05:3;
    p=zeros(1,length(x0)); 
    
    XN=-log(rand(1,N));
    ind1=logical(x0>0);
    p(ind1)=exp(-x0(ind1));
    
    p_=vkernel(x0,XN,h_N);
    err(t) = mean(abs(p(:) - p_(:)));
end

figure;
plot(RR, err); 



function p_=vkernel(x,XN,h_N)
    [n1,mx]=size(x);
    [n2,N]=size(XN);
    if n1==n2, n=n1;
        p_=zeros(1,mx);
        fit=zeros(N,mx); 
        for i=1:N,
            p_k=zeros(n,mx); mx_i=repmat(XN(:,i),[1,mx]);        
            ro=abs(x-mx_i)/h_N; 
            for k=1:n 
               ind=logical(ro(k,:)<1); p_k(k,ind)=(1-ro(k,ind));
            end
            fit(i,:)=prod(p_k,1)/h_N^n;         
        end
        if N>1
            p_=sum(fit)/N; 
        else
            p_=fit;
        end
    else %n1=~n2
        error('ðàçìåðíîñòè äàííûõ (n1 è n2) íå ñîâïàäàþò');
    end
end
