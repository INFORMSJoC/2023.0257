function [w_kcv,v_kcv] = kcv_lsM2(x,y,n,p,K,M,w0,options)

% K-fold Cross-Validation 
n_K=floor(n/K);
Ki=1:1:K;
idx=kron(Ki',ones(n_K,1));
n0=length(idx);
ee=zeros(n0,M);
for k=1:K
    idx_k=(idx~=k);
    x_k=x(idx_k,:);
    y_k=y(idx_k);
    x_k1=x_k(:,1:(p-1));
    b_k1=(x_k1'*x_k1)\(x_k1'*y_k);  
    b_k2=(x_k'*x_k)\(x_k'*y_k);     
    idx_k0=(idx==k);
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);
    x_k01=x_k0(:,1:(p-1));
    ee(idx_k0,1)=y_k0-x_k01*b_k1; 
    ee(idx_k0,2)=y_k0-x_k0*b_k2;
end
a1=ee'*ee;
w_kcv=quadprog(a1,zeros(M,1),[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   
v_kcv=(w_kcv'*a1*w_kcv)/n0;

end

