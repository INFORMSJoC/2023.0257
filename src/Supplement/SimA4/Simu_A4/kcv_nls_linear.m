function [w_kcv,v_kcv] = kcv_nls_linear(x,y,n,K,M,w0,options)

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
    idx_k0=(idx==k);
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);    
    for m=1:M
        x_km=x_k(:,1:m);
        x_k0m=x_k0(:,1:m);
        beta_km=(x_km'*x_km)\(x_km'*y_k); 
        y_km=x_k0m*beta_km;
        ee(idx_k0,m)=y_k0-y_km;  
    end
end
a1=ee'*ee;
w_kcv=quadprog(a1,zeros(M,1),[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   
v_kcv=(w_kcv'*a1*w_kcv)/n0;

end

