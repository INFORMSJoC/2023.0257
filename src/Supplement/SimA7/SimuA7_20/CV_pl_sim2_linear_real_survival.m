% RMA targeted PL cross-validation
function [L_n] = CV_pl_sim2_linear_real_survival(x,y,K,idx,n0,n_K,M,w)

L_k=zeros(K,1);
predict_all=zeros(n0,M);
for k=1:K
    idx_k=(idx~=k);
    x_k=x(idx_k,:);
    y_k=y(idx_k);
    idx_k0=(idx==k);
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);    
    for m=1:M
        x_km=x_k(:,20*(m-1)+1:20*m);
        x_k0m=x_k0(:,20*(m-1)+1:20*m);
        beta_km=(x_km'*x_km)\(x_km'*y_k); 
        predict_all(idx_k0,m)=x_k0m*beta_km;            
    end 
    L_k(k)=pairwise_loss_real(n_K,y_k0,predict_all(idx_k0,:),w); 
end

L_n=sum(L_k)/K;

end