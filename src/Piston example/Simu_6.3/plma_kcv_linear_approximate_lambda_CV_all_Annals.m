function [w_kcv,v_kcv,wplma,vplma,losspl] = plma_kcv_linear_approximate_lambda_CV_all_Annals(x,y,n,K,M,w0,options)

% 将两个程序（K折交叉验证和PLMA）整合

% 这个函数遍历候选的lambda选择近似函数sigmod的最优lamda值中最接近真实的loss function的lambda对应的权重

% K-fold Cross-Validation 

n_K=floor(n/K);
Ki=1:1:K;
idx=kron(Ki',ones(n_K,1));
n0=length(idx);

ee=zeros(n0,M);
predict_all=zeros(n0,M);
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
        predict_all(idx_k0,m)=x_k0m*beta_km;      
    end    
end
a1=ee'*ee;
w_kcv=quadprog(a1,zeros(M,1),[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   
v_kcv=(w_kcv'*a1*w_kcv)/n0;

opt=optimset('Display','off');

plfun=@(w)CV_pl_sim2_linear_Annals(x,y,K,idx,n0,n_K,M,w);
Aeq=ones(1,M);
beq=1;
lb=zeros(1,M);
ub=ones(1,M);
w_pl=fmincon(plfun,w0',[],[],Aeq,beq,lb,ub,[],opt)';
w_pl=w_pl.*(w_pl>0);
w_pl=w_pl/sum(w_pl);
wplma=w_pl;

% losspl=CV_pl_sim2_linear_Annals(x,y,K,idx,n0,n_K,M,w_pl');

losspl=CV_pl_sim2_linear_real(x,y,K,idx,n0,n_K,M,w_pl');

vplma=(wplma'*a1*wplma)/n0;

end

