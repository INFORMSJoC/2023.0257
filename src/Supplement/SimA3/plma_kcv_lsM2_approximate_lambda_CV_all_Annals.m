function [w_kcv,v_kcv,w_plma,vplma,losspl] = plma_kcv_lsM2_approximate_lambda_CV_all_Annals(x,y,n,p,K,M,w0,options)

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
    x_k1=x_k(:,1:(p-1));
    b_k1=(x_k1'*x_k1)\(x_k1'*y_k);  
    b_k2=(x_k'*x_k)\(x_k'*y_k);     
    idx_k0=(idx==k);
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);
    x_k01=x_k0(:,1:(p-1));
    ee(idx_k0,1)=y_k0-x_k01*b_k1;
    ee(idx_k0,2)=y_k0-x_k0*b_k2; 
    predict_all(idx_k0,1)=x_k01*b_k1;
    predict_all(idx_k0,2)=x_k0*b_k2;
end
a1=ee'*ee;
w_kcv=quadprog(a1,zeros(M,1),[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   
v_kcv=(w_kcv'*a1*w_kcv)/n0;

cv_plma_fun=@(w)CV_pl_Annals(x,y,p,K,idx,n0,n_K,M,w);
opt=optimset('Display','off');
Aeq=ones(1,M);
beq=1;
lb=zeros(1,M);
ub=ones(1,M);
w_pl=fmincon(cv_plma_fun,w0',[],[],Aeq,beq,lb,ub,[],opt)';
w_pl=w_pl.*(w_pl>0);
w_pl=w_pl/sum(w_pl);
w_plma=w_pl; 

% % % losspl=CV_pl_Annals(x,y,p,K,idx,n0,n_K,M,w_pl');

losspl=CV_pl_real(x,y,p,K,idx,n0,n_K,M,w_pl');

vplma=(w_plma'*a1*w_plma)/n0;

end

