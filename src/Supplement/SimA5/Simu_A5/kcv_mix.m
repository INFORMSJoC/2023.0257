function [w_kcv,v_kcv] = kcv_mix(x,y,n,K,M,t0,w0,myfun,options)

% K-fold Cross-Validation
n_K=floor(n/K);
Ki=1:1:K;
idx=kron(Ki',ones(n_K,1));
n0=length(idx);
e_1=zeros(n0,1);
e_2=zeros(n0,1);
for k=1:K
    idx_k=(idx~=k);
    x_k=x(idx_k,:);
    y_k=y(idx_k);
    b_k1=(x_k'*x_k)\(x_k'*y_k);                               %第一个模型拟合的结果
    b_k2=nlinfit(x_k,y_k,myfun,t0);                           %第二个模型拟合的结果
    idx_k0=(idx==k);
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);
    e_1(idx_k0)=y_k0-x_k0*b_k1;                               % 第一个模型的预测误差(K折CV的测试集)
    e_2(idx_k0)=y_k0-myfun(b_k2,x_k0);                        % 第二个模型的预测误差(K折CV的测试集)
end
ee=[e_1,e_2];
a1=ee'*ee;
w_kcv=quadprog(a1,zeros(M,1),[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   
v_kcv=(w_kcv'*a1*w_kcv)/n0;

end

