function [w_kcv,v_kcv,wplma,vplma,losspl] = plma_kcv_mix_approximate_lambda_CV_all_Annals_normal(x,y,n,K,M,w0,options)

% 将两个程序（K折交叉验证和PLMA）整合

% 这个函数遍历候选的lambda选择近似函数sigmod的最优lamda值中最接近真实的loss function的lambda对应的权重

% pairwise loss function of K-fold Cross-Validation 

n_K=floor(n/K);                                               %分为组后每组的样本个数
Ki=1:1:K;                                                     
idx=kron(Ki',ones(n_K,1));                                    %采用张量分组
n0=length(idx);
e_1=zeros(n0,1);
e_2=zeros(n0,1);
p_1=zeros(n0,1);
p_2=zeros(n0,1);
for k=1:K                                                     %进行折交叉验证(K组)
    idx_k=(idx~=k);                                           %训练集指标
    x_k=x(idx_k,:);
    y_k=y(idx_k);
    b_k1=(x_k(:,1:5)'*x_k(:,1:5))\(x_k(:,1:5)'*y_k);          %第一个模型拟合的结果
    b_k2=(x_k'*x_k)\(x_k'*y_k);                               %第二个模型拟合的结果
    idx_k0=(idx==k);                                          %测试集指标
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);
    e_1(idx_k0)=y_k0-x_k0(:,1:5)*b_k1;                               % 第一个模型的预测误差(K折CV的测试集)
    e_2(idx_k0)=y_k0-x_k0*b_k2;                                      % 第二个模型的预测误差(K折CV的测试集)
    p_1(idx_k0)=x_k0(:,1:5)*b_k1;                                    % 第一个模型的预测值
    p_2(idx_k0)=x_k0*b_k2;                                           % 第二个模型的预测值
end
ee=[e_1,e_2];
a1=ee'*ee;
w_kcv=quadprog(a1,zeros(M,1),[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   
v_kcv=(w_kcv'*a1*w_kcv)/n0;

opt=optimset('Display','off');

plfun=@(w)CV_pl_mix_Annals_normal(x,y,K,idx,n0,n_K,M,w);            % 连续情形Annals近似方案1和2
Aeq=ones(1,M);
beq=1;
lb=zeros(1,M);
ub=ones(1,M);
w_pl=fmincon(plfun,w0',[],[],Aeq,beq,lb,ub,[],opt)';
w_pl=w_pl.*(w_pl>0);
w_pl=w_pl/sum(w_pl);
wplma=w_pl;

% % % losspl=CV_pl_sim1_Annals(x,y,K,idx,n0,n_K,M,w_pl');           % 连续情形(近似CV) 

losspl=CV_pl_mix_real_normal(x,y,K,idx,n0,n_K,M,w_pl');             % 连续情形(真实CV)  

vplma=(wplma'*a1*wplma)/n0;

end

