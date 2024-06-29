% RMA targeted PL cross-validation
function [L_n] = CV_pl_mix_real_normal(x,y,K,idx,n0,n_K,M,w)

L_k=zeros(K,1);
predict_all=zeros(n0,M);
for k=1:K                                                     % 进行折交叉验证(K组)
    idx_k=(idx~=k);                                           % 训练集指标
    x_k=x(idx_k,:);
    y_k=y(idx_k);
    b_k1=(x_k(:,1:5)'*x_k(:,1:5))\(x_k(:,1:5)'*y_k);          %第一个模型拟合的结果
    b_k2=(x_k'*x_k)\(x_k'*y_k);                               %第二个模型拟合的结果
    idx_k0=(idx==k);                                          % 测试集指标
    x_k0=x(idx_k0,:);
    y_k0=y(idx_k0);
    predict_all(idx_k0,1)=x_k0(:,1:5)*b_k1;            % 第一个模型的预测值
    predict_all(idx_k0,2)=x_k0*b_k2;                   % 第二个模型的预测值  
    L_k(k)=pairwise_loss_real(n_K,y_k0,predict_all(idx_k0,:),w);   
end

L_n=sum(L_k)/K;

end