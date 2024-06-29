% RMA cross-validation
function [L_n] = CV_pl_Annals(x,y,p,K,idx,n0,n_K,M,w)

L_k=zeros(K,1);
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
    predict_all(idx_k0,1)=x_k01*b_k1;
    predict_all(idx_k0,2)=x_k0*b_k2;
% % %     L_k(k)=pairwise_loss_approximate_Annals_1(n_K,y_k0,predict_all(idx_k0,:),w);
    L_k(k)=pairwise_loss_approximate_Annals_2(n_K,y_k0,predict_all(idx_k0,:),w);
% % %     L_k(k)=pairwise_loss_approximate_Annals_3(n_K,y_k0,predict_all(idx_k0,:),w);
% % %     L_k(k)=pairwise_loss_approximate_Annals_4(n_K,y_k0,predict_all(idx_k0,:),w);
end

L_n=sum(L_k)/K;

end