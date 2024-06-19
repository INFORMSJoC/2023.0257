%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation A.6 
%  piston data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; tic;
rseed=RandStream('mt19937ar','Seed',202301);  %set seed for reproducing results
RandStream.setGlobalStream(rseed);
options=optimset('algorithm','interior-point-convex','Display','off','MaxIter',500);
warning('off')

% Load data and define variables
load('piston100.m')
data=piston100;
Y=data(:,8);             
X=zscore(data(:,2:7));
N=length(Y);
p=size(X,2);

% slection matrix
s0=zeros(1,p);
s0(1)=1;
s1=zeros(1,p);
Pi=zeros(2^p,p);
for i=2:2^p
    s1=s0+s1;
    for j=1:p
        if s1(j)==2
            s1(j+1)=s1(j+1)+1;
            s1(j)=0;
        end
    end
    Pi(i,:)=s1;
end
Pi=[ones(size(Pi,1),1) Pi];
Pi=(Pi==1);
M=length(Pi);        
w0=ones(p+1,1)/(p+1);       

% parameters
S=50;                             
n2s=40;                         % evaluation set  60 40 20 10
n1=60;                          % training set  40 60 80 90  
if n1==40; Ks=floor([2,5,10,n1/3,n1/2,n1]);  end
if n1==60; Ks=floor([2,5,10,n1/3,n1/2,n1]);  end
if n1==80; Ks=floor([2,5,10,n1/3,n1/2,n1]);  end
if n1==90; Ks=floor([2,5,10,n1/3,n1/2,n1]);  end
kn=length(Ks);
cons=ones(n1,1);

cn=size(n2s,2);
% r1_cn=zeros(cn,2*kn+9);
r2_cn=zeros(cn,2*kn+9);
% r1_cn_in=zeros(cn,2*kn+9);
r2_cn_in=zeros(cn,2*kn+9);
% w_cn=zeros(cn,kn+1);
% w_cn_pl=zeros(cn,kn);
% w_cn_pl_loss=zeros(cn,kn);
k_cn=zeros(cn,kn);
pl_cn=zeros(cn,kn-1);
pl_loss_cn=zeros(cn,kn-1);
rank_score_cn=zeros(cn,2*kn+9);
rank_score_cn_in=zeros(cn,2*kn+9);
ppl_cn=zeros(cn,2*kn+9);
ppl_cn_in=zeros(cn,2*kn+9);
ppl_Annals_1_cn=zeros(cn,2*kn+9);
ppl_Annals_1_cn_in=zeros(cn,2*kn+9);
% % % % % ppl_Annals_2_cn=zeros(cn,2*kn+9);
% % % % % ppl_Annals_2_cn_in=zeros(cn,2*kn+9);
% % % % % ppl_Annals_3_cn=zeros(cn,2*kn+9);
% % % % % ppl_Annals_3_cn_in=zeros(cn,2*kn+9);
% % % % % ppl_Annals_4_cn=zeros(cn,2*kn+9);
% % % % % ppl_Annals_4_cn_in=zeros(cn,2*kn+9);  

for c=1:cn 
n2=n2s(c);
% r1s=zeros(S,2*kn+9);
r2s=zeros(S,2*kn+9);
% r3s=zeros(S,p+1);
% r1s_in=zeros(S,2*kn+9);
r2s_in=zeros(S,2*kn+9);
% r3s_in=zeros(S,M);
rws=zeros(p+1,2*kn+8); 
% ws_cv=zeros(S,kn+1);
% wpl_cv=zeros(S,kn);
% wpl_loss_cv=zeros(S,kn);
ks_min=zeros(S,kn);
pl_min=zeros(S,kn-1);
pl_loss_min=zeros(S,kn-1);
ppl=zeros(S,2*kn+9);
ppl_in=zeros(S,2*kn+9);
% % % % % ppl_Annals_1=zeros(S,2*kn+9);
% % % % % ppl_Annals_1_in=zeros(S,2*kn+9);
ppl_Annals_2=zeros(S,2*kn+9);
ppl_Annals_2_in=zeros(S,2*kn+9);
% % % % % ppl_Annals_3=zeros(S,2*kn+9);
% % % % % ppl_Annals_3_in=zeros(S,2*kn+9);
% % % % % ppl_Annals_4=zeros(S,2*kn+9);
% % % % % ppl_Annals_4_in=zeros(S,2*kn+9);


parfor_progress(S);
parfor s=1:S
if floor(s/50)*50==s; disp([n1,c,s/50]); end

idx0=randperm(N,n1+n2)';
idx1=idx0(1:n1);
idx2=idx0(n1+1:n1+n2);
x=[ones(n1,1),X(idx1,:)];
y=Y(idx1);
x0=[ones(n2,1),X(idx2,:)];
y0=Y(idx2);


M=length(Pi);                          
aic=zeros(M,1);
bic=zeros(M,1);
yin=zeros(M,size(y,1));                
yy=zeros(M,n2);                      
eem=zeros(n1,M);
sigg=zeros(M,1);
pim=zeros(M,1);

for m=1:M
    xm=x(:,Pi(m,:)); 

    betam=(xm'*xm)\(xm'*y); 
    ym=xm*betam;    
    em=y-ym;
    sig2=(em'*em)/n1;
    km=sum(Pi(m,:));
    aic(m)=n1*log(sig2)+2*km;
    bic(m)=n1*log(sig2)+log(n1)*km;
    yy(m,:)=x0(:,Pi(m,:))*betam;
    yin(m,:)=x(:,Pi(m,:))*betam;
    
    % MMA  
    eem(:,m)=em;
    sigg(m)=sig2;
    pim(m)=km;
end

Pi_m=zeros(p+1,p+1);
aic_m=zeros(p+1,1);
bic_m=zeros(p+1,1);
yin_m=zeros(p+1,size(y,1));                      % 样本内训练集
yy_m=zeros(p+1,n2);                            % test number  
eem_m=zeros(n1,p+1);
sigg_m=zeros(p+1,1);
pim_m=zeros(p+1,1);

for mm=1:p+1
   m_loc=find(sum(Pi,2)==mm);
   [~,b_loc]=min(bic(m_loc));
   Pi_m(mm,:)=Pi(m_loc(b_loc),:);
   aic_m(mm,:)=aic(m_loc(b_loc));
   bic_m(mm,:)=bic(m_loc(b_loc));
   yy_m(mm,:)=yy(m_loc(b_loc),:);
   yin_m(mm,:)=yin(m_loc(b_loc),:);
   eem_m(:,mm)=eem(:,m_loc(b_loc));
   sigg_m(mm)=sigg(m_loc(b_loc));
   pim_m(mm)=pim(m_loc(b_loc));
end

M=size(Pi_m,1);

Pi_m=(Pi_m==1);

% AIC & BIC & MMA 
aic_m=aic_m-min(aic_m);    
bic_m=bic_m-min(bic_m);  
w_aic=0+logical(aic_m==min(aic_m));
w_bic=0+logical(bic_m==min(bic_m));  
w_saic=exp(-0.5*aic_m)./sum(exp(-0.5*aic_m));
w_sbic=exp(-0.5*bic_m)./sum(exp(-0.5*bic_m));
if sum(exp(-0.5*aic_m))==0; w_saic=w0; end
if sum(exp(-0.5*bic_m))==0; w_sbic=w0; end
sighat=min(sigg_m);
a1=eem_m'*eem_m;  
a2=pim_m*sighat;
w_nma=quadprog(a1,a2,[],[],ones(1,M),1,zeros(M,1),ones(M,1),w0,options);   

% Data-driven K 
w_ks=zeros(M,kn);
w_pl=zeros(M,kn-1);
v_ks=zeros(1,kn);
v_pl=zeros(1,kn-1);
loss_pl=zeros(1,kn-1);
y_ks=zeros(size(yy_m,2),kn);          
y_pl=zeros(size(yy_m,2),kn-1);          
y_ks_in=zeros(size(yin_m,2),kn);       
y_pl_in=zeros(size(yin_m,2),kn-1);     
ppl_ks=zeros(1,kn);
ppl_pl=zeros(1,kn-1);
ppl_ks_in=zeros(1,kn);
ppl_pl_in=zeros(1,kn-1);
% % % % % ppl_Annals_1_ks=zeros(1,kn);
% % % % % ppl_Annals_1_pl=zeros(1,kn-1);
% % % % % ppl_Annals_1_ks_in=zeros(1,kn);
% % % % % ppl_Annals_1_pl_in=zeros(1,kn-1);
ppl_Annals_2_ks=zeros(1,kn);
ppl_Annals_2_pl=zeros(1,kn-1);
ppl_Annals_2_ks_in=zeros(1,kn);
ppl_Annals_2_pl_in=zeros(1,kn-1);
% % % % % ppl_Annals_3_ks=zeros(1,kn);
% % % % % ppl_Annals_3_pl=zeros(1,kn-1);
% % % % % ppl_Annals_3_ks_in=zeros(1,kn);
% % % % % ppl_Annals_3_pl_in=zeros(1,kn-1);
% % % % % ppl_Annals_4_ks=zeros(1,kn);
% % % % % ppl_Annals_4_pl=zeros(1,kn-1);
% % % % % ppl_Annals_4_ks_in=zeros(1,kn);
% % % % % ppl_Annals_4_pl_in=zeros(1,kn-1);

for k=1:kn
    if k~=kn
        K=Ks(k);
       [w_k,v_k,wplma,vplma,losspl]=plma_kcv_linearM_approximate_lambda_CV_all_Annals(x,y,n1,K,M,Pi_m,w0,options);
        w_ks(:,k)=w_k;  
        v_ks(k)=v_k;  
        y_ks(:,k)=yy_m'*w_k;
        y_ks_in(:,k)=yin_m'*w_k;
        w_pl(:,k)=wplma;
        v_pl(k)=vplma;    
        y_pl(:,k)=yy_m'*wplma;
        y_pl_in(:,k)=yin_m'*wplma;
        loss_pl(k)=losspl;
        ppl_ks(k)=pairwise_loss(length(y0),y0,yy_m',w_k');      
        ppl_pl(k)=pairwise_loss(length(y0),y0,yy_m',wplma');
        ppl_ks_in(k)=pairwise_loss(length(y),y,yin_m',w_k');      
        ppl_pl_in(k)=pairwise_loss(length(y),y,yin_m',wplma');        

% % % % %         ppl_Annals_1_ks(k)=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_k'); 
% % % % %         ppl_Annals_1_pl(k)=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',wplma'); 
% % % % %         ppl_Annals_1_ks_in(k)=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_k'); 
% % % % %         ppl_Annals_1_pl_in(k)=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',wplma'); 

        ppl_Annals_2_ks(k)=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_k');   
        ppl_Annals_2_pl(k)=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',wplma');
        ppl_Annals_2_ks_in(k)=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_k');   
        ppl_Annals_2_pl_in(k)=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',wplma');

% % % % %         ppl_Annals_3_ks(k)=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_k'); 
% % % % %         ppl_Annals_3_pl(k)=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',wplma');
% % % % %         ppl_Annals_3_ks_in(k)=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_k'); 
% % % % %         ppl_Annals_3_pl_in(k)=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',wplma');   

% % % % %         ppl_Annals_4_ks(k)=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_k'); 
% % % % %         ppl_Annals_4_pl(k)=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',wplma'); 
% % % % %         ppl_Annals_4_ks_in(k)=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_k'); 
% % % % %         ppl_Annals_4_pl_in(k)=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',wplma'); 

    elseif k==kn
        K=Ks(k);
        [w_k,v_k]=kcv_nlsM_linear(x,y,n1,K,M,Pi_m,w0,options);
        w_ks(:,k)=w_k;  
        v_ks(k)=v_k;
        y_ks(:,k)=yy_m'*w_k;
        y_ks_in(:,k)=yin_m'*w_k;
        ppl_ks(k)=pairwise_loss(length(y0),y0,yy_m',w_k');
        ppl_ks_in(k)=pairwise_loss(length(y),y,yin_m',w_k');

% % % % %         ppl_Annals_1_ks(k)=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_k');    
% % % % %         ppl_Annals_1_ks_in(k)=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_k'); 

        ppl_Annals_2_ks(k)=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_k'); 
        ppl_Annals_2_ks_in(k)=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_k');   

% % % % %         ppl_Annals_3_ks(k)=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_k'); 
% % % % %         ppl_Annals_3_ks_in(k)=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_k'); 

% % % % %         ppl_Annals_4_ks(k)=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_k');
% % % % %         ppl_Annals_4_ks_in(k)=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_k'); 

    end
end

% % % kmin=logical(v_ks==min(v_ks));
[~,kmin]=min(v_ks);
w_kmin=w_ks(:,kmin);
y_kmin=yy_m'*w_kmin;
y_kmin_in=yin_m'*w_kmin;
ks_min(s,:)=0+kmin;

% % % plmin=logical(v_pl==min(v_pl));
[~,plmin]=min(v_pl);
w_plmin=w_pl(:,plmin);
y_plmin=yy_m'*w_plmin;
y_plmin_in=yin_m'*w_plmin;
pl_min(s,:)=0+plmin;

% 求解最小的loss下的权重
[~,plloss_min]=min(loss_pl);
w_pl_loss_min=w_pl(:,plloss_min);
y_pl_loss_min=yy_m'*w_pl_loss_min;
y_pl_loss_min_in=yin_m'*w_pl_loss_min;
pl_loss_min(s,:)=0+plloss_min;

% Prediciton Risk
y_max=yy(size(yy,1),:)';    
y1=yy_m'*w_aic;
y2=yy_m'*w_bic;
y3=yy_m'*w_saic;
y4=yy_m'*w_sbic;
y5=yy_m'*w_nma;
y6=yy_m'*w0;
ys=[y_max,y1,y2,y3,y4,y5,y6,y_kmin,y_plmin,y_pl_loss_min,y_ks,y_pl];
% r1s(s,:)=mean((ys-yT).^2);
r2s(s,:)=mean((ys-y0).^2);
% r3s(s,:)=mean((yy'-yT).^2);
[~,I0]=sort(y0,'descend');
[~,I]=sort(ys,'descend');
rank_score(s,:)=sum(I0==I,1)/length(y0);
rws=rws+[w_aic,w_bic,w_saic,w_sbic,w_nma,w0,w_kmin,w_plmin,w_pl_loss_min,w_ks,w_pl];

% Risk in sample
y_max_in=yin(size(yin,1),:)';    
y1_in=yin_m'*w_aic;
y2_in=yin_m'*w_bic;
y3_in=yin_m'*w_saic;
y4_in=yin_m'*w_sbic;
y5_in=yin_m'*w_nma;
y6_in=yin_m'*w0;
ys_in=[y_max_in,y1_in,y2_in,y3_in,y4_in,y5_in,y6_in,y_kmin_in,y_plmin_in,y_pl_loss_min_in,y_ks_in,y_pl_in];
% r1s_in(s,:)=mean((ys_in-yT_in).^2);
r2s_in(s,:)=mean((ys_in-y).^2);
% r3s_in(s,:)=mean((yin'-yT_in).^2);
[~,I0_in]=sort(y,'descend');
[~,I_in]=sort(ys_in,'descend');
rank_score_in(s,:)=sum(I0_in==I_in,1)/length(y);
% Prediciton Pairwise loss
ppl_y_max=pairwise_loss(length(y0),y0,y_max,1);
ppl_y1=pairwise_loss(length(y0),y0,yy_m',w_aic');
ppl_y2=pairwise_loss(length(y0),y0,yy_m',w_bic');
ppl_y3=pairwise_loss(length(y0),y0,yy_m',w_saic');
ppl_y4=pairwise_loss(length(y0),y0,yy_m',w_sbic');
ppl_y5=pairwise_loss(length(y0),y0,yy_m',w_nma');
ppl_y6=pairwise_loss(length(y0),y0,yy_m',w0');
ppl_ykmin=pairwise_loss(length(y0),y0,yy_m',w_kmin');
ppl_yplmin=pairwise_loss(length(y0),y0,yy_m',w_plmin');
ppl_loss_yplmin=pairwise_loss(length(y0),y0,yy_m',w_pl_loss_min');
ppl(s,:)=[ppl_y_max,ppl_y1,ppl_y2,ppl_y3,ppl_y4,ppl_y5,ppl_y6,ppl_ykmin,ppl_yplmin,ppl_loss_yplmin,ppl_ks,ppl_pl];

% Pairwise loss in sample
ppl_y_max_in=pairwise_loss(length(y),y,y_max_in,1);
ppl_y1_in=pairwise_loss(length(y),y,yin_m',w_aic');
ppl_y2_in=pairwise_loss(length(y),y,yin_m',w_bic');
ppl_y3_in=pairwise_loss(length(y),y,yin_m',w_saic');
ppl_y4_in=pairwise_loss(length(y),y,yin_m',w_sbic');
ppl_y5_in=pairwise_loss(length(y),y,yin_m',w_nma');
ppl_y6_in=pairwise_loss(length(y),y,yin_m',w0');
ppl_ykmin_in=pairwise_loss(length(y),y,yin_m',w_kmin');
ppl_yplmin_in=pairwise_loss(length(y),y,yin_m',w_plmin');
ppl_loss_yplmin_in=pairwise_loss(length(y),y,yin_m',w_pl_loss_min');
ppl_in(s,:)=[ppl_y_max_in,ppl_y1_in,ppl_y2_in,ppl_y3_in,ppl_y4_in,ppl_y5_in,ppl_y6_in,ppl_ykmin_in,ppl_yplmin_in,ppl_loss_yplmin_in,ppl_ks_in,ppl_pl_in];

ppl_std(s,:)=(ppl(s,:)-ppl_y_max_in)/ppl_y_max_in;

% % % % % % Prediciton Pairwise loss Annals 1
% % % % % ppl_Annals_1_y_max=pairwise_loss_approximate_Annals_1(length(y0),y0,y_max,1);
% % % % % ppl_Annals_1_y1=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_aic');
% % % % % ppl_Annals_1_y2=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_bic');
% % % % % ppl_Annals_1_y3=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_saic');
% % % % % ppl_Annals_1_y4=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_sbic');
% % % % % ppl_Annals_1_y5=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_nma');
% % % % % ppl_Annals_1_y6=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w0');
% % % % % ppl_Annals_1_ykmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_kmin');
% % % % % ppl_Annals_1_yplmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_plmin');
% % % % % ppl_loss_Annals_1_yplmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy_m',w_pl_loss_min');
% % % % % ppl_Annals_1(s,:)=[ppl_Annals_1_y_max,ppl_Annals_1_y1,ppl_Annals_1_y2,ppl_Annals_1_y3,ppl_Annals_1_y4,ppl_Annals_1_y5,ppl_Annals_1_y6,ppl_Annals_1_ykmin,ppl_Annals_1_yplmin,ppl_loss_Annals_1_yplmin,ppl_Annals_1_ks,ppl_Annals_1_pl];
% % % % % 
% % % % % % Pairwise loss in sample Annals 1
% % % % % ppl_Annals_1_y_max_in=pairwise_loss_approximate_Annals_1(length(y),y,y_max_in,1);
% % % % % ppl_Annals_1_y1_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_aic');
% % % % % ppl_Annals_1_y2_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_bic');
% % % % % ppl_Annals_1_y3_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_saic');
% % % % % ppl_Annals_1_y4_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_sbic');
% % % % % ppl_Annals_1_y5_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_nma');
% % % % % ppl_Annals_1_y6_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w0');
% % % % % ppl_Annals_1_ykmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_kmin');
% % % % % ppl_Annals_1_yplmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_plmin');
% % % % % ppl_loss_Annals_1_yplmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin_m',w_pl_loss_min');
% % % % % ppl_Annals_1_in(s,:)=[ppl_Annals_1_y_max_in,ppl_Annals_1_y1_in,ppl_Annals_1_y2_in,ppl_Annals_1_y3_in,ppl_Annals_1_y4_in,ppl_Annals_1_y5_in,ppl_Annals_1_y6_in,ppl_Annals_1_ykmin_in,ppl_Annals_1_yplmin_in,ppl_loss_Annals_1_yplmin_in,ppl_Annals_1_ks_in,ppl_Annals_1_pl_in];

% Prediciton Pairwise loss Annals 2
ppl_Annals_2_y_max=pairwise_loss_approximate_Annals_2(length(y0),y0,y_max,1);
ppl_Annals_2_y1=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_aic');
ppl_Annals_2_y2=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_bic');
ppl_Annals_2_y3=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_saic');
ppl_Annals_2_y4=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_sbic');
ppl_Annals_2_y5=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_nma');
ppl_Annals_2_y6=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w0');
ppl_Annals_2_ykmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_kmin');
ppl_Annals_2_yplmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_plmin');
ppl_loss_Annals_2_yplmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy_m',w_pl_loss_min');
ppl_Annals_2(s,:)=[ppl_Annals_2_y_max,ppl_Annals_2_y1,ppl_Annals_2_y2,ppl_Annals_2_y3,ppl_Annals_2_y4,ppl_Annals_2_y5,ppl_Annals_2_y6,ppl_Annals_2_ykmin,ppl_Annals_2_yplmin,ppl_loss_Annals_2_yplmin,ppl_Annals_2_ks,ppl_Annals_2_pl];

% Pairwise loss in sample Annals 2
ppl_Annals_2_y_max_in=pairwise_loss_approximate_Annals_2(length(y),y,y_max_in,1);
ppl_Annals_2_y1_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_aic');
ppl_Annals_2_y2_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_bic');
ppl_Annals_2_y3_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_saic');
ppl_Annals_2_y4_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_sbic');
ppl_Annals_2_y5_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_nma');
ppl_Annals_2_y6_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w0');
ppl_Annals_2_ykmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_kmin');
ppl_Annals_2_yplmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_plmin');
ppl_loss_Annals_2_yplmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin_m',w_pl_loss_min');
ppl_Annals_2_in(s,:)=[ppl_Annals_2_y_max_in,ppl_Annals_2_y1_in,ppl_Annals_2_y2_in,ppl_Annals_2_y3_in,ppl_Annals_2_y4_in,ppl_Annals_2_y5_in,ppl_Annals_2_y6_in,ppl_Annals_2_ykmin_in,ppl_Annals_2_yplmin_in,ppl_loss_Annals_2_yplmin_in,ppl_Annals_2_ks_in,ppl_Annals_2_pl_in];

ppl_Annals_2_std(s,:)=(ppl_Annals_2(s,:)-ppl_Annals_2_y_max_in)/ppl_Annals_2_y_max_in;

% % % % % % Prediciton Pairwise loss Annals 3
% % % % % ppl_Annals_3_y_max=pairwise_loss_approximate_Annals_3(length(y0),y0,y_max,1);
% % % % % ppl_Annals_3_y1=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_aic');
% % % % % ppl_Annals_3_y2=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_bic');
% % % % % ppl_Annals_3_y3=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_saic');
% % % % % ppl_Annals_3_y4=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_sbic');
% % % % % ppl_Annals_3_y5=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_nma');
% % % % % ppl_Annals_3_y6=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w0');
% % % % % ppl_Annals_3_ykmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_kmin');
% % % % % ppl_Annals_3_yplmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_plmin');
% % % % % ppl_loss_Annals_3_yplmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy_m',w_pl_loss_min');
% % % % % ppl_Annals_3(s,:)=[ppl_Annals_3_y_max,ppl_Annals_3_y1,ppl_Annals_3_y2,ppl_Annals_3_y3,ppl_Annals_3_y4,ppl_Annals_3_y5,ppl_Annals_3_y6,ppl_Annals_3_ykmin,ppl_Annals_3_yplmin,ppl_loss_Annals_3_yplmin,ppl_Annals_3_ks,ppl_Annals_3_pl];
% % % % % 
% % % % % % Pairwise loss in sample Annals 3
% % % % % ppl_Annals_3_y_max_in=pairwise_loss_approximate_Annals_3(length(y),y,y_max_in,1);
% % % % % ppl_Annals_3_y1_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_aic');
% % % % % ppl_Annals_3_y2_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_bic');
% % % % % ppl_Annals_3_y3_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_saic');
% % % % % ppl_Annals_3_y4_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_sbic');
% % % % % ppl_Annals_3_y5_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_nma');
% % % % % ppl_Annals_3_y6_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w0');
% % % % % ppl_Annals_3_ykmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_kmin');
% % % % % ppl_Annals_3_yplmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_plmin');
% % % % % ppl_loss_Annals_3_yplmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin_m',w_pl_loss_min');
% % % % % ppl_Annals_3_in(s,:)=[ppl_Annals_3_y_max_in,ppl_Annals_3_y1_in,ppl_Annals_3_y2_in,ppl_Annals_3_y3_in,ppl_Annals_3_y4_in,ppl_Annals_3_y5_in,ppl_Annals_3_y6_in,ppl_Annals_3_ykmin_in,ppl_Annals_3_yplmin_in,ppl_loss_Annals_3_yplmin_in,ppl_Annals_3_ks_in,ppl_Annals_3_pl_in];

% % % % % % Prediciton Pairwise loss Annals 4
% % % % % ppl_Annals_4_y_max=pairwise_loss_approximate_Annals_4(length(y0),y0,y_max,1);
% % % % % ppl_Annals_4_y1=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_aic');
% % % % % ppl_Annals_4_y2=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_bic');
% % % % % ppl_Annals_4_y3=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_saic');
% % % % % ppl_Annals_4_y4=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_sbic');
% % % % % ppl_Annals_4_y5=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_nma');
% % % % % ppl_Annals_4_y6=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w0');
% % % % % ppl_Annals_4_ykmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_kmin');
% % % % % ppl_Annals_4_yplmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_plmin');
% % % % % ppl_loss_Annals_4_yplmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy_m',w_pl_loss_min');
% % % % % ppl_Annals_4(s,:)=[ppl_Annals_4_y_max,ppl_Annals_4_y1,ppl_Annals_4_y2,ppl_Annals_4_y3,ppl_Annals_4_y4,ppl_Annals_4_y5,ppl_Annals_4_y6,ppl_Annals_4_ykmin,ppl_Annals_4_yplmin,ppl_loss_Annals_4_yplmin,ppl_Annals_4_ks,ppl_Annals_4_pl];
% % % % % 
% % % % % % Pairwise loss in sample Annals 4
% % % % % ppl_Annals_4_y_max_in=pairwise_loss_approximate_Annals_4(length(y),y,y_max_in,1);
% % % % % ppl_Annals_4_y1_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_aic');
% % % % % ppl_Annals_4_y2_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_bic');
% % % % % ppl_Annals_4_y3_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_saic');
% % % % % ppl_Annals_4_y4_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_sbic');
% % % % % ppl_Annals_4_y5_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_nma');
% % % % % ppl_Annals_4_y6_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w0');
% % % % % ppl_Annals_4_ykmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_kmin');
% % % % % ppl_Annals_4_yplmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_plmin');
% % % % % ppl_loss_Annals_4_yplmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin_m',w_pl_loss_min');
% % % % % ppl_Annals_4_in(s,:)=[ppl_Annals_4_y_max_in,ppl_Annals_4_y1_in,ppl_Annals_4_y2_in,ppl_Annals_4_y3_in,ppl_Annals_4_y4_in,ppl_Annals_4_y5_in,ppl_Annals_4_y6_in,ppl_Annals_4_ykmin_in,ppl_Annals_4_yplmin_in,ppl_loss_Annals_4_yplmin_in,ppl_Annals_4_ks_in,ppl_Annals_4_pl_in];


parfor_progress;


end

parfor_progress(0);

% r1_cn(c,:)=mean(r1s.^2);
r2_cn(c,:)=mean(r2s.^2);
% r1_cn_in(c,:)=mean(r1s_in.^2);
r2_cn_in(c,:)=mean(r2s_in.^2);
% w_cn(c,:)=mean(ws_cv);
% w_cn_pl(c,:)=mean(wpl_cv);
% w_cn_pl_loss(c,:)=mean(wpl_loss_cv);
k_cn(c,:)=mean(ks_min);
pl_cn(c,:)=mean(pl_min);
pl_loss_cn(c,:)=mean(pl_loss_min);
rank_score_cn(c,:)=mean(rank_score);
rank_score_cn_in(c,:)=mean(rank_score_in);
ppl_cn(c,:)=mean(ppl);
ppl_cn_std(c,:)=mean(ppl_std);
ppl_cn_in(c,:)=mean(ppl_in);
ppl_cn_sd(c,:)=std(ppl);
ppl_cn_in_sd(c,:)=std(ppl_in);
ppl_cn_med(c,:)=median(ppl);
ppl_cn_in_med(c,:)=median(ppl_in);
% % % ppl_Annals_1_cn(c,:)=mean(ppl_Annals_1);
% % % ppl_Annals_1_cn_in(c,:)=mean(ppl_Annals_1_in);
ppl_Annals_2_cn(c,:)=mean(ppl_Annals_2);
ppl_Annals_2_cn_std(c,:)=mean(ppl_Annals_2_std);
ppl_Annals_2_cn_in(c,:)=mean(ppl_Annals_2_in);
ppl_Annals_2_cn_sd(c,:)=std(ppl_Annals_2);
ppl_Annals_2_cn_in_sd(c,:)=std(ppl_Annals_2_in);
ppl_Annals_2_cn_med(c,:)=median(ppl_Annals_2);
ppl_Annals_2_cn_in_med(c,:)=median(ppl_Annals_2_in);
% % % ppl_Annals_3_cn(c,:)=mean(ppl_Annals_3);
% % % ppl_Annals_3_cn_in(c,:)=mean(ppl_Annals_3_in);
% % % ppl_Annals_4_cn(c,:)=mean(ppl_Annals_4);
% % % ppl_Annals_4_cn_in(c,:)=mean(ppl_Annals_4_in);

end

time=toc; disp(['elapsed time: ',num2str(time/60),' minutes']);

xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA6_Piston_exapmle_mean.xlsx',ppl_cn(:,[1:8,10,17:19,21])./ppl_cn(:,18),1,'A2'); 
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA6_Piston_exapmle_mean.xlsx',ppl_Annals_2_cn(:,[1:8,10,17:19,21])./ppl_Annals_2_cn(:,18),1,'A6'); 
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA6_Piston_exapmle_std.xlsx',ppl_cn_sd(:,[1:8,10,17:19,21])./ppl_cn_sd(:,18),1,'A2'); 
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA6_Piston_exapmle_std.xlsx',ppl_Annals_2_cn_sd(:,[1:8,10,17:19,21])./ppl_Annals_2_cn_sd(:,18),1,'A6');  
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_exapmle_com_mean.xlsx',ppl_cn(:,10),1,'B2'); 
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_exapmle_com_mean.xlsx',ppl_Annals_2_cn(:,10),1,'B6');
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_exapmle_com_std.xlsx',ppl_cn_sd(:,10),1,'B2'); 
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_exapmle_com_std.xlsx',ppl_Annals_2_cn_sd(:,10),1,'B6');
