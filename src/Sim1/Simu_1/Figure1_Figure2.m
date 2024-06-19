%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  Simulation 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%  This Matlab program executes the simulation reported in
%  "RMA: Ranking based on model averaging."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; tic;
warning('off','all')
rseed=RandStream('mt19937ar','Seed',202301);
RandStream.setGlobalStream(rseed);
options=optimset('algorithm','interior-point-convex','Display','off','MaxIter',500);

S=500;          
V=1;
p=4;
r=16;
q=p+r;
M=2;
w0=ones(M,1)/M;  

rho=0.5;
% sig_x=rho*(rho*(ones(r,r)-eye(r))+eye(r)); 
% sig_x=rho*(ones(q,q)-eye(q))+eye(q);          % correlation
sig_x=eye(q);                                   % uncorrelation
cs=[0.01,0.1:0.1:0.9,0.99];
cn=length(cs);

for jj=2:5 
if jj==1; n=50; Ks=floor([2,5,n/8,n/7,n/6,n/5,n/4,n/3,n/2,n]);  end 
if jj==2; n=100; Ks=floor([2,5,10,n/8,n/6,n/5,n/4,n/3,n/2,n]);  end 
if jj==3; n=200; Ks=floor([2,5,10,20,30,n/5,n/4,n/3,n/2,n]);  end 
if jj==4; n=500; Ks=floor([2,5,10,20,50,n/5,n/4,n/3,n/2,n]);  end
if jj==5; n=1000; Ks=floor([2,5,10,20,50,100,n/5,n/3,n/2,n]);  end

kn=length(Ks);
r1_cn=zeros(cn,2*kn+5);
r1_cn_in=zeros(cn,2*kn+5);
r2_cn=zeros(cn,2*kn+5);
r2_cn_in=zeros(cn,2*kn+5);
w_cn=zeros(cn,kn+1);
w_cn_pl=zeros(cn,kn);
w_cn_pl_loss=zeros(cn,kn);
k_cn=zeros(cn,kn);
pl_cn=zeros(cn,kn-1);
pl_loss_cn=zeros(cn,kn-1);
rank_score_cn=zeros(cn,2*kn+5);
rank_score_cn_in=zeros(cn,2*kn+5);
ppl_cn=zeros(cn,2*kn+5);
ppl_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_Annals_1_cn=zeros(cn,2*kn+5);
% % % % % ppl_Annals_1_cn_in=zeros(cn,2*kn+5);
ppl_Annals_2_cn=zeros(cn,2*kn+5);
ppl_Annals_2_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_Annals_3_cn=zeros(cn,2*kn+5);
% % % % % ppl_Annals_3_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_Annals_4_cn=zeros(cn,2*kn+5);
% % % % % ppl_Annals_4_cn_in=zeros(cn,2*kn+5);

for c=1:cn 
beta=[1;1;1;cs(c);zeros(r,1)];
r1s=zeros(S,2*kn+5);
r1s_in=zeros(S,2*kn+5);
r2s=zeros(S,2*kn+5);
r2s_in=zeros(S,2*kn+5);
ws_cv=zeros(S,kn+1);
wpl_cv=zeros(S,kn);
wpl_loss_cv=zeros(S,kn);
rws=zeros(M,2*kn+3);
ks_min=zeros(S,kn);
pl_min=zeros(S,kn-1);
pl_loss_min=zeros(S,kn-1);
ppl=zeros(S,2*kn+5);
ppl_in=zeros(S,2*kn+5);
% % % % % ppl_Annals_1=zeros(S,2*kn+5);
% % % % % ppl_Annals_1_in=zeros(S,2*kn+5);
ppl_Annals_2=zeros(S,2*kn+5);
ppl_Annals_2_in=zeros(S,2*kn+5);
% % % % % ppl_Annals_3=zeros(S,2*kn+5);
% % % % % ppl_Annals_3_in=zeros(S,2*kn+5);
% % % % % ppl_Annals_4=zeros(S,2*kn+5);
% % % % % ppl_Annals_4_in=zeros(S,2*kn+5);

parfor_progress(S);
parfor s=1:S
if floor(s/50)*50==s; disp([n,c,s/50]); end
% DGP
X=mvnrnd(zeros(1,q),sig_x,n+10);        
if V==1; sig_e=1; end            
if V==2; sig_e=0.5+0.5*(X(:,p).^2);  end 
e=sqrt(sig_e).*randn(n+10,1);            
Y=X*beta+e;

% LS Estimation
x=X(1:n,:);
y=Y(1:n);
x0=X(n+1:n+10,:);                        
y0=Y(n+1:n+10);                          
yT=x0*beta;
yT_in=x*beta;
x1=x(:,1:(p-1));
beta1=(x1'*x1)\(x1'*y);  
beta2=(x'*x)\(x'*y); 
y1=x0(:,1:(p-1))*beta1;
y1_in=x(:,1:(p-1))*beta1;
y2=x0*beta2;
y2_in=x*beta2;
yy=[y1 y2];
yin=[y1_in y2_in];
y3=yy*w0;
y3_in=yin*w0;

% Data-driven K 
w_ks=zeros(M,kn);
w_pl=zeros(M,kn-1);
v_ks=zeros(1,kn);
v_pl=zeros(1,kn-1);
loss_pl=zeros(1,kn-1);
y_ks=zeros(size(yy,1),kn);
y_pl=zeros(size(yy,1),kn-1);
y_ks_in=zeros(size(yin,1),kn);
y_pl_in=zeros(size(yin,1),kn-1);
ppl_ks=zeros(1,kn);
ppl_pl=zeros(1,kn-1);
ppl_ks_in=zeros(1,kn);
ppl_pl_in=zeros(1,kn-1);
% % % ppl_Annals_1_ks=zeros(1,kn);
% % % ppl_Annals_1_pl=zeros(1,kn-1);
% % % ppl_Annals_1_ks_in=zeros(1,kn);
% % % ppl_Annals_1_pl_in=zeros(1,kn-1);
ppl_Annals_2_ks=zeros(1,kn);
ppl_Annals_2_pl=zeros(1,kn-1);
ppl_Annals_2_ks_in=zeros(1,kn);
ppl_Annals_2_pl_in=zeros(1,kn-1);
% % % ppl_Annals_3_ks=zeros(1,kn);
% % % ppl_Annals_3_pl=zeros(1,kn-1);
% % % ppl_Annals_3_ks_in=zeros(1,kn);
% % % ppl_Annals_3_pl_in=zeros(1,kn-1);
% % % ppl_Annals_4_ks=zeros(1,kn);
% % % ppl_Annals_4_pl=zeros(1,kn-1);
% % % ppl_Annals_4_ks_in=zeros(1,kn);
% % % ppl_Annals_4_pl_in=zeros(1,kn-1);

for k=1:kn
    if k~=kn
        K=Ks(k);
        [w_k,v_k,wplma,vplma,losspl]=plma_kcv_lsM2_approximate_lambda_CV_all_Annals(x,y,n,p,K,M,w0,options);
        w_ks(:,k)=w_k;  
        v_ks(k)=v_k;
        y_ks(:,k)=yy*w_k; 
        y_ks_in(:,k)=yin*w_k; 
        w_pl(:,k)=wplma;
        v_pl(k)=vplma;    
        y_pl(:,k)=yy*wplma;
        y_pl_in(:,k)=yin*wplma;
        loss_pl(k)=losspl;
        ppl_ks(k)=pairwise_loss(length(y0),y0,yy,w_k');      
        ppl_pl(k)=pairwise_loss(length(y0),y0,yy,wplma');
        ppl_ks_in(k)=pairwise_loss(length(y),y,yin,w_k');      
        ppl_pl_in(k)=pairwise_loss(length(y),y,yin,wplma');
        
% % % % %         ppl_Annals_1_ks(k)=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_k'); 
% % % % %         ppl_Annals_1_pl(k)=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,wplma');  
% % % % %         ppl_Annals_1_ks_in(k)=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_k'); 
% % % % %         ppl_Annals_1_pl_in(k)=pairwise_loss_approximate_Annals_1(length(y),y,yin,wplma'); 

        ppl_Annals_2_ks(k)=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_k');  
        ppl_Annals_2_pl(k)=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,wplma'); 
        ppl_Annals_2_ks_in(k)=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_k');  
        ppl_Annals_2_pl_in(k)=pairwise_loss_approximate_Annals_2(length(y),y,yin,wplma'); 

% % % % %         ppl_Annals_3_ks(k)=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_k'); 
% % % % %         ppl_Annals_3_pl(k)=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,wplma');
% % % % %         ppl_Annals_3_ks_in(k)=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_k'); 
% % % % %         ppl_Annals_3_pl_in(k)=pairwise_loss_approximate_Annals_3(length(y),y,yin,wplma');

% % % % %         ppl_Annals_4_ks(k)=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_k'); 
% % % % %         ppl_Annals_4_pl(k)=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,wplma');  
% % % % %         ppl_Annals_4_ks_in(k)=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_k'); 
% % % % %         ppl_Annals_4_pl_in(k)=pairwise_loss_approximate_Annals_4(length(y),y,yin,wplma');  

    elseif k==kn
        K=Ks(k);
        [w_k,v_k]=kcv_lsM2(x,y,n,p,K,M,w0,options);
        w_ks(:,k)=w_k;  
        v_ks(k)=v_k;
        y_ks(:,k)=yy*w_k;
        y_ks_in(:,k)=yin*w_k;
        ppl_ks(k)=pairwise_loss(length(y0),y0,yy,w_k'); 
        ppl_ks_in(k)=pairwise_loss(length(y0),y0,yy,w_k'); 

% % % % %         ppl_Annals_1_ks(k)=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_k'); 
% % % % %         ppl_Annals_1_ks_in(k)=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_k'); 

        ppl_Annals_2_ks(k)=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_k'); 
        ppl_Annals_2_ks_in(k)=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_k'); 

% % % % %         ppl_Annals_3_ks(k)=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_k'); 
% % % % %         ppl_Annals_3_ks_in(k)=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_k'); 

% % % % %         ppl_Annals_4_ks(k)=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_k'); 
% % % % %         ppl_Annals_4_ks_in(k)=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_k'); 

    end
end

% % % kmin=logical(v_ks==min(v_ks));
[~,kmin]=min(v_ks);
w_kmin=w_ks(:,kmin);
y_kmin=yy*w_kmin;
y_kmin_in=yin*w_kmin;
ks_min(s,:)=0+kmin;

% % % plmin=logical(v_pl==min(v_pl));
[~,plmin]=min(v_pl);
w_plmin=w_pl(:,plmin);
y_plmin=yy*w_plmin;
y_plmin_in=yin*w_plmin;
pl_min(s,:)=0+plmin;

[~,plloss_min]=min(loss_pl);
w_pl_loss_min=w_pl(:,plloss_min);
y_pl_loss_min=yy*w_pl_loss_min;
y_pl_loss_min_in=yin*w_pl_loss_min;
pl_loss_min(s,:)=0+plloss_min;

% Prediciton Risk
ys=[y1,y2,y3,y_kmin,y_plmin,y_pl_loss_min,y_ks,y_pl];
r1s(s,:)=mean(ys-yT);
r2s(s,:)=mean(ys-y0);
[~,I0]=sort(y0,'descend');
[~,I]=sort(ys,'descend');
rank_score(s,:)=sum(I0==I,1)/length(y0);
ws_cv(s,:)=[w_kmin(1),w_ks(1,:)];
wpl_cv(s,:)=[w_plmin(1),w_pl(1,:)];
wpl_loss_cv(s,:)=[w_pl_loss_min(1),w_pl(1,:)];

rws=rws+[w0,w_kmin,w_plmin,w_pl_loss_min,w_ks,w_pl];

% Risk in sample
ys_in=[y1_in,y2_in,y3_in,y_kmin_in,y_plmin_in,y_pl_loss_min_in,y_ks_in,y_pl_in];
r1s_in(s,:)=mean(ys_in-yT_in);
r2s_in(s,:)=mean(ys_in-y);
[~,I0_in]=sort(y,'descend');
[~,I_in]=sort(ys_in,'descend');
rank_score_in(s,:)=sum(I0_in==I_in,1)/length(y);

% Prediciton Pairwise loss
ppl_y1=pairwise_loss(length(y0),y0,y1,1);
ppl_y2=pairwise_loss(length(y0),y0,y2,1);
ppl_y3=pairwise_loss(length(y0),y0,yy,w0');
ppl_ykmin=pairwise_loss(length(y0),y0,yy,w_kmin');
ppl_yplmin=pairwise_loss(length(y0),y0,yy,w_plmin');
ppl_loss_yplmin=pairwise_loss(length(y0),y0,yy,w_pl_loss_min');
ppl(s,:)=[ppl_y1,ppl_y2,ppl_y3,ppl_ykmin,ppl_yplmin,ppl_loss_yplmin,ppl_ks,ppl_pl];

% Pairwise loss in sample 
ppl_y1_in=pairwise_loss(length(y),y,y1_in,1);
ppl_y2_in=pairwise_loss(length(y),y,y2_in,1);
ppl_y3_in=pairwise_loss(length(y),y,yin,w0');
ppl_ykmin_in=pairwise_loss(length(y),y,yin,w_kmin');
ppl_yplmin_in=pairwise_loss(length(y),y,yin,w_plmin');
ppl_loss_yplmin_in=pairwise_loss(length(y),y,yin,w_pl_loss_min');
ppl_in(s,:)=[ppl_y1_in,ppl_y2_in,ppl_y3_in,ppl_ykmin_in,ppl_yplmin_in,ppl_loss_yplmin_in,ppl_ks_in,ppl_pl_in];

% % % % % % Prediciton Pairwise loss Annals 1
% % % % % ppl_Annals_1_y1=pairwise_loss_approximate_Annals_1(length(y0),y0,y1,1);
% % % % % ppl_Annals_1_y2=pairwise_loss_approximate_Annals_1(length(y0),y0,y2,1);
% % % % % ppl_Annals_1_y3=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w0');
% % % % % ppl_Annals_1_ykmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_kmin');
% % % % % ppl_Annals_1_yplmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_plmin');
% % % % % ppl_loss_Annals_1_yplmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_pl_loss_min');
% % % % % ppl_Annals_1(s,:)=[ppl_Annals_1_y1,ppl_Annals_1_y2,ppl_Annals_1_y3,ppl_Annals_1_ykmin,ppl_Annals_1_yplmin,ppl_loss_Annals_1_yplmin,ppl_Annals_1_ks,ppl_Annals_1_pl];

% % % % % % % % % Pairwise loss in sample Annals approximate 1
% % % % % ppl_Annals_1_y1_in=pairwise_loss_approximate_Annals_1(length(y),y,y1_in,1);
% % % % % ppl_Annals_1_y2_in=pairwise_loss_approximate_Annals_1(length(y),y,y2_in,1);
% % % % % ppl_Annals_1_y3_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w0');
% % % % % ppl_Annals_1_ykmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_kmin');
% % % % % ppl_Annals_1_yplmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_plmin');
% % % % % ppl_loss_Annals_1_yplmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_pl_loss_min');
% % % % % ppl_Annals_1_in(s,:)=[ppl_Annals_1_y1_in,ppl_Annals_1_y2_in,ppl_Annals_1_y3_in,ppl_Annals_1_ykmin_in,ppl_Annals_1_yplmin_in,ppl_loss_Annals_1_yplmin_in,ppl_Annals_1_ks_in,ppl_Annals_1_pl_in];

% Prediciton Pairwise loss Annals 2
ppl_Annals_2_y1=pairwise_loss_approximate_Annals_2(length(y0),y0,y1,1);
ppl_Annals_2_y2=pairwise_loss_approximate_Annals_2(length(y0),y0,y2,1);
ppl_Annals_2_y3=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w0');
ppl_Annals_2_ykmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_kmin');
ppl_Annals_2_yplmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_plmin');
ppl_loss_Annals_2_yplmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_pl_loss_min');
ppl_Annals_2(s,:)=[ppl_Annals_2_y1,ppl_Annals_2_y2,ppl_Annals_2_y3,ppl_Annals_2_ykmin,ppl_Annals_2_yplmin,ppl_loss_Annals_2_yplmin,ppl_Annals_2_ks,ppl_Annals_2_pl];

% % % % Pairwise loss in sample Annals approximate 2
ppl_Annals_2_y1_in=pairwise_loss_approximate_Annals_2(length(y),y,y1_in,1);
ppl_Annals_2_y2_in=pairwise_loss_approximate_Annals_2(length(y),y,y2_in,1);
ppl_Annals_2_y3_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w0');
ppl_Annals_2_ykmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_kmin');
ppl_Annals_2_yplmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_plmin');
ppl_loss_Annals_2_yplmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_pl_loss_min');
ppl_Annals_2_in(s,:)=[ppl_Annals_2_y1_in,ppl_Annals_2_y2_in,ppl_Annals_2_y3_in,ppl_Annals_2_ykmin_in,ppl_Annals_2_yplmin_in,ppl_loss_Annals_2_yplmin_in,ppl_Annals_2_ks_in,ppl_Annals_2_pl_in];

% % % % % % Prediciton Pairwise loss Annals 3
% % % % % ppl_Annals_3_y1=pairwise_loss_approximate_Annals_3(length(y0),y0,y1,1);
% % % % % ppl_Annals_3_y2=pairwise_loss_approximate_Annals_3(length(y0),y0,y2,1);
% % % % % ppl_Annals_3_y3=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w0');
% % % % % ppl_Annals_3_ykmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_kmin');
% % % % % ppl_Annals_3_yplmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_plmin');
% % % % % ppl_loss_Annals_3_yplmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_pl_loss_min');
% % % % % ppl_Annals_3(s,:)=[ppl_Annals_3_y1,ppl_Annals_3_y2,ppl_Annals_3_y3,ppl_Annals_3_ykmin,ppl_Annals_3_yplmin,ppl_loss_Annals_3_yplmin,ppl_Annals_3_ks,ppl_Annals_3_pl];

% % % % % % % % % Pairwise loss in sample Annals approximate 3
% % % % % ppl_Annals_3_y1_in=pairwise_loss_approximate_Annals_3(length(y),y,y1_in,1);
% % % % % ppl_Annals_3_y2_in=pairwise_loss_approximate_Annals_3(length(y),y,y2_in,1);
% % % % % ppl_Annals_3_y3_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w0');
% % % % % ppl_Annals_3_ykmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_kmin');
% % % % % ppl_Annals_3_yplmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_plmin');
% % % % % ppl_loss_Annals_3_yplmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_pl_loss_min');
% % % % % ppl_Annals_3_in(s,:)=[ppl_Annals_3_y1_in,ppl_Annals_3_y2_in,ppl_Annals_3_y3_in,ppl_Annals_3_ykmin_in,ppl_Annals_3_yplmin_in,ppl_loss_Annals_3_yplmin_in,ppl_Annals_3_ks_in,ppl_Annals_3_pl_in];

% % % % % % Prediciton Pairwise loss Annals 4
% % % % % ppl_Annals_4_y1=pairwise_loss_approximate_Annals_4(length(y0),y0,y1,1);
% % % % % ppl_Annals_4_y2=pairwise_loss_approximate_Annals_4(length(y0),y0,y2,1);
% % % % % ppl_Annals_4_y3=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w0');
% % % % % ppl_Annals_4_ykmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_kmin');
% % % % % ppl_Annals_4_yplmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_plmin');
% % % % % ppl_loss_Annals_4_yplmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_pl_loss_min');
% % % % % ppl_Annals_4(s,:)=[ppl_Annals_4_y1,ppl_Annals_4_y2,ppl_Annals_4_y3,ppl_Annals_4_ykmin,ppl_Annals_4_yplmin,ppl_loss_Annals_4_yplmin,ppl_Annals_4_ks,ppl_Annals_4_pl];
% % % % % 
% % % % % % % % % Pairwise loss in sample Annals approximate 4
% % % % % ppl_Annals_4_y1_in=pairwise_loss_approximate_Annals_4(length(y),y,y1_in,1);
% % % % % ppl_Annals_4_y2_in=pairwise_loss_approximate_Annals_4(length(y),y,y2_in,1);
% % % % % ppl_Annals_4_y3_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w0');
% % % % % ppl_Annals_4_ykmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_kmin');
% % % % % ppl_Annals_4_yplmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_plmin');
% % % % % ppl_loss_Annals_4_yplmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_pl_loss_min');
% % % % % ppl_Annals_4_in(s,:)=[ppl_Annals_4_y1_in,ppl_Annals_4_y2_in,ppl_Annals_4_y3_in,ppl_Annals_4_ykmin_in,ppl_Annals_4_yplmin_in,ppl_loss_Annals_4_yplmin_in,ppl_Annals_4_ks_in,ppl_Annals_4_pl_in];


parfor_progress;

end  % s %

parfor_progress(0);

r1_cn(c,:)=mean(r1s.^2);
r1_cn_in(c,:)=mean(r1s_in.^2);
r2_cn(c,:)=mean(r2s.^2);
r2_cn_in(c,:)=mean(r2s_in.^2);
w_cn(c,:)=mean(ws_cv);
w_cn_pl(c,:)=mean(wpl_cv);
w_cn_pl_loss(c,:)=mean(wpl_loss_cv);
w(c,:)=rws(1,:)/S;
k_cn(c,:)=mean(ks_min);
pl_cn(c,:)=mean(pl_min);
pl_loss_cn(c,:)=mean(pl_loss_min);
rank_score_cn(c,:)=mean(rank_score);
rank_score_cn_in(c,:)=mean(rank_score_in);
ppl_cn(c,:)=mean(ppl);
ppl_cn_in(c,:)=mean(ppl_in);
% % % ppl_Annals_1_cn(c,:)=mean(ppl_Annals_1);
% % % ppl_Annals_1_cn_in(c,:)=mean(ppl_Annals_1_in);
ppl_Annals_2_cn(c,:)=mean(ppl_Annals_2);
ppl_Annals_2_cn_in(c,:)=mean(ppl_Annals_2_in);
% % % ppl_Annals_3_cn(c,:)=mean(ppl_Annals_3);
% % % ppl_Annals_3_cn_in(c,:)=mean(ppl_Annals_3_in);
% % % ppl_Annals_4_cn(c,:)=mean(ppl_Annals_4);
% % % ppl_Annals_4_cn_in(c,:)=mean(ppl_Annals_4_in);

end  % c %

time=toc; disp(['elapsed time: ',num2str(time/60),' minutes']);
if jj==1; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_PL50.xlsx',[ppl_Annals_2_cn(:,[1:2,6,17:18,22,25])./ppl_Annals_2_cn(:,18);ppl_cn(:,[1:2,6,17:18,22,25])./ppl_cn(:,18)]); 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_w50.xlsx',1-w(:,[4,15,16,20,23])); 
end
if jj==2; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_PL.xlsx',[ppl_Annals_2_cn(:,[1:2,6,17:19,25])./ppl_Annals_2_cn(:,18);ppl_cn(:,[1:2,6,17:19,25])./ppl_cn(:,18)],1,'A1'); 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_w.xlsx',1-w(:,[4,15:17,23]),1,'A1'); 
end
if jj==3; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_PL.xlsx',[ppl_Annals_2_cn(:,[1:2,6,17:19,25])./ppl_Annals_2_cn(:,18);ppl_cn(:,[1:2,6,17:19,25])./ppl_cn(:,18)],1,'A25'); 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_w.xlsx',1-w(:,[4,15:17,23]),1,'A14'); 
end
if jj==4; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_PL.xlsx',[ppl_Annals_2_cn(:,[1:2,6,17:19,25])./ppl_Annals_2_cn(:,18);ppl_cn(:,[1:2,6,17:19,25])./ppl_cn(:,18)],1,'A49'); 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_w.xlsx',1-w(:,[4,15:17,23]),1,'A27'); 
end
if jj==5; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_PL.xlsx',[ppl_Annals_2_cn(:,[1:2,6,17:19,25])./ppl_Annals_2_cn(:,18);ppl_cn(:,[1:2,6,17:19,25])./ppl_cn(:,18)],1,'A73'); 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Sim1_w.xlsx',1-w(:,[4,15:17,23]),1,'A40'); 
end

end