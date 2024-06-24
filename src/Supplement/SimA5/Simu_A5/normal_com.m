%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation A.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; tic;
warning('off','all')
rseed=RandStream('mt19937ar','Seed',20222022);
RandStream.setGlobalStream(rseed);
options=optimset('algorithm','interior-point-convex','Display','off','MaxIter',500);

S=500;          
M=2;              
w0=ones(M,1)/M; 

%cs=[0.01,0.1:0.1:0.9,0.91:0.01:0.99];   
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
% % % % % ppl_new_cn=zeros(cn,2*kn+5);
% % % % % ppl_new_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_JASA_cn=zeros(cn,2*kn+5);
% % % % % ppl_JASA_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_Annals_1_cn=zeros(cn,2*kn+5);
% % % % % ppl_Annals_1_cn_in=zeros(cn,2*kn+5);
ppl_Annals_2_cn=zeros(cn,2*kn+5);
ppl_Annals_2_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_Annals_3_cn=zeros(cn,2*kn+5);
% % % % % ppl_Annals_3_cn_in=zeros(cn,2*kn+5);
% % % % % ppl_Annals_4_cn=zeros(cn,2*kn+5);
% % % % % ppl_Annals_4_cn_in=zeros(cn,2*kn+5);
% % % % % mspe_cn=zeros(cn,2*kn+5);
% % % % % mspe_cn_in=zeros(cn,2*kn+5);
% % % % % class_id_cn=zeros(cn,2*kn+5);
% % % % % class_id_cn_in=zeros(cn,2*kn+5);

for c=1:cn                                 
alpha=cs(c);
r1s=zeros(S,2*kn+5);
r1s_in=zeros(S,2*kn+5);
r2s=zeros(S,2*kn+5);
r2s_in=zeros(S,2*kn+5);
ws_cv=zeros(S,kn+1);
wpl_cv=zeros(S,kn);
wpl_loss_cv=zeros(S,kn);
ks_min=zeros(S,kn);
pl_min=zeros(S,kn-1);
pl_loss_min=zeros(S,kn-1);

parfor_progress(S);

parfor s=1:S                                
if floor(s/100)*100==s; disp([c,s/100]); end
% DGP
X_all=2*rand(10^5,10)-1;
% X_all=normrnd(0,1,10^5,10);
myfun=@(theta,X)exp(X*theta);
% e1=randn(10^(5),1);
% e2=randn(10^(5),1);
e1=normrnd(0,1,10^5,1);
e2=lognrnd(0,1,10^5,1);
% theta=[ones(1,5),0.5*ones(1,3),0.1,0]';
% theta=[ones(1,3),0.5*ones(1,2),0.1,zeros(1,4)]';
theta=[ones(1,3),0.5*ones(1,2),0.1,zeros(1,4)]';
Yc1=X_all*theta+e1;
Yc2=myfun(theta,X_all)+e2;
Xc=X_all(Yc1>0 & Yc2>0,:);
Y_p1=Yc1(Yc1>0 & Yc2>0);
Y_p2=Yc2(Yc1>0 & Yc2>0);
X=Xc(1:n+10,:);
loc=binornd(1,alpha*ones(n+10,1));
Y1=Y_p1(1:n+10,:);

% X=rand(n+10,3);
% myfun=@(theta,X)exp(X*theta);
% e1=randn(n+10,1);
% e2=lognrnd(0,1,n+10,1);
% theta=[1,1,1]';
% loc=binornd(1,alpha*ones(n+10,1));
% Y1=X*theta+e1;


Y2=Y_p2(1:n+10,:);
Y=Y1.*(1-loc)+Y2.*loc;

x=X(1:n,:);                                  
y=Y(1:n);                            
beta1=(x(:,1:5)'*x(:,1:5))\(x(:,1:5)'*y);
% t0=(x'*x)\(x'*y);
beta2=(x'*x)\(x'*y);
x0=X(n+1:n+10,:);                           
y0=Y(n+1:n+10);                             
y1=x0(:,1:5)*beta1;                                
y2=x0*beta2;                      
y1_in=x(:,1:5)*beta1;                               
y2_in=x*beta2;                        
yy=[y1,y2];                                 
yin=[y1_in,y2_in];
y3=yy*w0;                                    
y3_in=yin*w0;                               

% Data-driven K 
w_ks=zeros(kn,M);
w_pl=zeros(kn-1,M);
v_ks=zeros(kn,1);
v_pl=zeros(kn-1,1);
loss_pl=zeros(kn-1,1);
y_ks=zeros(kn,length(yy));
y_pl=zeros(kn-1,length(yy));
y_ks_in=zeros(kn,length(yin));
y_pl_in=zeros(kn-1,length(yin));
ppl_ks=zeros(1,kn);
ppl_pl=zeros(1,kn-1);
ppl_ks_in=zeros(1,kn);
ppl_pl_in=zeros(1,kn-1);
% % % % % ppl_new_ks=zeros(1,kn);
% % % % % ppl_new_pl=zeros(1,kn-1);
% % % % % ppl_new_ks_in=zeros(1,kn);
% % % % % ppl_new_pl_in=zeros(1,kn-1);
% % % % % ppl_JASA_ks=zeros(1,kn);
% % % % % ppl_JASA_pl=zeros(1,kn-1);
% % % % % ppl_JASA_ks_in=zeros(1,kn);
% % % % % ppl_JASA_pl_in=zeros(1,kn-1);
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
        [w_k,v_k,wplma,vplma,losspl]=plma_kcv_mix_approximate_lambda_CV_all_Annals_normal(x,y,n,K,M,w0,options);

        w_ks(k,:)=w_k';
        v_ks(k)=v_k;
        y_ks(k,:)=yy*w_k;
        y_ks_in(k,:)=yin*w_k;
        w_pl(k,:)=wplma';
        v_pl(k)=vplma;
        y_pl(k,:)=yy*wplma;
        y_pl_in(k,:)=yin*wplma;
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
        [w_k,v_k]=kcv_mix_normal(x,y,n,K,M,w0,options);
        w_ks(k,:)=w_k;  
        v_ks(k)=v_k;
        y_ks(k,:)=yy*w_k;
        y_ks_in(k,:)=yin*w_k;
        ppl_ks(k)=pairwise_loss(length(y0),y0,yy,w_k');                  
        ppl_ks_in(k)=pairwise_loss(length(y),y,yin,w_k');               

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
w_kmin=w_ks(kmin,:);
y_kmin=yy*w_kmin';
y_kmin_in=yin*w_kmin';
ks_min(s,:)=0+kmin';

% % % plmin=logical(v_pl==min(v_pl));
[~,plmin]=min(v_pl);
w_plmin=w_pl(plmin,:);
y_plmin=yy*w_plmin';
y_plmin_in=yin*w_plmin';
pl_min(s,:)=0+plmin';

[~,plloss_min]=min(loss_pl);
w_pl_loss_min=w_pl(plloss_min,:);
y_pl_loss_min=yy*w_pl_loss_min';
y_pl_loss_min_in=yin*w_pl_loss_min';
pl_loss_min(s,:)=0+plloss_min';

% Prediciton Risk
ys=[y1,y2,y3,y_kmin,y_plmin,y_pl_loss_min,y_ks',y_pl'];
% r1s(s,:)=mean(ys-yT);
r2s(s,:)=mean(ys-y0);

% Risk in sample
ys_in=[y1_in,y2_in,y3_in,y_kmin_in,y_plmin_in,y_pl_loss_min_in,y_ks_in',y_pl_in'];
% r1s_in(s,:)=mean(ys_in-yT_in);
r2s_in(s,:)=mean(ys_in-y);

ws_cv(s,:)=[w_kmin(1),w_ks(:,1)'];
wpl_cv(s,:)=[w_plmin(1),w_pl(:,1)'];
wpl_loss_cv(s,:)=[w_pl_loss_min(1),w_pl(:,1)'];

% Prediciton Pairwise loss 
ppl_y1=pairwise_loss(length(y0),y0,y1,1);
ppl_y2=pairwise_loss(length(y0),y0,y2,1);
ppl_y3=pairwise_loss(length(y0),y0,yy,w0');
ppl_ykmin=pairwise_loss(length(y0),y0,yy,w_kmin);
ppl_yplmin=pairwise_loss(length(y0),y0,yy,w_plmin);
ppl_loss_yplmin=pairwise_loss(length(y0),y0,yy,w_pl_loss_min);
ppl(s,:)=[ppl_y1,ppl_y2,ppl_y3,ppl_ykmin,ppl_yplmin,ppl_loss_yplmin,ppl_ks,ppl_pl];

% Pairwise loss in sample 
ppl_y1_in=pairwise_loss(length(y),y,y1_in,1);
ppl_y2_in=pairwise_loss(length(y),y,y2_in,1);
ppl_y3_in=pairwise_loss(length(y),y,yin,w0');
ppl_ykmin_in=pairwise_loss(length(y),y,yin,w_kmin);
ppl_yplmin_in=pairwise_loss(length(y),y,yin,w_plmin);
ppl_loss_yplmin_in=pairwise_loss(length(y),y,yin,w_pl_loss_min);
ppl_in(s,:)=[ppl_y1_in,ppl_y2_in,ppl_y3_in,ppl_ykmin_in,ppl_yplmin_in,ppl_loss_yplmin_in,ppl_ks_in,ppl_pl_in];

% % % % % % % % % Prediciton Pairwise loss Annals approximate 1
% % % % % ppl_Annals_1_y1=pairwise_loss_approximate_Annals_1(length(y0),y0,y1,1);
% % % % % ppl_Annals_1_y2=pairwise_loss_approximate_Annals_1(length(y0),y0,y2,1);
% % % % % ppl_Annals_1_y3=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w0');
% % % % % ppl_Annals_1_ykmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_kmin);
% % % % % ppl_Annals_1_yplmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_plmin);
% % % % % ppl_loss_Annals_1_yplmin=pairwise_loss_approximate_Annals_1(length(y0),y0,yy,w_pl_loss_min);
% % % % % ppl_Annals_1(s,:)=[ppl_Annals_1_y1,ppl_Annals_1_y2,ppl_Annals_1_y3,ppl_Annals_1_ykmin,ppl_Annals_1_yplmin,ppl_loss_Annals_1_yplmin,ppl_Annals_1_ks,ppl_Annals_1_pl];

% % % % % % % % % Pairwise loss in sample Annals approximate 1
% % % % % ppl_Annals_1_y1_in=pairwise_loss_approximate_Annals_1(length(y),y,y1_in,1);
% % % % % ppl_Annals_1_y2_in=pairwise_loss_approximate_Annals_1(length(y),y,y2_in,1);
% % % % % ppl_Annals_1_y3_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w0');
% % % % % ppl_Annals_1_ykmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_kmin);
% % % % % ppl_Annals_1_yplmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_plmin);
% % % % % ppl_loss_Annals_1_yplmin_in=pairwise_loss_approximate_Annals_1(length(y),y,yin,w_pl_loss_min);
% % % % % ppl_Annals_1_in(s,:)=[ppl_Annals_1_y1_in,ppl_Annals_1_y2_in,ppl_Annals_1_y3_in,ppl_Annals_1_ykmin_in,ppl_Annals_1_yplmin_in,ppl_loss_Annals_1_yplmin_in,ppl_Annals_1_ks_in,ppl_Annals_1_pl_in];

% % % % Prediciton Pairwise loss Annals approximate 2
ppl_Annals_2_y1=pairwise_loss_approximate_Annals_2(length(y0),y0,y1,1);
ppl_Annals_2_y2=pairwise_loss_approximate_Annals_2(length(y0),y0,y2,1);
ppl_Annals_2_y3=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w0');
ppl_Annals_2_ykmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_kmin);
ppl_Annals_2_yplmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_plmin);
ppl_loss_Annals_2_yplmin=pairwise_loss_approximate_Annals_2(length(y0),y0,yy,w_pl_loss_min);
ppl_Annals_2(s,:)=[ppl_Annals_2_y1,ppl_Annals_2_y2,ppl_Annals_2_y3,ppl_Annals_2_ykmin,ppl_Annals_2_yplmin,ppl_loss_Annals_2_yplmin,ppl_Annals_2_ks,ppl_Annals_2_pl];

% % % % Pairwise loss in sample Annals approximate 2
ppl_Annals_2_y1_in=pairwise_loss_approximate_Annals_2(length(y),y,y1_in,1);
ppl_Annals_2_y2_in=pairwise_loss_approximate_Annals_2(length(y),y,y2_in,1);
ppl_Annals_2_y3_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w0');
ppl_Annals_2_ykmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_kmin);
ppl_Annals_2_yplmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_plmin);
ppl_loss_Annals_2_yplmin_in=pairwise_loss_approximate_Annals_2(length(y),y,yin,w_pl_loss_min);
ppl_Annals_2_in(s,:)=[ppl_Annals_2_y1_in,ppl_Annals_2_y2_in,ppl_Annals_2_y3_in,ppl_Annals_2_ykmin_in,ppl_Annals_2_yplmin_in,ppl_loss_Annals_2_yplmin_in,ppl_Annals_2_ks_in,ppl_Annals_2_pl_in];

% % % % % % % % % Prediciton Pairwise loss Annals approximate 3
% % % % % ppl_Annals_3_y1=pairwise_loss_approximate_Annals_3(length(y0),y0,y1,1);
% % % % % ppl_Annals_3_y2=pairwise_loss_approximate_Annals_3(length(y0),y0,y2,1);
% % % % % ppl_Annals_3_y3=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w0');
% % % % % ppl_Annals_3_ykmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_kmin);
% % % % % ppl_Annals_3_yplmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_plmin);
% % % % % ppl_loss_Annals_3_yplmin=pairwise_loss_approximate_Annals_3(length(y0),y0,yy,w_pl_loss_min);
% % % % % ppl_Annals_3(s,:)=[ppl_Annals_3_y1,ppl_Annals_3_y2,ppl_Annals_3_y3,ppl_Annals_3_ykmin,ppl_Annals_3_yplmin,ppl_loss_Annals_3_yplmin,ppl_Annals_3_ks,ppl_Annals_3_pl];

% % % % % % % % % Pairwise loss in sample Annals approximate 3
% % % % % ppl_Annals_3_y1_in=pairwise_loss_approximate_Annals_3(length(y),y,y1_in,1);
% % % % % ppl_Annals_3_y2_in=pairwise_loss_approximate_Annals_3(length(y),y,y2_in,1);
% % % % % ppl_Annals_3_y3_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w0');
% % % % % ppl_Annals_3_ykmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_kmin);
% % % % % ppl_Annals_3_yplmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_plmin);
% % % % % ppl_loss_Annals_3_yplmin_in=pairwise_loss_approximate_Annals_3(length(y),y,yin,w_pl_loss_min);
% % % % % ppl_Annals_3_in(s,:)=[ppl_Annals_3_y1_in,ppl_Annals_3_y2_in,ppl_Annals_3_y3_in,ppl_Annals_3_ykmin_in,ppl_Annals_3_yplmin_in,ppl_loss_Annals_3_yplmin_in,ppl_Annals_3_ks_in,ppl_Annals_3_pl_in];

% % % % % % % % % Prediciton Pairwise loss Annals approximate 4
% % % % % ppl_Annals_4_y1=pairwise_loss_approximate_Annals_4(length(y0),y0,y1,1);
% % % % % ppl_Annals_4_y2=pairwise_loss_approximate_Annals_4(length(y0),y0,y2,1);
% % % % % ppl_Annals_4_y3=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w0');
% % % % % ppl_Annals_4_ykmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_kmin);
% % % % % ppl_Annals_4_yplmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_plmin);
% % % % % ppl_loss_Annals_4_yplmin=pairwise_loss_approximate_Annals_4(length(y0),y0,yy,w_pl_loss_min);
% % % % % ppl_Annals_4(s,:)=[ppl_Annals_4_y1,ppl_Annals_4_y2,ppl_Annals_4_y3,ppl_Annals_4_ykmin,ppl_Annals_4_yplmin,ppl_loss_Annals_4_yplmin,ppl_Annals_4_ks,ppl_Annals_4_pl];

% % % % % % % % % Pairwise loss in sample Annals approximate 4
% % % % % ppl_Annals_4_y1_in=pairwise_loss_approximate_Annals_4(length(y),y,y1_in,1);
% % % % % ppl_Annals_4_y2_in=pairwise_loss_approximate_Annals_4(length(y),y,y2_in,1);
% % % % % ppl_Annals_4_y3_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w0');
% % % % % ppl_Annals_4_ykmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_kmin);
% % % % % ppl_Annals_4_yplmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_plmin);
% % % % % ppl_loss_Annals_4_yplmin_in=pairwise_loss_approximate_Annals_4(length(y),y,yin,w_pl_loss_min);
% % % % % ppl_Annals_4_in(s,:)=[ppl_Annals_4_y1_in,ppl_Annals_4_y2_in,ppl_Annals_4_y3_in,ppl_Annals_4_ykmin_in,ppl_Annals_4_yplmin_in,ppl_loss_Annals_4_yplmin_in,ppl_Annals_4_ks_in,ppl_Annals_4_pl_in];



parfor_progress;

end  % s %

parfor_progress(0);

% r1_cn(c,:)=mean(r1s.^2);
% r1_cn_in(c,:)=mean(r1s_in.^2);
r2_cn(c,:)=mean(r2s.^2);
r2_cn_in(c,:)=mean(r2s_in.^2);
w_cn(c,:)=mean(ws_cv);
w_cn_pl(c,:)=mean(wpl_cv);
w_cn_pl_loss(c,:)=mean(wpl_loss_cv);
k_cn(c,:)=mean(ks_min);
pl_cn(c,:)=mean(pl_min);
pl_loss_cn(c,:)=mean(pl_loss_min);
% rank_score_cn(c,:)=mean(rank_score);
% rank_score_cn_in(c,:)=mean(rank_score_in);
ppl_cn(c,:)=mean(ppl);
ppl_cn_in(c,:)=mean(ppl_in);
% % % % % ppl_new_cn(c,:)=mean(ppl_new);
% % % % % ppl_new_cn_in(c,:)=mean(ppl_in_new);
% % % % % ppl_JASA_cn(c,:)=mean(ppl_JASA);
% % % % % ppl_JASA_cn_in(c,:)=mean(ppl_JASA_in);
% % % % % ppl_Annals_1_cn(c,:)=mean(ppl_Annals_1);
% % % % % ppl_Annals_1_cn_in(c,:)=mean(ppl_Annals_1_in);
ppl_Annals_2_cn(c,:)=mean(ppl_Annals_2);
ppl_Annals_2_cn_in(c,:)=mean(ppl_Annals_2_in);
% % % % % ppl_Annals_3_cn(c,:)=mean(ppl_Annals_3);
% % % % % ppl_Annals_3_cn_in(c,:)=mean(ppl_Annals_3_in);
% % % % % ppl_Annals_4_cn(c,:)=mean(ppl_Annals_4);
% % % % % ppl_Annals_4_cn_in(c,:)=mean(ppl_Annals_4_in);
% % % % % mspe_cn(c,:)=mean(MSPE);
% % % % % mspe_cn_in(c,:)=mean(MSPE_in);
% class_id_cn(c,:)=mean(class_id);
% class_id_cn_in(c,:)=mean(class_id_in);

end  % c %

time=toc; disp(['elapsed time: ',num2str(time/60),' minutes']);
if jj==1; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_Pl50.xlsx',[ppl_Annals_2_cn(:,[1:2,6,17:18,22,25])./ppl_Annals_2_cn(:,18);ppl_cn(:,[1:2,6,17:18,22,25])./ppl_cn(:,18)]); 
    
end
if jj==2;
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_comprae.xlsx',[ppl_Annals_2_cn(:,6);ppl_cn(:,6)],1,'A1'); 
end         
if jj==3; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_comprae.xlsx',[ppl_Annals_2_cn(:,6);ppl_cn(:,6)],1,'A25');     
end   
if jj==4; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_comprae.xlsx',[ppl_Annals_2_cn(:,6);ppl_cn(:,6)],1,'A49');     
end     
if jj==5; 
    xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_comprae.xlsx',[ppl_Annals_2_cn(:,6);ppl_cn(:,6)],1,'A73');         
end 

end  % jj %


