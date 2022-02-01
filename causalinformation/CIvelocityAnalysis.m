%Causality information rate experiments
%Adrian J Guel C
%May 2021
clear all;
close all;
clc

D= @(D0,a,b,t0,t) D0+b*(exp(-((t-t0).^2)/(a.^2))/(sqrt(pi)*abs(a)));

x0=-0.5;
v0=0.5;
d11=.01;
d12=0;
d21=0;
d22=.01;
a=0.5;
c=0.5;
t0=4;
t20=8;
b=0;
d=0;

tf1=40;
tf2=10;

%D0=0;
D0=[1e-4];
%D0=flip(D0);
w=1;
g=[0 1];
clrs={'k' ':b' '-.r' ':k' '--b' '--r' '-.k'};
%clrs=flip(clrs);

plt=0; %Plot T and E then plt=1, plt=0 otherwise. 

for j=1:length(g)
    fig = figure;
    left_color = [0 0 0];
    right_color = [0 0 1];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    set(fig,'color','w');
    if plt==1
        set(fig, 'Position',  [100, 100, 900, 450])
    else
        set(fig, 'Position',  [100, 100, 800, 250])
    end

    for i=1:length(D0)
         if j==1
            t=0:0.001:tf1;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0(i)==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);                        
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);                        
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                    end
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;                           
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));                  
            elseif g(j)==2*w
                S=SigmaCD(w,t,t0,a,b,d11,d12,d22,D0(i));
                Xm=meanscd(w,t,t20,c,d,x0,v0);
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;        
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));  
            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);                    
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;         
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));            
            end
         else
            t=0:0.001:tf2;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0(i)==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);                        
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);                        
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                    end
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;           
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*S(1,:).^2);
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));                  
            elseif g(j)==2*w
                S=SigmaCD(w,t,t0,a,b,d11,d12,d22,D0(i));
                Xm=meanscd(w,t,t20,c,d,x0,v0);
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;            
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*S(1,:).^2);
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));  
            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);                    
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;           
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*S(1,:).^2);
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));            
            end
         end
         if plt==1
                 if b==0 && d==0            
                      subplot(2,3,1)
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                  
                       subplot(2,3,2)
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(2,3,3)
                            plot(real(Txv),real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')               
                       subplot(2,3,4)
                            indx=find(abs(real(Evx))<4e2);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(2,3,5)
                            indx=find(abs(real(Exv))<4e2);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(2,3,6)
                            indx=find(abs(real(Exv))<4e2);
                            plot(real(Exv(indx)),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')               
                 elseif b==1 && d==0
                      subplot(2,3,1)
                            yyaxis left
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t,D(0,a,b,t0,t),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,2)
                            yyaxis left
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t,D(0,a,b,t0,t),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')
                       subplot(2,3,3)
                            plot(real(Txv),real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                      
                       subplot(2,3,4)
                            yyaxis left
                            indx=find(abs(real(Evx))<4e2);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,5)
                            yyaxis left
                            indx=find(abs(real(Exv))<4e2);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                             yyaxis right
                             plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                        subplot(2,3,6)
                            indx=find(abs(real(Exv))<4e2);
                            plot(real(Exv(indx)),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                      
                 elseif b==0 && d==1
                      subplot(2,3,1)
                            yyaxis left
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,2)
                            yyaxis left
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square')
                       subplot(2,3,3)
                            plot(real(Txv),real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                      
                       subplot(2,3,4)
                            yyaxis left
                            indx=find(abs(real(Evx))<1e13);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,5)
                            yyaxis left
                            indx=find(abs(real(Exv))<1e13);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                             yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                        subplot(2,3,6)
                            indx=find(abs(real(Exv))<1e13);
                            plot(real(Exv(indx)),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')              
                 else
                      subplot(2,3,1)
                            yyaxis left
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                             hold on
                             plot(t,D(0,a,b,t0,t),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,2)
                            yyaxis left
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                             hold on
                             plot(t,D(0,a,b,t0,t),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')  
                       subplot(2,3,3)
                            plot(real(Txv),real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                      
                       subplot(2,3,4)
                            yyaxis left
                            indx=find(abs(real(Evx))<1e13);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                             hold on
                              plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,5)
                            yyaxis left
                            indx=find(abs(real(Exv))<1e13);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                             hold on
                              plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')  
                        subplot(2,3,6)
                            indx=find(abs(real(Exv))<1e13);
                            plot(real(Exv(indx)),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')             
                 end
         else             
                       subplot(1,3,1)
                            indx=find(abs(real(Evx))<4e2);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(1,3,2)
                            indx=find(abs(real(Exv))<4e2);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(1,3,3)
                            indx=find(abs(real(Exv))<4e2);
                            plot(real(Exv(indx)),real(Evx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                 
         end
    end
    str = sprintf('\\gamma=%.2f, \\omega=%0.2f, D_0=%0.2f', g(j),w,D0(i));
    sgtitle(str)
end
%AddLabels(D0,clrs)