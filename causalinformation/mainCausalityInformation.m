%Causal Information of the Kramers equation
%Plots of Exv, Evx, Txv, Tvx 
%Adrian J Guel C
%April 2021
clear all;
close all;
clc

%impulse input
D= @(D0,a,b,t0,t) D0+b*(exp(-((t-t0).^2)/(a.^2))/(sqrt(pi)*abs(a)));

%Kramers equation parameters
x0=-0.5;
v0=0.5;
d11=.01;
d12=0;
d21=0;
d22=.01;

%Settings of impulse functions
a=0.5;%1/sqrt(pi);
c=0.5;

%setting impulses activation time
t0=4;
t20=8;

%Activate impulse on D(t)
b=0;
%Activate impulse on u(t)
d=0;

D0=0;

tf1=15; %time for the first element in g
tf2=10; %time for the rest of elements in g

w=0;
g=[0 1];
clrs={'k' ':b' '-.r'};

for j=1:length(g)
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
set(fig,'color','w');
set(fig, 'Position',  [100, 100, 800, 400])
    for i=1:length(a)
         if j==1
            t=0:0.001:tf1;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a(i),b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);                        
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a(i),b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0D0(w,t,t0,a(i),b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);                        
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0);
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                    end                   
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0,a(i),b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0,a(i),b,t0,t)+Q2.*(D(D0,a(i),b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;                           
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0,a(i),b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));                  
            elseif g(j)==2*w
                S=SigmaCD(w,t,t0,a(i),b,d11,d12,d22,D0);
                Xm=meanscd(w,t,t20,c,d,x0,v0);
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0,a(i),b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0,a(i),b,t0,t)+Q2.*(D(D0,a(i),b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;        
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0,a(i),b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));  
            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0);
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);                    
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0);
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0,a(i),b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0,a(i),b,t0,t)+Q2.*(D(D0,a(i),b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;         
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0,a(i),b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));            
            end
         else
            t=0:0.001:tf2;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0==0
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
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0);
                        Xm=meansg0w0(t,t20,c,d,x0,v0);                        
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                    end
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0,a(i),b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0,a(i),b,t0,t)+Q2.*(D(D0,a(i),b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;           
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*S(1,:).^2);
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0,a(i),b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));                  
            elseif g(j)==2*w
                S=SigmaCD(w,t,t0,a(i),b,d11,d12,d22,D0);
                Xm=meanscd(w,t,t20,c,d,x0,v0);
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0,a(i),b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0,a(i),b,t0,t)+Q2.*(D(D0,a(i),b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;            
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*S(1,:).^2);
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0,a(i),b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));  
            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0);
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);                    
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0);
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                end
                detS=(S(1,:).*S(4,:)-S(2,:).^2);
                deltaxv=-2*g(j)*detS+2*D(D0,a(i),b,t0,t).*S(1,:);
                dvdt=(gradient(Xm(2,:)) ./ gradient(t(:)'));
                Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0,a(i),b,t0,t)+Q2.*(D(D0,a(i),b,t0,t).^2));
                Exv=((S(2,:).*dvdt).^2)./(detS.*S(4,:))+Q;           
                Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*S(1,:).^2);
                Tvx=S(2,:)./S(1,:);                
                Txv=-(w^2)*S(2,:)./S(4,:)-D(D0,a(i),b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));            
            end
        end
        if d==0
            if D0==0 && w==1 && g(j)==0 && b==0
                subplot(2,2,1)
                    yyaxis left
                    indx=find(abs(real(Exv))<1e10);
                    plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t(indx),D(D0,a(i),b,t0,t(indx)),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$\mathcal{E}_{x\to v}$','$D(t)$');
                    set(leg1,'Interpreter','latex');                              
                 subplot(2,2,2)
                    yyaxis left
                    indx=find(abs(real(Evx))<1e10);
                    plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t(indx),D(D0,a(i),b,t0,t(indx)),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$\mathcal{E}_{v\to x}$','$D(t)$');
                    set(leg1,'Interpreter','latex');
                 subplot(2,2,3)
                    yyaxis left
                    plot(t,real(Txv),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$T_{x\to v}$','$D(t)$');
                    set(leg1,'Interpreter','latex');
                 subplot(2,2,4)
                    yyaxis left
                    plot(t,real(Tvx),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$T_{v\to x}$','$D(t)$');
                    set(leg1,'Interpreter','latex');                      
            else
                subplot(2,3,1)
                    yyaxis left
                    indx=find(abs(real(Exv))<1e10);
                    plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t(indx),D(D0,a(i),b,t0,t(indx)),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$\mathcal{E}_{x\to v}$','$D(t)$');
                    set(leg1,'Interpreter','latex');                              
                 subplot(2,3,2)
                    yyaxis left
                    indx=find(abs(real(Evx))<1e10);
                    plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t(indx),D(D0,a(i),b,t0,t(indx)),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$\mathcal{E}_{v\to x}$','$D(t)$');
                    set(leg1,'Interpreter','latex');
                 subplot(2,3,4)
                    yyaxis left
                    plot(t,real(Txv),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$T_{x\to v}$','$D(t)$');
                    set(leg1,'Interpreter','latex');
                 subplot(2,3,5)
                    yyaxis left
                    plot(t,real(Tvx),clrs{i},'LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                    yyaxis right
                     plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
                    xlabel('$t$','Interpreter','Latex','FontSize', 14)
                    ylabel('$D(t)$','Interpreter','Latex','FontSize', 14)
                    axis('square') 
                    leg1 = legend('$T_{v\to x}$','$D(t)$');
                    set(leg1,'Interpreter','latex');                      
                subplot(2,3,3)
                    indx=find(abs(real(Exv))<1e10);
                    plot(real(Exv(indx)),real(Evx(indx)),'b','LineWidth',1)
                    xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                    ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                    axis('square')
                subplot(2,3,6)
                    plot(real(Txv),real(Tvx),'b','LineWidth',1)
                    xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                    ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                    axis('square')
            end
        else
            if t0==t20
                        subplot(2,3,1)
                            yyaxis left
                            indx=find(abs(real(Exv))<1e13);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right 
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$\mathcal{E}_{x\to v}$','$u(t)$');
                            set(leg1,'Interpreter','latex');                              
                         subplot(2,3,2)
                            yyaxis left
                            indx=find(abs(real(Evx))<1e13);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right           
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$\mathcal{E}_{v\to x}$','$u(t)$');
                            set(leg1,'Interpreter','latex');
                         subplot(2,3,4)
                           yyaxis left
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$T_{x\to v}$','$u(t)$');
                            set(leg1,'Interpreter','latex');
                         subplot(2,3,5)
                            yyaxis left
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$T_{v\to x}$','$u(t)$');
                            set(leg1,'Interpreter','latex');                       
                        subplot(2,3,3)
                            indx=find(abs(real(Exv))<1e13);
                            plot(real(Exv(indx)),real(Evx(indx)),'b','LineWidth',1)
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                        subplot(2,3,6)
                            plot(real(Txv),real(Tvx),'b','LineWidth',1)
                            xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square') 
            else
                        subplot(2,3,1)
                            yyaxis left
                            indx=find(abs(real(Exv))<1e13);
                            plot(t(indx),real(Exv(indx)),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right 
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                             hold on
                             plot(t(indx),D(D0,a(i),b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$\mathcal{E}_{x\to v}$','$u(t)$','$D(t)$');
                            set(leg1,'Interpreter','latex');                              
                         subplot(2,3,2)
                            yyaxis left
                            indx=find(abs(real(Evx))<1e13);
                            plot(t(indx),real(Evx(indx)),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right           
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                             hold on
                             plot(t(indx),D(D0,a(i),b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$\mathcal{E}_{v\to x}$','$u(t)$','$D(t)$');
                            set(leg1,'Interpreter','latex');
                         subplot(2,3,4)
                           yyaxis left
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                             hold on
                             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                            leg1 = legend('$T_{x\to v}$','$u(t)$','$D(t)$');
                            set(leg1,'Interpreter','latex');
                         subplot(2,3,5)
                            yyaxis left
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            yyaxis right
                             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
                             hold on
                             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)$','Interpreter','Latex','FontSize', 14)
                            
                            axis('square') 
                            leg1 = legend('$T_{v\to x}$','$u(t)$','$D(t)$');
                            set(leg1,'Interpreter','latex');                       
                        subplot(2,3,3)
                            indx=find(abs(real(Exv))<1e13);
                            plot(real(Exv(indx)),real(Evx(indx)),'b','LineWidth',1)
                            xlabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                        subplot(2,3,6)
                            plot(real(Txv),real(Tvx),'b','LineWidth',1)
                            xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                 
            end           
        end
    end
    str = sprintf('\\gamma=%.2f', g(j));
    sgtitle(str)
end
