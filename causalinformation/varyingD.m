%Causality information rate for harmonically bound particle
%Plots of Txv, Tvx, Exv, Evx for different D0
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

a=0.5; % width impulse in D
c=0.1; % width impulse in u
t0=4;  % time impulse in D
t20=4; %time impulse in u
b=0; %activation of impulse in D
d=0; %activation of impulse in u

tf1=20; %time for the first element in g
tf2=10; %time for the rest of elements in g
h=1e-3;

D0=[1e-2]; %you can add more or less values for D0
w=1;
g=[1];
clrs={'k' ':b' '-.r' ':k' '--b' '--r' '-.k'};

plt=1; %Plot T and E then plt=1, plt=0 otherwise. 
tplot=0; %plot only T then tplot=1, tplot=0 otherwise

for j=1:length(g)
    fig = figure;
    left_color = [0 0 0];
    right_color = [0 0 1];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    set(fig,'color','w');
    if plt==1
        set(fig, 'Position',  [100, 100, 1400, 250])
    else
        set(fig, 'Position',  [100, 100, 1400, 250])
    end

    for i=1:length(D0)
         if j==1
            t=0:h:tf1;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0(i)==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                        
                    end
                end
                    detS=(S(1,:).*S(4,:)-S(2,:).^2);
                    deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                    dvdt=dMdt(g(j),w,t,t20,c,d,x0,v0);
                    Evs=(S(1,:)./detS).*(dvdt(2,:).^2)+(1./(2*detS.^2)).*(2*detS.*((-g(j)*S(2,:)-(w^2).*S(1,:)).^2)+deltaxv.^2);
                    Exs=(S(4,:)./detS).*(Xm(2,:).^2)+(S(4,:).^2)./detS;
                    Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                    Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                    Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                    Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                    Exv=((S(2,:).*dvdt(2,:)).^2)./(detS.*S(4,:))+Q;                           
                    Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                    Gvx=sqrt(S(4,:).*(Xm(2,:).^2)./detS+(S(4,:).^2)./detS)-sqrt((Xm(2,:).^2)./S(1,:)+2*(S(2,:)./S(1,:)).^2);
                    Gxv=sqrt(Evs)-sqrt(velocityv);
                    Tvx=S(2,:)./S(1,:);                
                    Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));                  

            elseif g(j)==2*w
                    S=SigmaCD(w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=meanscd(w,t,t20,c,d,x0,v0);
                    if D0==0
                        if b==0
                            velocity=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,'g');
                            velocityx=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,'x');
                            velocityv=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,'v');
                            %velocity(isnan(velocity))=0;                             
                        else
                            velocity=EpsilonCDD0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'g');
                            velocityx=EpsilonCDD0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'x');
                            velocityv=EpsilonCDD0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'v');
                            %velocity(isnan(velocity))=0; 
                        end
                    else
                        velocity=EpsilonCD(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=EpsilonCD(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=EpsilonCD(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;
                    end
                    detS=(S(1,:).*S(4,:)-S(2,:).^2);
                    deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                    dvdt=dMdt(g(j),w,t,t20,c,d,x0,v0);
                    Evs=(S(1,:)./detS).*(dvdt(2,:).^2)+(1./(2*detS.^2)).*(2*detS.*((-g(j)*S(2,:)-(w^2).*S(1,:)).^2)+deltaxv.^2);
                    Exs=(S(4,:)./detS).*(Xm(2,:).^2)+(S(4,:).^2)./detS;
                    Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                    Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                    Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                    Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                    Exv=((S(2,:).*dvdt(2,:)).^2)./(detS.*S(4,:))+Q;        
                    Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                    Gvx=sqrt(S(4,:).*(Xm(2,:).^2)./detS+(S(4,:).^2)./detS)-sqrt((Xm(2,:).^2)./S(1,:)+2*(S(2,:)./S(1,:)).^2);
                    Gxv=sqrt(Evs)-sqrt(velocityv);                    
                    Tvx=S(2,:)./S(1,:);                
                    Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));  

            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);
                    velocity=Epsilonw0(g(j),a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'g');
                    velocityx=Epsilonw0(g(j),a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'x');
                    velocityv=Epsilonw0(g(j),a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'v');
                    %velocity(isnan(velocity))=0;                     
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                    velocity=Epsilon(g(j),w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                    velocityx=Epsilon(g(j),w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                    velocityv=Epsilon(g(j),w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                    %velocity(isnan(velocity))=0;                    
                end
                    detS=(S(1,:).*S(4,:)-S(2,:).^2);
                    deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                    dvdt=dMdt(g(j),w,t,t20,c,d,x0,v0);
                    Evs=(S(1,:)./detS).*(dvdt(2,:).^2)+(1./(2*detS.^2)).*(2*detS.*((-g(j)*S(2,:)-(w^2).*S(1,:)).^2)+deltaxv.^2);
                    Exs=(S(4,:)./detS).*(Xm(2,:).^2)+(S(4,:).^2)./detS;
                    Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                    Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                    Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                    Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                    Exv=((S(2,:).*dvdt(2,:)).^2)./(detS.*S(4,:))+Q;         
                    Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                    Gvx=sqrt(S(4,:).*(Xm(2,:).^2)./detS+(S(4,:).^2)./detS)-sqrt((Xm(2,:).^2)./S(1,:)+2*(S(2,:)./S(1,:)).^2);
                    Gxv=sqrt(Evs)-sqrt(velocityv);                    
                    Tvx=S(2,:)./S(1,:);                
                    Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));            

            end
        else
            t=0:h:tf2;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0(i)==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilonw0g0(a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=Epsilong0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;                        
                    end
                end
                    detS=(S(1,:).*S(4,:)-S(2,:).^2);
                    deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                    dvdt=dMdt(g(j),w,t,t20,c,d,x0,v0);
                    Evs=(S(1,:)./detS).*(dvdt(2,:).^2)+(1./(2*detS.^2)).*(2*detS.*((-g(j)*S(2,:)-(w^2).*S(1,:)).^2)+deltaxv.^2);
                    Exs=(S(4,:)./detS).*(Xm(2,:).^2)+(S(4,:).^2)./detS;
                    Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                    Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                    Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                    Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                    Exv=((S(2,:).*dvdt(2,:)).^2)./(detS.*S(4,:))+Q;                           
                    Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                    Gvx=sqrt(S(4,:).*(Xm(2,:).^2)./detS+(S(4,:).^2)./detS)-sqrt((Xm(2,:).^2)./S(1,:)+2*(S(2,:)./S(1,:)).^2);
                    Gxv=sqrt(Evs)-sqrt(velocityv);                    
                    Tvx=S(2,:)./S(1,:);                
                    Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));                  

            elseif g(j)==2*w
                    S=SigmaCD(w,t,t0,a(i),b,d11,d12,d22,D0);
                    Xm=meanscd(w,t,t20,c,d,x0,v0);
                    if D0==0
                        if b==0
                            velocity=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,'g');
                            velocityx=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,'x');
                            velocityv=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,'v');
                            %velocity(isnan(velocity))=0;                             
                        else
                            velocity=EpsilonCDD0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'g');
                            velocityx=EpsilonCDD0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'x');
                            velocityv=EpsilonCDD0(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'v');
                            %velocity(isnan(velocity))=0; 
                        end
                    else
                        velocity=EpsilonCD(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                        velocityx=EpsilonCD(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                        velocityv=EpsilonCD(w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                        %velocity(isnan(velocity))=0;
                    end
                    detS=(S(1,:).*S(4,:)-S(2,:).^2);
                    deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                    dvdt=dMdt(g(j),w,t,t20,c,d,x0,v0);
                    Evs=(S(1,:)./detS).*(dvdt(2,:).^2)+(1./(2*detS.^2)).*(2*detS.*((-g(j)*S(2,:)-(w^2).*S(1,:)).^2)+deltaxv.^2);
                    Exs=(S(4,:)./detS).*(Xm(2,:).^2)+(S(4,:).^2)./detS;
                    Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                    Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                    Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                    Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                    Exv=((S(2,:).*dvdt(2,:)).^2)./(detS.*S(4,:))+Q;        
                    Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                    Gvx=sqrt(S(4,:).*(Xm(2,:).^2)./detS+(S(4,:).^2)./detS)-sqrt((Xm(2,:).^2)./S(1,:)+2*(S(2,:)./S(1,:)).^2);
                    Gxv=sqrt(Evs)-sqrt(velocityv);                    
                    Tvx=S(2,:)./S(1,:);                
                    Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));  

            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);
                    velocity=Epsilonw0(g(j),a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'g');
                    velocityx=Epsilonw0(g(j),a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'x');
                    velocityv=Epsilonw0(g(j),a,b,c,d,t,t0,t20,v0,d11,d12,d22,D0(i),'v');
                    %velocity(isnan(velocity))=0;                     
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                    velocity=Epsilon(g(j),w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'g');
                    velocityx=Epsilon(g(j),w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'x');
                    velocityv=Epsilon(g(j),w,a,b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0(i),'v');
                    %velocity(isnan(velocity))=0;                    
                end
                    detS=(S(1,:).*S(4,:)-S(2,:).^2);
                    deltaxv=-2*g(j)*detS+2*D(D0(i),a,b,t0,t).*S(1,:);
                    dvdt=dMdt(g(j),w,t,t20,c,d,x0,v0);
                    Evs=(S(1,:)./detS).*(dvdt(2,:).^2)+(1./(2*detS.^2)).*(2*detS.*((-g(j)*S(2,:)-(w^2).*S(1,:)).^2)+deltaxv.^2);
                    Exs=(S(4,:)./detS).*(Xm(2,:).^2)+(S(4,:).^2)./detS;
                    Q0=(g(j)*S(2,:).*S(4,:)+(2*S(2,:).^2-S(1,:).*S(4,:))*(w^2)).^2+2*(w^4)*(S(2,:).^2).*detS;
                    Q1=4*S(2,:).*(-g(j)*S(2,:).*S(4,:)+detS*(w^2));
                    Q2=(2*S(2,:).^4)./detS+4*S(2,:).^2;                
                    Q=(1./(detS.*(S(4,:).^2))).*(Q0+Q1.*D(D0(i),a,b,t0,t)+Q2.*(D(D0(i),a,b,t0,t).^2));
                    Exv=((S(2,:).*dvdt(2,:)).^2)./(detS.*S(4,:))+Q;         
                    Evx=((S(2,:).*Xm(2,:)).^2)./(detS.*S(1,:))+(detS.^2+S(2,:).^4)./(detS.*(S(1,:).^2));
                    Gvx=sqrt(S(4,:).*(Xm(2,:).^2)./detS+(S(4,:).^2)./detS)-sqrt((Xm(2,:).^2)./S(1,:)+2*(S(2,:)./S(1,:)).^2);
                    Gxv=sqrt(Evs)-sqrt(velocityv);                    
                    Tvx=S(2,:)./S(1,:);                
                    Txv=-(w^2)*S(2,:)./S(4,:)-D(D0(i),a,b,t0,t).*(S(2,:).^2)./((S(1,:).*S(4,:)-S(2,:).^2).*S(4,:));            

            end
         end
         if plt==1
                 if b==0 && d==0            
                      subplot(1,4,1)
                            plot(t,real(Tvx),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                  
                       subplot(1,4,2)
                            plot(t,real(Txv),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
%                        subplot(2,3,3)
%                             plot(real(Txv),real(Tvx),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
%                             ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
%                             axis('square')               
                       subplot(1,4,3)
                            indx=find(abs(real(Gvx))<4e2);
                            plot(t(indx),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(1,4,4)
                            indx=find(abs(real(Gxv))<4e2);
                            plot(t(indx),real(Gxv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
%                        subplot(2,3,6)
%                             indx=find(abs(real(Gxv))<4e2);
%                             plot(real(Gxv(indx)),real(Gvx(indx)),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 12)
%                             ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
%                             axis('square') 
%                        subplot(1,5,5)
%                             indx=find(abs(real(Gxv))<4e2);
%                             plot(real(Gxv./sqrt(Evs)),real(Gvx./sqrt(Exs)),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$\Gamma_{x\to v}/\Gamma_v^*$','Interpreter','Latex','FontSize', 15)
%                             ylabel('$\Gamma_{v\to x}/\Gamma_x^*$','Interpreter','Latex','FontSize', 15)
%                             axis('square')
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
                            indx=find(abs(real(Gvx))<4e2);
                            plot(t(indx),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(2,3,5)
                            yyaxis left
                            indx=find(abs(real(Gxv))<4e2);
                            plot(t(indx),real(Gxv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                             yyaxis right
                             plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
                        subplot(2,3,6)
                            indx=find(abs(real(Gxv))<4e2);
                            plot(real(Gxv(indx)),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')                      
                 elseif b==0 && d==1
                      subplot(1,4,1)
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
                       subplot(1,4,2)
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
%                        subplot(1,5,3)
%                             plot(real(Txv),real(Tvx),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
%                             ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
%                             axis('square')                      
                       subplot(1,4,3)
                            yyaxis left
                            indx=find(abs(real(Gvx))<1e13);
                            plot(t(indx),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square')                     
                       subplot(1,4,4)
                            yyaxis left
                            indx=find(abs(real(Gxv))<1e13);
                            plot(t(indx),real(Gxv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                             yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
                            axis('square') 
%                        subplot(1,5,5)
%                             indx=find(abs(real(Gxv))<4e2);
%                             plot(real(Gxv./sqrt(Evs)),real(Gvx./sqrt(Exs)),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$\Gamma_{x\to v}/\Gamma_v^*$','Interpreter','Latex','FontSize', 15)
%                             ylabel('$\Gamma_{v\to x}/\Gamma_x^*$','Interpreter','Latex','FontSize', 15)
%                             axis('square')              
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
                            indx=find(abs(real(Gvx))<1e13);
                            plot(t(indx),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
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
                            indx=find(abs(real(Gxv))<1e13);
                            plot(t(indx),real(Gxv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                            yyaxis right
                             plot(t(indx),D(0,c,d,t20,t(indx)),'r:','LineWidth',1)
                             hold on
                              plot(t(indx),D(0,a,b,t0,t(indx)),'b:','LineWidth',1)
                            xlabel('$t$','Interpreter','Latex','FontSize', 14)
                            ylabel('$u(t),D(t)-D_0$','Interpreter','Latex','FontSize', 14)
                            axis('square')  
                        subplot(2,3,6)
                            indx=find(abs(real(Gxv))<1e13);
                            plot(real(Gxv(indx)),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')             
                 end
         elseif tplot==1
                       subplot(2,1,1)
                            indx=find(abs(real(Tvx))<4e2);
                            plot(t(indx),real(Tvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(2,1,2)
                            indx=find(abs(real(Txv))<4e2);
                            plot(t(indx),real(Txv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
%                        subplot(1,3,3)
%                             indx=find(abs(real(Txv))<4e2);
%                             plot(real(Txv(indx)),real(Tvx(indx)),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 12)
%                             ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 15)
%                             axis('square')                
         else
                       subplot(2,1,1)
                            indx=find(abs(real(Gvx))<4e2);
                            plot(t(indx),real(Gvx(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
                       subplot(2,1,2)
                            indx=find(abs(real(Gxv))<4e2);
                            plot(t(indx),real(Gxv(indx)),clrs{i},'LineWidth',1)
                            hold on
                            xlabel('$t$','Interpreter','Latex','FontSize', 12)
                            ylabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 15)
                            axis('square')
%                        subplot(1,3,3)
%                             indx=find(abs(real(Gxv))<4e2);
%                             plot(real(Gxv(indx)),real(Gvx(indx)),clrs{i},'LineWidth',1)
%                             hold on
%                             xlabel('$\Gamma_{x\to v}$','Interpreter','Latex','FontSize', 12)
%                             ylabel('$\Gamma_{v\to x}$','Interpreter','Latex','FontSize', 15)
%                             axis('square')                 
         end
    end
    %str = sprintf('\\gamma=%.2f, \\omega=%0.2f, D_0=%0.2f', g(j),w,D0(i));
    %sgtitle(str)
end
%AddLabels(D0,clrs)