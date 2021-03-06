%Causal Information of the Kramers equation
%Plots for Ex*, Ev*, Ex, Ev, E
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
c=0.1;
t0=4;
t20=4;

%Activate impulse on D(t)
b=0;
%Activate impulse on u(t)
d=1;
D0=1e-2;

tf1=20; %time for the first element in g
tf2=20; %time for the rest of elements in g
h=1e-3; %timestep

w=1;
g=[1];
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
            t=0:h:tf1;
            x1=zeros(1,length(t));
            x2=zeros(1,length(t));
            if g(j)==0
                if D0(i)==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
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
                            velocity=EpsilonCDD0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'g');
                            velocityx=EpsilonCDD0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'x');
                            velocityv=EpsilonCDD0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'v');
                            %velocity(isnan(velocity))=0; 
                        end
                    else
                        velocity=EpsilonCD(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=EpsilonCD(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=EpsilonCD(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
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
                    velocity=Epsilonw0(g(j),a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'g');
                    velocityx=Epsilonw0(g(j),a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'x');
                    velocityv=Epsilonw0(g(j),a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'v');
                    %velocity(isnan(velocity))=0;                     
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                    velocity=Epsilon(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                    velocityx=Epsilon(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                    velocityv=Epsilon(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
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
                        velocity=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                        velocity=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilonw0g0(a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'v');
                        %velocity(isnan(velocity))=0;                         
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                        velocity=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=Epsilong0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
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
                            velocity=EpsilonCDD0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'g');
                            velocityx=EpsilonCDD0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'x');
                            velocityv=EpsilonCDD0(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,'v');
                            %velocity(isnan(velocity))=0; 
                        end
                    else
                        velocity=EpsilonCD(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                        velocityx=EpsilonCD(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                        velocityv=EpsilonCD(w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
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
                    velocity=Epsilonw0(g(j),a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'g');
                    velocityx=Epsilonw0(g(j),a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'x');
                    velocityv=Epsilonw0(g(j),a(i),b,c,d,t,t0,t20,v0,d11,d12,d22,D0,'v');
                    %velocity(isnan(velocity))=0;                     
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                    velocity=Epsilon(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'g');
                    velocityx=Epsilon(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'x');
                    velocityv=Epsilon(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0,'v');
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
         if d==1
        subplot(2,3,1)
           yyaxis left
            plot(t,sqrt(real(velocityx)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{x}$','Interpreter','Latex','FontSize', 15)
            axis('square')

           yyaxis right
             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 14)
            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
            axis('square')     
        subplot(2,3,2)
        yyaxis left
            plot(t,sqrt(real(velocityv)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{v}$','Interpreter','Latex','FontSize', 15)
            axis('square')
           yyaxis right
             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 14)
            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
            axis('square')             
%         subplot(2,4,3)
%             plot(t,real(Exv),clrs{i},'LineWidth',1)
%             hold on
%             xlabel('$t$','Interpreter','Latex','FontSize', 12)
%             ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
%             axis('square')
        subplot(2,3,3)
        yyaxis left
            plot(t,real(Gxv./sqrt(Evs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{x\to v}/\Gamma_v^*$','Interpreter','Latex','FontSize', 15)
            axis('square')
           yyaxis right
             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 14)
            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
            axis('square')             
        subplot(2,3,4)
        yyaxis left
             plot(t,sqrt(real(Exs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_x^*$','Interpreter','Latex','FontSize', 15)
            axis('square') 
           yyaxis right
             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 14)
            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
            axis('square')             
        subplot(2,3,5)
        yyaxis left
             plot(t,sqrt(real(Evs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_v^*$','Interpreter','Latex','FontSize', 15)
            axis('square')
           yyaxis right
             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 14)
            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
            axis('square')             
%         subplot(2,4,7)
%              plot(t,real(Evx),clrs{i},'LineWidth',1)
%             hold on
%             xlabel('$t$','Interpreter','Latex','FontSize', 12)
%             ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
%             axis('square')       
        subplot(2,3,6)
        yyaxis left
             plot(t,real(Gvx./sqrt(Exs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{v\to x}/\Gamma_x^*$','Interpreter','Latex','FontSize', 15)
            axis('square')
           yyaxis right
             plot(t,D(0,c,d,t20,t),'r:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 14)
            ylabel('$u(t)$','Interpreter','Latex','FontSize', 14)
            axis('square')             
         else
        subplot(2,3,1)
            plot(t,sqrt(real(velocityx)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{x}$','Interpreter','Latex','FontSize', 15)
            axis('square')
        subplot(2,3,2)
            plot(t,sqrt(real(velocityv)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{v}$','Interpreter','Latex','FontSize', 15)
            axis('square')
%         subplot(2,4,3)
%             plot(t,real(Exv),clrs{i},'LineWidth',1)
%             hold on
%             xlabel('$t$','Interpreter','Latex','FontSize', 12)
%             ylabel('$\mathcal{E}_{x\to v}$','Interpreter','Latex','FontSize', 15)
%             axis('square')
        subplot(2,3,3)
            plot(t,real(Gxv./sqrt(Evs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{x\to v}/\Gamma_v^*$','Interpreter','Latex','FontSize', 15)
            axis('square')
        subplot(2,3,4)
             plot(t,sqrt(real(Exs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_x^*$','Interpreter','Latex','FontSize', 15)
            axis('square')       
        subplot(2,3,5)
             plot(t,sqrt(real(Evs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_v^*$','Interpreter','Latex','FontSize', 15)
            axis('square')             
%         subplot(2,4,7)
%              plot(t,real(Evx),clrs{i},'LineWidth',1)
%             hold on
%             xlabel('$t$','Interpreter','Latex','FontSize', 12)
%             ylabel('$\mathcal{E}_{v\to x}$','Interpreter','Latex','FontSize', 15)
%             axis('square')       
        subplot(2,3,6)
             plot(t,real(Gvx./sqrt(Exs)),clrs{i},'LineWidth',1)
            hold on
            xlabel('$t$','Interpreter','Latex','FontSize', 12)
            ylabel('$\Gamma_{v\to x}/\Gamma_x^*$','Interpreter','Latex','FontSize', 15)
            axis('square')   
         end
    end
   % str = sprintf('\\gamma=%.2f, \\omega=%0.2f, D_0=%0.2f', g(j),w,D0);
   % sgtitle(str)
end
