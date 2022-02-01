%Information length for harmonically bound particle
%Snapshots and values of \Sigma
%Adrian J Guel C
%March 2021

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
c=0.5; % width impulse in u
t0=4;  % time impulse in D
t20=4; %time impulse in u
b=0; %activation of impulse in D
d=0; %activation of impulse in u

D0=0;

tf1=20; %time for the first element in g
tf2=10; %time for the rest of elements in g

w=0;
g=[0];
clrs={'k' ':b' '-.r'};


for j=1:length(g)
    if g(j)==0 && w==0
        x1 = -5:0.01:5;
        x2 = -5:0.01:5;
        [X1,X2] = meshgrid(x1,x2);
        X = [X1(:) X2(:)];    
    else
        x1 = -5:0.01:5;
        x2 = -5:0.01:5;
        [X1,X2] = meshgrid(x1,x2);
        X = [X1(:) X2(:)];
    end
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
set(fig,'color','w');
set(fig, 'Position',  [100, 100, 800, 400])
    for i=1:length(a)
           t=0:0.01:tf1;
           subplot(3,2,[1 3 5])
           if g(j)==0
                if D0(i)==0
                    if w==1
                        S=Sigmag0D0w1(t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                         for r=1:50:length(t)
                                y = mvnpdf(X,meansg0(w,t(r),t20,c,d,x0,v0)',nearestSPD(Sigmag0D0w1plot(t(r),t0,a,b,d11,d12,d22)));
                                y = reshape(y,length(x2),length(x1));
                                contour(real(x1),real(x2),real(y),1)
                                hold on              
                         end                          
                    elseif w==0
                        S=Sigmag0D0w0(t,t0,a,b,d11,d12,d22); 
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                         for r=1:50:length(t)
                                y = mvnpdf(X,meansg0w0(t(r),t20,c,d,x0,v0)',nearestSPD(Sigmag0D0w0plot(t(r),t0,a,b,d11,d12,d22)));
                                y = reshape(y,length(x2),length(x1));
                                contour(real(x1),real(x2),real(y),1)
                                hold on              
                         end                          
                    else
                        S=Sigmag0D0(w,t,t0,a,b,d11,d12,d22);
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                         for r=1:50:length(t)
                                y = mvnpdf(X,meansg0(w,t(r),t20,c,d,x0,v0)',nearestSPD(Sigmag0D0plot(w,t(r),t0,a,b,d11,d12,d22)));
                                y = reshape(y,length(x2),length(x1));
                                contour(real(x1),real(x2),real(y),1)
                                hold on              
                         end                          
                    end
                else
                    if w==0
                        S=Sigmag0w0(t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0w0(t,t20,c,d,x0,v0);
                         for r=1:50:length(t)
                                y = mvnpdf(X,meansg0w0(t(r),t20,c,d,x0,v0)',nearestSPD(Sigmag0w0plot(t(r),t0,a,b,d11,d12,d22,D0(i))));
                                y = reshape(y,length(x2),length(x1));
                                contour(real(x1),real(x2),real(y),1)
                                hold on              
                         end                          
                    else
                        S=Sigmag0(w,t,t0,a,b,d11,d12,d22,D0(i));
                        Xm=meansg0(w,t,t20,c,d,x0,v0);
                         for r=1:50:length(t)
                                y = mvnpdf(X,meansg0(w,t(r),t20,c,d,x0,v0)',nearestSPD(Sigmag0plot(w,t(r),t0,a,b,d11,d12,d22,D0(i))));
                                y = reshape(y,length(x2),length(x1));
                                contour(real(x1),real(x2),real(y),1)
                                hold on              
                         end                          
                    end
                end
            elseif g(j)==2*w
                S=SigmaCD(w,t,t0,a,b,d11,d12,d22,D0(i));
                Xm=meanscd(w,t,t20,c,d,x0,v0);
                 for r=1:50:length(t)
                        y = mvnpdf(X,meanscd(w,t(r),t20,c,d,x0,v0)',nearestSPD(SigmaCDplot(w,t(r),t0,a,b,d11,d12,d22,D0(i))));
                        y = reshape(y,length(x2),length(x1));
                        contour(real(x1),real(x2),real(y),1)
                        hold on              
                 end                  
            else
                if w==0
                    S=Sigmaw0(g(j),t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=meansw0(g(j),t,t20,c,d,x0,v0);
                     for r=1:50:length(t)
                            y = mvnpdf(X,meansw0(g(j),t(r),t20,c,d,x0,v0)',nearestSPD(Sigmaw0plot(g(j),t(r),t0,a,b,d11,d12,d22,D0(i))));
                            y = reshape(y,length(x2),length(x1));
                            contour(real(x1),real(x2),real(y),1)
                            hold on              
                     end                      
                else
                    S=sigmaDv2(g(j),w,t,t0,a,b,d11,d12,d22,D0(i));
                    Xm=means(g(j),w,t,t20,c,d,x0,v0);
                     for r=1:50:length(t)
                            y = mvnpdf(X,means(g(j),w,t(r),t20,c,d,x0,v0)',nearestSPD(sigmaDv2plot(g(j),w,t(r),t0,a,b,d11,d12,d22,D0(i))));
                            y = reshape(y,length(x2),length(x1));
                            contour(real(x1),real(x2),real(y),1)
                            hold on              
                     end                      
                end
            end                        
            hold on
            %aux=means(g(j),w,t,t20,c,d,x0,v0);
            plot(real(Xm(1,:)),real(Xm(2,:)),'k:','LineWidth',1)
            axis('square')
            grid on
            xlabel('$x_1$','Interpreter','Latex','FontSize', 15)
            ylabel('$x_2$','Interpreter','Latex','FontSize', 15)            
            leg1 = legend('$\mathbf{x}\sim N(\langle\mathbf{x}\rangle,\Sigma)$');
            set(leg1,'Interpreter','latex');
            subplot(3,2,2)
            plot(t,real(S(1,:)),'k','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 15)
            ylabel('$\Sigma_{xx}$','Interpreter','Latex','FontSize', 15)            
            subplot(3,2,4)
             plot(t,real(S(4,:)),'k','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 15)
            ylabel('$\Sigma_{vv}$','Interpreter','Latex','FontSize', 15)            
            subplot(3,2,6)
             plot(t,real(S(2,:)),'k','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 15)
            ylabel('$\Sigma_{xv}$','Interpreter','Latex','FontSize', 15)
             
    end
   str = sprintf('\\gamma=%.2f, \\omega=%0.2f and D_0=%0.2f', g(j),w,D0);
   sgtitle(str)
end
