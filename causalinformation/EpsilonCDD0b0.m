function ECD=EpsilonCDD0b0(w,c,d,t,t20,x0,v0,d11,d12,d22,type)
if type=='g'
ECD=(1/2).*((((-1).*d12.^2.*exp(1).^((-4).*t.*w)+d11.*d22.*exp(1).^((-4).* ...
  t.*w)).^(-1).*(exp(1).^((-2).*t.*w).*t.*(d22+d12.*w)+exp(1).^((-2).*t.* ...
  w).*(d12+d11.*w).*(1+t.*w)+exp(1).^((-2).*t.*w).*w.*(d11+d12.*t+d11.*t.* ...
  w)+(-2).*exp(1).^((-2).*t.*w).*w.*(1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp( ...
  1).^((-2).*t.*w).*(d12+d22.*t+d12.*t.*w)+(-2).*exp(1).^((-2).*t.*w).*t.* ...
  w.*(d12+d22.*t+d12.*t.*w)).*(exp(1).^((-2).*t.*w).*t.*w.^2.*((-1).*d12+ ...
  d12.*t.*w+d11.*t.*w.^2)+exp(1).^((-2).*t.*w).*((-1)+t.*w).*((-1).*d22+ ...
  d22.*t.*w+d12.*t.*w.^2))+((-1).*d12.^2.*exp(1).^((-4).*t.*w)+d11.*d22.* ...
  exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-1).*t.*w).*w.*((-1).*d11.*exp( ...
  1).^((-1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))+( ...
  -1).*exp(1).^((-1).*t.*w).*w.*(1+t.*w).*((-1).*d11.*exp(1).^((-1).*t.*w) ...
  .*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))+exp(1).^((-1).*t.* ...
  w).*((-1).*d12.*exp(1).^((-1).*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w) ...
  .*(1+(-1).*t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.*w.*((-1).*d12.*exp(1) ...
  .^((-1).*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))+exp( ...
  1).^((-1).*t.*w).*(1+t.*w).*((-1).*d12.*exp(1).^((-1).*t.*w).*w+(-1).* ...
  d11.*exp(1).^((-1).*t.*w).*w.^2+d11.*exp(1).^((-1).*t.*w).*t.*w.^3+(-1) ...
  .*d12.*exp(1).^((-1).*t.*w).*w.*(1+(-1).*t.*w))+exp(1).^((-1).*t.*w).* ...
  t.*((-1).*d22.*exp(1).^((-1).*t.*w).*w+(-1).*d12.*exp(1).^((-1).*t.*w).* ...
  w.^2+d12.*exp(1).^((-1).*t.*w).*t.*w.^3+(-1).*d22.*exp(1).^((-1).*t.*w) ...
  .*w.*(1+(-1).*t.*w))).*(exp(1).^((-1).*t.*w).*t.*w.^2.*(d12.*exp(1).^(( ...
  -1).*t.*w).*t+d11.*exp(1).^((-1).*t.*w).*(1+t.*w))+exp(1).^((-1).*t.*w) ...
  .*((-1)+t.*w).*(d22.*exp(1).^((-1).*t.*w).*t+d12.*exp(1).^((-1).*t.*w).* ...
  (1+t.*w)))).^2+2.*(((-1).*d12.^2.*exp(1).^((-4).*t.*w)+d11.*d22.*exp(1) ...
  .^((-4).*t.*w)).^(-1).*(exp(1).^((-2).*t.*w).*t.*(d22+d12.*w)+exp(1).^(( ...
  -2).*t.*w).*(d12+d11.*w).*(1+t.*w)+exp(1).^((-2).*t.*w).*w.*(d11+d12.*t+ ...
  d11.*t.*w)+(-2).*exp(1).^((-2).*t.*w).*w.*(1+t.*w).*(d11+d12.*t+d11.*t.* ...
  w)+exp(1).^((-2).*t.*w).*(d12+d22.*t+d12.*t.*w)+(-2).*exp(1).^((-2).*t.* ...
  w).*t.*w.*(d12+d22.*t+d12.*t.*w)).*((-1).*exp(1).^((-1).*t.*w).*(1+t.*w) ...
  .*((-1).*d11.*exp(1).^((-1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*( ...
  1+(-1).*t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.*((-1).*d12.*exp(1).^((-1) ...
  .*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w)))+((-1).* ...
  d12.^2.*exp(1).^((-4).*t.*w)+d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*( ...
  exp(1).^((-2).*t.*w).*(1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp(1).^((-2).* ...
  t.*w).*t.*(d12+d22.*t+d12.*t.*w)).*(exp(1).^((-1).*t.*w).*w.*((-1).* ...
  d11.*exp(1).^((-1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*(1+(-1).* ...
  t.*w))+(-1).*exp(1).^((-1).*t.*w).*w.*(1+t.*w).*((-1).*d11.*exp(1).^(( ...
  -1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))+exp(1).^( ...
  (-1).*t.*w).*((-1).*d12.*exp(1).^((-1).*t.*w).*t.*w.^2+d22.*exp(1).^(( ...
  -1).*t.*w).*(1+(-1).*t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.*w.*((-1).* ...
  d12.*exp(1).^((-1).*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).*(1+(-1).* ...
  t.*w))+exp(1).^((-1).*t.*w).*(1+t.*w).*((-1).*d12.*exp(1).^((-1).*t.*w) ...
  .*w+(-1).*d11.*exp(1).^((-1).*t.*w).*w.^2+d11.*exp(1).^((-1).*t.*w).*t.* ...
  w.^3+(-1).*d12.*exp(1).^((-1).*t.*w).*w.*(1+(-1).*t.*w))+exp(1).^((-1).* ...
  t.*w).*t.*((-1).*d22.*exp(1).^((-1).*t.*w).*w+(-1).*d12.*exp(1).^((-1).* ...
  t.*w).*w.^2+d12.*exp(1).^((-1).*t.*w).*t.*w.^3+(-1).*d22.*exp(1).^((-1) ...
  .*t.*w).*w.*(1+(-1).*t.*w)))).*(((-1).*d12.^2.*exp(1).^((-4).*t.*w)+ ...
  d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-2).*t.*w).*t.*w.^2.*( ...
  d12.*w+d11.*w.^2)+exp(1).^((-2).*t.*w).*((-1)+t.*w).*(d22.*w+d12.*w.^2)+ ...
  exp(1).^((-2).*t.*w).*w.^2.*((-1).*d12+d12.*t.*w+d11.*t.*w.^2)+(-2).* ...
  exp(1).^((-2).*t.*w).*t.*w.^3.*((-1).*d12+d12.*t.*w+d11.*t.*w.^2)+exp(1) ...
  .^((-2).*t.*w).*w.*((-1).*d22+d22.*t.*w+d12.*t.*w.^2)+(-2).*exp(1).^(( ...
  -2).*t.*w).*w.*((-1)+t.*w).*((-1).*d22+d22.*t.*w+d12.*t.*w.^2)).*(exp(1) ...
  .^((-1).*t.*w).*t.*w.^2.*(d12.*exp(1).^((-1).*t.*w).*t+d11.*exp(1).^(( ...
  -1).*t.*w).*(1+t.*w))+exp(1).^((-1).*t.*w).*((-1)+t.*w).*(d22.*exp(1).^( ...
  (-1).*t.*w).*t+d12.*exp(1).^((-1).*t.*w).*(1+t.*w)))+((-1).*d12.^2.*exp( ...
  1).^((-4).*t.*w)+d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-2).* ...
  t.*w).*t.*w.^2.*((-1).*d12+d12.*t.*w+d11.*t.*w.^2)+exp(1).^((-2).*t.*w) ...
  .*((-1)+t.*w).*((-1).*d22+d22.*t.*w+d12.*t.*w.^2)).*((-1).*exp(1).^((-1) ...
  .*t.*w).*w.^2.*(d12.*exp(1).^((-1).*t.*w).*t+d11.*exp(1).^((-1).*t.*w).* ...
  (1+t.*w))+exp(1).^((-1).*t.*w).*t.*w.^3.*(d12.*exp(1).^((-1).*t.*w).*t+ ...
  d11.*exp(1).^((-1).*t.*w).*(1+t.*w))+(-1).*exp(1).^((-1).*t.*w).*w.*( ...
  d22.*exp(1).^((-1).*t.*w).*t+d12.*exp(1).^((-1).*t.*w).*(1+t.*w))+exp(1) ...
  .^((-1).*t.*w).*w.*((-1)+t.*w).*(d22.*exp(1).^((-1).*t.*w).*t+d12.*exp( ...
  1).^((-1).*t.*w).*(1+t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.*w.^2.*(d12.* ...
  exp(1).^((-1).*t.*w)+d11.*exp(1).^((-1).*t.*w).*w+(-1).*d12.*exp(1).^(( ...
  -1).*t.*w).*t.*w+(-1).*d11.*exp(1).^((-1).*t.*w).*w.*(1+t.*w))+(-1).* ...
  exp(1).^((-1).*t.*w).*((-1)+t.*w).*(d22.*exp(1).^((-1).*t.*w)+d12.*exp( ...
  1).^((-1).*t.*w).*w+(-1).*d22.*exp(1).^((-1).*t.*w).*t.*w+(-1).*d12.* ...
  exp(1).^((-1).*t.*w).*w.*(1+t.*w))))+(((-1).*d12.^2.*exp(1).^((-4).*t.* ...
  w)+d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-2).*t.*w).*(1+t.* ...
  w).*(d11+d12.*t+d11.*t.*w)+exp(1).^((-2).*t.*w).*t.*(d12+d22.*t+d12.*t.* ...
  w)).*(exp(1).^((-2).*t.*w).*t.*w.^2.*(d12.*w+d11.*w.^2)+exp(1).^((-2).* ...
  t.*w).*((-1)+t.*w).*(d22.*w+d12.*w.^2)+exp(1).^((-2).*t.*w).*w.^2.*((-1) ...
  .*d12+d12.*t.*w+d11.*t.*w.^2)+(-2).*exp(1).^((-2).*t.*w).*t.*w.^3.*((-1) ...
  .*d12+d12.*t.*w+d11.*t.*w.^2)+exp(1).^((-2).*t.*w).*w.*((-1).*d22+d22.* ...
  t.*w+d12.*t.*w.^2)+(-2).*exp(1).^((-2).*t.*w).*w.*((-1)+t.*w).*((-1).* ...
  d22+d22.*t.*w+d12.*t.*w.^2))+((-1).*d12.^2.*exp(1).^((-4).*t.*w)+d11.* ...
  d22.*exp(1).^((-4).*t.*w)).^(-1).*((-1).*exp(1).^((-1).*t.*w).*(1+t.*w) ...
  .*((-1).*d11.*exp(1).^((-1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*( ...
  1+(-1).*t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.*((-1).*d12.*exp(1).^((-1) ...
  .*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))).*((-1).* ...
  exp(1).^((-1).*t.*w).*w.^2.*(d12.*exp(1).^((-1).*t.*w).*t+d11.*exp(1).^( ...
  (-1).*t.*w).*(1+t.*w))+exp(1).^((-1).*t.*w).*t.*w.^3.*(d12.*exp(1).^(( ...
  -1).*t.*w).*t+d11.*exp(1).^((-1).*t.*w).*(1+t.*w))+(-1).*exp(1).^((-1).* ...
  t.*w).*w.*(d22.*exp(1).^((-1).*t.*w).*t+d12.*exp(1).^((-1).*t.*w).*(1+ ...
  t.*w))+exp(1).^((-1).*t.*w).*w.*((-1)+t.*w).*(d22.*exp(1).^((-1).*t.*w) ...
  .*t+d12.*exp(1).^((-1).*t.*w).*(1+t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.* ...
  w.^2.*(d12.*exp(1).^((-1).*t.*w)+d11.*exp(1).^((-1).*t.*w).*w+(-1).* ...
  d12.*exp(1).^((-1).*t.*w).*t.*w+(-1).*d11.*exp(1).^((-1).*t.*w).*w.*(1+ ...
  t.*w))+(-1).*exp(1).^((-1).*t.*w).*((-1)+t.*w).*(d22.*exp(1).^((-1).*t.* ...
  w)+d12.*exp(1).^((-1).*t.*w).*w+(-1).*d22.*exp(1).^((-1).*t.*w).*t.*w+( ...
  -1).*d12.*exp(1).^((-1).*t.*w).*w.*(1+t.*w)))).^2+2.*(((-2).*exp(1).^(( ...
  -1).*t.*w).*v0.*w+exp(1).^((-1).*t.*w).*t.*v0.*w.^2+(-1).*exp(1).^((-1) ...
  .*t.*w).*w.^2.*x0+exp(1).^((-1).*t.*w).*t.*w.^3.*x0+d.*exp(1).^((-1).* ...
  t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+( ...
  -1/2).*c.*w).^2).*pi.^(-1/2).*abs(c).^(-1)+(-1).*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).* ...
  c.*w).^2).*pi.^(-1/2).*t.*w.*abs(c).^(-1)+d.*exp(1).^((-1).*t.*w+t20.*w+ ...
  (1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2) ...
  .*pi.^(-1/2).*t20.*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).* ...
  c.*w).^2).*pi.^(-1/2).*w.^2.*abs(c).^(-1)+(-1/2).*c.^2.*d.*exp(1).^((-1) ...
  .*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t20+(1/2).*c.*w).^2).* ...
  pi.^(-1/2).*w.^2.*abs(c).^(-1)+(-1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).*t+2.*t20+c.^2.*w).^2).* ...
  pi.^(-1/2).*w.*((-1).*w+c.^(-2).*((-2).*t+2.*t20+c.^2.*w)).*abs(c).^(-1) ...
  +(-1).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.*abs(c).^( ...
  -1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.^2.*abs(c).^(-1).*erfz(c.^( ...
  -1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1).*t.* ...
  w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+( ...
  -1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2).*w.^3.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^( ...
  -1).*t20+(-1/2).*c.*w)+(-1).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2).*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.^2.*abs(c).^(-1) ...
  .*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.* ...
  w+(1/4).*c.^2.*w.^2).*t20.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).* ...
  c.*w)+(-1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).* ...
  w.^3.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)).*(((-1).*d12.^2.*exp( ...
  1).^((-4).*t.*w)+d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-1).* ...
  t.*w).*t.*w.^2.*(d12.*exp(1).^((-1).*t.*w).*t+d11.*exp(1).^((-1).*t.*w) ...
  .*(1+t.*w))+exp(1).^((-1).*t.*w).*((-1)+t.*w).*(d22.*exp(1).^((-1).*t.* ...
  w).*t+d12.*exp(1).^((-1).*t.*w).*(1+t.*w))).*(exp(1).^((-1).*t.*w).*v0+( ...
  -1).*exp(1).^((-1).*t.*w).*t.*v0.*w+(-1).*exp(1).^((-1).*t.*w).*t.* ...
  w.^2.*x0+d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1) ...
  .*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t.*abs(c).^(-1)+( ...
  -1).*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+( ...
  -1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t20.*abs(c).^(-1)+( ...
  -1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^( ...
  -1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1) ...
  +(1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*( ...
  c.^(-1).*t20+(1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2).*c.^2.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).* ...
  t+2.*t20+c.^2.*w).^2).*pi.^(-1/2).*((-1).*w+c.^(-2).*((-2).*t+2.*t20+ ...
  c.^2.*w)).*abs(c).^(-1)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2).*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.* ...
  w)+(-1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.* ...
  abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).* ...
  erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/4).*c.^3.*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs(c).^(-1).*erfz(c.^(-1) ...
  .*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2).*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+( ...
  -1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c) ...
  .^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2) ...
  .*c.*w)+(1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).* ...
  w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w))+((-1).*d12.^2.*exp(1) ...
  .^((-4).*t.*w)+d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-2).* ...
  t.*w).*(1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp(1).^((-2).*t.*w).*t.*(d12+ ...
  d22.*t+d12.*t.*w)).*((-2).*exp(1).^((-1).*t.*w).*v0.*w+exp(1).^((-1).* ...
  t.*w).*t.*v0.*w.^2+(-1).*exp(1).^((-1).*t.*w).*w.^2.*x0+exp(1).^((-1).* ...
  t.*w).*t.*w.^3.*x0+d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1) ...
  .*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*abs(c) ...
  .^(-1)+(-1).*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^( ...
  -1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t.*w.*abs(c).^( ...
  -1)+d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+( ...
  -1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t20.*w.*abs(c).^(-1)+( ...
  1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^( ...
  -1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*w.^2.*abs(c).^( ...
  -1)+(-1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).* ...
  (c.^(-1).*t20+(1/2).*c.*w).^2).*pi.^(-1/2).*w.^2.*abs(c).^(-1)+(-1/2).* ...
  c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*( ...
  (-2).*t+2.*t20+c.^2.*w).^2).*pi.^(-1/2).*w.*((-1).*w+c.^(-2).*((-2).*t+ ...
  2.*t20+c.^2.*w)).*abs(c).^(-1)+(-1).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+( ...
  1/4).*c.^2.*w.^2).*w.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+( ...
  -1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2) ...
  .*t.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w) ...
  +(-1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.* ...
  w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+( ...
  -1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^3.* ...
  abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1).*c.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.*abs(c).^(-1).*erfz( ...
  c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2).*t.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+( ...
  -1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.^2.* ...
  abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1).^(( ...
  -1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^3.*abs(c).^(-1).*erfz(c.^(-1).* ...
  t20+(1/2).*c.*w)))+(exp(1).^((-1).*t.*w).*v0+(-1).*exp(1).^((-1).*t.*w) ...
  .*t.*v0.*w+(-1).*exp(1).^((-1).*t.*w).*t.*w.^2.*x0+d.*exp(1).^((-1).*t.* ...
  w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).* ...
  c.*w).^2).*pi.^(-1/2).*t.*abs(c).^(-1)+(-1).*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).* ...
  c.*w).^2).*pi.^(-1/2).*t20.*abs(c).^(-1)+(-1/2).*c.^2.*d.*exp(1).^((-1) ...
  .*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+( ...
  -1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp(1).^(( ...
  -1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t20+(1/2).*c.*w).^2) ...
  .*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.* ...
  w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).*t+2.*t20+c.^2.*w).^2).*pi.^( ...
  -1/2).*((-1).*w+c.^(-2).*((-2).*t+2.*t20+c.^2.*w)).*abs(c).^(-1)+(1/2).* ...
  c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*abs(c).^(-1).*erfz( ...
  c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1) ...
  .*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1).*erfz(c.^(-1).*t+( ...
  -1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+ ...
  (1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).* ...
  t20+(-1/2).*c.*w)+(1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2).*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+( ...
  -1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2) ...
  .*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^(( ...
  -1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1).*erfz(c.^(-1).* ...
  t20+(1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.* ...
  w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/4).*c.^3.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs(c).^(-1).* ...
  erfz(c.^(-1).*t20+(1/2).*c.*w)).*(((-1).*d12.^2.*exp(1).^((-4).*t.*w)+ ...
  d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*(exp(1).^((-2).*t.*w).*t.*w.^2.*( ...
  (-1).*d12+d12.*t.*w+d11.*t.*w.^2)+exp(1).^((-2).*t.*w).*((-1)+t.*w).*(( ...
  -1).*d22+d22.*t.*w+d12.*t.*w.^2)).*(exp(1).^((-1).*t.*w).*v0+(-1).*exp( ...
  1).^((-1).*t.*w).*t.*v0.*w+(-1).*exp(1).^((-1).*t.*w).*t.*w.^2.*x0+d.* ...
  exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).* ...
  c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t.*abs(c).^(-1)+(-1).*d.* ...
  exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).* ...
  c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t20.*abs(c).^(-1)+(-1/2).* ...
  c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+ ...
  (-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2) ...
  .*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).* ...
  t20+(1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).*t+2.*t20+ ...
  c.^2.*w).^2).*pi.^(-1/2).*((-1).*w+c.^(-2).*((-2).*t+2.*t20+c.^2.*w)).* ...
  abs(c).^(-1)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2) ...
  .*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).* ...
  c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1) ...
  .*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^( ...
  (-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1) ...
  .*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/4).*c.^3.*d.*exp(1).^((-1).*t.* ...
  w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).* ...
  c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2).*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1).* ...
  erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+( ...
  1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+( ...
  1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs( ...
  c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w))+((-1).*d12.^2.*exp(1).^((-4).* ...
  t.*w)+d11.*d22.*exp(1).^((-4).*t.*w)).^(-1).*((-1).*exp(1).^((-1).*t.*w) ...
  .*(1+t.*w).*((-1).*d11.*exp(1).^((-1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1) ...
  .*t.*w).*(1+(-1).*t.*w))+(-1).*exp(1).^((-1).*t.*w).*t.*((-1).*d12.*exp( ...
  1).^((-1).*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))).*( ...
  (-2).*exp(1).^((-1).*t.*w).*v0.*w+exp(1).^((-1).*t.*w).*t.*v0.*w.^2+(-1) ...
  .*exp(1).^((-1).*t.*w).*w.^2.*x0+exp(1).^((-1).*t.*w).*t.*w.^3.*x0+d.* ...
  exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).* ...
  c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*abs(c).^(-1)+(-1).*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).* ...
  t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t.*w.*abs(c).^(-1)+d.*exp(1).^((-1).* ...
  t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+( ...
  -1/2).*c.*w).^2).*pi.^(-1/2).*t20.*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp( ...
  1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1) ...
  .*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*w.^2.*abs(c).^(-1)+(-1/2).*c.^2.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t20+( ...
  1/2).*c.*w).^2).*pi.^(-1/2).*w.^2.*abs(c).^(-1)+(-1/2).*c.^2.*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).*t+2.*t20+ ...
  c.^2.*w).^2).*pi.^(-1/2).*w.*((-1).*w+c.^(-2).*((-2).*t+2.*t20+c.^2.*w)) ...
  .*abs(c).^(-1)+(-1).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2) ...
  .*w.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2) ...
  .*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.^2.*abs(c) ...
  .^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).*c.*d.* ...
  exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.^2.*abs(c).^(-1).* ...
  erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^3.*abs(c).^(-1).*erfz(c.^(-1) ...
  .*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1).*c.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2).*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.* ...
  w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.^2.* ...
  abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1) ...
  .*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).* ...
  t20+(1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2).*w.^3.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)))));
elseif type=='x'
ECD=(1/2).*((exp(1).^((-2).*t.*w).*(1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp(1) ...
  .^((-2).*t.*w).*t.*(d12+d22.*t+d12.*t.*w)).^(-2).*(exp(1).^((-2).*t.*w) ...
  .*t.*(d22+d12.*w)+exp(1).^((-2).*t.*w).*(d12+d11.*w).*(1+t.*w)+exp(1).^( ...
  (-2).*t.*w).*w.*(d11+d12.*t+d11.*t.*w)+(-2).*exp(1).^((-2).*t.*w).*w.*( ...
  1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp(1).^((-2).*t.*w).*(d12+d22.*t+d12.* ...
  t.*w)+(-2).*exp(1).^((-2).*t.*w).*t.*w.*(d12+d22.*t+d12.*t.*w)).^2+2.*( ...
  exp(1).^((-2).*t.*w).*(1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp(1).^((-2).* ...
  t.*w).*t.*(d12+d22.*t+d12.*t.*w)).^(-1).*(exp(1).^((-1).*t.*w).*v0+(-1) ...
  .*exp(1).^((-1).*t.*w).*t.*v0.*w+(-1).*exp(1).^((-1).*t.*w).*t.*w.^2.* ...
  x0+d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+( ...
  -1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t.*abs(c).^(-1)+(-1).* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).* ...
  c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*t20.*abs(c).^(-1)+(-1/2).* ...
  c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+ ...
  (-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2) ...
  .*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).* ...
  t20+(1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp(1) ...
  .^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).*t+2.*t20+ ...
  c.^2.*w).^2).*pi.^(-1/2).*((-1).*w+c.^(-2).*((-2).*t+2.*t20+c.^2.*w)).* ...
  abs(c).^(-1)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2) ...
  .*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).* ...
  c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1) ...
  .*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^( ...
  (-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1) ...
  .*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/4).*c.^3.*d.*exp(1).^((-1).*t.* ...
  w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).* ...
  c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2).*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1).* ...
  erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+( ...
  1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+( ...
  1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs( ...
  c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)).^2);
elseif type=='v'
ECD=(1/2).*((exp(1).^((-2).*t.*w).*t.*w.^2.*((-1).*d12+d12.*t.*w+d11.*t.* ...
  w.^2)+exp(1).^((-2).*t.*w).*((-1)+t.*w).*((-1).*d22+d22.*t.*w+d12.*t.* ...
  w.^2)).^(-2).*(exp(1).^((-2).*t.*w).*t.*w.^2.*(d12.*w+d11.*w.^2)+exp(1) ...
  .^((-2).*t.*w).*((-1)+t.*w).*(d22.*w+d12.*w.^2)+exp(1).^((-2).*t.*w).* ...
  w.^2.*((-1).*d12+d12.*t.*w+d11.*t.*w.^2)+(-2).*exp(1).^((-2).*t.*w).*t.* ...
  w.^3.*((-1).*d12+d12.*t.*w+d11.*t.*w.^2)+exp(1).^((-2).*t.*w).*w.*((-1) ...
  .*d22+d22.*t.*w+d12.*t.*w.^2)+(-2).*exp(1).^((-2).*t.*w).*w.*((-1)+t.*w) ...
  .*((-1).*d22+d22.*t.*w+d12.*t.*w.^2)).^2+2.*(exp(1).^((-2).*t.*w).*t.* ...
  w.^2.*((-1).*d12+d12.*t.*w+d11.*t.*w.^2)+exp(1).^((-2).*t.*w).*((-1)+t.* ...
  w).*((-1).*d22+d22.*t.*w+d12.*t.*w.^2)).^(-1).*((-2).*exp(1).^((-1).*t.* ...
  w).*v0.*w+exp(1).^((-1).*t.*w).*t.*v0.*w.^2+(-1).*exp(1).^((-1).*t.*w).* ...
  w.^2.*x0+exp(1).^((-1).*t.*w).*t.*w.^3.*x0+d.*exp(1).^((-1).*t.*w+t20.* ...
  w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w) ...
  .^2).*pi.^(-1/2).*abs(c).^(-1)+(-1).*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).* ...
  pi.^(-1/2).*t.*w.*abs(c).^(-1)+d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2).*pi.^( ...
  -1/2).*t20.*w.*abs(c).^(-1)+(1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+( ...
  1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w).^2) ...
  .*pi.^(-1/2).*w.^2.*abs(c).^(-1)+(-1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2+(-1).*(c.^(-1).*t20+(1/2).*c.*w).^2).*pi.^( ...
  -1/2).*w.^2.*abs(c).^(-1)+(-1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+( ...
  1/4).*c.^2.*w.^2+(-1/4).*c.^(-2).*((-2).*t+2.*t20+c.^2.*w).^2).*pi.^( ...
  -1/2).*w.*((-1).*w+c.^(-2).*((-2).*t+2.*t20+c.^2.*w)).*abs(c).^(-1)+(-1) ...
  .*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.*abs(c).^(-1).* ...
  erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^(( ...
  -1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.^2.*abs(c).^(-1).*erfz(c.^(-1) ...
  .*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1) ...
  .*c.^(-1).*t20+(-1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.* ...
  w+(1/4).*c.^2.*w.^2).*w.^3.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).* ...
  t20+(-1/2).*c.*w)+(-1).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.* ...
  w.^2).*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.*d.*exp( ...
  1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.^2.*abs(c).^(-1).*erfz( ...
  c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2).*t20.*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+( ...
  -1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^3.* ...
  abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)).^2);
else
       DISP('not enough arguments');    
end
end