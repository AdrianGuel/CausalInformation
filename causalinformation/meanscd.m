function MCD=meanscd(w,t,t20,c,d,x0,v0)
MCD=[exp(1).^((-1).*t.*w).*t.*v0+exp(1).^((-1).*t.*w).*x0+exp(1).^((-1).*t.* ...
  w).*t.*w.*x0+(-1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.* ...
  w.^2+(-1).*(c.^(-1).*t20+(1/2).*c.*w).^2).*pi.^(-1/2).*abs(c).^(-1)+( ...
  1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+(-1/4).*c.^( ...
  -2).*((-2).*t+2.*t20+c.^2.*w).^2).*pi.^(-1/2).*abs(c).^(-1)+(1/2).*c.* ...
  d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*abs(c).^(-1).*erfz( ...
  c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).*c.*d.*exp(1).^((-1) ...
  .*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1) ...
  .*c.^(-1).*t20+(-1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.* ...
  w+(1/4).*c.^2.*w.^2).*w.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).* ...
  t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.* ...
  w.^2).*t.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.*d.*exp( ...
  1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*abs(c).^(-1).*erfz(c.^( ...
  -1).*t20+(1/2).*c.*w)+(-1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2).*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w);exp(1).^(( ...
  -1).*t.*w).*v0+(-1).*exp(1).^((-1).*t.*w).*t.*v0.*w+(-1).*exp(1).^((-1) ...
  .*t.*w).*t.*w.^2.*x0+(1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2+(-1).*(c.^(-1).*t20+(1/2).*c.*w).^2).*pi.^(-1/2).*w.*abs(c) ...
  .^(-1)+(-1/2).*c.^2.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2+( ...
  -1/4).*c.^(-2).*((-2).*t+2.*t20+c.^2.*w).^2).*pi.^(-1/2).*w.*abs(c).^( ...
  -1)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*abs(c) ...
  .^(-1).*erfz(c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(-1/2).*c.*d.* ...
  exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1).*erfz( ...
  c.^(-1).*t+(-1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).* ...
  t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t+( ...
  -1).*c.^(-1).*t20+(-1/2).*c.*w)+(1/4).*c.^3.*d.*exp(1).^((-1).*t.*w+ ...
  t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs(c).^(-1).*erfz(c.^(-1).*t+(-1).*c.^( ...
  -1).*t20+(-1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).* ...
  c.^2.*w.^2).*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(-1/2).*c.*d.* ...
  exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*t.*w.*abs(c).^(-1).*erfz( ...
  c.^(-1).*t20+(1/2).*c.*w)+(1/2).*c.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4) ...
  .*c.^2.*w.^2).*t20.*w.*abs(c).^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)+(1/4) ...
  .*c.^3.*d.*exp(1).^((-1).*t.*w+t20.*w+(1/4).*c.^2.*w.^2).*w.^2.*abs(c) ...
  .^(-1).*erfz(c.^(-1).*t20+(1/2).*c.*w)];
end