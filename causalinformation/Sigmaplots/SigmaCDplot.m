function SCD=SigmaCDplot(w,t,t0,a,b,d11,d12,d22,D0)
SCD=[exp(1).^((-2).*t.*w).*(1+t.*w).*(d11+d12.*t+d11.*t.*w)+exp(1).^((-2).* ...
  t.*w).*t.*(d12+d22.*t+d12.*t.*w)+(1/2).*pi.^(-1/2).*w.^(-1).*abs(a).^( ...
  -1).*((1/2).*D0.*exp(1).^((-2).*t.*w).*pi.^(1/2).*w.^(-2).*((-1)+exp(1) ...
  .^(2.*t.*w)+(-2).*t.*w).*abs(a)+(1/2).*D0.*exp(1).^((-2).*t.*w).*pi.^( ...
  1/2).*w.^(-2).*((-1)+exp(1).^(2.*t.*w)+(-2).*t.*w+(-4).*t.^2.*w.^2).* ...
  abs(a)+(a.^(-2)).^(-1/2).*b.*(2.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.* ...
  w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).* ...
  a.^2.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).* ...
  (t0+a.^2.*w).^2).*w.^(1/2)).*w.^(1/2).*(t+(-1).*t0+(-1).*a.^2.*w)+exp(1) ...
  .^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((-1/2).*(a.^(-2)).^(-1/2).*exp(1) ...
  .^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+(1/2) ...
  .*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w) ...
  .^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2))+(1/4).*(a.^(-2)).^(-1/2) ...
  .*a.^(-1).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*(2.*((-1).*t+t0+2.* ...
  a.^2.*w)+(2.*t.*w.^(1/2)+(-2).*t0.*w.^(1/2)+(-2).*a.^2.*w.^(3/2)).^2).*( ...
  pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+pi.^(1/2).*erfz(a.^(-1).*(t+(-1).*t0+( ...
  -1).*a.^2.*w))))+b.*pi.^(1/2).*(exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2) ...
  .*((-1).*a.^2.*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*pi.^(-1/2)+a.^2.* ...
  exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w).^2).*pi.^(-1/2))+exp(1).^(( ...
  -2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((-1).*t+t0+a.^2.*w).*((-1).*a.*erfz(a.^( ...
  -1).*t0+a.*w)+a.*erfz(a.^(-1).*((-1).*t+t0+a.^2.*w))))+(a.^(-2)).^(-1/2) ...
  .*b.*((-2).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^(-1/2) ...
  .*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w).^2).*w.^(1/2)).*w.^( ...
  1/2).*(t+(-1).*t0+(-1).*a.^2.*w)+exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.* ...
  w.^2).*((-1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2) ...
  .*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+(1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1) ...
  .*a.^(-2).*((-1).*t+t0+a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.* ...
  w.^2))+(1/4).*(a.^(-2)).^(-1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*(2.*((-1).*t+t0+2.*a.^2.*w)+((-2).*t.*w.^(1/2)+2.*t0.*w.^( ...
  1/2)+2.*a.^2.*w.^(3/2)).^2).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+(-1).* ...
  pi.^(1/2).*erfz(a.^(-1).*((-1).*t+t0+a.^2.*w))))),(-1).*exp(1).^((-1).* ...
  t.*w).*t.*w.^2.*(d12.*exp(1).^((-1).*t.*w).*t+d11.*exp(1).^((-1).*t.*w) ...
  .*(1+t.*w))+(-1).*exp(1).^((-1).*t.*w).*((-1)+t.*w).*(d22.*exp(1).^((-1) ...
  .*t.*w).*t+d12.*exp(1).^((-1).*t.*w).*(1+t.*w))+(1/4).*pi.^(-1/2).*w.^( ...
  -1).*abs(a).^(-1).*(2.*D0.*exp(1).^((-2).*t.*w).*pi.^(1/2).*t.*abs(a)+ ...
  2.*D0.*exp(1).^((-2).*t.*w).*pi.^(1/2).*w.^(-1).*((-1)+exp(1).^(2.*t.*w) ...
  +(-2).*t.*w).*abs(a)+(-2).*D0.*exp(1).^((-2).*t.*w).*pi.^(1/2).*w.^(-1) ...
  .*((-1)+exp(1).^(2.*t.*w)+(-1).*t.*w+(-2).*t.^2.*w.^2).*abs(a)+(a.^(-2)) ...
  .^(-1/2).*b.*(2.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w).^2).*w.^(1/2)+ ...
  (-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*w.^( ...
  1/2)).*w.^(1/2)+2.*(a.^(-2)).^(-1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.* ...
  t0.*w+a.^2.*w.^2).*w.*(t+(-1).*t0+(-1).*a.^2.*w).*(pi.^(1/2).*erfz(a.^( ...
  -1).*t0+a.*w)+pi.^(1/2).*erfz(a.^(-1).*(t+(-1).*t0+(-1).*a.^2.*w))))+( ...
  a.^(-2)).^(-1/2).*b.*((-4).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*(( ...
  a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w).^2) ...
  .*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w) ...
  .^2).*w.^(1/2)).*w.^(3/2).*(t+(-1).*t0+(-1).*a.^2.*w)+(-2).*exp(1).^(( ...
  -2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*((-1/2).*(a.^(-2)).^(-1/2).*exp(1).^( ...
  (-1).*a.^(-2).*(t0+a.^2.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+(1/2).*( ...
  a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w).^2) ...
  .*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2))+(-1/2).*(a.^(-2)).^(-1/2).* ...
  a.^(-1).*(exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2)+exp(1).^((-2).*t.*w+ ...
  2.*t0.*w+a.^2.*w.^2).*w.*(2.*((-1).*t+t0+2.*a.^2.*w)+(2.*t.*w.^(1/2)+( ...
  -2).*t0.*w.^(1/2)+(-2).*a.^2.*w.^(3/2)).^2)).*(pi.^(1/2).*erfz(a.^(-1).* ...
  t0+a.*w)+pi.^(1/2).*erfz(a.^(-1).*(t+(-1).*t0+(-1).*a.^2.*w))))+b.*pi.^( ...
  1/2).*(2.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*(a.^2.*exp(1).^((-1) ...
  .*(a.^(-1).*t0+a.*w).^2).*pi.^(-1/2)+(-1).*a.^2.*exp(1).^((-1).*a.^(-2) ...
  .*((-1).*t+t0+a.^2.*w).^2).*pi.^(-1/2)).*w+(exp(1).^((-2).*t.*w+2.*t0.* ...
  w+a.^2.*w.^2)+2.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*((-1).*t+ ...
  t0+a.^2.*w)).*(a.*erfz(a.^(-1).*t0+a.*w)+(-1).*a.*erfz(a.^(-1).*((-1).*t+ ...
  t0+a.^2.*w))))+(a.^(-2)).^(-1/2).*b.*((-2).*exp(1).^((-2).*t.*w+2.*t0.* ...
  w+a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2) ...
  .*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+ ...
  a.^2.*w).^2).*w.^(1/2)).*w.^(1/2)+2.*(a.^(-2)).^(-1/2).*a.^(-1).*exp(1) ...
  .^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*(t+(-1).*t0+(-1).*a.^2.*w).*( ...
  pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+(-1).*pi.^(1/2).*erfz(a.^(-1).*((-1).*t+ ...
  t0+a.^2.*w))))+(a.^(-2)).^(-1/2).*b.*(4.*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).* ...
  w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+ ...
  a.^2.*w).^2).*w.^(1/2)).*w.^(3/2).*(t+(-1).*t0+(-1).*a.^2.*w)+(-2).*exp( ...
  1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*((-1/2).*(a.^(-2)).^(-1/2).* ...
  exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+( ...
  1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w) ...
  .^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2))+(-1/2).*(a.^(-2)).^( ...
  -1/2).*a.^(-1).*(exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2)+exp(1).^((-2) ...
  .*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*(2.*((-1).*t+t0+2.*a.^2.*w)+((-2).*t.* ...
  w.^(1/2)+2.*t0.*w.^(1/2)+2.*a.^2.*w.^(3/2)).^2)).*(pi.^(1/2).*erfz(a.^( ...
  -1).*t0+a.*w)+(-1).*pi.^(1/2).*erfz(a.^(-1).*((-1).*t+t0+a.^2.*w))))); ...
  exp(1).^((-1).*t.*w).*(1+t.*w).*((-1).*d11.*exp(1).^((-1).*t.*w).*t.* ...
  w.^2+d12.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))+exp(1).^((-1).*t.*w).* ...
  t.*((-1).*d12.*exp(1).^((-1).*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).* ...
  (1+(-1).*t.*w))+(1/4).*pi.^(-1/2).*w.^(-1).*abs(a).^(-1).*(2.*D0.*exp(1) ...
  .^((-2).*t.*w).*pi.^(1/2).*t.*abs(a)+2.*D0.*exp(1).^((-2).*t.*w).*pi.^( ...
  1/2).*w.^(-1).*((-1)+exp(1).^(2.*t.*w)+(-2).*t.*w).*abs(a)+(-2).*D0.* ...
  exp(1).^((-2).*t.*w).*pi.^(1/2).*w.^(-1).*((-1)+exp(1).^(2.*t.*w)+(-1).* ...
  t.*w+(-2).*t.^2.*w.^2).*abs(a)+(a.^(-2)).^(-1/2).*b.*(2.*exp(1).^((-2).* ...
  t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*( ...
  t+(-1).*t0+(-1).*a.^2.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1) ...
  .^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*w.^(1/2)).*w.^(1/2)+2.*(a.^(-2)).^( ...
  -1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*(t+(-1).* ...
  t0+(-1).*a.^2.*w).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+pi.^(1/2).*erfz(a.^( ...
  -1).*(t+(-1).*t0+(-1).*a.^2.*w))))+(a.^(-2)).^(-1/2).*b.*((-4).*exp(1) ...
  .^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).* ...
  a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*w.^(1/2)).*w.^(3/2).*( ...
  t+(-1).*t0+(-1).*a.^2.*w)+(-2).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2) ...
  .*w.*((-1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w) ...
  .^2).*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+(1/2).*(a.^(-2)).^(-1/2).*exp(1).^(( ...
  -1).*a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+ ...
  2.*a.^2.*w.^2))+(-1/2).*(a.^(-2)).^(-1/2).*a.^(-1).*(exp(1).^((-2).*t.* ...
  w+2.*t0.*w+a.^2.*w.^2)+exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*(2.* ...
  ((-1).*t+t0+2.*a.^2.*w)+(2.*t.*w.^(1/2)+(-2).*t0.*w.^(1/2)+(-2).*a.^2.* ...
  w.^(3/2)).^2)).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+pi.^(1/2).*erfz(a.^(-1) ...
  .*(t+(-1).*t0+(-1).*a.^2.*w))))+b.*pi.^(1/2).*(2.*exp(1).^((-2).*t.*w+ ...
  2.*t0.*w+a.^2.*w.^2).*(a.^2.*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).* ...
  pi.^(-1/2)+(-1).*a.^2.*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w).^2) ...
  .*pi.^(-1/2)).*w+(exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2)+2.*exp(1).^(( ...
  -2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*((-1).*t+t0+a.^2.*w)).*(a.*erfz(a.^( ...
  -1).*t0+a.*w)+(-1).*a.*erfz(a.^(-1).*((-1).*t+t0+a.^2.*w))))+(a.^(-2)).^( ...
  -1/2).*b.*((-2).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)) ...
  .^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w).^2).*w.^(1/2)).* ...
  w.^(1/2)+2.*(a.^(-2)).^(-1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*w.*(t+(-1).*t0+(-1).*a.^2.*w).*(pi.^(1/2).*erfz(a.^(-1).*t0+ ...
  a.*w)+(-1).*pi.^(1/2).*erfz(a.^(-1).*((-1).*t+t0+a.^2.*w))))+(a.^(-2)).^( ...
  -1/2).*b.*(4.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)) ...
  .^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w).^2).*w.^(1/2)).* ...
  w.^(3/2).*(t+(-1).*t0+(-1).*a.^2.*w)+(-2).*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*w.*((-1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1).*t0+ ...
  a.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+(1/2).*(a.^(-2)).^(-1/2).*exp( ...
  1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+ ...
  2.*a.^2.*w.^2))+(-1/2).*(a.^(-2)).^(-1/2).*a.^(-1).*(exp(1).^((-2).*t.* ...
  w+2.*t0.*w+a.^2.*w.^2)+exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.*(2.* ...
  ((-1).*t+t0+2.*a.^2.*w)+((-2).*t.*w.^(1/2)+2.*t0.*w.^(1/2)+2.*a.^2.*w.^( ...
  3/2)).^2)).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+(-1).*pi.^(1/2).*erfz(a.^( ...
  -1).*((-1).*t+t0+a.^2.*w))))),(-1).*exp(1).^((-1).*t.*w).*t.*w.^2.*((-1) ...
  .*d11.*exp(1).^((-1).*t.*w).*t.*w.^2+d12.*exp(1).^((-1).*t.*w).*(1+(-1) ...
  .*t.*w))+(-1).*exp(1).^((-1).*t.*w).*((-1)+t.*w).*((-1).*d12.*exp(1).^(( ...
  -1).*t.*w).*t.*w.^2+d22.*exp(1).^((-1).*t.*w).*(1+(-1).*t.*w))+(1/8).* ...
  pi.^(-1/2).*w.^(-1).*abs(a).^(-1).*((-6).*D0.*exp(1).^((-2).*t.*w).* ...
  pi.^(1/2).*((-1)+exp(1).^(2.*t.*w)+(-2).*t.*w).*abs(a)+(-2).*D0.*exp(1) ...
  .^((-2).*t.*w).*pi.^(1/2).*((-1)+exp(1).^(2.*t.*w)+(-2).*t.*w+(-4).* ...
  t.^2.*w.^2).*abs(a)+4.*D0.*exp(1).^((-2).*t.*w).*pi.^(1/2).*((-3)+3.* ...
  exp(1).^(2.*t.*w)+(-2).*t.*w+(-4).*t.^2.*w.^2).*abs(a)+(-2).*(a.^(-2)) ...
  .^(-1/2).*b.*(4.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).*a.^2.*w).^2).*w.^(1/2)+ ...
  (-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*w.^( ...
  1/2)).*w.^(3/2)+4.*(a.^(-2)).^(-1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.* ...
  t0.*w+a.^2.*w.^2).*w.^2.*(t+(-1).*t0+(-1).*a.^2.*w).*(pi.^(1/2).*erfz( ...
  a.^(-1).*t0+a.*w)+pi.^(1/2).*erfz(a.^(-1).*(t+(-1).*t0+(-1).*a.^2.*w))))+ ...
  (-4).*(a.^(-2)).^(-1/2).*b.*w.^2.*(2.*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+( ...
  -1).*a.^2.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^( ...
  -2).*(t0+a.^2.*w).^2).*w.^(1/2)).*w.^(1/2).*(t+(-1).*t0+(-1).*a.^2.*w)+ ...
  exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((-1/2).*(a.^(-2)).^(-1/2).* ...
  exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.*w.^2)+ ...
  (1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+(-1).*t0+(-1).* ...
  a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2))+(1/4).*(a.^(-2)) ...
  .^(-1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*(2.*((-1).* ...
  t+t0+2.*a.^2.*w)+(2.*t.*w.^(1/2)+(-2).*t0.*w.^(1/2)+(-2).*a.^2.*w.^(3/2) ...
  ).^2).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+pi.^(1/2).*erfz(a.^(-1).*(t+(-1) ...
  .*t0+(-1).*a.^2.*w))))+2.*(a.^(-2)).^(-1/2).*b.*(8.*exp(1).^((-2).*t.*w+ ...
  2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t+( ...
  -1).*t0+(-1).*a.^2.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^(( ...
  -1).*a.^(-2).*(t0+a.^2.*w).^2).*w.^(1/2)).*w.^(5/2).*(t+(-1).*t0+(-1).* ...
  a.^2.*w)+4.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.^2.*((-1/2).*( ...
  a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(t0+a.^2.*w).^2).*((-1)+2.* ...
  t0.*w+2.*a.^2.*w.^2)+(1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*( ...
  t+(-1).*t0+(-1).*a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2)) ...
  +(1/2).*(a.^(-2)).^(-1/2).*a.^(-1).*(4.*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*w+2.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.^2.*(2.*(( ...
  -1).*t+t0+2.*a.^2.*w)+(2.*t.*w.^(1/2)+(-2).*t0.*w.^(1/2)+(-2).*a.^2.* ...
  w.^(3/2)).^2)).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+pi.^(1/2).*erfz(a.^(-1) ...
  .*(t+(-1).*t0+(-1).*a.^2.*w))))+4.*b.*pi.^(1/2).*w.^2.*(exp(1).^((-2).* ...
  t.*w+2.*t0.*w+a.^2.*w.^2).*((-1).*a.^2.*exp(1).^((-1).*(a.^(-1).*t0+a.* ...
  w).^2).*pi.^(-1/2)+a.^2.*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+a.^2.*w) ...
  .^2).*pi.^(-1/2))+exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((-1).*t+t0+ ...
  a.^2.*w).*((-1).*a.*erfz(a.^(-1).*t0+a.*w)+a.*erfz(a.^(-1).*((-1).*t+t0+ ...
  a.^2.*w))))+2.*(a.^(-2)).^(-1/2).*b.*(4.*exp(1).^((-2).*t.*w+2.*t0.*w+ ...
  a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).* ...
  w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+ ...
  a.^2.*w).^2).*w.^(1/2)).*w.^(3/2)+(-4).*(a.^(-2)).^(-1/2).*a.^(-1).*exp( ...
  1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.^2.*(t+(-1).*t0+(-1).*a.^2.*w).* ...
  (pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+(-1).*pi.^(1/2).*erfz(a.^(-1).*((-1).* ...
  t+t0+a.^2.*w))))+(-4).*(a.^(-2)).^(-1/2).*b.*w.^2.*((-2).*exp(1).^((-2) ...
  .*t.*w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1) ...
  .*t0+a.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2) ...
  .*((-1).*t+t0+a.^2.*w).^2).*w.^(1/2)).*w.^(1/2).*(t+(-1).*t0+(-1).* ...
  a.^2.*w)+exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*((-1/2).*(a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.* ...
  w.^2)+(1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+ ...
  a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2))+(1/4).*(a.^(-2)) ...
  .^(-1/2).*a.^(-1).*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*(2.*((-1).* ...
  t+t0+2.*a.^2.*w)+((-2).*t.*w.^(1/2)+2.*t0.*w.^(1/2)+2.*a.^2.*w.^(3/2)) ...
  .^2).*(pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+(-1).*pi.^(1/2).*erfz(a.^(-1).*(( ...
  -1).*t+t0+a.^2.*w))))+2.*(a.^(-2)).^(-1/2).*b.*((-8).*exp(1).^((-2).*t.* ...
  w+2.*t0.*w+a.^2.*w.^2).*((a.^(-2)).^(-1/2).*exp(1).^((-1).*(a.^(-1).*t0+ ...
  a.*w).^2).*w.^(1/2)+(-1).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*(( ...
  -1).*t+t0+a.^2.*w).^2).*w.^(1/2)).*w.^(5/2).*(t+(-1).*t0+(-1).*a.^2.*w)+ ...
  4.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.^2.*((-1/2).*(a.^(-2)).^( ...
  -1/2).*exp(1).^((-1).*(a.^(-1).*t0+a.*w).^2).*((-1)+2.*t0.*w+2.*a.^2.* ...
  w.^2)+(1/2).*(a.^(-2)).^(-1/2).*exp(1).^((-1).*a.^(-2).*((-1).*t+t0+ ...
  a.^2.*w).^2).*((-1)+(-2).*t.*w+2.*t0.*w+2.*a.^2.*w.^2))+(1/2).*(a.^(-2)) ...
  .^(-1/2).*a.^(-1).*(4.*exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w+2.* ...
  exp(1).^((-2).*t.*w+2.*t0.*w+a.^2.*w.^2).*w.^2.*(2.*((-1).*t+t0+2.* ...
  a.^2.*w)+((-2).*t.*w.^(1/2)+2.*t0.*w.^(1/2)+2.*a.^2.*w.^(3/2)).^2)).*( ...
  pi.^(1/2).*erfz(a.^(-1).*t0+a.*w)+(-1).*pi.^(1/2).*erfz(a.^(-1).*((-1).*t+ ...
  t0+a.^2.*w)))))];
end