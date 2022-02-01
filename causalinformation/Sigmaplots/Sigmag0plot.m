function Sg0=Sigmag0plot(w,t,t0,a,b,d11,d12,d22,D0)
Sg0=[(1/4).*exp(1).^((-2).*t.*((-1).*w.^2).^(1/2)).*(1+exp(1).^(2.*t.*((-1) ...
  .*w.^2).^(1/2))).*((-1).*w.^2).^(-1/2).*((-1).*d12+d12.*exp(1).^(2.*t.*( ...
  (-1).*w.^2).^(1/2))+d11.*((-1).*w.^2).^(1/2)+d11.*exp(1).^(2.*t.*((-1).* ...
  w.^2).^(1/2)).*((-1).*w.^2).^(1/2))+(-1/4).*exp(1).^((-2).*t.*((-1).* ...
  w.^2).^(1/2)).*((-1)+exp(1).^(2.*t.*((-1).*w.^2).^(1/2))).*w.^(-2).*(( ...
  -1).*d22+d22.*exp(1).^(2.*t.*((-1).*w.^2).^(1/2))+d12.*((-1).*w.^2).^( ...
  1/2)+d12.*exp(1).^(2.*t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2))+( ...
  -1/4).*exp(1).^((-1).*a.^2.*w.^2+(-2).*t.*((-1).*w.^2).^(1/2)+(-2).*(t+( ...
  -1).*t0).*((-1).*w.^2).^(1/2)).*w.^(-2).*((-1).*w.^2).^(-1/2).*abs(a).^( ...
  -1).*((-1).*D0.*exp(1).^(a.^2.*w.^2+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2) ...
  ).*abs(a)+D0.*exp(1).^(a.^2.*w.^2+4.*t.*((-1).*w.^2).^(1/2)+2.*(t+(-1).* ...
  t0).*((-1).*w.^2).^(1/2)).*abs(a)+(-4).*D0.*exp(1).^(a.^2.*w.^2+2.*t.*(( ...
  -1).*w.^2).^(1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*t.*((-1).*w.^2) ...
  .^(1/2).*abs(a)+(-2).*a.*b.*exp(1).^(a.^2.*w.^2+2.*t.*((-1).*w.^2).^( ...
  1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2).*erfz( ...
  a.^(-1).*t0)+2.*a.*b.*exp(1).^(a.^2.*w.^2+2.*t.*((-1).*w.^2).^(1/2)+2.*( ...
  t+(-1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2).*erfz((1/2).*a.^( ...
  -1).*((-2).*t+2.*t0))+a.*b.*exp(1).^(2.*t.*((-1).*w.^2).^(1/2)+4.*(t+( ...
  -1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2).*erfz(a.^(-1).*t0+( ...
  -1).*a.*((-1).*w.^2).^(1/2))+a.*b.*exp(1).^(2.*t.*((-1).*w.^2).^(1/2)).* ...
  ((-1).*w.^2).^(1/2).*erfz(a.^(-1).*t0+a.*((-1).*w.^2).^(1/2))+a.*b.*exp( ...
  1).^(2.*t.*((-1).*w.^2).^(1/2)+4.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*(( ...
  -1).*w.^2).^(1/2).*erfz(a.^(-1).*(t+(-1).*t0+a.^2.*((-1).*w.^2).^(1/2)))+ ...
  (-1).*a.*b.*exp(1).^(2.*t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2).* ...
  erfz(a.^(-1).*((-1).*t+t0+a.^2.*((-1).*w.^2).^(1/2)))),(1/2).*(exp(1).^(( ...
  -1).*t.*((-1).*w.^2).^(1/2))+(-1).*exp(1).^(t.*((-1).*w.^2).^(1/2))).* ...
  w.^2.*((-1).*w.^2).^(-1/2).*((-1/2).*d12.*(exp(1).^((-1).*t.*((-1).* ...
  w.^2).^(1/2))+(-1).*exp(1).^(t.*((-1).*w.^2).^(1/2))).*((-1).*w.^2).^( ...
  -1/2)+(-1/2).*d11.*((-1).*w.^2).^(-1/2).*((-1).*exp(1).^((-1).*t.*((-1) ...
  .*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)+(-1).*exp(1).^(t.*((-1).*w.^2).^( ...
  1/2)).*((-1).*w.^2).^(1/2)))+(-1/2).*((-1).*w.^2).^(-1/2).*((-1).*exp(1) ...
  .^((-1).*t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)+(-1).*exp(1).^(t.* ...
  ((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)).*((-1/2).*d22.*(exp(1).^((-1) ...
  .*t.*((-1).*w.^2).^(1/2))+(-1).*exp(1).^(t.*((-1).*w.^2).^(1/2))).*((-1) ...
  .*w.^2).^(-1/2)+(-1/2).*d12.*((-1).*w.^2).^(-1/2).*((-1).*exp(1).^((-1) ...
  .*t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)+(-1).*exp(1).^(t.*((-1).* ...
  w.^2).^(1/2)).*((-1).*w.^2).^(1/2)))+(-1/4).*pi.^(-1/2).*w.^(-2).*abs(a) ...
  .^(-1).*(D0.*((-1)+exp(1).^(2.*t.*((-1).*w.^2).^(1/2))).*pi.^(1/2).*abs( ...
  a)+(-1).*D0.*exp(1).^((-2).*t.*((-1).*w.^2).^(1/2)).*((-1)+exp(1).^(2.* ...
  t.*((-1).*w.^2).^(1/2))).*pi.^(1/2).*abs(a)+2.*(a.^(-2)).^(-1/2).*b.* ...
  exp(1).^((-1).*a.^2.*w.^2+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*((-1).* ...
  w.^2).^(1/2).*((1/2).*(a.^(-2)).^(1/2).*pi.^(1/2).*(2.*(t+(-1).*t0)+2.* ...
  a.^2.*((-1).*w.^2).^(1/2)).*(a.^(-2).*(2.*(t+(-1).*t0)+2.*a.^2.*((-1).* ...
  w.^2).^(1/2)).^2).^(-1/2).*erfz((1/2).*(a.^(-2).*(2.*(t+(-1).*t0)+2.* ...
  a.^2.*((-1).*w.^2).^(1/2)).^2).^(1/2))+(1/2).*(a.^(-2)).^(-1/2).*pi.^( ...
  1/2).*(2.*a.^(-2).*t0+(-2).*((-1).*w.^2).^(1/2)).*(a.^2.*(2.*a.^(-2).* ...
  t0+(-2).*((-1).*w.^2).^(1/2)).^2).^(-1/2).*erfz((1/2).*(a.^(-2).*((-2).* ...
  t0+2.*a.^2.*((-1).*w.^2).^(1/2)).^2).^(1/2)))+(-2).*(a.^(-2)).^(-1/2).* ...
  b.*exp(1).^((-1).*a.^2.*w.^2+(-2).*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*( ...
  (-1).*w.^2).^(1/2).*((1/2).*(a.^(-2)).^(-1/2).*pi.^(1/2).*(2.*a.^(-2).* ...
  t0+2.*((-1).*w.^2).^(1/2)).*(a.^2.*(2.*a.^(-2).*t0+2.*((-1).*w.^2).^( ...
  1/2)).^2).^(-1/2).*erfz((1/2).*(a.^(-2).*(2.*t0+2.*a.^2.*((-1).*w.^2).^( ...
  1/2)).^2).^(1/2))+(-1/2).*(a.^(-2)).^(1/2).*pi.^(1/2).*((-2).*t+2.*t0+ ...
  2.*a.^2.*((-1).*w.^2).^(1/2)).*(a.^(-2).*((-2).*t+2.*t0+2.*a.^2.*((-1).* ...
  w.^2).^(1/2)).^2).^(-1/2).*erfz((1/2).*(a.^(-2).*((-2).*t+2.*t0+2.*a.^2.* ...
  ((-1).*w.^2).^(1/2)).^2).^(1/2))));(-1/2).*((-1).*w.^2).^(-1/2).*((-1).* ...
  exp(1).^((-1).*t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)+(-1).*exp(1) ...
  .^(t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)).*((1/2).*d11.*(exp(1) ...
  .^((-1).*t.*((-1).*w.^2).^(1/2))+(-1).*exp(1).^(t.*((-1).*w.^2).^(1/2))) ...
  .*w.^2.*((-1).*w.^2).^(-1/2)+(-1/2).*d12.*((-1).*w.^2).^(-1/2).*((-1).* ...
  exp(1).^((-1).*t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)+(-1).*exp(1) ...
  .^(t.*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2)))+(-1/2).*(exp(1).^((-1) ...
  .*t.*((-1).*w.^2).^(1/2))+(-1).*exp(1).^(t.*((-1).*w.^2).^(1/2))).*((-1) ...
  .*w.^2).^(-1/2).*((1/2).*d12.*(exp(1).^((-1).*t.*((-1).*w.^2).^(1/2))+( ...
  -1).*exp(1).^(t.*((-1).*w.^2).^(1/2))).*w.^2.*((-1).*w.^2).^(-1/2)+( ...
  -1/2).*d22.*((-1).*w.^2).^(-1/2).*((-1).*exp(1).^((-1).*t.*((-1).*w.^2) ...
  .^(1/2)).*((-1).*w.^2).^(1/2)+(-1).*exp(1).^(t.*((-1).*w.^2).^(1/2)).*(( ...
  -1).*w.^2).^(1/2)))+(-1/4).*pi.^(-1/2).*w.^(-2).*abs(a).^(-1).*(D0.*(( ...
  -1)+exp(1).^(2.*t.*((-1).*w.^2).^(1/2))).*pi.^(1/2).*abs(a)+(-1).*D0.* ...
  exp(1).^((-2).*t.*((-1).*w.^2).^(1/2)).*((-1)+exp(1).^(2.*t.*((-1).* ...
  w.^2).^(1/2))).*pi.^(1/2).*abs(a)+2.*(a.^(-2)).^(-1/2).*b.*exp(1).^((-1) ...
  .*a.^2.*w.^2+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2) ...
  .*((1/2).*(a.^(-2)).^(1/2).*pi.^(1/2).*(2.*(t+(-1).*t0)+2.*a.^2.*((-1).* ...
  w.^2).^(1/2)).*(a.^(-2).*(2.*(t+(-1).*t0)+2.*a.^2.*((-1).*w.^2).^(1/2)) ...
  .^2).^(-1/2).*erfz((1/2).*(a.^(-2).*(2.*(t+(-1).*t0)+2.*a.^2.*((-1).* ...
  w.^2).^(1/2)).^2).^(1/2))+(1/2).*(a.^(-2)).^(-1/2).*pi.^(1/2).*(2.*a.^( ...
  -2).*t0+(-2).*((-1).*w.^2).^(1/2)).*(a.^2.*(2.*a.^(-2).*t0+(-2).*((-1).* ...
  w.^2).^(1/2)).^2).^(-1/2).*erfz((1/2).*(a.^(-2).*((-2).*t0+2.*a.^2.*((-1) ...
  .*w.^2).^(1/2)).^2).^(1/2)))+(-2).*(a.^(-2)).^(-1/2).*b.*exp(1).^((-1).* ...
  a.^2.*w.^2+(-2).*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2) ...
  .*((1/2).*(a.^(-2)).^(-1/2).*pi.^(1/2).*(2.*a.^(-2).*t0+2.*((-1).*w.^2) ...
  .^(1/2)).*(a.^2.*(2.*a.^(-2).*t0+2.*((-1).*w.^2).^(1/2)).^2).^(-1/2).* ...
  erfz((1/2).*(a.^(-2).*(2.*t0+2.*a.^2.*((-1).*w.^2).^(1/2)).^2).^(1/2))+( ...
  -1/2).*(a.^(-2)).^(1/2).*pi.^(1/2).*((-2).*t+2.*t0+2.*a.^2.*((-1).*w.^2) ...
  .^(1/2)).*(a.^(-2).*((-2).*t+2.*t0+2.*a.^2.*((-1).*w.^2).^(1/2)).^2).^( ...
  -1/2).*erfz((1/2).*(a.^(-2).*((-2).*t+2.*t0+2.*a.^2.*((-1).*w.^2).^(1/2)) ...
  .^2).^(1/2)))),(1/4).*exp(1).^((-1).*a.^2.*w.^2+(-2).*t.*((-1).*w.^2).^( ...
  1/2)+(-2).*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^(-2).*abs(a).^(-1).*( ...
  d22.*exp(1).^(a.^2.*w.^2+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^2.* ...
  abs(a)+2.*d22.*exp(1).^(a.^2.*w.^2+2.*t.*((-1).*w.^2).^(1/2)+2.*(t+(-1) ...
  .*t0).*((-1).*w.^2).^(1/2)).*w.^2.*abs(a)+d22.*exp(1).^(a.^2.*w.^2+4.* ...
  t.*((-1).*w.^2).^(1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^2.*abs( ...
  a)+4.*D0.*exp(1).^(a.^2.*w.^2+2.*t.*((-1).*w.^2).^(1/2)+2.*(t+(-1).*t0) ...
  .*((-1).*w.^2).^(1/2)).*t.*w.^2.*abs(a)+(-1).*d11.*exp(1).^(a.^2.*w.^2+ ...
  2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^4.*abs(a)+2.*d11.*exp(1).^( ...
  a.^2.*w.^2+2.*t.*((-1).*w.^2).^(1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^( ...
  1/2)).*w.^4.*abs(a)+(-1).*d11.*exp(1).^(a.^2.*w.^2+4.*t.*((-1).*w.^2).^( ...
  1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^4.*abs(a)+D0.*exp(1).^( ...
  a.^2.*w.^2+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2).* ...
  abs(a)+(-1).*D0.*exp(1).^(a.^2.*w.^2+4.*t.*((-1).*w.^2).^(1/2)+2.*(t+( ...
  -1).*t0).*((-1).*w.^2).^(1/2)).*((-1).*w.^2).^(1/2).*abs(a)+(-2).*d12.* ...
  exp(1).^(a.^2.*w.^2+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^2.*((-1).* ...
  w.^2).^(1/2).*abs(a)+2.*d12.*exp(1).^(a.^2.*w.^2+4.*t.*((-1).*w.^2).^( ...
  1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^2.*((-1).*w.^2).^(1/2).* ...
  abs(a)+2.*a.*b.*exp(1).^(a.^2.*w.^2+2.*t.*((-1).*w.^2).^(1/2)+2.*(t+(-1) ...
  .*t0).*((-1).*w.^2).^(1/2)).*w.^2.*erfz(a.^(-1).*t0)+(-2).*a.*b.*exp(1) ...
  .^(a.^2.*w.^2+2.*t.*((-1).*w.^2).^(1/2)+2.*(t+(-1).*t0).*((-1).*w.^2).^( ...
  1/2)).*w.^2.*erfz((1/2).*a.^(-1).*((-2).*t+2.*t0))+a.*b.*exp(1).^(2.*t.*( ...
  (-1).*w.^2).^(1/2)+4.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^2.*erfz(a.^( ...
  -1).*t0+(-1).*a.*((-1).*w.^2).^(1/2))+a.*b.*exp(1).^(2.*t.*((-1).*w.^2) ...
  .^(1/2)).*w.^2.*erfz(a.^(-1).*t0+a.*((-1).*w.^2).^(1/2))+a.*b.*exp(1).^( ...
  2.*t.*((-1).*w.^2).^(1/2)+4.*(t+(-1).*t0).*((-1).*w.^2).^(1/2)).*w.^2.* ...
  erfz(a.^(-1).*(t+(-1).*t0+a.^2.*((-1).*w.^2).^(1/2)))+(-1).*a.*b.*exp(1) ...
  .^(2.*t.*((-1).*w.^2).^(1/2)).*w.^2.*erfz(a.^(-1).*((-1).*t+t0+a.^2.*(( ...
  -1).*w.^2).^(1/2))))];
end