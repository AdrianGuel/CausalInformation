function Mw0=meansw0(g,t,t20,c,d,x0,v0)
Mw0=[g.^(-1).*(v0+(-1).*exp(1).^((-1).*g.*t).*v0)+x0+(1/2).*d.*g.^(-1).*( ...
  erfz(c.^(-1).*(t+(-1).*t20))+erfz(c.^(-1).*t20)+exp(1).^((1/4).*g.*(c.^2.* ...
  g+(-4).*t+4.*t20)).*(erfz((1/2).*c.^(-1).*(c.^2.*g+(-2).*t+2.*t20))+(-1) ...
  .*erfz((1/2).*c.*g+c.^(-1).*t20))).*sign(c);exp(1).^((-1).*g.*t).*v0+( ...
  1/2).*d.*exp(1).^((1/4).*g.*(c.^2.*g+(-4).*t+4.*t20)).*((-1).*erfz((1/2) ...
  .*c.^(-1).*(c.^2.*g+(-2).*t+2.*t20))+erfz((1/2).*c.*g+c.^(-1).*t20)).* ...
  sign(c)];
end