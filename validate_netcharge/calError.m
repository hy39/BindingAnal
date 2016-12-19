function [ plb pub ] = calError( sizeage, pbin )
%CALERROR Summary of this function goes here
%   Detailed explanation goes here
  y0 = pbin;
  options = optimoptions('fsolve','Display','off'); % Option to display output
  x1 = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
  plb = x1;
  x2 = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
  pub = x2;
end

