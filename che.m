T = [145.10, 134.40, 129.30, 123.20, 118.30, 116.20, 113.10, 108.20, 104.40, 102.10, 101.20];
x1_0ri = [0.0850, 0.1760, 0.2520, 0.3510, 0.4390, 0.4820, 0.5460, 0.6760, 0.8000, 0.8900, 0.9290];
y1_0ri = [0.4940, 0.6860, 0.7560, 0.8230, 0.8650, 0.8750, 0.8980, 0.9350, 0.9630, 0.9810, 0.9870];
Y1= @(x2) exp((x2)^2*(0.18*((0.95/((1-x2)+(x2)*0.95))^2)+0.32*0.91/(((x2)+(1-x2)*0.91))^2));
Y2= @(x2) exp((1-x2)^2*(0.32*((0.91/((1-x2)+(x2)*0.91))^2)+0.18*0.95/(((1-x2)+(x2)*0.95)^2)));
Antoine= @(T,A,B,C) 10^(A-B/(T+C));
x2 = (11); % mole fraction N,N Dimethylacetamide in liquid
res = (11); % mole fraction of water in liquid
pes = (11); % mole fraction of water in vapor
%solving for x2 and x1 by applying modified rault's law
for i = 1:11
P1=Antoine(T(i),8.01770,1715.700,234.268); % vapour pressure of pure water as function of temperature
P2=Antoine(T(i),7.76228,1889.100,221.000); % vapour pressure of pure N,N Dimethylacetamide
NRTL = @(x2) 760 - P1*Y1(x2)*(1-x2) - P2*Y2(x2)*(x2);
x2(i) = fsolve(NRTL,0.5);
res(i) = 1-x2(i);
%disp(res(i));
end
%solving for y1(pes:- mole fraction of water in liquid phase) 
for i = 1:11
p1 = Antoine(T(i),8.01770,1715.700,234.268); 
EQY = @(y1) y1 - ((p1*Y1(x2(i))*res(i))/760);
pes(i) = fsolve(EQY,0.5);
%disp(pes(i));
end
tiledlayout(3,1)
nexttile
plot(res,T);
hold on
plot(x1_0ri,T);
nexttile
plot(pes,T);
hold on
plot(y1_0ri,T);
nexttile
plot(res,pes);
hold on
plot(x1_0ri,y1_0ri);
% fsolve(@F,)
% x1(i)=;