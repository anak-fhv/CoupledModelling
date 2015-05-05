clear all
close all
clc

currDir = pwd;
fname = 'lmod_scrPolBins.res';
fpath = strcat(pwd,'/',fname);
fid = fopen(fpath,'r');
nbin = 20;
nrays = 1e6;
% pray = 1.0/nrays;
pray = 0.15/1e6;
p = fscanf(fid,'%f',[nbin,2*nbin]);
p = p';
%--------power values------------%
p = p.*pray;
%--------power values------------%
% s = p./sum(sum(p));
s = p;
bY = s(1:2:2*nbin-1,:);
bB = s(2:2:2*nbin,:);
bYov = sum(bY,2);
bBov = sum(bB,2);
pov = bYov+bBov;
ybr = bYov./bBov;
theta = linspace(0,pi/2,nbin+1);
sa = zeros(nbin,1);
cs = zeros(nbin,1);
for i=1:nbin
    sa(i) = (cos(theta(i)) - cos(theta(i+1)))*2*pi;
    cs(i) = cos((theta(i+1)+theta(i))/2);
end
th = linspace(0,pi/2,nbin)';
p1 = bYov./sa;
p1Int = bYov./(sa.*cs);
p2 = bBov./sa;
p2Int = bBov./(sa.*cs);
pC = pov./sa;
pCInt = pov./(sa.*cs);

p1Norm = p1./max(pC);
p2Norm = p2./max(pC);

pCNorm = pC./max(pC);

rYB = p1./p2;

figure(1);
polar(th,pCNorm,'-k');
hold on
polar(th,p1Norm,'--k');
polar(th,p2Norm,'-.k');
hold off
view([90,-90]);
xlim([0,Inf]);
ylim([0,Inf]);

figure(2);
plot(th.*180/pi,pCNorm,'-k');
hold on
plot(th.*180/pi,p1Norm,'--k');
plot(th.*180/pi,p2Norm,'-.k');
plot(th.*180/pi,rYB,':k');
xlim([0,90]);
xlabel('$\theta$ (-)','Interpreter','Latex');
ylabel('$q" \: (-),\: r_\mathrm{y/b} (-)$','Interpreter','Latex');

set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|:|-.|-')
set(0,'DefaulttextFontName','Times','DefaulttextFontSize',10)
set(0,'DefaultAxesUnits','normalized')
figHandles = findall(0,'Type','figure');
set(figHandles,'PaperUnits','centimeters','PaperSize',[10 10],'PaperPosition',[0 0 10 10]);
allAxes = findall(0,'type','axes');
set(allAxes,'FontName','Times','FontSize',10,'Units','centimeters','Position',[1.5 1.5 8, 8]);
allText = findall(figHandles,'Type','text');
set(allText,'FontName','Times','FontSize',10);

figure(1);
legend({'$q"_\mathrm{tot}$ (-)','$q"_\mathrm{y}$ (-)','$q"_\mathrm{b}$ (-)'},...
    'Interpreter','Latex','Location','best');
legend('boxoff');

figure(2);
legend({'$q"_\mathrm{tot}$ (-)','$q"_\mathrm{y}$ (-)','$q"_\mathrm{b}$ (-)',...
    '$r_\mathrm{y/b}$'},'Interpreter','Latex','Location','best');
legend('boxoff');

figure(1);
print -deps qPolar.eps

figure(2);
print -deps qLinear.eps