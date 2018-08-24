% Find sparse sensor locations for reconstruction of eigenfaces,
% incorporating a cost function.

clear all; close all; clc

Z = load('YaleB_32x32.mat');
Z = Z.fea';

[N,L] = size(Z);
n = 32;
r = round(L*.80); % The size of the test set

% Construct the cost function (uncomment your choice)
% %Gaussian
[X1,X2] = meshgrid(-n/2:n/2-1,-n/2:n/2-1);
f2 = exp(-0.05*X1.^2 - 0.05*X2.^2);
f = reshape(f2,N,1);
Gamma = 10^4*[0:0.1:1,1.25,1.5:0.5:4,5:2:19]; % Cost function weightings
Gplot = [1,17,26]; Cmin = -20; Contours = [1 25];

% %Center blocked
% f2 = zeros(n,n);
% f2(11:22,11:22) = 1;
% f = reshape(f2,N,1);
% Gamma = 10^4*[0:0.05:0.5,0.6:0.1:1,1.5:0.5:6];
% Gplot = [1,8,26]; Cmin = -40; Contours = [1 50];

% %Left blocked
% f2 = zeros(n,n);
% f2(:,1:11) = 1;
% f = reshape(f2,N,1);
% Gamma = 10^4*[0:0.05:0.5,0.6:0.1:1,1.5:0.5:6];
% Gplot = [1,17,26]; Cmin = -40; Contours = [1 50];


p = [5,10,15,25,50,75,100:50:400]; % The number of sensors
LCV = 20; % The number of cross validations

% Obtain sensor locations, reconstruction errors, and costs
[Cost,Error,A,stdE,stdC] = CostError(Z,f,Gamma,p,r,LCV);


% Plotting
% Plot cost function and sensor locations, as in Figure 4
[U,Sig,V] = svd(Z,'econ');
EigFace = U(:,1);

% Cost function
figure(1)
subplot(1,5,1)
pcolor(f2), shading interp, colormap hot
axis square off
set(gca,'FontSize',12)
set(gcf,'Position',get(0,'Screensize'))
title('Cost Function','FontName','Times')

% Sensor locations
for k = 1:3
    subplot(1,5,k+1)
    imagesc(reshape(-EigFace,32,32))
    axis square off
    unfreezeColors
    colormap gray
    freezeColors
    hold on
    
    A2 = reshape(A(:,Gplot(k),9),32,32);
    x = [];
    y = [];
    msize = [];
    clearvars msize2 msize3 msize4
    for i = 1:n
        for j = 1:n
            if A2(i,j) ~= 0
                x = [x; j];
                y = [y; i];
                msize = [msize; A2(i,j)];
            end
        end
    end
    
    hold on
    scatter(x,y,20*msize,msize/50,'filled','MarkerEdgeColor','k',...
        'LineWidth',1)
    colormap hot
    
    title(['\gamma = ',num2str(Gamma(Gplot(k)))],'FontName','Times')
    if k > 1
        title(['\gamma = ',num2str(Gamma(Gplot(k)),'%0.1E')],'FontName',...
            'Times')
    end
    set(gca,'FontSize',12)
end

% Cost and error versus Gamma
E200 = Error(:,9);
C200 = Cost(:,9);
DevE200 = stdE(:,9);
DevC200 = stdC(:,9);

EU = E200 + DevE200;
EL = E200 - DevE200;
CU = C200 + DevC200;
CL = C200 - DevC200;
Eerr = [EU; flipud(EL)];
Cerr = [CU; flipud(CL)];
Gamma2 = [Gamma, fliplr(Gamma)];

xx = (Gplot(2))*ones(100,1);
yy = linspace(0,max(C200)+10,100);

subplot(1,5,5),hold on
yyaxis left
fill(Gamma2,Cerr,[0.5,0.5,0.5])
hold on
p1 = plot(Gamma,C200,'k-','LineWidth',1.5);
h = gca; h.YColor = 'k';
set(gca,'FontSize',12)
xlabel('\gamma / 10^4','FontName','Times')
ylabel('Cost','FontName','Times')
xlim([0 max(Gamma)])
ylim([0 max(C200)+10])
hold on
plot(xx,yy,'-','Color',[0.5,0.5,0.5],'LineWidth',1.5)
alpha(0.5)
axis square

yyaxis right
fill(Gamma2,Eerr,[1,0.5,0.5])
p2 = plot(Gamma,E200,'r-','LineWidth',1.5);
h = gca; h.YColor = 'r';
ylabel('Error','FontName','Times','Color','r')
box on


% Cost colormap plot, as in Figure 5
p2 = repmat(p,size(Error,1),1);

figure(2)
pcolor(p2,Error,Cost)
shading interp, colormap jet
caxis([Cmin max(max(Cost))])
xlim([0,400])
set(gca,'FontSize',12);
xlabel('# of Sensors','FontName','Times')
ylabel('Error','FontName','Times')
axis square
hold on
[Cont,Height] = contour(p2,Error,Cost,Contours,...
    'k','LineWidth',1.5);
clabel(Cont,Height,'FontSize',12,'FontWeight','bold')

clrbar = colorbar;
clrbar.Limits = [0 max(max(Cost))];
clrbar.Label.String = 'Cost'; clrbar.Label.FontName = 'Times';
clrbar.Label.FontSize = 12;


% Plot cost versus error for several numbers of sensors, as in Figure 7
figure(3)
plot(Error(:,4),Cost(:,4),'.','Color',[0,0.8 0],'MarkerSize',12)
hold on
plot(Error(:,7),Cost(:,7),'.','Color',[1 0.8 0],'MarkerSize',12)
plot(Error(:,9),Cost(:,9),'.','Color',[1 0.5 0],'MarkerSize',12)
plot(Error(:,11),Cost(:,11),'.','Color',[0.9,0,0],'MarkerSize',12)
set(gca,'FontSize',12)
xlabel('Error','FontName','Times')
ylabel('Cost','FontName','Times')

Leg = legend('25 Sensors','100 Sensors','200 Sensors','300 Sensors');
Leg.FontSize = 12; Leg.FontName = 'Times';














