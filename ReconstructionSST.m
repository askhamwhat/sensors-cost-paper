% Find sparse sensor locations for reconstruction of sea surface temps,
% incorporating a cost function.

clear all; close all; clc

sst = ncread('sst.wkmean.1990-present.nc','sst');
mask = ncread('lsmask.nc','mask');

[m,n] = size(mask);
N = m*n;
L = 1400;
r = 1100; % The size of the training set

for j = 1:L
    sst_mask(:,:,j) = sst(:,:,j).*mask;
    Z(:,j) = reshape(sst_mask(:,:,j),N,1); % The data matrix
end

% Construct a cost function
f2 = ones(m,n);
for j = 3:357
   for jj = 3:177
      if mask(j,jj) == 0
          f2(j-2:j+2,jj-2:jj+2) = 0;
      end
   end
end
f = reshape(f2,N,1);

Gamma = [0:10:100,125:25:225]; % Cost function weightings
p = [5,10,25,50,75,100:50:300]; % The number of sensors
LCV = 10; % The number of cross validations

% Obtain sensor locations and reconstructions
[Cost,Error,A,stdE,stdC] = CostError(Z,f,Gamma,p,r,LCV);


% Plotting
% Plot cost function and sensor locations, as in Figure 8
cmap = [0,0,0;
    0,0,1;
    0,0,1;
    0,0.9,0;
    0,0.9,0;
    1,1,0;
    1,1,0;
    1,0.5,0;
    1,0.5,0;
    1,0.2,0;
    1,0.2,0;
    1,1,1];

num = 0;
figure(1)
set(gcf,'Position',get(0,'Screensize'));
% Sensor locations
for k = [1,11,16]
    num = num+1;
    A2 = reshape(A(:,k,8),360,180)';
    x = [];
    y = [];
    msize = [];
    mcol = [];
    for i = 1:180
        for j = 1:360
            if A2(i,j) ~= 0
                x = [x; j];
                y = [y; i];
                msize = [msize; A2(i,j)];
                mcol = [mcol; cmap(10*A2(i,j)+1,:)];
            end
        end
    end
    
    subplot(2,2,num)
    imagesc(mask');
    shading interp
    colormap(cmap)
    hold on
    scatter(x,y,10*msize,mcol,'filled','MarkerEdgeColor',[0,0,0.5],...
        'LineWidth',0.9);
    set(gca,'YDir','reverse','FontSize',12)
    axis equal off
    axis off
    title(['\gamma = ',num2str(Gamma(k))],'FontName','Times')
end

% Cost and error versus gamma
E200 = Error(:,8);
C200 = Cost(:,8);
stdE200 = stdE(:,8);
stdC200 = stdC(:,8);

EU = E200 + stdE200; EL = E200 - stdE200;
CU = C200 + stdC200; CL = C200 - stdC200;
Eerr = [EU; flipud(EL)]; Cerr = [CU; flipud(CL)];
Gamma2 = [Gamma, fliplr(Gamma)];

xx = 100*ones(1,101);
yy = 0:100;

subplot(2,2,4)
yyaxis right
fill(Gamma2,Eerr,[1,0.5,0.5])
hold on
plot(Gamma,E200,'r-','LineWidth',1.5)
h = gca; h.YColor = 'r';
ylim([0 max(EU)])
alpha(0.5)
ylabel('Error','Color','r','FontName','Times')

yyaxis left
fill(Gamma2,Cerr,[0.5,0.5,0.5])
hold on
plot(Gamma,C200,'k-','LineWidth',1.5)
h = gca; h.YColor = 'k';
ylim([0 max(CU)])
alpha(0.5)
set(gca,'FontSize',12)
xlabel('\gamma','FontName','Times')
ylabel('Cost','FontName','Times')
hold on
plot(xx,yy,'-','Color',[0.5,0.5,0.5],'LineWidth',1.5)


% Cost colormap plot, as in Figure 9
p2 = repmat(p,size(Error,1),1);

figure(2)
pcolor(p2,Error,Cost)
shading interp, colormap jet
caxis([-50 max(max(Cost))])
xlim([0,300])
set(gca,'FontSize',12);
xlabel('# of Sensors','FontName','Times')
ylabel('Error','FontName','Times')
axis square
hold on
[Cont,Height] = contour(p2,Error,Cost,[1,50],...
    'k','LineWidth',1.5);
clabel(Cont,Height,'FontSize',12,'FontWeight','bold','LabelSpacing',150)

clrbar = colorbar;
clrbar.Limits = [0 max(max(Cost))];
clrbar.Label.String = 'Cost'; clrbar.Label.FontName = 'Times';
clrbar.Label.FontSize = 12;


% Cost versus error, as in Figure 10
figure(3)
plot(Error(:,3),Cost(:,3),'.','Color',[0,0.8 0],'MarkerSize',12)
hold on
plot(Error(:,6),Cost(:,6),'.','Color',[1 0.8 0],'MarkerSize',12)
plot(Error(:,8),Cost(:,8),'.','Color',[1 0.5 0],'MarkerSize',12)
% plot(Error(:,10),Cost(:,10),'.','Color',[0.9,0,0],'MarkerSize',12)
set(gca,'FontSize',12)
xlabel('Error','FontName','Times')
ylabel('Cost','FontName','Times')
Leg = legend('25 Sensors','100 Sensors','200 Sensors');%,'300 Sensors');
Leg.FontSize = 12; Leg.FontName = 'Times';









