% Find sparse sensor locations for reconstruction of cylinder flow,
% incorporating a cost function.

clear all; close all; clc

cyl = load('CYLINDER_ALL.mat');

Z = cyl.VORTALL;
m = cyl.m;
n = cyl.n;
N = m*n;
L = size(Z,2);
r = 120; % The size of the training set

% Construct a cost function
f2 = zeros(m,n);
f2(1:99,:) = 1;
f = reshape(f2,N,1);

p = 2:2:20; % The number of sensors
Gamma = [0:0.1:1,1.5:0.5:5,6:1:15]; % The cost function weightings
LCV = 30; % The number of cross validations

% Obtain sensor locations and reconstructions
[Cost,Error,A,stdE,stdC] = CostError(Z,f,Gamma,p,r,LCV);


% Plotting
% Sensor locations, as in Figure 11
vortmin = -5;  % only plot what is in -5 to 5 range
vortmax = 5;
Z(Z>vortmax) = vortmax;  % cutoff at vortmax
Z(Z<vortmin) = vortmin;

theta = (1:100)/100'*2*pi;
xx = 49+25*sin(theta);
yy = 99+25*cos(theta);

vort12 = reshape(Z(:,12),m,n);
kk = [1,11,29];
% Sensor locations
figure(1)
set(gcf,'Position',get(0,'Screensize'));
for k = 1:3
    subplot(2,2,k)
    pcolor(vort12(25:175,1:340))
    shading interp
    unfreezeColors
    colormap gray
    freezeColors
    hold on
    fill(xx,(yy-25),[.3 .3 .3])
    plot(xx,(yy-25),'k','LineWidth',1.5)
    axis equal off
    
    AQR = reshape(A(:,kk(k),7),m,n);
    x = [];
    y = [];
    msize = [];
    mcol = [];
    clearvars msize2 msize3 msize4
    for i = 1:m
        for j = 1:n
            if AQR(i,j) ~= 0
                x = [x; j];
                y = [y; i];
                msize = [msize; AQR(i,j)];
            end
        end
    end
    ind = find(msize == 1);
    msize2 = [msize(1:ind-1); msize(ind+1:end)];
    val1 = min(msize2);
    val2 = max(msize2);
    msize3 = [msize(1:ind-1); val2; msize(ind+1:end)];
    
    for j = 1:length(msize3)
        if msize3(j) == val1
            msize4(j) = vortmin;
        elseif msize3(j) == val2
            msize4(j) = vortmax;
        else
            msize4(j) = msize3(j)*10*(val2-val1);
        end
        if msize3(j) > 0.3
            msize3(j) = msize3(j)/2;
        end
    end
    
    hold on
    scatter(x,(y-25),30*msize3,msize4,'filled','MarkerEdgeColor','k')
    colormap hot
    
    title(['\gamma = ',num2str(Gamma(kk(k)))],'FontName','Times',...
        'Position',[170.5,155,0])
end

% Cost and error versus gamma
C14 = Cost(:,7);
E14 = Error(:,7);
stdE14 = stdE(:,7);
stdC14 = stdC(:,7);

EU = E14 + stdE14; EL = E14 - stdE14;
CU = C14 + stdC14; CL = C14 - stdC14;
Eerr = [EU; flipud(EL)]; Cerr = [CU; flipud(CL)];
Gamma2 = [Gamma, fliplr(Gamma)];

xvert = Gamma(11)*ones(1,100);
yvert = linspace(-0.1,10,100);

subplot(2,2,4)
yyaxis right
fill(Gamma2,Eerr,[1,0.5,0.5])
hold on
plot(Gamma,E14,'r-','LineWidth',1.5)
ylabel('Error','FontName','Times','Color','r')
h = gca; h.YColor = 'r';
ylim([0 max(EU)])
set(gca,'FontSize',12)
alpha(0.5)

yyaxis left
fill(Gamma2,Cerr,[0.5,0.5,0.5])
hold on
plot(Gamma,C14,'k-','LineWidth',1.5)
h = gca; h.YColor = 'k';
set(gca,'FontSize',12)
xlabel('\gamma','FontName','Times')
ylabel('Cost','FontName','Times')
hold on
plot(xvert,yvert,'-','Color',[0.5,0.5,0.5],'LineWidth',1.5)
ylim([0 max(CU)])
alpha(0.5)


% Cost landscape, as in Figure 12
p2 = repmat(p,size(Error,1),1);
figure(2)
pcolor(p2,Error,Cost)
shading interp, colormap gray
hold on
p1 = plot(p,Error(1,:),'r--','LineWidth',1.5);
p2 = plot(p,Error(29,:),'b-.','LineWidth',1.5);
caxis([-15 max(max(Cost))])
set(gca,'FontSize',12);
xlabel('# of Sensors','FontName','Times')
ylabel('Error','FontName','Times')
axis square
Leg = legend([p1,p2],'\gamma = 0','\gamma = 15');
Leg.FontSize = 12; Leg.FontName = 'Times';


% Cost versus Error, as in Figure 13
figure(3)
plot(Error(:,2),Cost(:,2),'.','Color',[0,0.8 0],'MarkerSize',12)
hold on
plot(Error(:,4),Cost(:,4),'.','Color',[1 0.8 0],'MarkerSize',12)
plot(Error(:,6),Cost(:,6),'.','Color',[1 0.5 0],'MarkerSize',12)
plot(Error(:,10),Cost(:,10),'.','Color',[0.9 0 0],'MarkerSize',12)
set(gca,'FontSize',12)
xlabel('Error','FontName','Times')
ylabel('Cost','FontName','Times')
Leg = legend('4 Sensors','8 Sensors','12 Sensors','20 Sensors');
Leg.FontSize = 12; Leg.FontName = 'Times';





