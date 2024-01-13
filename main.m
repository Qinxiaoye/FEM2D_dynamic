% function  main()
clear;
% 2d dynamic problems, Newmark method
sumTime = 1;
dt = 0.01; % fixed time step

% read the element and boundary condation 
fprintf(1,'read the modal\n')
node = load('NLIST.DAT');
sumNode = size(node,1);
elem = load('ELIST.DAT');
sumElem = size(elem,1);
fixNode = load('fixNode.dat');
nodeForce = load('nodeForce.dat'); % 节点力
nodeForce(1,:) = [];

% -------------------------------------------------------------------------
mat = [200000,0.3,1;]; % 弹性模量、泊松比、密度, 多种材料需要换行
ndim = 2;isStress = 1; % 1-plane stress 0-plane strain

% dynamic parameters-------
alpha1 = 0;alpha2 = 0;  % 瑞利阻尼系数
alpha = 0.25250625;  % alpha and delta
delta = 0.505;

h = 1;   % 应用于二维问题，默认是1

matID  = elem(:,2); % material num
EX = mat(matID,1); % E
mu = mat(matID,2); % nu
ro = mat(matID,3);

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:,1:ndim);

mnode = size(elem,2);  % 单元类型

if mnode == 4
    reduce = 0;  % 是否采用减缩积分，为1表示就采用减缩积分，0表示全积分。
else
    reduce = 1;  
end

% ---------------------------转化压强---------------------------------------
if exist('press.dat','file')
    fprintf(1,'trans the pressure\n')
    nodeForceNew = transPres(node,elem,ndim,mnode);
    nodeForce = [nodeForce;nodeForceNew];
end

% ----------------------------形成总体刚度阵--------------------------------
fprintf(1,'sort the global K\n')

[GK_u,GK_v,GK_a] = globalK2D(EX,mu,h,elem,node,isStress,reduce,mnode);
M = globalM(ro,ndim,mnode,node,elem);
% ----------------------------第一类边界条件-----------------------------------
fprintf(1,'boundary condition\n')

% 对角元改1法施加第一类边界条件
[GK,force] = boundary_chang_one(GK_u,GK_v,GK_a,fixNode,nodeForce,sumNode,ndim);

% ----------------------------求解方程--------------------------------------
fprintf(1,'solve the dynamic equation\n')
% 生成阻尼阵
C = alpha1*M+alpha2*GK; % 瑞利阻尼
x = NewMark(alpha,delta,GK,M,C,force,sumTime,dt,ndim);% x 为非稀疏矩阵

% post
% show the contour at step 20
step = 100;
u = reshape(x(:,step),ndim,[])';

% ux
figure('color',[1 1 1])
axis equal;
showContour(node,elem,u(:,1));
axis off;


% uy
figure('color',[1 1 1])
axis equal;
showContour(node,elem,u(:,2));
axis off;


figure;
plot(dt:dt:sumTime,x(4,:),'LineWidth',2)
