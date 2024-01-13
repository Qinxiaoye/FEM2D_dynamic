function elemM = elemMass(ndim,mnode,ro,coor)
% 得到单元质量阵
% ndim - 维数，mnode-单元类型，ro-密度，coor-单元节点的坐标
% coor - [x,y,z]
% 单维方向积分点个数为3
x = coor(:,1);
y = coor(:,2);

nip = 3; % 单维方向的积分点数
ks = [-0.774596669241483,0,0.774596669241483]; % 积分点
w = [0.55555555555556,0.888888888888888889,0.55555555555556]; % 权系数

elemM = zeros(mnode*ndim,mnode*ndim);
if ndim == 2
    for m = 1:nip
        for n = 1:nip
            N = fun(ks(m),ks(n),0,ndim,mnode);
            if mnode == 4
                N = [N(1),0,N(2),0,N(3),0,N(4),0;0,N(1),0,N(2),0,N(3),0,N(4)];
            else
                N = [N(1),0,N(2),0,N(3),0,N(4),0,N(5),0,N(6),0,N(7),0,N(8),0;...
                    0,N(1),0,N(2),0,N(3),0,N(4),0,N(5),0,N(6),0,N(7),0,N(8)];
            end
            [N_ks,N_yt] = dfun2D(ks(m),ks(n),mnode);
            J = jacobi2D(x,y,N_ks,N_yt);
            elemM = elemM + w(m)*w(n)*(N'*N)*det(J)*ro;
        end
    end
end

