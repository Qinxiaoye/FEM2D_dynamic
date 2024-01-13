function M = globalM(ro,ndim,mnode,node,elem)
% 形成总体质量阵
% 调用elemMass
% 适用于现有的所有情况
sumNode = size(node,1);
sumElem = size(elem,1);

u = zeros(sumElem*(mnode*ndim)*(mnode*ndim),1); % 稀疏矩阵的行-列-值
v = zeros(sumElem*(mnode*ndim)*(mnode*ndim),1);
a = zeros(sumElem*(mnode*ndim)*(mnode*ndim),1);

sumndf = (ndim*mnode)^2;

for n = 1:sumElem
    coor = node(elem(n,:),:);
    
    elemM = elemMass(ndim,mnode,ro(n),coor); % 所得到的单元质量阵
    nodeID = elem(n,:);
    if ndim == 2
        ndfID(1:ndim:mnode*ndim-1) = nodeID*ndim-1;
        ndfID(2:ndim:mnode*ndim) = nodeID*ndim;
    else
        ndfID(1:ndim:mnode*ndim-2) = nodeID*3-2;
        ndfID(2:ndim:mnode*ndim-1) = nodeID*3-1;
        ndfID(3:ndim:mnode*ndim) = nodeID*3;
    end
    ddd = repmat(ndfID,ndim*mnode,1); % ddd 为临时变量，无实际意义
    u(sumndf*(n-1)+1:sumndf*n) = ddd(:);
    v(sumndf*(n-1)+1:sumndf*n) = repmat(ndfID,1,ndim*mnode);
    a(sumndf*(n-1)+1:sumndf*n) = elemM(:);
end

M = sparse(u,v,a,sumNode*ndim,sumNode*ndim);