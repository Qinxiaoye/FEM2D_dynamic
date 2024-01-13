function x = NewMark(alpha,delta,K,M,C,F,sumTime,dt,ndim)
% New-Mark 逐步积分解方程
% F为外力，不能随时间变化


sumStep = fix(sumTime/dt);
sumNode = size(K,1)/ndim;

%--------------------------初始计算----------------------------------------
x = zeros(sumNode*ndim,sumStep);
xn = sparse(zeros(sumNode*ndim,1));
vn = xn;
% 计算a(0)----a为加速度
an = M\(F-C*vn-K*xn);

% 计算不变的参数
a0 = 1/(alpha*dt^2);
a1 = delta/(alpha*dt);
a2 = 1/(alpha*dt);
a3 = 1/(2*alpha)-1;
a4 = delta/alpha-1;
a5 = dt/2*(delta/alpha-2);
a6 = dt*(1-delta);
a7 = delta*dt;


% 计算等效刚度阵K
K = K + a0*M + a1*C;
% ---------------------每一个时间步计算------------------------------------

for n = 1:sumStep
    f = F + M*(a0*xn+a2*vn+a3*an)+C*(a1*xn+a4*vn+a5*an);  % 等效力

    xnp = K\f;
    anp = a0*(xnp-xn)-a2*vn-a3*an;
    vnp = vn+a6*an+a7*anp;
    an = anp;
    vn = vnp;
    xn = xnp;
    x(:,n) = full(xn);

end



