
function     y=AdaptFunc(x)
%适应度函数（罚函数法） 
%输入参数：x 第i个粒子的位置参数 ParSwarm(i,1:ParticleSize)；
%输出参数：y 第i个粒子的适应度值，放入ParSwarm(i,2*ParticleSize+1)；
z=100000*((max((x(1)+x(2)+x(4)+x(5)-400),0)).^2+max((2*x(1)+x(2)+6*x(3)-100),0).^2+max((x(3)+x(4)+5*x(5)-200),0).^2);
y=x(1).^2+x(2).^2-3*x(3).^2+4*x(4).^2+2*x(5).^2-z;
end