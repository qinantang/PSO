function [Particle,Result]=PsoProcess(SwarmSize,ParticleSize,ParticleScope,Maxgen,Integer)
%[ParSwarm,OptSwarm,BadSwarm]=InitSwarm(SwarmSize,ParticleSize,ParticleScope,AdaptFunc)
 %
 %输入参数：SwarmSize:种群大小的个数
 %输入参数：ParticleSize：一个粒子的维数,其中前integer个是整数变量，其它为非整数变量；
 %输入参数：ParticleScope:一个粒子在运算中各维的范围；
 %　　ParticleScope格式:
 %　　3维粒子的ParticleScope格式:
 %                               [x1Min,x1Max]
 %　　　　　　　　　　　　　　　　  x2Min,x2Max
 %                                x3Min,x3Max]
 %
 %输入参数：AdaptFunc：适应度函数
 %输出：Result 迭代过程中的全局最优解
 %输出：Particle 全局最优解的位置

 Result=zeros(1,Maxgen);
 Particle=zeros(1,ParticleSize);
[ParSwarm,OptSwarm,MaxV,MinV]=InitSwarm(SwarmSize,ParticleSize,ParticleScope,Integer);
for i=1:Maxgen
    [ParSwarm,OptSwarm,maxValue]=BaseStepPso(ParSwarm,OptSwarm,ParticleScope,Integer,MaxV,MinV,i,Maxgen);
    Result(1,i)=maxValue;
end
Particle=OptSwarm(SwarmSize+1,:);
plot(Result);
disp(Result(1,Maxgen));
end















