function [Position,Result]=PsoProcess(SwarmSize,ParticleSize,ParticleScope,Maxgen,Count,Var,Int,Filter)
%[ParSwarm,OptSwarm,BadSwarm]=InitSwarm(SwarmSize,ParticleSize,ParticleScope,AdaptFunc)
 %
 %输入参数：SwarmSize:种群大小的个数
 %输入参数：ParticleSize：一个粒子的维数,变量按01变量，整形变量，非整型变量顺序排列；
 %输入参数：ParticleScope:一个粒子在运算中各维的范围；
 %　　ParticleScope格式:
 %　　3维粒子的ParticleScope格式:
 %                               [x1Min,x1Max]
 %　　　　　　　　　　　　　　　　  x2Min,x2Max
 %                                x3Min,x3Max]
 %
 %输入参数：AdaptFunc：适应度函数
 %输入参数：Maxgen：迭代次数
 %输入参数：Num：重复计算次数
 %输入参数：Var：01变量的个数
 %输入参数：Int：整形变量的个数    
 %输入参数：Filter：筛选因子，根据适应度函数确定
 %输出：Result 迭代过程中的全局最优解
 %输出：Particle 全局最优解的位置

 Result=zeros(Count,Maxgen);
 Position=zeros(Count,ParticleSize);
for i=1:Count
   [ParSwarm,OptSwarm,MaxV]=InitSwarm(SwarmSize,ParticleSize,ParticleScope,Var,Int,Filter);
   for j=1:Maxgen
       [ParSwarm,OptSwarm,MaxValue]=BaseStepPso(ParSwarm,OptSwarm,ParticleScope,Int,Var,MaxV,j,Maxgen);
       Result(i,j)=MaxValue;
       %{
       %如果所有粒子都聚集在一起或者适应度的差小于一个优化过程可以接受的小数，结束迭代
       %提高运算效率,但也容易陷入局部最优解
       Delta=0.01;
       Logic=1;
       for k=1:SwarmSize
           if (MaxValue-ParSwarm(k,2*ParticleSize+1))<=Delta
               Logic=Logic && 1;
           else
               Logic=Logic && 0;
               break
           end
       end
       if Logic
           Result(i,j+1:Maxgen)=MaxValue;
           break
       end
       %}
   end
   %每次重复计算的最优位置
   Position(i,:)=OptSwarm(SwarmSize+1,:);  
end

%输出结果
%disp(Result);
%disp(Particle)
end















