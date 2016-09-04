function [ParSwarm,OptSwarm,maxValue]=BaseStepPso(ParSwarm,OptSwarm,ParticleScope,Integer,MaxV,MinV,CurCount,LoopCount)
%功能描述：全局版本：基本的粒子群算法的单步更新位置,速度的算法
%
%[ParSwarm,OptSwarm]=BaseStepPso(ParSwarm,OptSwarm,AdaptFunc,ParticleScope,MaxW,MinW,LoopCount,CurCount)
%
%输入参数：ParSwarm:粒子群矩阵，包含粒子的位置，速度与当前的目标函数值，其中前Integer个是整数变量，其它为非整数变量；
%输入参数：OptSwarm：包含粒子群个体最优解与全局最优解的矩阵
%输入参数：ParticleScope:一个粒子在运算中各维的范围；
%输入参数：AdaptFunc：适应度函数
%输入参数：Integer: 变量（粒子维数）中整数变量的个数
%输入参数:MaxV最大速度,MinV最小速度,
%输入参数:CurCount当前迭代次数，LoopCount迭代总数
%返回值：ParSwarm OptSwarm
%

%开始单步更新的操作
%得到粒子群群体大小以及一个粒子维数的信息
[ParRow,ParCol]=size(ParSwarm);
%得到粒子的维数
ParCol=(ParCol-1)/2;


%*********************************************
%*****更改下面的代码，可以更改惯性因子的变化*****
%****不同问题，好的惯性因子也不同***************
%{
%线形递减策略
w=zeros(1,ParCol);
for i=1:ParCol
   w(1,i)=MaxV(i)-CurCount*((MaxV(i)-MinV(i))/LoopCount);
end
%}


%w固定不变策略
w=zeros(1,ParCol);
for i=1:ParCol
   w(1,i)=0.7;
end


%{
%w非线形递减，以凹函数递减
w=zeros(1,ParCol);
for i=1:ParCol
    w(1,i)=(MaxV(i)-MinV(i))*(CurCount/LoopCount)^2+(MinV(i)-MaxV(i))*(2*CurCount/LoopCount)+MaxV(i);
end
%}

%{
%w非线形递减，以凹函数递减
w=zeros(1,ParCol);
for i=1:ParCol
    w(1,i)=MinV(i)*(MaxV(i)/MinV(i))^(1/(1+10*CurCount/LoopCount));
end
%}
%*********************************************


SubTract1=OptSwarm(1:ParRow,:)-ParSwarm(:,1:ParCol);

%*********************************************
%*****更改下面的代码，可以更改c1,c2的变化*****
c1=2;
c2=2;
%
%con=1;
%c1=4-exp(-con*abs(mean(ParSwarm(:,2*ParCol+1))-AdaptFunc(OptSwarm(ParRow+1,:))));
%c2=4-c1;
%
%
%*********************************************

TempV=zeros(ParRow,ParCol);
TempPos=zeros(ParRow,ParCol);
for row=1:ParRow
    SubTract2=OptSwarm(ParRow+1,:)-ParSwarm(row,1:ParCol);
    for col=1:ParCol
        TempV(row,col)=w(1,col)*ParSwarm(row,ParCol+col)+c1*unifrnd(0,1).*SubTract1(row,col)+c2*unifrnd(0,1).*SubTract2(col);
    end
    for h=1:ParCol
        if TempV(row,h)>0.2*(ParticleScope(h,2)-ParticleScope(h,1));
            TempV(row,h)=0.2*(ParticleScope(h,2)-ParticleScope(h,1));
        end
        if TempV(row,h)<-0.2*(ParticleScope(h,2)-ParticleScope(h,1));
            TempV(row,h)=-0.2*(ParticleScope(h,2)-ParticleScope(h,1));
            %加1e-10防止适应度函数被零除
        end
    end
    %限制速度的代码;
    ParSwarm(row,ParCol+1:2*ParCol)=TempV(row,:);
    %更新速度
     
     %*********************************************
     %*****更改下面的代码，可以更改约束因子的变化*****
     %
     %a=1;
     %
     a=0.729;
     %*********************************************
     for col=1:Integer
         TempPos(row,col)=floor(ParSwarm(row,col)+a*TempV(row,col));
     end
     for col=Integer+1:ParCol
         TempPos(row,col)=ParSwarm(row,col)+a*TempV(row,col);
     end
     %位置更新；
      for h=1:ParCol
                if TempPos(row,h)>ParticleScope(h,2)
                     TempPos(row,h)=ParticleScope(h,2);
                end
                if TempPos(row,h)<=ParticleScope(h,1)
                    TempPos(row,h)=ParticleScope(h,1);
                end
      end
      
      %不等式约束，用罚函数法代替；
      %{
         n=1;
         while(n<100)
            %位置的常数限制
             if (TempPos(row,1)+TempPos(row,2)+TempPos(row,3)+TempPos(row,4)+TempPos(row,5)<=400)&&((TempPos(row,1)+2*TempPos(row,2)+2*TempPos(row,3)+TempPos(row,4)+6*TempPos(row,5)<=800))&&((2*TempPos(row,1)+TempPos(row,2)+6*TempPos(row,3)<=200))&&((TempPos(row,3)+TempPos(row,4)+5*TempPos(row,5)<=200))
             %不等式约束
             break;
             else
                 random=rand;
                 for i=1:ParCol
                      TempPos(row,i)=TempPos(row,i)*random;
                     %位置衰减，个人认为可以看作一个变异；
                 end
             end
             n=n+1;
         end
     %}
      %{
      while(~((TempPos(row,1)+TempPos(row,2)+TempPos(row,3)+TempPos(row,4)+TempPos(row,5)<=400) && ((TempPos(row,1)+2*TempPos(row,2)+2*TempPos(row,3)+TempPos(row,4)+6*TempPos(row,5)<=800)) && ((2*TempPos(row,1)+TempPos(row,2)+6*TempPos(row,3)<=200)) && ((TempPos(row,3)+TempPos(row,4)+5*TempPos(row,5)<=200))))
           random=rand;
                 for i=1:ParCol
                     TempPos(row,i)=TempPos(row,i)*random;
                     %位置衰减，个人认为可以看作一个变异；
                 end
       end
      %} 
      
     ParSwarm(row,1:ParCol)=TempPos(row,:);
     %更新位置

     ParSwarm(row,2*ParCol+1)=AdaptFunc(ParSwarm(row,1:ParCol));
     %计算每个粒子的新的适应度值
     if ParSwarm(row,2*ParCol+1)>AdaptFunc(OptSwarm(row,1:ParCol))
         OptSwarm(row,1:ParCol)=ParSwarm(row,1:ParCol);
     end
     %每个粒子适应度值得更新；
end
%for循环结束

%寻找适应度函数值最大的解在矩阵中的位置(行数)，进行全局最优的改变 
[maxValue,row]=max(ParSwarm(:,2*ParCol+1));
if  maxValue>AdaptFunc(OptSwarm(ParRow+1,:))
    OptSwarm(ParRow+1,:)=ParSwarm(row,1:ParCol);
end
end