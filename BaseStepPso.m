function [ParSwarm,OptSwarm,MaxValue]=BaseStepPso(ParSwarm,OptSwarm,ParticleScope,Int,Var,MaxV,CurCount,LoopCount)
%功能描述：全局版本：基本的粒子群算法的单步更新位置,速度的算法
%
%[ParSwarm,OptSwarm]=BaseStepPso(ParSwarm,OptSwarm,AdaptFunc,ParticleScope,MaxV)
%
%输入参数：ParSwarm:粒子群矩阵，包含粒子的位置，速度与当前的目标函数值
%输入参数：OptSwarm：包含粒子群个体最优解与全局最优解的位置
%输入参数：ParticleScope:一个粒子在运算中各维的范围；
%输入参数：AdaptFunc：适应度函数
%输入参数：Int: 变量（粒子维数）中整数变量的个数
%输入参数:MaxV最大速度,-MaxV最小速度,设正无穷为速度正方向
%输入参数:LoopCount迭代总数


%输出参数：ParSwarm：更新后的粒子群矩阵
%输出参数：OptSwarm：更新后的粒子最优解与群最优解位置矩阵
%输出参数：MaxValue：群最优解的适应度值（若采用罚函数法MaxValue就是最优解）


%开始单步更新的操作
%得到粒子群群体大小以及一个粒子维数的信息
[ParRow,ParCol]=size(ParSwarm);
%得到粒子的列数
ParCol=(ParCol-1)/2;


%*********************************************
%*****惯性因子参数*****
%****不同问题，好的惯性因子也不同***************


%线形递减策略
w=zeros(1,ParCol);
for i=1:ParCol
   w(1,i)=MaxV(i)-CurCount*((MaxV(i)+MaxV(i))/LoopCount);
end



%{
%w固定不变策略
w=zeros(1,ParCol);
for i=1:ParCol
   w(1,i)=0.7;
end
%}

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



%*********************************************
%*****速度更新参数*****
c1=2;
c2=2;
%
%con=1;
%c1=4-exp(-con*abs(mean(ParSwarm(:,2*ParCol+1))-AdaptFunc(OptSwarm(ParRow+1,:))));
%c2=4-c1;
%
%
%*********************************************



%*********************************************
%*****约束因子*****
%
%a=1;
%
a=0.729;
%*********************************************





%速度位置矩阵；
TempV=zeros(ParRow,ParCol);
TempPos=zeros(ParRow,ParCol);
%实际位置与粒子最优位置之间的差；
SubTract1=OptSwarm(1:ParRow,:)-ParSwarm(:,1:ParCol);

%初始种群最优解
MaxValue=AdaptFunc(OptSwarm(ParRow+1,:));

for row=1:ParRow
    %实际位置与群最优解之间的差；
    SubTract2=OptSwarm(ParRow+1,:)-ParSwarm(row,1:ParCol);
    %更新速度
    for col=1:ParCol
        TempV(row,col)=w(1,col)*ParSwarm(row,ParCol+col)+c1*random('unif',0,1).*SubTract1(row,col)+c2*random('unif',0,1).*SubTract2(col);
    end
    %限制粒子速度；
    for i=1:ParCol
        if TempV(row,i)>MaxV(i);
            TempV(row,i)=MaxV(i);
        end
        if TempV(row,i)<-MaxV(i);
            TempV(row,i)=-MaxV(i);
        end
    end
    ParSwarm(row,ParCol+1:2*ParCol)=TempV(row,:);
  
     %位置更新；
     for i=1:col
         TempPos(row,i)=ParSwarm(row,i)+a*TempV(row,i);
     end
     
      %整形变量和01变量的处理
      if Var>0 && Int>0
         for i=1:Var
             TempPos(row,i)=round(TempPos(row,i));     %四舍五入;
         end 
         for i=Var+1:Int+Var
             %TempPos(row,i)=floor(TempPos(row,i));    %负方向取整;
             %TempPos(row,i)=ceil(TempPos(row,i));     %正方向取整;
             TempPos(row,i)=round(TempPos(row,i));    %四舍五入;
             %TempPos(row,i)=fix(TempPos(row,i));      %取离零近的整数;
             %如果整形变量的定义域不是整数
             if TempPos(row,i)>ParticleScope(i,2)
                  TempPos(row,i)=TempPos(row,i)-1;
             elseif TempPos(row,i)<ParticleScope(i,1)
                  TempPos(row,i)=TempPos(row,i)+1;
             end
         end
      elseif Var>0
              for i=1:Var
                 TempPos(row,i)=round(TempPos(row,i));     %四舍五入;
              end 
      elseif Int>0
              for i=1:Int
                 %TempPos(row,i)=floor(TempPos(row,i));    %负方向取整;
                 %TempPos(row,i)=ceil(TempPos(row,i));     %正方向取整;
                 TempPos(row,i)=round(TempPos(row,i));    %四舍五入;
                 %TempPos(row,i)=fix(TempPos(row,i));       %取离零近的整数;
                 %如果整形变量的定义域不是整数
                 if TempPos(row,i)>ParticleScope(i,2)
                      TempPos(row,i)=TempPos(row,i)-1;
                 elseif ParSwarm(row,i)<ParticleScope(i,1)
                      TempPos(row,i)=TempPos(row,i)+1;
                 end
              end
      end
      %限制位置；
       for i=1:ParCol
           if TempPos(row,i)>=ParticleScope(i,2)
                 TempPos(row,i)=ParticleScope(i,2);
           end
           if TempPos(row,i)<=ParticleScope(i,1)
                 TempPos(row,i)=ParticleScope(i,1);
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
      
     
     %更新位置 
     ParSwarm(row,1:ParCol)=TempPos(row,:);
     %计算每个粒子的新的适应度值
     ParSwarm(row,2*ParCol+1)=AdaptFunc(ParSwarm(row,1:ParCol));
     if ParSwarm(row,2*ParCol+1)>AdaptFunc(OptSwarm(row,1:ParCol))
         OptSwarm(row,1:ParCol)=ParSwarm(row,1:ParCol);
     end
     %每个粒子适应度值得更新；
end


%全局最优解位置更新
[Max,row]=max(ParSwarm(:,2*ParCol+1));
if  Max>MaxValue
    %最优位置更新
    OptSwarm(ParRow+1,:)=ParSwarm(row,1:ParCol);
    %最优解更新
    MaxValue=Max;
end
end