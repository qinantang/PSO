function [ParSwarm,OptSwarm,MaxV]=InitSwarm(SwarmSize,ParticleSize,ParticleScope,Var,Int,Filter)
 %[ParSwarm,OptSwarm,BadSwarm]=InitSwarm(SwarmSize,ParticleSize,ParticleScope,AdaptFunc)
 %
 %输入参数：SwarmSize:种群大小的个数
 %输入参数：ParticleSize：一个粒子的维数
 %输入参数：ParticleScope:一个粒子在运算中各维的范围；
 %　　ParticleScope格式:
 %　　3维粒子的ParticleScope格式:
 %                               [x1Min,x1Max]
 %　　　　　　　　　　　　　　　　  x2Min,x2Max
 %                                x3Min,x3Max]
 %
 %输入参数：Var：01变量的个数
 %输入参数：Int：整形变量的个数    
 %输入参数：Filter：筛选因子
 
 %输出：ParSwarm初始化的粒子群
 %输出：OptSwarm粒子群最优解与群最优解位置矩阵
 %ParSwarm前ParticleSize列为位置参数，ParticleSize~2*ParticleSize列为速度参数，最后为适应度函数
  
  %速度系数；
  coefV=0.2;
  
  ParSwarm=zeros(SwarmSize,2*ParticleSize+1);
  %随机生成粒子,构建种群，有下种两方法:
  for i=1:ParticleSize
      ParSwarm(i,:)=random('unif',ParticleScope(i,1),ParticleScope(i,2),1,2*ParticleSize+1);
      %unif按平均分布生成随机位置，norm按正态分布生成随机位置，poiss按泊西分布生成随机位置。
  end
  %1.用while语句筛选出符合不等式约束的粒子建立种群，但是筛选因数受主观因素影响较大
  for n=1:SwarmSize
    %选取符合要求的粒子；
    while (1)
      %整形变量和01变量的处理
      if Var>0 && Int>0
         for i=1:Var
             ParSwarm(n,i)=round(ParSwarm(n,i));    %四舍五入;
         end 
         for i=Var+1:Int+Var
             %ParSwarm(n,i)=floor(ParSwarm(n,i));    %负方向取整;
             %ParSwarm(n,i)=ceil(ParSwarm(n,i));     %正方向取整;
             ParSwarm(n,i)=round(ParSwarm(n,i));     %四舍五入;
             %ParSwarm(n,i)=fix(ParSwarm(n,i));      %取离零近的整数;             
             %如果整形变量的定义域不是整数
             if ParSwarm(n,i)>ParticleScope(i,2)
                 ParSwarm(n,i)=ParSwarm(n,i)-1;
             elseif ParSwarm(n,i)<ParticleScope(i,1)
                 ParSwarm(n,i)=ParSwarm(n,i)+1;
             end
         end
      elseif Int>0
              for i=1:Int
                 %ParSwarm(n,i)=floor(ParSwarm(n,i));    %负方向取整;
                 %ParSwarm(n,i)=ceil(ParSwarm(n,i));     %正方向取整;
                 ParSwarm(n,i)=round(ParSwarm(n,i));     %四舍五入;
                 %ParSwarm(n,i)=fix(ParSwarm(n,i));      %取离零近的整数;
                 %如果整形变量的定义域不是整数
                 if ParSwarm(n,i)>ParticleScope(i,2)
                     ParSwarm(n,i)=ParSwarm(n,i)-1;
                 elseif ParSwarm(n,i)<ParticleScope(i,1)
                     ParSwarm(n,i)=ParSwarm(n,i)+1;
                 end
              end
             
      elseif Var>0
              for i=1:Var
                 ParSwarm(n,i)=round(ParSwarm(n,i));    %四舍五入;
              end 
      else
      end
      %利用筛选因子评价出事粒子，
      if AdaptFunc(ParSwarm(n,:))>-Filter 
           %初始速度设置，系数0.2可调
           ParSwarm(n,ParticleSize:2*ParticleSize)=coefV*random('unif',0,1)*ParSwarm(n,ParticleSize:2*ParticleSize);
           break;
      else
           ParSwarm(n,:)=random('unif',ParticleScope(n,1),ParticleScope(n,2),1,2*ParticleSize+1);
           %unif按平均分布生成随机位置，norm按正态分布生成随机位置，poiss按泊西分布生成随机位置。
      end
    end
 end
  
  %{
  %2.用while语句筛选出符合不等式约束的一个“好的”粒子，以此为中心，可行域为直径，向任意方向生长，形成种群，
     %缺点在于，如果筛选出的粒子不够“好”或者可行域离散分布时，容易陷入局部最优解!!!，但在某些问题上效率是否更高？？？
     Filter;       %筛选因数,根据罚函数确定；
     Swarm=rand(1,ParticleSize);
      for i=1:ParticleSize
                Swarm(1,i)=Swarm(1,i)*(ParticleScope(i,2)-ParticleScope(i,1))+ParticleScope(i,1);
      end
    while (1)
      %整形变量和01变量的处理
      if Var>0 && Int>0
         for i=1:Var
             round(Swarm(1,ParticleSize+i));    %四舍五入;
         end 
         for i=Var+1:Int+Var
             %floor(Swarm(1,ParticleSize+i));    %负方向取整;
             %ceil(Swarm(1,ParticleSize+i));     %正方向取整;
             round(Swarm(1,ParticleSize+i));    %四舍五入;
             %fix(Swarm(1,ParticleSize+i));      %取离零近的整数;
         end
       else if Var>0
              for i=1:Var
                 round(Swarm(1,ParticleSize+i));    %四舍五入;
              end 
          else  Int>0
              for i=1:Int
                 %floor(Swarm(1,ParticleSize+i));    %负方向取整;
                 %ceil(Swarm(1,ParticleSize+i));     %正方向取整;
                 round(Swarm(1,ParticleSize+i));    %四舍五入;
                 %fix(Swarm(1,ParticleSize+i));      %取离零近的整数;
              end
          end
      end
      if AdaptFunc(Swarm(1,:))>Filter   
           break;
      else
         Swarm=rand(1,ParticleSize);
         for i=1:ParticleSize
            Swarm(1,i)=Swarm(1,i)*(ParticleScope(i,2)-ParticleScope(i,1))+ParticleScope(i,1);
         end
      end
    end 
      ParSwarm(1,1:ParticleSize)=Swarm(1,:);   %筛选出的“好的”粒子，整个种群的“母代”粒子；
      for n=2:SwarmSize
          for i=1:ParticleSize
               ParSwarm(n,i)=ParSwarm(1,i)+0.5*(2*rand-1)*(ParticleScope(i,2)-ParticleScope(i,1));
              if ParSwarm(n,i)>ParticleScope(i,2)
                  ParSwarm(n,i)=ParticleScope(i,2);
              end
              if ParSwarm(n,i)<ParticleScope(i,1)
              ParSwarm(n,i)=ParticleScope(i,1);
              end
          end
      end
 %初始速度设置
 for i=1:Integer
    %调节速度，使速度与位置的范围一致
    ParSwarm(:,ParticleSize+i)=0.2*(ParSwarm(:,ParticleSize+i)*(ParticleScope(i,2)-ParticleScope(i,1))+ParticleScope(i,1));
 end
 for i=Integer+1:ParticleSize
    %调节速度，使速度与位置的范围一致
    ParSwarm(:,ParticleSize+i)=0.2*(ParSwarm(:,ParticleSize+i)*(ParticleScope(i,2)-ParticleScope(i,1))+ParticleScope(i,1));
 end
%} 
MaxV=zeros(1,ParticleSize);
for i=1:ParticleSize
    MaxV(i)=coefV*(ParticleScope(i,2)-ParticleScope(i,1));
end
 
%对每一个粒子计算其适应度函数的值
for i=1:SwarmSize
    ParSwarm(i,2*ParticleSize+1)=AdaptFunc(ParSwarm(i,1:ParticleSize));
end

 %初始化粒子群最优解矩阵
 %粒子群最优解矩阵全部设为零
 OptSwarm=zeros(SwarmSize+1,ParticleSize);
 %寻找适应度函数值最大的解在矩阵中的位置(行数)
 [~,row]=max(ParSwarm(:,2*ParticleSize+1));
 %粒子最优解位置参数
 OptSwarm(1:SwarmSize,1:ParticleSize)=ParSwarm(1:SwarmSize,1:ParticleSize);
 %种群最优解位置参数
 OptSwarm(SwarmSize+1,:)=ParSwarm(row,1:ParticleSize);
 end




