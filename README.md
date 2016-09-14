PSO

Particle Swarm Optimization

Working on Matlab

基本函数：

PSO          主函数，包括参数的输入，结果输出

PsoProcess   种群的迭代计算，多次重复计算

InitSwarm    生成初始种群

BaseStepPso  种群单步更新

AdaptFunc    适应度函数，约束条件用罚函数来体现（效率比较低）


输入参数基本要求：

粒子位置变量按01变量，整形变量，非整型变量排列

需要输入所有变量的可行域（不能为无穷大）
  
