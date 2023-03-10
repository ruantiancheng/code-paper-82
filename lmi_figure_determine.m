function [MM,NN,YY11,YY12,YY22,LL1,LL2] = lmi_figure_determine(formatSpec,num, tau,hm,hij,kij)

n = 3*num; % vector number

[Psi,Psid] = model_form(num,tau,kij,hij); % system modeling

setlmis([]);
M = lmivar(1,[n,1]);
N = lmivar(1,[n,1]);
Y11 = lmivar(1,[n,1]);
Y12 = lmivar(1,[n,1]);
Y22 = lmivar(1,[n,1]);
L1 = lmivar(1,[n,1]);
L2 = lmivar(1,[n,1]);

lmiterm([1 1 1 M],Psi,1,'s');
lmiterm([1 1 1 L1],1,1,'s');
lmiterm([1 1 1 Y11],hm,1);

lmiterm([1 1 2 M],1,Psid);
lmiterm([1 1 2 L1],-1,1);
lmiterm([1 1 2 -L2],1,1);
lmiterm([1 1 2 Y12],hm,1);

lmiterm([1 2 2 L2],-1,1,'s');
lmiterm([1 2 2 Y22],hm,1);
 
lmiterm([1 1 3 N],hm*Psi',1);

lmiterm([1 2 3 N],hm*Psid',1);

lmiterm([1 3 3 N],-hm,1);



lmiterm([2 1 1 Y11],-1,1);

lmiterm([2 1 2 Y12],-1,1);

lmiterm([2 1 3 L1],-1,1);

lmiterm([2 2 2 Y22],-1,1);

lmiterm([2 2 3 L2],-1,1);

lmiterm([2 3 3 N],-1,1);

% 正定矩阵限制

lmiterm([3 1 1 M],-1,1);
lmiterm([4 1 1 N],-1,1);
lmiterm([5 1 1 Y11],-1,1);
lmiterm([5 1 2 Y12],-1,1);
lmiterm([5 2 2 Y22],-1,1);


% 求解LMI系统
lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
tmin;
% % 
MM=dec2mat(lmisys,xfeas,M);
NN=dec2mat(lmisys,xfeas,N);
YY11=dec2mat(lmisys,xfeas,Y11);
YY12=dec2mat(lmisys,xfeas,Y12);
YY22=dec2mat(lmisys,xfeas,Y22);
LL1=dec2mat(lmisys,xfeas,L1);
LL2=dec2mat(lmisys,xfeas,L2);

str = sprintf(formatSpec,i);
csvwrite([str,'M.csv'],MM);
csvwrite([str,'N.csv'],NN);
csvwrite([str,'Y11.csv'],YY11);
csvwrite([str,'Y12.csv'],YY12);
csvwrite([str,'Y22.csv'],YY22);
csvwrite([str,'L1.csv'],LL1);
csvwrite([str,'L2.csv'],LL2);
