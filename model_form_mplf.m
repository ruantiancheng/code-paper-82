function [Etau1,Etau2] = model_form_mplf(n,tau,kij,hij)
%MODEL_FORM 此处显示有关此函数的摘要
%   此处显示详细说明
A = [0 1 0;
     0 0 1;
     0 0 -1/tau];
B = [0;0;1/tau]; 
B_ = kron(eye(n),B);
An = zeros(n,n);

infor_num = 2;
for i =1:n
    for j = 1:n
        if j<i
            An(i,j) = 1/(i-1)
        end
    end
end

F = [];
J = [];
for i = 1:n
    D = [];
    Hi = [];
    fi = [1 hij(i) 0;
            0 1 0;
            0 0 1];
    all_fi = [];
    for j =1:n
        temp = An(i,j) * kij(i,:);
        all_fi = [all_fi fi];
        D = [D temp];
    end    
    Hi = temp * all_fi;
    F = blkdiag(F,Hi);
    J = blkdiag(J,D);
end
I3 = eye(3);
I2 = [];
I1 = [];
for i = 1:n
   I2 = blkdiag(I2,I3);
   I1 = [I1,transpose(I3)];
end
I1 = transpose(I1);

E1 = [];
E2 = [];
for i = 1:n
    E1 = blkdiag(E1,I1);
    E2 = [E2,transpose(I2)];
end
E2 = transpose(E2);
A_ = kron(eye(n),A);
Etau2 = B_ * J*-E2;
Etau1 =A_ -B_ * F*E1;
end

