function lmi_figure_hinfty(formatSpec,num, hm,tau,hij,kij_1)

cvx_begin sdp

n = 3*num;


[Psi,Psi_d] = model_form(num,tau,kij_1,hij); % system modeling
% Define the variables (adjust n for appropriate dimensions)
variable S(n, n) symmetric
variable T(n, n) symmetric
variable U(n, n) symmetric
variable V(n, n) symmetric
variable W(n, n)
variable X(n, n)

% Define the system matrices Psi and Psi_d and delay bound h
% Replace these with actual values

h = hm; % Define h

% Define matrices Y1, Y2, Y3, Y4
Y1 = S*Psi + Psi'*S + T;
Y2 = W*X - W - X;
Y3 = Psi_d*S + Psi_d*S; % Adjusted to include Psi_d
Y4 = S;

% Define the LMIs
LMI1 = [Y1 + Y2 + Y2', sqrt(h)*Y3'*U; sqrt(h)*Y3*U, -U];
LMI2 = V*W'*U > 0;

% Define the constraints
subject to
    S > 0
    T > 0
    U > 0
    V > 0
    LMI1 < 0
    LMI2 > 0

% Solve the problem
cvx_end

% Check feasibility and define the controller
if strcmp(cvx_status, 'Solved')
    % System is asymptotically stable
    disp('System is asymptotically stable.')
else
    % System may not be stable under given conditions
    disp('System may not be stable under the given conditions.')
end