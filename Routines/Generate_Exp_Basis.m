% this is a program for generating exponentail type basis function by 
% using Gram-Schmidt process
% created on Oct. 14, 2015; successful;

% function Generate_Exp_Basis

clear;

Num_basis = 10; % at most 10, otherwise accuracy issue
%a = -4;   % positive is growth, negative is decay
T = 1001; % time grid points 

t = linspace(0,1,T)';
my_power = linspace(-5,-1,Num_basis);
%temp_power = linspace(-1,-10,Num_basis-1);
%my_power = [0,temp_power];
Exp_Basis = zeros(T,Num_basis);
Ortho_Exp_Basis = zeros(T,Num_basis);

% create original exponential function
for k = 1:Num_basis
    %if mod(k,2)==0
        Exp_Basis(:,k) = exp(my_power(k)*t);
    %else
        %Exp_Basis(:,k) = exp(my_power(k)*t);
   % end
end
plot(t,Exp_Basis);

% do Gram-Schmidt Process
Ortho_Exp_Basis(:,1) = Exp_Basis(:,1);
for k = 2:Num_basis
    temp_poly_sum = zeros(T,1); 
    for j = 1:k-1
        coef_top = L2_InnerProduct(Exp_Basis(:,k),Ortho_Exp_Basis(:,j));
        coef_bot = L2_InnerProduct(Ortho_Exp_Basis(:,j),Ortho_Exp_Basis(:,j));
        temp_poly_sum = temp_poly_sum + coef_top/coef_bot * Ortho_Exp_Basis(:,j);
    end
    Ortho_Exp_Basis(:,k) = Exp_Basis(:,k) - temp_poly_sum;
end

%Check_Orthogonality(Ortho_Exp_Poly(:,5),Ortho_Exp_Poly(:,10),t)

% from orthogonal into orthonormal
for k = 1:Num_basis
    Ortho_Exp_Basis(:,k) = Ortho_Exp_Basis(:,k)/L2norm(t,Ortho_Exp_Basis(:,k));
end

%plot(t,Ortho_Exp_Poly(:,6:10));
plot(t,Ortho_Exp_Basis);

