%----------------------------------------------------------------------------------------------------------------
% Authors: Hassan Talebi, Ute Mueller, Raimon Tolosana-Delgado
% Paper: "Joint simulation of compositional and categorical data via direct sampling technique – Application to
%         improve mineral resource confidence"
% Journal: Computer & Geosciences
%----------------------------------------------------------------------------------------------------------------

% The ilr transformation in this MATLAB function is based on the method proposed by Egozcue et al, 2003.
% Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and Barceló-Vidal, C., 2003. Isometric Logratio Transformations for Compositional Data Analysis. 
% Mathematical Geology, 35(3):279-300



function [compositions]=b_ilr (ilr,contrast_matrix,constant)

if nargin <2
    error('Matrix of ilr data and contrast matrix are needed')
end

[n,balances]=size(ilr);
D=balances+1;
compositions=zeros(n,D);
delta=exp(ilr*contrast_matrix);

if nargin ==2
    for i=1:n
        for j=1:D
            compositions(i,j)=delta(i,j)/sum(delta(i,:));
        end
    end
end

if nargin ==3
    for i=1:n
        for j=1:D
            compositions(i,j)=delta(i,j)*constant/sum(delta(i,:));
        end
    end
end
