%----------------------------------------------------------------------------------------------------------------
% Authors: Hassan Talebi, Ute Mueller, Raimon Tolosana-Delgado
% Paper: "Joint simulation of compositional and categorical data via direct sampling technique – Application to
%         improve mineral resource confidence"
% Journal: Computer & Geosciences
%----------------------------------------------------------------------------------------------------------------

% The ilr transformation in this MATLAB function is based on the method proposed by Egozcue et al, 2003.
% Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and Barceló-Vidal, C., 2003. Isometric Logratio Transformations for Compositional Data Analysis. 
% Mathematical Geology, 35(3):279-300





function [ilr,clr,contrast_matrix] = ilr (data)

if nargin < 1
    error('Matrix of compositional data is needed')
end
    
[n,D]=size(data);
contrast_matrix=zeros(D-1,D);
clr=zeros(n,D);

for i=1:D-1
    for j=1:D
        if j<=D-i
            contrast_matrix(i,j)=sqrt(1/((D-i)*(D-i+1)));
        elseif j==(D-i+1)
            contrast_matrix(i,j)=-sqrt((D-i)/(D-i+1));
        else
            contrast_matrix(i,j)=0;
        end
    end
end

for i=1:n
    for j=1:D
        clr(i,j)=log(data(i,j)/exp(1/size(data(i,:),2)*sum(log(data(i,:)),2)));
    end
end

ilr=clr*contrast_matrix';



            