function X = get_columVec(X)
% it also works with a matrix, in whihc case, it transposes is the nCol is
% larger than nRow

if size(X,1)< size(X,2)
   
    X        = X';
    
end



end