function [whitenedM, whiteningMatrix, dewhiteningMatrix, means] = preprocessing(M)
    % find length
    l = max(size(M));
    
    %Centering M - subtract out means (add back in later)
    means = [mean(M(1,:)); mean(M(2,:))];
    Mc = M - repmat(means,1,l);

    %Whiten mixed signals
    %EVD of covariance matrix of M
    covM=cov(Mc',1);
    [E,D]=svd(covM,'econ');
    D2=sqrt(D);
    %whiteningMatrix = E*((D2)^-1)*E';
    whiteningMatrix = inv(D2)*E';
    dewhiteningMatrix = E*D2;
    %whitenedM=(E*D2*E.'*Mc')';
    % covM
    % E
    % size(whiteningMatrix)
    % size(Mc)
    whitenedM = whiteningMatrix*Mc;
end
