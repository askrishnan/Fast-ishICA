function [A,W] = icaRounds(whiteningMatrix, dewhiteningMatrix, Mw, means)
    a1 = 1;
    
    iterations=5;
    
    nic = min(size(whiteningMatrix)); % number of independent components
    nsamp = max(size(Mw));
    
    savedW = ones(nic,max(size(whiteningMatrix)));

    for k=1:iterations
    
        w = randn(nic,1);
        B = zeros(nic);
        round = 1;
        while round <= nic
            % wnew = mean(M*g(w'*M))-mean(dg(w'*M))*w;
            wnew = w-B*B'*w;
            w = wnew/norm(wnew);

            i = 1;
            while i <= 1000  % you might want to try how many times this loops throug
                
                wnew = w-B*B'*w;
                if wnew(1)==0 && wnew(2)==0
                    wnew=[.001;.001];
                end
                w = wnew/norm(wnew);
                wold = zeros(size(w));
                wold2 = zeros(size(w));
                % check for convergence
                c = dot(wold, w);
                if abs(c-1) <= .0001
                    break;
                else
                    B(:, round) = w;
                   
                    A(:, round) = dewhiteningMatrix*w;
                    W(round, :) = w'*whiteningMatrix;
                    wold2 = wold;
                    wold = w;
                    %hypTan = tanh(a1 * Mw' * w);
                    %w = (Mw * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / nsamp;
                    w = (Mw*((Mw'*w).^3))/nsamp-3*w;
                    w = w/norm(w);
                    i = i + 1; 
                end
            end
          round = round + 1;
        end  
       
       newW=W; 
        %Test non-gaussianity! 
        sig1=savedW*Mw;
        sig2=newW*Mw;
        
        kurt1=mean(sig1.^4)-3;
        kurt2=mean(sig2.^4)-3;
        
        if abs(kurt1) < abs(kurt2)
            savedW = savedW;
        else
            savedW=newW;
        end
        
    end
    
    W=savedW;
    
%     if ~isreal(A)
%         A=real(A);
%         W=real(W);
%     end

%icasig = W * Mw + (W * means) * ones(1, max(size(whiteningMatrix)));
end