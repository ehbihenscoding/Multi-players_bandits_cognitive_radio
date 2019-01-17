function j=subroutineMC(A,N,T0,Nstar,C)
%This subroutine is describe in the article about Musical Chairs
    for k=1:N-T0
        i = floor(rand*min(Nstar,size(C,2))+1);
        if C(i) == false
            j = A(i);
            return
        end
    end
    j = i;
end