clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Implementation of the Musical Chairs            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 60;  %range of time
K = 9;    %number of arm
M = 6;     %number of player

%%% special parameters for Musical Chairs
T0 = 10;

%%% paramètres de mesure
sumcolisions = zeros(1,N);
nbexp = 2000;

%%% Generation of the probability for mu
%mu = rand(K,1);
%mu = ones(K,1)/K;
mu = [.1;.2;.3;.4;.5;.6;.7;.8;.9]/4.5;

for iter = 1:nbexp

%%% Generation of all the Y(k,t)
for k = 1:K
    for t = 1:N
        if rand()>mu(k)
            Y{t}(k,1) = 0;
        else
            Y{t}(k,1) = 1;
        end
    end
end

%%% Initialisation of the choise of arm
A{1} = floor(rand(M,1)*K+1);
C{1} = zeros(1,M);

%detection of colisions and sensing?
for j1 = 1:M
    for j2 =1:M
        if A{1}(j1) == A{1}(j2) & not(j1 == j2)
            C{1}(j1) = true;
        end
    end
end

CT0 = zeros(1,M);
muest = zeros(K,M);
o = zeros(1,K);
s{1} = zeros(M,K);

%%% time iteration
for t=1:N
    if t<T0
        for j = 1:M
            A{t+1}(j) = floor(rand*K+1);
            if C{t}(j) == false
                o(j) = o(j)+1;
                s{t+1}(j,A{t+1}(j)) = s{t}(j,A{t}(j)) + not(Y{t}(A{t}(j)));
            else
                CT0(j) = CT0(j) + 1;
                s{t+1}(j,A{t+1}(j)) = s{t}(j,A{t}(j));
            end
        end
    else
        for j = 1:M
            for i=1:K
                if o(i) == 0
                    muest(i,j) = 0;
                else
                    muest(i,j) = s{T0}(j,i)/o(i);
                end
            end
            %%% Ordering of mu
            [sortV,Ordre] = sort(muest(:,j));
            AMC = [];
            for i=0:K-1
                AMC(i+1) = Ordre(end-i);
            end
            if 1 == T0
                Nstar = M;
            else
                Nstar = floor(log(abs(T0-CT0(j))/T0)/log(1-1/M)+1);
            end
            A{t+1}(j) = subroutineMC(AMC,N,T0,Nstar,C{t});
        end
    end
    
    %%% detection of colision at the time t+1
    for j1 = 1:M
        C{t+1}(1,j1) = false;
        for j2 =1:M
            if A{t+1}(j1) == A{t+1}(j2) & not(j1 == j2)
                C{t+1}(j1) = true;
            end
        end
    end
end

%%% computation of the number of colision
for t = 1:N
    sumcolisions(t) = sumcolisions(t) + sum(C{t});
end

end

%plot of A
figure();
for j=1:M
    for t=1:N
        a(t) = A{t}(j);
    end
    hold on;
    plot(1:N,a);
end

% plot colisions
figure();
plot(1:N,sumcolisions/nbexp);