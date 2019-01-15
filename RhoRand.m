clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Implementation of the RhoRand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 60;  %range of time
K = 9;    %number of arm
M = 6;     %number of player

%%% paramètres de mesure
sumcolisions = zeros(1,N);
nbexp = 2;

%%% Generation of the probability for mu
%mu = rand(K,1);
mu = ones(K,1)/K;
%mu = [.1;.2;.3;.4;.5;.6;.7;.8;.9]/4.5;

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

%%% Ordering of mu
[sortV,Ordre] = sort(mu);
best = [];
for i=0:M-1
    best(i+1) = Ordre(end-i);
end

%%% Initialisation of the choise of arm
A{1} = floor(rand(M,1)*K+1);
C{1} = zeros(1,M);
s{1} = zeros(1,M);

%detection of colisions and sensing?
for j1 = 1:M
    for j2 =1:M
        if A{1}(j1) == A{1}(j2) & not(j1 == j2)
            C{1}(j1) = true;
        end
    end
end

%%% time iteration
for t=1:N
    for j = 1:M
        if C{t}(j) == 1
            A{t+1}(j) = best(floor(rand*M+1));
        else
            A{t+1}(j) = A{t}(j);
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