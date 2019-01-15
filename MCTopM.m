clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Implementation of the MCTopM                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 60;  %range of time
K = 9;    %number of arm
M = 6;     %number of player

%%% paramètres de mesure
sumcolisions = zeros(1,N);
nbexp = 2000;

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

%%% Computation of g for t=1
for k=1:K   %computation of T and  S and mubar also
    for j=1:M
        if A{1}(j)==k
            T{1}(k,j) = 1;
            S{1}(k,j) = Y{1}(k);
            mubar{1}(k,j) = Y{1}(k);
        else
            T{1}(k,j) = 0;
            S{1}(k,j) = 0;
            mubar{1}(k,j) = 0;
        end
    end
end
g{1} = mubar{1}; %because log(1)=0

%%% Computation of Mhat
for j=1:M
    [sortV,Ordre] = sort(g{1}(:,j));    %tri des valeurs de g
    for m=0:M-1
        Mhat{1}(j,m+1) = Ordre(end-m); %on récupère indices des M plus grandes valeurs
    end
end

%%% time iteration
for t=1:N
   %%% player iteration
    for j = 1:M
       if min(A{t}(j)~=Mhat{t}(:,j))%si A n'appartient pas à Mhat
           % computation of the set of the possible values
           set=[];
           for i=1:M
               if t == 1
                   set = Mhat{t}(i,j);
               elseif g{t-1}(Mhat{t}(i,j),j)<=g{t-1}(A{t}(j),j)
                   set(end+1) = Mhat{t}(i,j);
               end
           end
           if not(size(set,1))
               set = Mhat{t}(i,j);
           end
           A{t+1}(j) = set(floor(rand*size(set,2)+1));% random 
           s{t+1}(j) = false;
       elseif C{t}(j) & not(s{t}(j))
           alea = floor(rand()*M +1); %compute a random number in 1:M (size of Mhat)
           A{t+1}(j) = Mhat{t}(alea,j); %a vérifier
           s{t+1}(j) = false;
       else
           A{t+1}(j) = A{t}(j);
           s{t+1}(j) = true;
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
    %%% Computation of g 
    T{t+1} = zeros(K,M);
    S{t+1} = zeros(K,M);
    mubar{t+1} = zeros(K,M);
    g{t+1} = zeros(K,M);
    for k=1:K   %computation of T and  S and mubar also
        for j=1:M
            if A{t}(j)==k
                T{t+1}(k,j) = T{t}(k,j)+1;
                S{t+1}(k,j) = S{t}(k,j)+Y{t}(k);
                mubar{t+1}(k,j) = S{t+1}(k,j)/T{t+1}(k,j);
                g{t+1}(k,t) = mubar{t+1}(k,j)+sqrt(log(t)/(2*T{t+1}(k,j)));
            else
                T{t+1}(k,j) = T{t}(k,j);
                S{t+1}(k,j) = S{t}(k,j);
                mubar{t+1}(k,j) = mubar{t}(k,j);
                if T{t+1}(k,j)==0
                    g{t+1}(k,t) = 0;
                else
                    g{t+1}(k,j) = mubar{t+1}(k,j)+sqrt(log(t)/(2*T{t+1}(k,j)));
                end
            end
        end
    end
    %%% Computation of Mhat
    for j=1:M
        [sortV,Ordre] = sort(g{t+1}(:,j));    %tri des valeurs de g
        for m=0:M-1
            Mhat{t+1}(j,m+1) = Ordre(end-m); %on récupère indices des M plus grandes valeurs
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