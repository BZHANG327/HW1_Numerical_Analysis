%===============Iterative methods for Sparse linear system=================
tic
clear;
m=8;
error_vec = zeros(m-1,1);
steps_vec = zeros(m-1,1);
N_vec = zeros(m-1,1);
for k = 2:m
    N = 2^k; %Initialize the number grid points in one axes
    x = 0:(1/N):1;
    y = 1:(-1/N):0;
    [X,Y] = meshgrid(x,y);
    U = zeros(N+1);
    Unew = zeros(N+1);
    %Unew = U;
    Unew(:,1) = cos(2*pi*y); % The function f
    %Unew(:,1) = sign(cos(2*pi*y)); % The function f_hat
    maxiter = 1000;
    steps = 0;
    %=========================Jacobi Iteration=================================
    while norm(Unew-U) >1/N^2%Set the error tolerance
        U = Unew;
        for j = 2:N %index of y axes
            for i = 2:N %index of x axes
                if j==2
                    Unew(j,i) = 1/3*(U(j,i-1)+U(j,i+1)+U(j+1,i));
                elseif j==N
                    Unew(j,i) = 1/3*(U(j,i-1)+U(j,i+1)+U(j-1,i));
                else
                    Unew(j,i) = 1/4*(U(j,i-1)+U(j,i+1)+U(j+1,i)+U(j-1,i));
                end
            end
        end
        Unew(1,2:N) = Unew(2,2:N);
        Unew(N+1,2:N) = Unew(N,2:N);
        steps = steps +1;
        %norm(Unew-U)
    end
    for i = 1:N+1
        for j = 1:N+1
            U(i,j) = cos(2*pi*Y(i,j))*(-exp(2*pi*X(i,j))/(exp(4*pi)-1)...
                +(exp(4*pi)/(exp(4*pi)-1))*exp(-2*pi*X(i,j)));
        end
    end
    error_vec(k-1) = norm(reshape(Unew,[],1)-reshape(U,[],1),inf);
    steps_vec(k-1) = steps;
    
    h_vec = (2:k);
    N_vec(k-1) = N;
end
handle_figure = figure(1);
hold on
plot(h_vec,-log2(error_vec),'*')
line(h_vec,-log2(error_vec),'LineStyle','-')
line([0,9],[0,18]);
hold off
handle_figure = figure(2);
hold on 
plot(N_vec,steps_vec,'*')
line(N_vec,steps_vec,'Linestyle','-')
plot((1:N),(1:N).^2)
toc