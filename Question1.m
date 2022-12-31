close all;
clear all;
clc;
options = input("Choose the method you want to use: \n Gauss Elimination (without pivoting)-1, \n Gauss Elimination (with partial pivoting)-2, \n Dolittle method (without pivoting)-3, \n Crout Method (without pivoting)-4, \n Cholesky Decomposition (without pivoting)-5 \n" ,'s');

%% Gauss Elimination (without pivoting)
if(options == '1')
   
    fileID = fopen("matrix.txt","rt");
    readsize1 = 1;
    n = fscanf(fileID, "%f", readsize1);
    readsize2 = [n+1 n];
    amatrix = fscanf(fileID, "%f", readsize2);
    amatrix = amatrix';
    % disp(amatrix);
    x = zeros(n,1);
    for j=1:n-1      % traversing through the column
        for i=j+1:n   % traversing through the rows
            factor = amatrix(i,j)/amatrix(j,j);
            amatrix(i,:) = amatrix(i,:) - factor*amatrix(j,:); 
        end
    end

    % back substitution to calculate the x vector
    x(n) = amatrix(n,n+1)/amatrix(n,n);
    for k=n-1:-1:1
        x(k) = (amatrix(k,n+1) - amatrix(k,k+1:n)*x(k+1:n))/amatrix(k,k);
    end

    disp(x)

    % Writing solutions to output file.
    fid = fopen("GuassElimination_withoutPivoting.txt", "w");
    fprintf(fid, "Gauss Elimination without Pivoting :\n");
    fprintf(fid, "%f\n", x);
    fclose(fid);
end

%% Gauss Elimination (with partial pivoting)
if(options == '2')
    fileID = fopen("matrix.txt","rt");
    readsize1 = 1;
    n = fscanf(fileID, "%f", readsize1);
    readsize2 = [n+1 n];
    amatrix = fscanf(fileID, "%f", readsize2);
    amatrix = amatrix';
    x = zeros(n,1);
    for j=1:n-1
        [max_, position] = max(abs(amatrix(j:n, j)));
        temp = amatrix(j, :);
        amatrix(j, :) = amatrix(position+j-1, :);
        amatrix(position+j-1, :) = temp;
        for i=j+1:n
            factor = amatrix(i,j)/amatrix(j,j);
            amatrix(i,:) = amatrix(i,:) - factor*amatrix(j,:);
        end
    end

    % back substitution to calculate the x values
    x(n) = amatrix(n,n+1)/amatrix(n,n);
    for k=n-1:-1:1
        x(k) = (amatrix(k,n+1) - amatrix(k,k+1:n)*x(k+1:n))/amatrix(k,k);
    end

    disp(x);

    % Writing solutions to output file.
    fid = fopen("GuassElimination_withPartialPivoting.txt", "w");
    fprintf(fid, "Gauss Elimination with Partial Pivoting :\n");
    fprintf(fid, "%f\n", x);
    fclose(fid);
end

%% LU Decomposition (Dolittle Method)
if(options == '3')

    fileID = fopen("matrix.txt","rt");
    readsize1 = 1;
    n = fscanf(fileID, "%f", readsize1);
    readsize2 = [n+1 n];
    amatrix = fscanf(fileID, "%f", readsize2);
    amatrix = amatrix';
    b = amatrix(:,4);
    amatrix = amatrix(:,1:3);
    % disp(amatrix);
    x = zeros(n,1);
    L = zeros(n,n);
    U = zeros(n,n);

    % Assigning 1 to the diagonal entries of the lower triangular matrix
    for a=1:n
        L(a,a) = 1;
    end

    % Calculating the coefficients
    U(1,:) = amatrix(1,:);
    L(:,1) = amatrix(:,1)/U(1,1);
    for i=2:n
        for j=i:n
            U(i,j) = amatrix(i,j) - L(i,1:i-1)*U(1:i-1,j);
        end
        for k=i+1:n
            L(k,i) = (amatrix(k,i) - L(k,1:i-1)*U(1:i-1,i))/U(i,i);
        end
    end

    disp(L);
    disp(U);

    Y = zeros(n,1);
    Y(1) = b(1)/L(1,1);
    for k=2:n
        Y(k) = (b(k) - L(k,1:k-1)*Y(1:k-1))/L(k,k);
    end
    x(n) = Y(n)/U(n,n);
    for k=n-1:-1:1
        x(k) = (Y(k) - U(k,k+1:n)*x(k+1:n))/U(k,k);
    end
 
    disp(x);
    % Writing solutions to output file.
    fid = fopen("DoLittle.txt", "w");
    fprintf(fid, "Dolittle Method :\n\n");
    fprintf(fid, "%f\n", x);
    fprintf(fid, "\nL:\n");
    for i=1:n
        fprintf(fid, "%f\t", L(i,:));
        fprintf(fid,"\n");
    end
    fprintf(fid, "\nU:\n");
    for i=1:n
        fprintf(fid, "%f\t", U(i,:));
        fprintf(fid, "\n");
    end
    fclose(fid);
    
end

%% LU Decomposition (Crout Method)
if(options == '4')
    
    fileID = fopen("matrix.txt","rt");
    readsize1 = 1;
    n = fscanf(fileID, "%f", readsize1);
    readsize2 = [n+1 n];
    amatrix = fscanf(fileID, "%f", readsize2);
    amatrix = amatrix';
    b = amatrix(:,4);
    amatrix = amatrix(:,1:3);
    % disp(amatrix);
    x = zeros(n,1);
    L = zeros(n,n);
    U = zeros(n,n);

    % Assigning 1 to the diagonal entries of the lower triangular matrix
    for a=1:n
        U(a,a) = 1;
    end

    % Calculating the coefficients
    U(1,:) = amatrix(1,:);
    L(:,1) = amatrix(:,1)/U(1,1);
    for i=2:n
        for k=i:n
            L(k,i) = (amatrix(k,i) - L(k,1:i-1)*U(1:i-1,i))/U(i,i);
        end
        for j=i:n
            U(i,j) = (amatrix(i,j) - L(i,1:i-1)*U(1:i-1,j))/L(i,i);
        end  
    end

    disp(L);
    disp(U);
    
    Y = zeros(n,1);
    Y(1) = b(1)/L(1,1);
    for k=2:n
        Y(k) = (b(k) - L(k,1:k-1)*Y(1:k-1))/L(k,k);
    end
    x(n) = Y(n)/U(n,n);
    for k=n-1:-1:1
        x(k) = (Y(k) - U(k,k+1:n)*x(k+1:n))/U(k,k);
    end
 
    disp(x);
    % Writing solutions to output file.
    fid = fopen("Crout.txt", "w");
    fprintf(fid, "Crout Method :\n\n");
    fprintf(fid, "%f\n", x);
    fprintf(fid, "\nL:\n");
    for i=1:n
        fprintf(fid, "%f\t", L(i,:));
        fprintf(fid,"\n");
    end
    fprintf(fid, "\nU:\n");
    for i=1:n
        fprintf(fid, "%f\t", U(i,:));
        fprintf(fid, "\n");
    end
    fclose(fid);

end

%% LU Decomposition (Cholesky Decomposition)
if(options == '5')
    
    fileID = fopen("matrix.txt","rt");
    readsize1 = 1;
    n = fscanf(fileID, "%f", readsize1);
    readsize2 = [n+1 n];
    amatrix = fscanf(fileID, "%f", readsize2);
    amatrix = amatrix';
    b = amatrix(:,4);
    amatrix = amatrix(:,1:3);
    % disp(amatrix);
    x = zeros(n,1);
    L = zeros(n,n);
    U = zeros(n,n);
    L(1,1) = sqrt(amatrix(1,1));
    U(1,1) = L(1,1);

    for a=2:n
        L(a,1) = amatrix(a,1)/L(1,1);
        U(1,a) = amatrix(1,a)/L(1,1);
    end

    for i=2:n
        for j=i:n
            if i==j
                L(j,i) = sqrt(amatrix(j,i) - L(j,1:i-1)*U(1:i-1,i));
                U(j,i) = L(j,i);
            else
                L(j,i) = (amatrix(j,i) - L(j,1:i-1)*U(1:i-1,i))/L(i,i);
            end
        end
        for k=i+1:n
            U(i,k) = (amatrix(i,k) - L(i,1:i-1)*U(1:i-1,k))/L(i,i);
        end  
    end

    Y = zeros(n,1);
    Y(1) = b(1)/L(1,1);
    for k=2:n
        Y(k) = (b(k) - L(k,1:k-1)*Y(1:k-1))/L(k,k);
    end
    x(n) = Y(n)/U(n,n);
    for k=n-1:-1:1
        x(k) = (Y(k) - U(k,k+1:n)*x(k+1:n))/U(k,k);
    end
    
    disp(x);
    % Writing solutions to output file.
    fid = fopen("Cholesky.txt", "w");
    fprintf(fid, "Cholesky Method :\n\n");
    fprintf(fid, "%f\n", x);
    fprintf(fid, "\nL:\n");
    for i=1:n
        fprintf(fid, "%f\t", L(i,:));
        fprintf(fid,"\n");
    end
    fprintf(fid, "\nU:\n");
    for i=1:n
        fprintf(fid, "%f\t", U(i,:));
        fprintf(fid, "\n");
    end
    fclose(fid);

end