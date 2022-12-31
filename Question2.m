clear all;
close all;
clc;

options = input("Choose the method you want to use: \n Power Method-1, \n Inverse Power Method-2, \n Inverse Power Method with Shift-3, \n QR Decomposition-4 \n" ,'s');

%% Power Method
if(options == '1')

    fid = fopen("question2.txt", "rt");
    readsize1 = 1;
    n = fscanf(fid, "%f", readsize1);
    readsize2 = [n n];
    amatrix = fscanf(fid, "%f", readsize2);
    readsize3 = 1;
    maxiterations = fscanf(fid, "%f", readsize3);
    readsize4 = 1;
    maxerror = fscanf(fid, "%f", readsize4);
    readsize5 = 1;
    shift = fscanf(fid, "%f", readsize5);
    xold = ones(n,1);
    error = 100;
    iterations = 0;
    lambda = 0;

    while((error > maxerror) && (iterations < maxiterations))
        v = amatrix*xold;
        maxeig = max(v);
        eigs(iterations+1) = maxeig;
        error = abs((maxeig - lambda)/maxeig)*100;
        v = v/maxeig;
        lambda = maxeig;
        xold = v;
        iterations = iterations + 1;
    end

    fileID = fopen("PowerMethod.txt", "w");
    fprintf(fileID, "Power Method: \n\n");
    fprintf(fileID, "Eigenvalue: \n");
    fprintf(fileID, "%f", maxeig);
    fprintf(fileID, "\n\n");
    fprintf(fileID, "Eigenvector: \n");
    fprintf(fileID, "%f\t", v);
    fprintf(fileID, "\n\n");
    fprintf(fileID, "Iterations: \n");
    fprintf(fileID, "%f\n\n", iterations);

    for i=1:iterations
        fprintf(fileID, "Iteration %f : %f\n", i, eigs(i));
    end


end

%% Inverse Power Method
if(options == '2')

    fid = fopen("question2.txt", "rt");
    readsize1 = 1;
    n = fscanf(fid, "%f", readsize1);
    readsize2 = [n n];
    amatrix = fscanf(fid, "%f", readsize2);
    readsize3 = 1;
    maxiterations = fscanf(fid, "%f", readsize3);
    readsize4 = 1;
    maxerror = fscanf(fid, "%f", readsize4);
    readsize5 = 1;
    shift = fscanf(fid, "%f", readsize5);
    xold = ones(n,1);
    error = 100;
    iterations = 0;
    lambda = 0;

    % calculating the inverse of the given matrix
    amatrix = inv(amatrix);

    while((error > maxerror) && (iterations < maxiterations))
        v = amatrix*xold;
        maxeig = max(v);
        eigs(iterations+1) = 1/maxeig;
        error = abs((maxeig - lambda)/maxeig)*100;
        v = v/maxeig;
        lambda = maxeig;
        xold = v;
        iterations = iterations + 1;
    end

    fileID = fopen("InversePowerMethod.txt", "w");
    fprintf(fileID, "Inverse Power Method: \n\n");
    fprintf(fileID, "Eigenvalue: \n");
    fprintf(fileID, "%f", 1/maxeig);
    fprintf(fileID, "\n\n");
    fprintf(fileID, "Eigenvector: \n");
    fprintf(fileID, "%f\t", v);
    fprintf(fileID, "\n\n");
    fprintf(fileID, "Iterations: \n");
    fprintf(fileID, "%f\n\n", iterations);
    
    for i=1:iterations
        fprintf(fileID, "Iteration %f : %f\n", i, eigs(i));
    end
end

%% Inverse Power Method with Shift
if(options == '3')
    
    fid = fopen("question2.txt", "rt");
    readsize1 = 1;
    n = fscanf(fid, "%f", readsize1);
    readsize2 = [n n];
    amatrix = fscanf(fid, "%f", readsize2);
    readsize3 = 1;
    maxiterations = fscanf(fid, "%f", readsize3);
    readsize4 = 1;
    maxerror = fscanf(fid, "%f", readsize4);
    readsize5 = 1;
    shift = fscanf(fid, "%f", readsize5);
    xold = ones(n,1);
    error = 100;
    iterations = 0;
    lambda = 0;

    amatrix = amatrix - shift*eye(n,n);

    % calculating the inverse of the given matrix
    amatrix = inv(amatrix);

    while((error > maxerror) && (iterations < maxiterations))
        v = amatrix*xold;
        maxeig = norm(abs(v));
        eigs(iterations + 1) = 1/maxeig+shift;
        error = abs((maxeig - lambda)/maxeig)*100;
        v = v/maxeig;
        lambda = maxeig;
        xold = v;
        iterations = iterations + 1;
    end

    fileID = fopen("InversePowerMethodwithShift.txt", "w");
    fprintf(fileID, "Inverse Power Method: \n\n");
    fprintf(fileID, "Eigenvalue: \n");
    fprintf(fileID, "%f",1/maxeig + shift);
    fprintf(fileID, "\n\n");
    fprintf(fileID, "Eigenvector: \n");
    fprintf(fileID, "%f\t", v);
    fprintf(fileID, "\n\n");
    fprintf(fileID, "Iterations: \n");
    fprintf(fileID, "%f\n\n", iterations);

    for i=1:iterations
        fprintf(fileID, "Iteration %f : %f\n", i, eigs(i));
    end

end

%% QR Decomposition
if(options == '4')
    
  fid = fopen("question2.txt", "rt");
  readsize1 = 1;
  n = fscanf(fid, "%f", readsize1);
  readsize2 = [n n];
  A = fscanf(fid, "%f", readsize2);
  readsize3 = 1;
  maxiterations = fscanf(fid, "%f", readsize3);
  readsize4 = 1;
  maxerror = fscanf(fid, "%f", readsize4);
  readsize5 = 1;
  shift = fscanf(fid, "%f", readsize5);
  q=zeros(n);
  r=zeros(n);
  m=zeros(n,1);
  value2=zeros(1,n);
  value1=[];
  i=0;
  j=0;
  k=0;
  for i= 1:100
    for j= 1:n
      m=zeros(n,1);
      value1= [];
      for k= 1:j-1
          m=m+(transpose(A(:,j))*q(:,k))*q(:,k);
      end
      q(:,j)=A(:,j)-m;    
      q(:,j)=q(:,j)/sqrt(sum(q(:,j).*q(:,j))); 
    end
    for index= 1:n
       for j= index:n
          r(index,j)=transpose(q(:,index))*A(:,j); 
       end
       value1=[value1,r(index,index)];
    end
    A=r*q;
    if i>1
       if max(abs((value1-value2)./value2)) < maxerror/100
         break;
       end
    end
    value2=value1;
  end

  fileID = fopen('QRDecomposition.txt','w');
  fprintf(fileID, "QR Decomposition Method: \n\n");
  fprintf(fileID,"Eigen Values: \n\n");
  fprintf(fileID, "%f\n", value2);
  fprintf(fileID, "\n");
  fprintf(fileID,'Iterations: %d\n',i);
  fclose(fileID);

end