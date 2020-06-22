%Gauss-Seidel implementation by Jainish Mehta
% The implementation was done as a function, where you can input matrix
% inputs.The screenshots for each circuit matrix are in the report.
%We need the following inputs:
% A is the coefficient matrix
% b is the constants vector
% initial_guess is the initial guess vector as a row vector (will get
% converted to column vector in this program)
%err_tolerance is the maximum tolerable error after which the program will
%stop
% max_iter is the maximum iteration
% lambda is the relaxation factor
% An alternative, similar approach would be like in the tutorial question
%where you find the strictly lower-triangular form of matrix A, strictly upper-triangular form of matrix A
%and the inverse (D+lambda*L) to plug into x = d*(lambda*b - (lambda*U + D*(lambda-1))*x_0)
function x = GaussSeidel(A,b,initial_guess,err_tolerance,max_iter, lambda)

[m,n] = size(A);

% Confirms matrix A is a square matrix, i.e the matrix has the same number
% of rows as columns
if m~=n
    error('Must be a square square matrix');
end
%Lambda by definition must be 0<x<2.
if (lambda >= 2)
     error('Lambda value too large');
end
if (lambda <=0)
     error('Lambda value too small');
end

%Checking if the matrix is diagonally dominant as per the definition,
%all rows non-diagonal elements must be equal to greater to diagonal
%elements with one row strictly greater than. The if statements check if 
%it can be possible, check and if not already 
%then the else tries to switch rows to make it happen
if (all((2*abs(diag(A))) >= sum(abs(A),2)))
    for L = 1:n
        if (2*abs(diag(A(L))) > sum(abs(A(L)),2))
            disp ('This matrix is diagonally dominant');
            break;
        end
    end
else
%This find the maximum value of each row and the column it is located in,
%which is stored in collocation
[maxrow,collocation] = max(abs(A),[],2);
  %This works by checking if for all the rows, the diagonal element is larger than 
  %the sum of the other elements. The column location of each largest
  %element of each row is ordered, taking into account the entire row. 
  %Numel returns the the number of elements, n, in the matrix
  if all(2*maxrow > sum(abs(A),2)) && isequal(sort(collocation),(1:numel(collocation))')
    A(collocation,:) = A;
  else
      %Either the matrix is already diagonally dominant from statement
      %above, or else it simply cannot be made diagonally dominant.
    disp('This matrix cannot altered to be made to be diagonally dominant. This does not mean it will not converge though!')
  end
end
% Stores the value of A to maintain its elements and initializes the
% iteration to increment it
temp_A = A;
iteration = 0;

% sets diagonal of temp_A to 0. The for loops iterates through the number
% of rows
for P = 1:n
    temp_A(P,P) = 0;
end

% Changes x from an array to a column matrix
x = initial_guess';
% divides matrix elements by each element on the diagonal
for L = 1:n
    temp_A(L,1:n) = temp_A(L,1:n)/A(L,L);   
end
% Divides vector "b" by each element on the diagonal
for k = 1:n
    b(k) = (b(k) / A(k,k));    
end
%  Initializes percent error
Percent_err = zeros(n);

for j = 1:n
    Percent_err(j) = 100;    
end

%Used for printing out iteration of the error values
iterationcounter = 1;

% Iterates until either the desired error or max iterations is reached
 while iteration < max_iter
    
    % Saves the current value of x at the start of the loop or
    % alternatively at the end
    x_old = x;
    %  Gauss Seidel Method/SOR where lambda is multiplied by new value and 
    % 1- lambda is multiplied by the old value i.e. over-relaxation means
    % the new x-value is given more weight, while the old is given less
    %This follows the equation of SOR that the new value is (1-w)old plus
    %(w/Aii) multiplied by (Bi - summation of Aij times new value) minus
    %summation of Aij times the old value). This is the formula. X is the
    %output column vector/matrix.
    for i = 1:n 
        x(i) = (lambda)*(b(i)-temp_A(i,:)*x) +(1-lambda)*x_old(i);
        
        % Calculates percent error when x is not zero as that would not be 
        %possible. This is calculated using (new-old)/new times 100
        if x(i) ~= 0
            Percent_err(i) = 100* abs((x(i) - x_old(i)) /x(i));    
        end
    end
      %Print out the error values for reference for each variable
    disp('ITERATION NUMBER:')
      disp(iterationcounter)
     for i = 1:n
      fprintf('The error of the %i th variable is %.3f percent\n',i,Percent_err(i));
     end
     fprintf('\n')
      iterationcounter = iterationcounter+1;
   % Terminates program if ALL the variable percent errors are less than
   % then the maximum_tolerable error.
    if(max(Percent_err, [], 'all') < err_tolerance)
        disp("Maximum percent error reached!")
        break;
    end
    % Increments iteration
    iteration = iteration + 1;
 end
  %Prints out maximum percent error at each iteration
     fprintf('The max percentage error is %.3f percent.', max(Percent_err, [], 'all'));
 %To implement the voltage calculations, you need to include Ohm's Law,
 %V=IR. It is omitted as it is not required and is highly dependent on the
 %circuit itself and would not make sense in a function as is my implementation.
 %However it is an easy calculation.
end
