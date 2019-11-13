function [T, Y, RCNTRL, ICNTRL, RSTAT, ISTAT] = ...
    Rosenbrock(Function, Tspan, Y0, Options, RCNTRL, ICNTRL)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Rosenbrock - Implementation of several Rosenbrock methods:
%               * Ros2                                                
%               * Ros3                                                
%               * Ros4                                                
%               * Rodas3                                              
%               * Rodas4                                              
%                                                                     
%    Solves the system y'=F(t,y) using a Rosenbrock method defined by:
%                                                                     
%     G = 1/(H*gamma(1)) - Jac(t0,Y0)                                 
%     T_i = t0 + Alpha(i)*H                                           
%     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
%     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
%         gamma(i)*dF/dT(t0, Y0)
%     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j
%
%    For details on Rosenbrock methods and their implementation consult:
%      E. Hairer and G. Wanner
%      "Solving ODEs II. Stiff and differential-algebraic problems".
%      Springer series in computational mathematics, Springer-Verlag, 1996.
%    The codes contained in the book inspired this implementation.
%
%    MATLAB implementation (C) John C. Linford (jlinford@vt.edu).      
%    Virginia Polytechnic Institute and State University             
%    November, 2009    
%
%    Based on the Fortran90 implementation (C) Adrian Sandu, August 2004
%    and revised by Philipp Miehe and Adrian Sandu, May 2006.
%    Virginia Polytechnic Institute and State University
%    Contact: sandu@cs.vt.edu
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Input Arguments :
%  The first four arguments are similar to the input arguments of 
%  MATLAB's ODE solvers
%      Function  - A function handle for the ODE function
%      Tspan     - The time space to integrate
%      Y0        - Initial value
%      Options   - ODE solver options created by odeset():           
%                  AbsTol, InitialStep, Jacobian, MaxStep, and RelTol
%                  are considered.  Other options are ignored.       
%                  'Jacobian' must be set.                           
%      RCNTRL    - real value input parameters (explained below)
%      ICNTRL    - integer input parameters (explained below)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Output Arguments:
%  The first two arguments are similar to the output arguments of
%  MATLAB's ODE solvers
%      T         - A vector of final integration times.
%      Y         - A matrix of function values.  Y(T(i),:) is the value of
%                  the function at the ith output time.
%      RCNTRL    - real value input parameters (explained below)
%      ICNTRL    - integer input parameters (explained below)
%      RSTAT     - real output parameters (explained below)
%      ISTAT     - integer output parameters (explained below)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    RCNTRL and ICNTRL on input and output:
%
%    Note: For input parameters equal to zero the default values of the
%       corresponding variables are used.
%
%    ICNTRL(1) = 1: F = F(y)   Independent of T (AUTONOMOUS)
%              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
%
%    ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
%              = 1: AbsTol, RelTol are scalars
%
%    ICNTRL(3)  -> selection of a particular Rosenbrock method
%        = 0 :    Rodas3 (default)
%        = 1 :    Ros2
%        = 2 :    Ros3
%        = 3 :    Ros4
%        = 4 :    Rodas3
%        = 5 :    Rodas4
%
%    ICNTRL(4)  -> maximum number of integration steps
%        For ICNTRL(4)=0) the default value of 100000 is used
%
%    RCNTRL(1)  -> Hmin, lower bound for the integration step size
%          It is strongly recommended to keep Hmin = ZERO
%    RCNTRL(2)  -> Hmax, upper bound for the integration step size
%    RCNTRL(3)  -> Hstart, starting value for the integration step size
%
%    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
%    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
%    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
%                          (default=0.1)
%    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
%         than the predicted value  (default=0.9)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    RSTAT and ISTAT on output:
%
%    ISTAT(1)  -> No. of function calls
%    ISTAT(2)  -> No. of jacobian calls
%    ISTAT(3)  -> No. of steps
%    ISTAT(4)  -> No. of accepted steps
%    ISTAT(5)  -> No. of rejected steps (except at very beginning)
%    ISTAT(6)  -> No. of LU decompositions
%    ISTAT(7)  -> No. of forward/backward substitutions
%    ISTAT(8)  -> No. of singular matrix decompositions
%
%    RSTAT(1)  -> Texit, the time corresponding to the
%                     computed Y upon return
%    RSTAT(2)  -> Hexit, last accepted step before exit
%    RSTAT(3)  -> Hnew, last predicted step (not yet taken)
%                   For multiple restarts, use Hnew as Hstart
%                     in the subsequent run
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~>  Statistics on the work performed by the Rosenbrock method
global Nfun Njac Nstp Nacc Nrej Ndec Nsol Nsng Ntexit Nhexit Nhnew
global ISTATUS RSTATUS

% Parse ODE options
AbsTol = odeget(Options, 'AbsTol');
if isempty(AbsTol)
    AbsTol = 1.0e-3;
end
Hstart = odeget(Options, 'InitialStep');
if isempty(Hstart)
    Hstart = 0;
end
Jacobian = odeget(Options, 'Jacobian');
if isempty(Jacobian)
    error('A Jacobian function is required.');
end
Hmax = odeget(Options, 'MaxStep');
if isempty(Hmax)
    Hmax = 0;
end
RelTol = odeget(Options, 'RelTol');
if isempty(RelTol)
    RelTol = 1.0e-4;
end

% Initialize statistics
Nfun=1; Njac=2; Nstp=3; Nacc=4; Nrej=5; Ndec=6; Nsol=7; Nsng=8;
Ntexit=1; Nhexit=2; Nhnew=3;
RSTATUS = zeros(20,1);
ISTATUS = zeros(20,1);

% Get problem size
steps = length(Tspan);
N = length(Y0);

% Initialize tolerances
ATOL = ones(N,1)*AbsTol;
RTOL = ones(N,1)*RelTol;

% Initialize step
RCNTRL(2) = max(0, Hmax);
RCNTRL(3) = max(0, Hstart);

% Integrate
Y = zeros(steps,N);
T = zeros(steps,1);
Y(1,:) = Y0;
T(1) = Tspan(1);
for i=2:steps
    [T(i), Y(i,:), IERR] = ...
        RosenbrockIntegrator(N, Y(i-1,:), T(i-1), Tspan(i), ...
            Function, Jacobian, ATOL, RTOL, RCNTRL, ICNTRL);

    if IERR < 0
        error(['Rosenbrock exited with IERR=',num2str(IERR)]);
    end

end

ISTAT = ISTATUS;
RSTAT = RSTATUS;

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ T, Y, IERR ] = RosenbrockIntegrator(N, Y, Tstart, Tend, ...
    OdeFunction, OdeJacobian, AbsTol, RelTol, RCNTRL, ICNTRL)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Advances the ODE system defined by OdeFunction, OdeJacobian and Y
% from time=Tstart to time=Tend
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~>  Statistics on the work performed by the Rosenbrock method
global Nfun Njac Nstp Nacc Nrej Ndec Nsol Nsng Ntexit Nhexit Nhnew Ntotal
global ISTATUS RSTATUS

%~~~>   Parameters
global ZERO ONE DeltaMin
ZERO = 0.0; ONE = 1.0; DeltaMin = 1.0E-5;

%~~~>  Initialize statistics
ISTATUS(1:8) = 0;
RSTATUS(1:3) = ZERO;

%~~~>  Autonomous or time dependent ODE. Default is time dependent.
Autonomous = ~(ICNTRL(1) == 0);

%~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
%   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
if ICNTRL(2) == 0
    VectorTol = true;
    UplimTol  = N;
else
    VectorTol = false;
    UplimTol  = 1;
end

%~~~>   Initialize the particular Rosenbrock method selected
switch (ICNTRL(3))
    case (1)
        ros_Param = Ros2;
    case (2)
        ros_Param = Ros3;
    case (3)
        ros_Param = Ros4;
    case {0,4}
        ros_Param = Rodas3;
    case (5)
        ros_Param = Rodas4;
    otherwise
        disp(['Unknown Rosenbrock method: ICNTRL(3)=',num2str(ICNTRL(3))]);
        IERR = ros_ErrorMsg(-2,Tstart,ZERO);
        return
end
ros_S     = ros_Param{1};
rosMethod = ros_Param{2};
ros_A     = ros_Param{3};
ros_C     = ros_Param{4};
ros_M     = ros_Param{5};
ros_E     = ros_Param{6};
ros_Alpha = ros_Param{7};
ros_Gamma = ros_Param{8};
ros_ELO   = ros_Param{9};
ros_NewF  = ros_Param{10};
ros_Name  = ros_Param{11};

%~~~>   The maximum number of steps admitted
if ICNTRL(4) == 0
    Max_no_steps = 200000;
elseif ICNTRL(4) > 0
    Max_no_steps=ICNTRL(4);
else
    disp(['User-selected max no. of steps: ICNTRL(4)=',num2str(ICNTRL(4))]);
    IERR = ros_ErrorMsg(-1,Tstart,ZERO);
    return
end

%~~~>  Unit roundoff (1+Roundoff>1)
Roundoff = eps;

%~~~>  Lower bound on the step size: (positive value)
if RCNTRL(1) == ZERO
    Hmin = ZERO;
elseif RCNTRL(1) > ZERO
    Hmin = RCNTRL(1);
else
    disp(['User-selected Hmin: RCNTRL(1)=', num2str(RCNTRL(1))]);
    IERR = ros_ErrorMsg(-3,Tstart,ZERO);
    return
end

%~~~>  Upper bound on the step size: (positive value)
if RCNTRL(2) == ZERO
    Hmax = abs(Tend-Tstart);
elseif RCNTRL(2) > ZERO
    Hmax = min(abs(RCNTRL(2)),abs(Tend-Tstart));
else
    disp(['User-selected Hmax: RCNTRL(2)=', num2str(RCNTRL(2))]);
    IERR = ros_ErrorMsg(-3,Tstart,ZERO);
    return
end

%~~~>  Starting step size: (positive value)
if RCNTRL(3) == ZERO
    Hstart = max(Hmin,DeltaMin);
elseif RCNTRL(3) > ZERO
    Hstart = min(abs(RCNTRL(3)),abs(Tend-Tstart));
else
    disp(['User-selected Hstart: RCNTRL(3)=', num2str(RCNTRL(3))]);
    IERR = ros_ErrorMsg(-3,Tstart,ZERO);
    return
end

%~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
if RCNTRL(4) == ZERO
    FacMin = 0.2;
elseif RCNTRL(4) > ZERO
    FacMin = RCNTRL(4);
else
    disp(['User-selected FacMin: RCNTRL(4)=', num2str(RCNTRL(4))]);
    IERR = ros_ErrorMsg(-4,Tstart,ZERO);
    return
end
if RCNTRL(5) == ZERO
    FacMax = 6.0;
elseif RCNTRL(5) > ZERO
    FacMax = RCNTRL(5);
else
    disp(['User-selected FacMax: RCNTRL(5)=', num2str(RCNTRL(5))]);
    IERR = ros_ErrorMsg(-4,Tstart,ZERO);
    return
end

%~~~>   FacRej: Factor to decrease step after 2 succesive rejections
if RCNTRL(6) == ZERO
    FacRej = 0.1;
elseif RCNTRL(6) > ZERO
    FacRej = RCNTRL(6);
else
    disp(['User-selected FacRej: RCNTRL(6)=', num2str(RCNTRL(6))]);
    IERR = ros_ErrorMsg(-4,Tstart,ZERO);
    return
end

%~~~>   FacSafe: Safety Factor in the computation of new step size
if RCNTRL(7) == ZERO
    FacSafe = 0.9;
elseif RCNTRL(7) > ZERO
    FacSafe = RCNTRL(7);
else
    disp(['User-selected FacSafe: RCNTRL(7)=', num2str(RCNTRL(7))]);
    IERR = ros_ErrorMsg(-4,Tstart,ZERO);
    return
end

%~~~>  Check if tolerances are reasonable
for i=1:UplimTol
    if  (AbsTol(i) <= ZERO) || (RelTol(i) <= 10.0*Roundoff) || (RelTol(i) >= 1.0)
        disp([' AbsTol(',i,') = ',num2str(AbsTol(i))]);
        disp([' RelTol(',i,') = ',num2str(RelTol(i))]);
        IERR = ros_ErrorMsg(-5,Tstart,ZERO);
        return
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Template for the implementation of a generic Rosenbrock method
%      defined by ros_S (no of stages)
%      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~ Local variables
Ynew  = zeros(N,1);
Fcn0  = zeros(N,1);
Fcn   = zeros(N,1);
K     = zeros(N,ros_S);
dFdT  = zeros(N,1);
Jac0  = sparse(N,N);
Ghimj = sparse(N,N);
Yerr  = zeros(N,1);
Pivot = zeros(N,1);

%~~~>  Initial preparations
T = Tstart;
RSTATUS(Nhexit) = ZERO;
H = min(max(abs(Hmin), abs(Hstart)), abs(Hmax));
if abs(H) <= 10.0*Roundoff
    H = DeltaMin;
end

if Tend >= Tstart
    Direction = +1;
else
    Direction = -1;
end
H = Direction*H;

RejectLastH = false;
RejectMoreH = false;

%~~~> Time loop begins below

while ( (Direction > 0) && ((T-Tend)+Roundoff <= ZERO) || ...
        (Direction < 0) && ((Tend-T)+Roundoff <= ZERO) )

    if ISTATUS(Nstp) > Max_no_steps   % Too many steps
        IERR = ros_ErrorMsg(-6,T,H);
        return
    end
    if ((T+0.1*H) == T) || (H <= Roundoff)   % Step size too small
        IERR = ros_ErrorMsg(-7,T,H);
        return
    end

    %~~~>  Limit H if necessary to avoid going beyond Tend
    H = min(H, abs(Tend-T));

    %~~~>   Compute the function at current time
    Fcn0 = OdeFunction(T, Y);
    ISTATUS(Nfun) = ISTATUS(Nfun) + 1;

    %~~~>  Compute the function derivative with respect to T
    if ~Autonomous
        dFdT = ros_FunTimeDerivative (T, OdeFunction, Roundoff, Y, Fcn0);
    end

    %~~~>   Compute the Jacobian at current time
    Jac0 = OdeJacobian(T,Y);
    ISTATUS(Njac) = ISTATUS(Njac) + 1;

    %~~~>  Repeat step calculation until current step accepted
    while(true)

        [H, Ghimj, Pivot, Singular] = ...
            ros_PrepareMatrix(N, H, Direction, ros_Gamma(1), Jac0);

        % Not calculating LU decomposition in ros_PrepareMatrix anymore 
        % so don't need to check if the matrix is singular
        %
        %if Singular % More than 5 consecutive failed decompositions
        %    IERR = ros_ErrorMsg(-8,T,H);
        %    return
        %end

        % For the 1st istage the function has been computed previously
        Fcn = Fcn0;
        K(:,1) = Fcn;
        if (~Autonomous) && (ros_Gamma(1) ~= ZERO)
            HG = Direction*H*ros_Gamma(1);
            K(:,1) = K(:,1) + HG * dFdT;
        end
        K(:,1) = ros_Solve(Ghimj, Pivot, K(:,1));
        
        %~~~>   Compute the remaining stages
        for istage=2:ros_S

            % istage>1 and a new function evaluation is needed
            if ros_NewF(istage)
                Ynew = Y;
                for j=1:istage-1
                    Ynew = Ynew + ros_A((istage-1)*(istage-2)/2+j)*K(:,j)';
                end
                Tau = T + ros_Alpha(istage)*Direction*H;
                Fcn = OdeFunction(Tau,Ynew);
                ISTATUS(Nfun) = ISTATUS(Nfun) + 1;
            end

            K(:,istage) = Fcn;
            for j=1:istage-1
                HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H);
                K(:,istage) = K(:,istage) + HC * K(:,j);
            end
            if (~Autonomous) && (ros_Gamma(istage) ~= ZERO)
                HG = Direction*H*ros_Gamma(istage);
                K(:,istage) = K(:,istage) + HG * dFdT;
            end

            % Linear system is solved here with MATLAB '\' operator
            % instead of LU decomposition.
            K(:,istage) = ros_Solve(Ghimj, Pivot, K(:,istage));

        end

        %~~~>  Compute the new solution
        Ynew = Y;
        for j=1:ros_S
            Ynew = Ynew + ros_M(j) * K(:,j)';
        end

        %~~~>  Compute the error estimation
        Yerr = zeros(N,1);
        for j=1:ros_S
            Yerr = Yerr + ros_E(j) * K(:,j);
        end
        Err = ros_ErrorNorm(N, Y, Ynew, Yerr, AbsTol, RelTol, VectorTol);

        %~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
        Fac  = min(FacMax, max(FacMin, FacSafe/(Err^(ONE/ros_ELO))));
        Hnew = H*Fac;

        %~~~>  Check the error magnitude and adjust step size
        ISTATUS(Nstp) = ISTATUS(Nstp) + 1;
        if  (Err <= ONE) || (H <= Hmin)   %~~~> Accept step
            ISTATUS(Nacc) = ISTATUS(Nacc) + 1;
            Y = Ynew;
            T = T + Direction*H;
            Hnew = max(Hmin, min(Hnew, Hmax));
            if RejectLastH  % No step size increase after a rejected step
                Hnew = min(Hnew, H);
            end
            RSTATUS(Nhexit) = H;
            RSTATUS(Nhnew)  = Hnew;
            RSTATUS(Ntexit) = T;
            RejectLastH = false;
            RejectMoreH = false;
            H = Hnew;
            break % EXIT THE LOOP: WHILE STEP NOT ACCEPTED
        else           %~~~> Reject step
            if RejectMoreH
                Hnew = H*FacRej;
            end
            RejectMoreH = RejectLastH;
            RejectLastH = true;
            H = Hnew;
            if ISTATUS(Nacc) >= 1
                ISTATUS(Nrej) = ISTATUS(Nrej) + 1;
            end
        end % Err <= 1

    end % UntilAccepted

end % TimeLoop

%~~~> Succesful exit
IERR = 1;  %~~~> The integration was successful

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ ErrNorm ] = ros_ErrorNorm(N, Y, Ynew, Yerr, AbsTol, RelTol, VectorTol)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~> Computes the "scaled norm" of the error vector Yerr
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Err = 0;
for i=1:N
    Ymax = max(abs(Y(i)), abs(Ynew(i)));
    if VectorTol
        Scale = AbsTol(i) + RelTol(i)*Ymax;
    else
        Scale = AbsTol(1) + RelTol(1)*Ymax;
    end
    Err = Err + (Yerr(i)/Scale)^2;
end
Err = sqrt(Err/N);
ErrNorm = max(Err,1.0e-10);

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ dFdT ] = ros_FunTimeDerivative (T, OdeFunction, Roundoff, Y, Fcn0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~> The time partial derivative of the function by finite differences
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~>  Statistics on the work performed by the Rosenbrock method
global Nfun ISTATUS

%~~~>   Parameters
global ZERO ONE DeltaMin

Delta = sqrt(Roundoff) * max(DeltaMin, abs(T));
dFdT = OdeFunction(T+Delta,Y);
ISTATUS(Nfun) = ISTATUS(Nfun) + 1;
dFdT = (dFdT - Fcn0) / Delta;

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [H, Ghimj, Pivot, Singular] = ...
    ros_PrepareMatrix(N, H, Direction, gam, Jac0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Prepares the LHS matrix for stage calculations
%  1.  Construct Ghimj = 1/(H*ham) - Jac0
%      "(Gamma H) Inverse Minus Jacobian"
%  2.  LU decomposition not performed here because
%      MATLAB solves the linear system with '\'.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ghimj = -Jac0;
ghinv = 1.0/(Direction*H*gam);
for i=1:N
    Ghimj(i,i) = Ghimj(i,i)+ghinv;
end

Pivot(1) = 1;
Singular = false;

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ IERR ] = ros_ErrorMsg(Code, T, H)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Handles all error messages
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IERR = Code;

disp('Forced exit from Rosenbrock due to the following error:');

switch (Code)
    case (-1)
        disp('--> Improper value for maximal no of steps');
    case (-2)
        disp('--> Selected Rosenbrock method not implemented');
    case (-3)
        disp('--> Hmin/Hmax/Hstart must be positive');
    case (-4)
        disp('--> FacMin/FacMax/FacRej must be positive');
    case (-5)
        disp('--> Improper tolerance values');
    case (-6)
        disp('--> No of steps exceeds maximum bound');
    case (-7)
        disp('--> Step size too small: T + 10*H = T or H < Roundoff');
    case (-8)
        disp('--> Matrix is repeatedly singular');
    otherwise
        disp(['Unknown Error code: ', num2str(Code)]);
end

disp(['T=',num2str(T),' and H=',num2str(H)]);

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [X] = ros_Solve(JVS, Pivot, X)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Template for the forward/backward substitution 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global Nsol ISTATUS

X = JVS\X;
Pivot(1) = 1;
ISTATUS(Nsol) = ISTATUS(Nsol) + 1;

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Ros2()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- AN L-STABLE METHOD, 2 stages, order 2
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

g = 1.0 + 1.0/sqrt(2.0);
rosMethod = 1;
%~~~> Name of the method
ros_Name = 'ROS-2';
%~~~> Number of stages
ros_S = 2;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = (1.0)/g;
ros_C(1) = (-2.0)/g;
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1) = true;
ros_NewF(2) = true;
%~~~> M_i = Coefficients for new step solution
ros_M(1)= (3.0)/(2.0*g);
ros_M(2)= (1.0)/(2.0*g);
% E_i = Coefficients for error estimator
ros_E(1) = 1.0/(2.0*g);
ros_E(2) = 1.0/(2.0*g);
%~~~> ros_ELO = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus one
ros_ELO = 2.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.0;
ros_Alpha(2) = 1.0;
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = g;
ros_Gamma(2) =-g;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Ros3()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 2;
%~~~> Name of the method
ros_Name = 'ROS-3';
%~~~> Number of stages
ros_S = 3;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1)= 1.0;
ros_A(2)= 1.0;
ros_A(3)= 0.0;

ros_C(1) = -0.10156171083877702091975600115545E+01;
ros_C(2) =  0.40759956452537699824805835358067E+01;
ros_C(3) =  0.92076794298330791242156818474003E+01;
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1) = true;
ros_NewF(2) = true;
ros_NewF(3) = false;
%~~~> M_i = Coefficients for new step solution
ros_M(1) =  0.1E+01;
ros_M(2) =  0.61697947043828245592553615689730E+01;
ros_M(3) = -0.42772256543218573326238373806514;
% E_i = Coefficients for error estimator
ros_E(1) =  0.5;
ros_E(2) = -0.29079558716805469821718236208017E+01;
ros_E(3) =  0.22354069897811569627360909276199;
%~~~> ros_ELO = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
ros_ELO = 3.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1)= 0.0;
ros_Alpha(2)= 0.43586652150845899941601945119356;
ros_Alpha(3)= 0.43586652150845899941601945119356;
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1)= 0.43586652150845899941601945119356;
ros_Gamma(2)= 0.24291996454816804366592249683314;
ros_Gamma(3)= 0.21851380027664058511513169485832E+01;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Ros4()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
%     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
%
%      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
%      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
%      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
%      SPRINGER-VERLAG (1990)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 3;
%~~~> Name of the method
ros_Name = 'ROS-4';
%~~~> Number of stages
ros_S = 4;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = 0.2000000000000000E+01;
ros_A(2) = 0.1867943637803922E+01;
ros_A(3) = 0.2344449711399156;
ros_A(4) = ros_A(2);
ros_A(5) = ros_A(3);
ros_A(6) = 0.0;

ros_C(1) =-0.7137615036412310E+01;
ros_C(2) = 0.2580708087951457E+01;
ros_C(3) = 0.6515950076447975;
ros_C(4) =-0.2137148994382534E+01;
ros_C(5) =-0.3214669691237626;
ros_C(6) =-0.6949742501781779;
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1)  = true;
ros_NewF(2)  = true;
ros_NewF(3)  = true;
ros_NewF(4)  = false;
%~~~> M_i = Coefficients for new step solution
ros_M(1) = 0.2255570073418735E+01;
ros_M(2) = 0.2870493262186792;
ros_M(3) = 0.4353179431840180;
ros_M(4) = 0.1093502252409163E+01;
%~~~> E_i  = Coefficients for error estimator
ros_E(1) =-0.2815431932141155;
ros_E(2) =-0.7276199124938920E-01;
ros_E(3) =-0.1082196201495311;
ros_E(4) =-0.1093502252409163E+01;
%~~~> ros_ELO  = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
ros_ELO  = 4.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.0;
ros_Alpha(2) = 0.1145640000000000E+01;
ros_Alpha(3) = 0.6552168638155900;
ros_Alpha(4) = ros_Alpha(3);
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = 0.5728200000000000;
ros_Gamma(2) =-0.1769193891319233E+01;
ros_Gamma(3) = 0.7592633437920482;
ros_Gamma(4) =-0.1049021087100450;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Rodas3()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- A STIFFLY-STABLE METHOD, 4 stages, order 3
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 4;
%~~~> Name of the method
ros_Name = 'RODAS-3';
%~~~> Number of stages
ros_S = 4;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = 0.0;
ros_A(2) = 2.0;
ros_A(3) = 0.0;
ros_A(4) = 2.0;
ros_A(5) = 0.0;
ros_A(6) = 1.0;

ros_C(1) = 4.0;
ros_C(2) = 1.0;
ros_C(3) =-1.0;
ros_C(4) = 1.0;
ros_C(5) =-1.0;
ros_C(6) =-(8.0/3.0);

%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1)  = true;
ros_NewF(2)  = false;
ros_NewF(3)  = true;
ros_NewF(4)  = true;
%~~~> M_i = Coefficients for new step solution
ros_M(1) = 2.0;
ros_M(2) = 0.0;
ros_M(3) = 1.0;
ros_M(4) = 1.0;
%~~~> E_i  = Coefficients for error estimator
ros_E(1) = 0.0;
ros_E(2) = 0.0;
ros_E(3) = 0.0;
ros_E(4) = 1.0;
%~~~> ros_ELO  = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
ros_ELO  = 3.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.0;
ros_Alpha(2) = 0.0;
ros_Alpha(3) = 1.0;
ros_Alpha(4) = 1.0;
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = 0.5;
ros_Gamma(2) = 1.5;
ros_Gamma(3) = 0.0;
ros_Gamma(4) = 0.0;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Rodas4()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
%
%      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
%      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
%      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
%      SPRINGER-VERLAG (1996)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 5;
%~~~> Name of the method
ros_Name = 'RODAS-4';
%~~~> Number of stages
ros_S = 6;

%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.000;
ros_Alpha(2) = 0.386;
ros_Alpha(3) = 0.210;
ros_Alpha(4) = 0.630;
ros_Alpha(5) = 1.000;
ros_Alpha(6) = 1.000;

%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = 0.2500000000000000;
ros_Gamma(2) =-0.1043000000000000;
ros_Gamma(3) = 0.1035000000000000;
ros_Gamma(4) =-0.3620000000000023E-01;
ros_Gamma(5) = 0.0;
ros_Gamma(6) = 0.0;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = 0.1544000000000000E+01;
ros_A(2) = 0.9466785280815826;
ros_A(3) = 0.2557011698983284;
ros_A(4) = 0.3314825187068521E+01;
ros_A(5) = 0.2896124015972201E+01;
ros_A(6) = 0.9986419139977817;
ros_A(7) = 0.1221224509226641E+01;
ros_A(8) = 0.6019134481288629E+01;
ros_A(9) = 0.1253708332932087E+02;
ros_A(10) =-0.6878860361058950;
ros_A(11) = ros_A(7);
ros_A(12) = ros_A(8);
ros_A(13) = ros_A(9);
ros_A(14) = ros_A(10);
ros_A(15) = 1.0;

ros_C(1) =-0.5668800000000000E+01;
ros_C(2) =-0.2430093356833875E+01;
ros_C(3) =-0.2063599157091915;
ros_C(4) =-0.1073529058151375;
ros_C(5) =-0.9594562251023355E+01;
ros_C(6) =-0.2047028614809616E+02;
ros_C(7) = 0.7496443313967647E+01;
ros_C(8) =-0.1024680431464352E+02;
ros_C(9) =-0.3399990352819905E+02;
ros_C(10) = 0.1170890893206160E+02;
ros_C(11) = 0.8083246795921522E+01;
ros_C(12) =-0.7981132988064893E+01;
ros_C(13) =-0.3152159432874371E+02;
ros_C(14) = 0.1631930543123136E+02;
ros_C(15) =-0.6058818238834054E+01;

%~~~> M_i = Coefficients for new step solution
ros_M(1) = ros_A(7);
ros_M(2) = ros_A(8);
ros_M(3) = ros_A(9);
ros_M(4) = ros_A(10);
ros_M(5) = 1.0;
ros_M(6) = 1.0;

%~~~> E_i  = Coefficients for error estimator
ros_E(1) = 0.0;
ros_E(2) = 0.0;
ros_E(3) = 0.0;
ros_E(4) = 0.0;
ros_E(5) = 0.0;
ros_E(6) = 1.0;

%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1) = true;
ros_NewF(2) = true;
ros_NewF(3) = true;
ros_NewF(4) = true;
ros_NewF(5) = true;
ros_NewF(6) = true;

%~~~> ros_ELO  = estimator of local order - the minimum between the
%        main and the embedded scheme orders plus 1
ros_ELO = 4.0;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

% End of INTEGRATE function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


