function [T,Y,RCNTRL,ICNTRL,RSTATUS,ISTATUS]=...
    SdirkInt(Function,Tspan,Y0,Options,RCNTRL,ICNTRL)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  SDIRK - Implementation of several SDIRK methods:
%               * SDIRK4a                                                
%               * SDIRK4b                                                
%               * SDIRK3a                                                
%               * SDIRK2b                                              
%               * SDIRK2a
%                                                                     
%    Solves the system y'=F(t,y) using a SDIRK Method 
%                                                                     
%
%    For details on SDIRK methods and their implementation consult:
%      E. Hairer and G. Wanner
%      "Solving ODEs II. Stiff and differential-algebraic problems".
%      Springer series in computational mathematics, Springer-Verlag, 1996.
%    The codes contained in the book inspired this implementation.
%
%    MATLAB implementation (C) Vishwas Rao (visrao@vt.edu).      
%    Virginia Polytechnic Institute and State University             
%    March, 2011    
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
%      RSTATUS     - real output parameters (explained below)
%      ISTATUS     - integer output parameters (explained below)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    RCNTRL and ICNTRL on input and output:
%
%    Note: For input parameters equal to zero the default values of the
%       corresponding variables are used.
%
%    
%
%    ICNTRL(1) = 0: AbsTol, RelTol are N-dimensional vectors
%              = 1: AbsTol, RelTol are scalars
%
%    ICNTRL(2)  -> selection of a particular SDIRK method
%        = 0 :    Sdirk2a(Default)
%        = 1 :    Sdirk2a
%        = 2 :    Sdirk2b
%        = 3 :    Sdirk3a
%        = 4 :    Sdirk4a
%        = 5 :    Sdirk4b
%        
%        
%
%    ICNTRL(3)  -> maximum number of integration steps
%        For ICNTRL(3)=0) the default value of 100000 is used
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
%    ISTATUS(1)  -> No. of function calls
%    ISTATUS(2)  -> No. of jacobian calls
%    ISTATUS(3)  -> No. of steps
%    ISTATUS(4)  -> No. of accepted steps
%    ISTATUS(5)  -> No. of rejected steps (except at very beginning)
%    ISTATUS(6)  -> No. of LU decompositions
%    ISTATUS(7)  -> No. of forward/backward substitutions
%    ISTATUS(8)  -> No. of singular matrix decompositions
%
%    RSTATUS(1)  -> Texit, the time corresponding to the
%                     computed Y upon return
%    RSTATUS(2)  -> Hexit, last accepted step before exit
%    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
%                   For multiple restarts, use Hnew as Hstart
%                     in the subsequent run
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
global Nfun Njac Nstp Nacc Nrej Ndec Nsol Nsng Ntexit Nhexit Nhnew
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
N = max(size(Y0));

% Initialize tolerances
ATOL = ones(N,1)*AbsTol;
RTOL = ones(N,1)*RelTol;

% Initialize step
RCNTRL(2) = max(0, Hmax);
RCNTRL(3) = max(0, Hstart);

% Integrate
Y = zeros(N,steps);
T = zeros(steps,1);
Y(:,1) = Y0;
T(1) = Tspan(1);
t=cputime;
for i=2:steps
    [T(i), Y(:,i), IERR,ISTATUS,RSTATUS] = ...
        Sdirk_Integrate(N, Y(:,i-1), T(i-1), Tspan(i), ...
            Function, Jacobian, ATOL, RTOL, RCNTRL, ICNTRL,ISTATUS,RSTATUS);

    if IERR < 0
        error(['SDIRK exited with IERR=',num2str(IERR)]);
    end

end
cputime-t



function [T,Y,Ierr,ISTATUS,RSTATUS] = Sdirk_Integrate(N,Y0,TIN, TOUT,OdeFunction1,...
    OdeJacobian1,ATOL,RTOL,RCNTRL,ICNTRL,ISTATUS,RSTATUS)
global ZERO ONE Nfun Njac Nstp Nacc Nrej Ndec Nsol Nsng Ntexit Nhexit Nhnew
global Smax S2A S2B S3A S4A S4B 
global rkTheta rkAlpha
ZERO = 0; ONE=1.0; Nfun = 1; Njac = 2; Nstp = 3; Nacc = 4; Nrej = 5;Ndec=6;
Nsol = 7; Nsng = 8; Ntexit = 1; Nhexit = 2; Nhnew = 3; Smax = 5; S2A = 1;
S2B = 2; S3A = 3; S4A = 4; S4B = 5;
rkTheta = zeros(Smax);
rkAlpha = zeros(Smax);

%Local Variables
ISTATUS = zeros(20,1);
RSTATUS = zeros(20,1);

% Fine-Tune the integrator
T1 = TIN;
T2 = TOUT;
RTOL=ones(size(Y0));
RTOL=1e-3*RTOL;
ATOL=RTOL;

    [T,Y,Ierr,RSTATUS,ISTATUS] = SDIRK(N, T1,T2,Y0,RTOL,ATOL,RCNTRL,ICNTRL,RSTATUS,ISTATUS,OdeFunction1,OdeJacobian1);

if Ierr < 0
    disp(['SDIRK:Unsuccessful exit at ',num2str(TIN), ' ',num2str(Ierr)]);
end
return;




function [T, Y, Ierr,RSTATUS,ISTATUS] = SDIRK(N,T1,T2,Y,AbsTol,RelTol,...
                        RCNTRL,ICNTRL,RSTATUS,ISTATUS,OdeFunction1,OdeJacobian1)

global ZERO ONE 


Max_no_steps = 0;Hmin = 0;Hmax = 0; Hstart = 0; Roundoff = eps; FacMin = 0;
FacMax = 0; FacSafe = 0; FacRej = 0; 
if(ICNTRL(1)==0)
    ITOL = 1;
else
    ITOL = 0;
end
Ierr = 0;
%~~~> Method Selection
switch (ICNTRL(2))
    
    case (0)
        Sdirk2a();
    case (1)
        Sdirk2a();
    case (2)
        Sdirk2b();
    case(3)
        Sdirk3a();
    case(4)
        Sdirk4a();
    case(5)
        Sdirk4b();
    otherwise
        Sdirk2a();
end

if (ICNTRL(3)==ZERO)
    Max_no_steps = 200000;
end
if(ICNTRL(3) > ZERO)
        Max_no_steps = ICNTRL(3);
end

if(ICNTRL(3) < ZERO)
    disp(['User-Selected ICNTRL(3)=',num2str(ICNTRL(3))]);
    SDIRK_ErrorMsg(-1,T1,ZERO,Ierr);
end;
%~~>StartNewton Extrapolate for starting values of Newton iterations
if(ICNTRL(5)==ZERO)
    StartNewton = 1;
else
    StartNewton = 0;
end;

if(ICNTRL(4)>=ZERO)
    if(ICNTRL(4)==ZERO)
        NewtonMaxit = 8;
    else
        NewtonMaxit = ICNTRL(4);
    end
else
    disp(['User-Selected NewtonMaxit ICNTRL(4)=',num2str(ICNTRL(4))]);
    SDIRK_ErrorMsg(-2,T1,ZERO,Ierr);
end;

%~~> Lower bound on the step size (Positive value)

if(RCNTRL(1)==ZERO)
    Hmin = ZERO;
end
if (RCNTRL(1)>ZERO)
        Hmin = RCNTRL(1);
end    
if(RCNTRL(1)<ZERO)
    disp(['User-Selected Hmin, RCNTRL(1)=',num2str(RCNTRL(1))]);
    SDIRK_ErrorMsg(-3,T1,ZERO,Ierr);
end;
% Upper bound on the step size : (Positive value)

if(RCNTRL(2)==ZERO)
    Hmax = abs(T2-T1);
end
if(RCNTRL(2)>ZERO)
    Hmax = min(abs(T2-T1),RCNTRL(2));
end;
if(RCNTRL(2)<ZERO)
    disp(['User-Selected Hmax, RCNTRL(2)=',num2str(RCNTRL(2))]);
    SDIRK_ErrorMsg(-3,T1,ZERO,Ierr);
end;
%~~> Starting step size

if(RCNTRL(3)>=ZERO)
    if(RCNTRL(3)==ZERO)
        Hstart = max(Hmin,Roundoff);
    else
        Hstart = min(abs(RCNTRL(3)),abs(T2-T1));
    end
else
    disp(['User-Selected Hstart, RCNTRL(3)=',num2str(RCNTRL(3))]);
    SDIRK_ErrorMsg(-3,T1,ZERO,Ierr);
end;
%~~> Stepsize can be changed such that FacMin < Hnew/Hexit < FacMax
if(RCNTRL(4) >= ZERO)
    if(RCNTRL(4) == ZERO)
        FacMin = 0.2;
    else
        FacMin = RCNTRL(4);
    end
else
    disp(['User-Selected FacMin, RCNTRL(4)=',num2str(RCNTRL(4))]);
    SDIRK_ErrorMsg(-4,T1,ZERO,Ierr);
end


if (RCNTRL(5)>=ZERO)
    if(RCNTRL(5)==ZERO)
        FacMax = 10.0;
    else
        FacMax = RCNTRL(5);
    end
else
    disp(['User-Selected FacMax, RCNTRL(5)=',num2str(RCNTRL(5))]);
    SDIRK_ErrorMsg(-4,T1,ZERO,Ierr);
end

%~~> FacRej Factor to decrease step after 2 successive rejections

if(RCNTRL(6)>=ZERO)
    if(RCNTRL(6) == ZERO)
        FacRej = 0.1;
    else
        FacRej = RCNTRL(6);
    end
else
    disp(['User-Selected FacRej, RCNTRL(6)=',num2str(RCNTRL(6))]);
    SDIRK_ErrorMsg(-4,T1,ZERO,Ierr);
end;

%~~> FacSafe Safety Factor in the computation of new step size

if(RCNTRL(7)>=ZERO)
    if(RCNTRL(7)==ZERO)
        FacSafe = 0.9;
    else
        FacSafe = RCNTRL(7);
    end
else
    disp(['User-Selected FacSafe, RCNTRL(7)=',num2str(RCNTRL(7))]);
    SDIRK_ErrorMsg(-4,T1,ZERO,Ierr);
end

%~~> ThetaMin decides whether the Jacobian should be recomputed

if(RCNTRL(8) == ZERO)
    ThetaMin = 1.0E-3;
else
    ThetaMin = RCNTRL(8);
end;

%~~> Stopping Criteria for Newtons method

if(RCNTRL(9) == ZERO)
    NewtonTol = 3.0e-4;
else
    NewtonTol = RCNTRL(9);
end;

%~~> Qmin, Qmax : IF Qmin <Hnew/Hold<Qmax STEP SIZE = constant

if(RCNTRL(10) == ZERO)
    Qmin = ONE;
else
    Qmin = RCNTRL(10);
end

if(RCNTRL(11) == ZERO)
    Qmax = 1.2;
else
    Qmax = RCNTRL(11);
end

%~~> Check if tolerances are reasonable
if(ITOL == ZERO)
    if((AbsTol(1)<=ZERO) || RelTol(1) <=10*Roundoff)
        SDIRK_ErrorMsg(-5,T1,ZERO,Ierr);
    end;
else
    for i=1:N
        if((AbsTol(i)<=ZERO)||RelTol(i)<=10*Roundoff)
            SDIRK_ErrorMsg(-5,T1,ZERO,Ierr);
        end
    end
end
if(Ierr < 0)
    return;
end



[T,Y,Ierr,ISTATUS,RSTATUS] = SDIRK_Inetgrator(N,T1,T2,Y,Hstart,Hmin,...
                            Hmax,Roundoff,AbsTol,RelTol,...
                            ITOL,Max_no_steps,StartNewton,NewtonTol,ThetaMin,...
                            FacSafe,FacMin,FacMax,FacRej,Qmin,Qmax,NewtonMaxit,ISTATUS,RSTATUS,Ierr,OdeFunction1,OdeJacobian1);
 


return;

function[T,Y,Ierr,ISTATUS,RSTATUS] = SDIRK_Inetgrator(N,T1,T2,Y,Hstart,Hmin,...
                            Hmax,Roundoff,AbsTol,RelTol,...
                            ITOL,Max_no_steps,StartNewton,NewtonTol,ThetaMin,...
                            FacSafe,FacMin,FacMax,FacRej,Qmin,Qmax,NewtonMaxit,ISTATUS,RSTATUS,Ierr,OdeFunction1,OdeJacobian1)
global ZERO ONE Nfun Nstp Nacc Nrej Ntexit Nhexit Nhnew
global rkS rkGamma rkELO
global rkC rkD rkE rkTheta rkAlpha
T =T1;
Tdirection = sign(T2-T1);
H = max(abs(Hmin),abs(Hstart));
if(H<=10.0*Roundoff)
    H = 1.0e-6;
end;
H = min(H,Hmax);
H=sign(Tdirection)*H;
SkipLU = 0;
SkipJac = 0;
Reject = 0;
FirstStep = 1;
CycleTloop = 0;
SCAL = SDIRK_ErrorScale(N,ITOL,AbsTol,RelTol,Y);
FJAC = zeros(N);
while((T2-T)*Tdirection - Roundoff > ZERO)
    if(SkipLU == 0)
        [E1L,E1U,ISTATUS,ISING]=PrepareMatrix(N,H,T,Y,FJAC,SkipJac, SkipLU,Reject,ISTATUS,OdeJacobian1);
        if(ISING ~= 0)
            SDIRK_ErrorMsg(-8,T,H,Ierr);
            return;
        end;
    end;
    if(ISTATUS(Nstp)>Max_no_steps)
        SDIRK_ErrorMsg(-6, T, H, Ierr);
        return;
    end;
    
    if(T+0.1*H==T || abs(H)<=Roundoff)
        SDIRK_ErrorMsg(-7,T,H,Ierr);
        return;
    end
    %~~> Simplified Newton Iterations
    Z=zeros(rkS,N);
    for istage=1:rkS
        G=zeros(N,1);
        if(istage>1)
            aa=rkTheta(1:istage,istage);
            G=Z(1:istage,:)'*aa;
            if(StartNewton==1)
                for j=1:istage
                    Z(istage,:)=Z(istage,:)+rkAlpha(j,istage)*Z(j,:);
                end
                
            end
                
		end
		NewtonDone = 0;
		Fac = 0.5;
		for	NewtonIter=0:NewtonMaxit-1
			TMP=Z(istage,:)'+Y;
			RHS=OdeFunction1(T+rkC(istage)*H,TMP);
			ISTATUS(Nfun)=ISTATUS(Nfun)+1;
			RHS=H*rkGamma*RHS;
			RHS=RHS-Z(istage,:)';
            RHS=RHS+G;
			[RHS,ISTATUS]=SDIRK_Solve(H,N,E1L,E1U,RHS,ISTATUS);
			NewtonIncrement = SDIRK_ErrorNorm(N,RHS,SCAL);
			if(NewtonIter==0)
				Theta = abs(ThetaMin);
				NewtonRate = 2.0;
			else
				Theta = NewtonIncrement/NewtonIncrementOld;
				if (Theta < 0.99)
					NewtonRate = Theta/(1-Theta);
					NewtonPredictedErr = NewtonIncrement*Theta^((NewtonMaxit - (NewtonIter+1))/(ONE-Theta));
					if(NewtonPredictedErr>=NewtonTol)
						Qnewton = min(10,NewtonPredictedErr/NewtonTol);
						Fac = 0.8*Qnewton ^(-ONE/(NewtonMaxit-NewtonIter));
						break;
					end
				else
					break;
				end
			end;
			NewtonIncrementOld = NewtonIncrement;
			Z(istage,:)=Z(istage,:)+RHS';
			NewtonDone=(NewtonRate*NewtonIncrement<=NewtonTol);
			if(NewtonDone == 1)
				break;
			end;
		end;

		if(NewtonDone == 0)
			H=Fac*H;
			Reject =1;
			SkipJac = 1;
			SkipLU = 0;
			CycleTloop = 1;
		end;
		if(CycleTloop == 1)
			CycleTloop = 0;
			break;
		end
	end;
	if(CycleTloop==0)
	%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
	%~~~~ Error Estimation ~~~~~~~~~~~~~~%
	%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
	ISTATUS(Nstp) = ISTATUS(Nstp)+1;
	TMP = zeros(N,1);
    rr=rkE(1:rkS);
    TMP = (rr*Z)';
	
	[TMP,ISTATUS]=	SDIRK_Solve(H,N,E1L,E1U,TMP,ISTATUS);
	Err = SDIRK_ErrorNorm (N,TMP,SCAL);
	%/*~~~~> Computation of new step size Hnew */

	Fac = FacSafe *Err^(-ONE/rkELO);
	Fac = max (FacMin, min(FacMax,Fac));
    Hnew = H*Fac;
	%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	%/*~~~>  Accept/Reject step    */
	%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	if(Err < ONE)
		FirstStep = 0;
		ISTATUS(Nacc)=ISTATUS(Nacc)+1;
		T=T+H;
        Y=Y+(rkD(1:rkS)*Z)';
		SCAL = SDIRK_ErrorScale(N,ITOL,AbsTol,RelTol,Y);
		Hnew = Tdirection*min(abs(Hnew),Hmax);
		RSTATUS(Ntexit)=T;
		RSTATUS(Nhexit)=H;
		RSTATUS(Nhnew)=Hnew;
		if(Reject == 1)
			Hnew = Tdirection*min(abs(Hnew),abs(H));
		end;
		Reject = 0;
		if((T+Hnew/Qmin - T2)*Tdirection > ZERO)
			H=T2-T;
		else
			Hratio = Hnew/H;
			SkipLU = ((Theta<=ThetaMin) && (Hratio >=Qmin) && (Hratio<=Qmax));
			if(SkipLU==0)
				H=Hnew;
			end
		end
		SkipJac = (Theta<=ThetaMin);
		SkipJac = 0;
	else
		if(FirstStep==1 || Reject==1)
			H=FacRej*H;
		else
			H=Hnew;
		end;
		Reject=1;
		SkipJac =1;
		SkipLU=0;
		if(ISTATUS(Nacc)>=1)
			ISTATUS(Nrej)=ISTATUS(Nrej)+1;
		end
	end
	end
end

Ierr =1;
return;




function [SCAL] = SDIRK_ErrorScale(N,ITOL,AbsTol,RelTol,Y)
ZERO =0;
if(ITOL==ZERO)
    temp = RelTol(1)*Y;
    temp =(AbsTol(1)+abs(temp));
    SCAL=1./temp;
else
    SCAL=1./(AbsTol+RelTol.*abs(Y));
end
return;


function [Err] = SDIRK_ErrorNorm(N,Y,SCAL)
ZERO=0;
temp=Y.*Y.*SCAL.*SCAL;
Err=sum(temp);
Err = max(sqrt(Err/N),1e-10);
return


function SDIRK_ErrorMsg(code,T,H,Ierr)
Ierr = code;

disp('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('\nForced exit from Sdirk due to the following error:\n');

switch (code) 

case -1
	disp('--> Improper value for maximal no of steps');
	
case -2
	disp('--> Selected SDIRK method not implemented');
	
case -3
	printf('--> Hmin/Hmax/Hstart must be positive');
	
case -4
	printf('--> FacMin/FacMax/FacRej must be positive');
	
case -5
	printf('-> Improper tolerance values');
case -6
	disp('--> No of steps exceeds maximum bound');
case -7
	disp('--> Step size too small (T + H/10 = T) or H < Roundoff');
case -8
	printf('--> Matrix is repeatedly singular');

end;
return


function [E1L,E1U,ISTATUS,ISING]=PrepareMatrix(N,H,T,Y,FJAC,SkipJac, SkipLU,Reject,ISTATUS,OdeJacobian1)
global ZERO Njac  Ndec  Nsng 
global rkGamma 


ConsecutiveSng = 0;
ISING =1;
while (ISING ~= ZERO)
    HGammaInv = 1/(H*rkGamma);
    if(SkipJac == 0)
        FJAC = OdeJacobian1(T,Y);
        ISTATUS(Njac)=ISTATUS(Njac)+1;
    end;
    E=eye(N);
    E=HGammaInv*E-FJAC;
    
    [E1L,E1U] = lu(E);
    if(min(rank(E1L),rank(E1U))<N)
        ISING = 1;
    else
        ISING = 0;
    end;
    ISTATUS(Ndec) =ISTATUS(Ndec)+1;
    if(ISING ~= 0 )
        ISTATUS(Nsng)=ISTATUS(Nsng)+1;
        ConsecutiveSng=ConsecutiveSng+1;
        
        if(ConsecutiveSng >= 6)
            return;
        end
        H = 0.5*H;
        SkipJac =1;
        SkipLU = 0;
        Reject =1 ;
    end;
end
return;


function [RHS,ISTATUS] = SDIRK_Solve(H,N,E1L,E1U,RHS,ISTATUS)
global  Nsol 
global  rkGamma 

HGammaInv = 1/(H*rkGamma);
RHS = HGammaInv*RHS;
size(RHS);
size(E1L);
size(E1U);
temp = E1L\RHS;
RHS = E1U\temp;
ISTATUS(Nsol)=ISTATUS(Nsol)+1;
return

function Sdirk4a()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
global ZERO ONE 
global  S4A  sdMethod rkS rkGamma rkA rkB rkELO rkBhat
global rkC rkD rkE rkTheta rkAlpha
sdMethod = S4A;

%/* Number of stages */
rkS = 5;

%/* Method Coefficients */
rkGamma = 0.2666666666666666666666666666666667;

rkA(1,1) = 0.2666666666666666666666666666666667;
rkA(1,2) = 0.5000000000000000000000000000000000;
rkA(2,2) = 0.2666666666666666666666666666666667;
rkA(1,3) = 0.3541539528432732316227461858529820;
rkA(2,3) = -0.5415395284327323162274618585298197e-01;
rkA(3,3) = 0.2666666666666666666666666666666667;
rkA(1,4) = 0.8515494131138652076337791881433756e-01;
rkA(2,4) = -0.6484332287891555171683963466229754e-01;
rkA(3,4) = 0.7915325296404206392428857585141242e-01;
rkA(4,4) = 0.2666666666666666666666666666666667;
rkA(1,5) = 2.100115700566932777970612055999074;
rkA(2,5) = -0.7677800284445976813343102185062276;
rkA(3,5) = 2.399816361080026398094746205273880;
rkA(4,5) = -2.998818699869028161397714709433394;
rkA(5,5) = 0.2666666666666666666666666666666667;
rkB(1) = 2.100115700566932777970612055999074;
rkB(2) = -0.7677800284445976813343102185062276;
rkB(3) = 2.399816361080026398094746205273880;
rkB(4) = -2.998818699869028161397714709433394;
rkB(5) = 0.2666666666666666666666666666666667;

rkBhat(1) = 2.885264204387193942183851612883390;
rkBhat(2) = -0.1458793482962771337341223443218041;
rkBhat(3) = 2.390008682465139866479830743628554;
rkBhat(4) = -4.129393538556056674929560012190140;
rkBhat(5) = ZERO;

rkC(1) = 0.2666666666666666666666666666666667;
rkC(2) = 0.7666666666666666666666666666666667;
rkC(3) = 0.5666666666666666666666666666666667;
rkC(4) = 0.3661315380631796996374935266701191;
rkC(5) = ONE;

%/* Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
rkD(1) = ZERO;
rkD(2) = ZERO;
rkD(3) = ZERO;
rkD(4) = ZERO;
rkD(5) = ONE;

%/* Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
rkE(1) = -0.6804000050475287124787034884002302;
rkE(2) = 1.558961944525217193393931795738823;
rkE(3) = -13.55893003128907927748632408763868;
rkE(4) = 15.48522576958521253098585004571302;
rkE(5) = ONE;

%/* Local order of Err estimate */
rkELO = 4;

%/* h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
rkTheta(1,2) = 1.875000000000000000000000000000000;
rkTheta(1,3) = 1.708847304091539528432732316227462;
rkTheta(2,3) = -0.2030773231622746185852981969486824;
rkTheta(1,4) = 0.2680325578937783958847157206823118;
rkTheta(2,4) = -0.1828840955527181631794050728644549;
rkTheta(3,4) = 0.2968246986151577397160821594427966;
rkTheta(1,5) = 0.9096171815241460655379433581446771;
rkTheta(2,5) = -3.108254967778352416114774430509465;
rkTheta(3,5) = 12.33727431701306195581826123274001;
rkTheta(4,5) = -11.24557012450885560524143016037523;

%/* Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j} */
rkAlpha(1,2) = 2.875000000000000000000000000000000;
rkAlpha(1,3) = 0.8500000000000000000000000000000000;
rkAlpha(2,3) = 0.4434782608695652173913043478260870;
rkAlpha(1,4) = 0.7352046091658870564637910527807370;
rkAlpha(2,4) = -0.9525565003057343527941920657462074e-01;
rkAlpha(3,4) = 0.4290111305453813852259481840631738;
rkAlpha(1,5) = -16.10898993405067684831655675112808;
rkAlpha(2,5) = 6.559571569643355712998131800797873;
rkAlpha(3,5) = -15.90772144271326504260996815012482;
rkAlpha(4,5) = 25.34908987169226073668861694892683;

rkELO = 4.0;

% /* end Sdirk4a */

%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
function Sdirk4b()
%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
global ZERO ONE 
global  S4B sdMethod rkS rkGamma rkA rkB rkELO rkBhat
global rkC rkD rkE rkTheta rkAlpha

sdMethod = S4B;

%/* Number of stages */
rkS = 5;

%/* Method coefficients */
rkGamma = 0.25;

rkA(1,1) = 0.25;
rkA(1,2) = 0.5;
rkA(2,2) = 0.25;
rkA(1,3) = 0.34;
rkA(2,3) = -0.40e-01;
rkA(3,3) = 0.25;
rkA(1,4) = 0.2727941176470588235294117647058824;
rkA(2,4) = -0.5036764705882352941176470588235294e-01;
rkA(3,4) = 0.2757352941176470588235294117647059e-01;
rkA(4,4) = 0.25;
rkA(1,5) = 1.041666666666666666666666666666667;
rkA(2,5) = -1.020833333333333333333333333333333;
rkA(3,5) = 7.812500000000000000000000000000000;
rkA(4,5) = -7.083333333333333333333333333333333;
rkA(5,5) = 0.25;

rkB(1) = 1.041666666666666666666666666666667;
rkB(2) = -1.020833333333333333333333333333333;
rkB(3) = 7.812500000000000000000000000000000;
rkB(4) = -7.083333333333333333333333333333333;
rkB(5) = 0.250000000000000000000000000000000;

rkBhat(1) = 1.069791666666666666666666666666667;
rkBhat(2) = -0.894270833333333333333333333333333;
rkBhat(3) = 7.695312500000000000000000000000000;
rkBhat(4) = -7.083333333333333333333333333333333;
rkBhat(5) = 0.212500000000000000000000000000000;

rkC(1) = 0.25;
rkC(2) = 0.75;
rkC(3) = 0.55;
rkC(4) = 0.5;
rkC(5) = ONE;

%/* Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
rkD(1) = ZERO;
rkD(2) = ZERO;
rkD(3) = ZERO;
rkD(4) = ZERO;
rkD(5) = ONE;

%/* Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
rkE(1) = 0.5750;
rkE(2) = 0.2125;
rkE(3) = -4.6875;
rkE(4) = 4.2500;
rkE(5) = 0.1500;

%/* Local order of Err estimate */
rkELO = 4;

%/* h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
rkTheta(1,2) = 2.0;
rkTheta(1,3) = 1.680000000000000000000000000000000;
rkTheta(2,3) = -0.1600000000000000000000000000000000;
rkTheta(1,4) = 1.308823529411764705882352941176471;
rkTheta(2,4) = -0.1838235294117647058823529411764706;
rkTheta(3,4) = 0.1102941176470588235294117647058824;
rkTheta(1,5) = -3.083333333333333333333333333333333;
rkTheta(2,5) = -4.291666666666666666666666666666667;
rkTheta(3,5) = 34.37500000000000000000000000000000;
rkTheta(4,5) = -28.3333333333333333333333333333;

%/* Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j} */
rkAlpha(1,2) = 3.0;
rkAlpha(1,3) = 0.8800000000000000000000000000000000;
rkAlpha(2,3) = 0.4400000000000000000000000000000000;
rkAlpha(1,4) = 0.1666666666666666666666666666666667;
rkAlpha(2,4) = -0.8333333333333333333333333333333333e-01;
rkAlpha(3,4) = 0.9469696969696969696969696969696970;
rkAlpha(1,5) = -6.0;
rkAlpha(2,5) = 9.0;
rkAlpha(3,5) = -56.81818181818181818181818181818182;
rkAlpha(4,5) = 54.0;

 %/* end Sdirk4b */

%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
function Sdirk2a()
%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
global ZERO ONE 
global  S2A  sdMethod rkS rkGamma rkA rkB rkELO rkBhat
global rkC rkD rkE rkTheta rkAlpha

sdMethod = S2A;

%/* ~~~> Number of stages    */
rkS = 2;

%/* ~~~> Method coefficients    */
rkGamma = 0.2928932188134524755991556378951510;
rkA(1,1) = 0.2928932188134524755991556378951510;
rkA(1,2) = 0.7071067811865475244008443621048490;
rkA(2,2) = 0.2928932188134524755991556378951510;
rkB(1) = 0.7071067811865475244008443621048490;
rkB(2) = 0.2928932188134524755991556378951510;
rkBhat(1) = 0.6666666666666666666666666666666667;
rkBhat(2) = 0.3333333333333333333333333333333333;
rkC(1) = 0.292893218813452475599155637895151;
rkC(2) = ONE;

%/* ~~~> Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}    */
rkD(1) = ZERO;
rkD(2) = ONE;

%/* ~~~> Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}    */
rkE(1) = 0.4714045207910316829338962414032326;
rkE(2) = -0.1380711874576983496005629080698993;

%/* ~~~> Local order of Err estimate    */
rkELO = 2;

%/* ~~~> h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}    */
rkTheta(1,2) = 2.414213562373095048801688724209698;

%/* ~~~> Starting value for Newton iterations */
rkAlpha(1,2) = 3.414213562373095048801688724209698;

% /* end Sdirk2a */

%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
function Sdirk2b()
%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
global ZERO ONE 
global  S2B  sdMethod rkS rkGamma rkA rkB rkELO rkBhat
global rkC rkD rkE rkTheta rkAlpha

sdMethod = S2B;

%/* ~~~> Number of stages    */
rkS      = 2;

%/* ~~~> Method coefficients    */
rkGamma = 1.707106781186547524400844362104849;
rkA(1,1) = 1.707106781186547524400844362104849;
rkA(1,2) = -0.707106781186547524400844362104849;
rkA(2,2) = 1.707106781186547524400844362104849;
rkB(1) = -0.707106781186547524400844362104849;
rkB(2) = 1.707106781186547524400844362104849;
rkBhat(1) = 0.6666666666666666666666666666666667;
rkBhat(2) = 0.3333333333333333333333333333333333;
rkC(1) = 1.707106781186547524400844362104849;
rkC(2) = ONE;

%/* ~~~> Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}    */
rkD(1) = ZERO;
rkD(2) = ONE;

%/* ~~~> Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}    */
rkE(1) = -0.4714045207910316829338962414032326;
rkE(2) = 0.8047378541243650162672295747365659;

%/* ~~~> Local order of Err estimate    */
rkELO = 2;

%/* ~~~> h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}    */
rkTheta(1,2) = -0.414213562373095048801688724209698;

%/* ~~~> Starting value for Newton iterations */
rkAlpha(1,2) = 0.5857864376269049511983112757903019;

% /* end Sdirk2b */

%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
function Sdirk3a()
%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
global ZERO ONE 
global  S3A  sdMethod rkS rkGamma rkA rkB rkELO rkBhat
global rkC rkD rkE rkTheta rkAlpha

sdMethod = S3A;

%/* ~~~> Number of stages    */
rkS = 3;

%/* ~~~> Method coefficients    */
rkGamma = 0.2113248654051871177454256097490213;
rkA(1,1) = 0.2113248654051871177454256097490213;
rkA(1,2) = 0.2113248654051871177454256097490213;
rkA(2,2) = 0.2113248654051871177454256097490213;
rkA(1,3) = 0.2113248654051871177454256097490213;
rkA(2,3) = 0.5773502691896257645091487805019573;
rkA(3,3) = 0.2113248654051871177454256097490213;
rkB(1) = 0.2113248654051871177454256097490213;
rkB(2) = 0.5773502691896257645091487805019573;
rkB(3) = 0.2113248654051871177454256097490213;
rkBhat(1)= 0.2113248654051871177454256097490213;
rkBhat(2)= 0.6477918909913548037576239837516312;
rkBhat(3)= 0.1408832436034580784969504064993475;
rkC(1) = 0.2113248654051871177454256097490213;
rkC(2) = 0.4226497308103742354908512194980427;
rkC(3) = ONE;

%/* ~~~> Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}    */
rkD(1) = ZERO;
rkD(2) = ZERO;
rkD(3) = ONE;

%/* ~~~> Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}    */
rkE(1) = 0.9106836025229590978424821138352906;
rkE(2) = -1.244016935856292431175815447168624;
rkE(3) = 0.3333333333333333333333333333333333;

%/* ~~~> Local order of Err estimate    */
rkELO    = 2;

%/* ~~~> h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}    */
rkTheta(1,2) =  ONE;
rkTheta(1,3) = -1.732050807568877293527446341505872;
rkTheta(2,3) = 2.732050807568877293527446341505872;

%/* ~~~> Starting value for Newton iterations */
rkAlpha(1,2) = 2.0;
rkAlpha(1,3) = -12.92820323027550917410978536602349;
rkAlpha(2,3) = 8.83012701892219323381861585376468;

% /* end Sdirk3a */    


      
