function [T, Y, RCNTRL, ICNTRL, RSTAT, ISTAT] = ...
    RK_Int(Function, Tspan, Y0, Options, RCNTRL, ICNTRL)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Implementation of Fully Implicit RK methods with the
%  Coefficients:
%               * Radau2A                                               
%               * Lobatto3C                                             
%               * Radau1A                                                
%               * Gauss                                             
%               * Lobatto3A                                              
%                                                                     
%    Solves the system y'=F(t,y) using a Fully Implicit RK method.
%
%    For details on Fully Implicit RK methods and their implementation consult:
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
%      RSTAT     - real output parameters (explained below)
%      ISTAT     - integer output parameters (explained below)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    RCNTRL and ICNTRL on input and output:
%
%    Note: For input parameters equal to zero the default values of the
%       corresponding variables are used.
%
%    
%    ICNTRL(1) = 0: AbsTol, RelTol are N-dimensional vectors
%              = 1: AbsTol, RelTol are scalars
%
%    ICNTRL(2)  -> selection of coefficients
%        = 0 :    Radau2A
%        = 1 :    Radau2A
%        = 2 :    Lobatto3C
%        = 3 :    Gauss
%        = 4 :    Radau1A
%        = 5 :    Lobatto3A
%       
%       
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
        RK_Integrate(N, Y(:,i-1), T(i-1), Tspan(i), ...
            Function, Jacobian, ATOL, RTOL, RCNTRL, ICNTRL,ISTATUS,RSTATUS);

    if IERR < 0
        error(['IRK exited with IERR=',num2str(IERR)]);
    end

end
cputime-t
RSTAT=RSTATUS;
ISTAT=ISTATUS;
%~~~>  Statistics on the work performed by the RK method
function [T,Y,Ierr,ISTATUS,RSTATUS] = RK_Integrate(N,Y0,TIN, TOUT,OdeFunction1,...
    OdeJacobian1,ATOL,RTOL,RCNTRL,ICNTRL,ISTATUS,RSTATUS)
global ZERO ONE Nfun Njac Nstp Nacc Nrej Ndec Nsol Nsng Ntexit Nhacc Nhnew
global R2A R1A L3C GAU L3A RKmax
ZERO =0;ONE=1;Nfun=1;Njac=2;Nstp=3;Nacc=4;Nrej=5;Ndec=6;Nsol=7;Nsng=8;
Ntexit=1;Nhacc=2;Nhnew=3;RKmax=3;R2A=1;R1A=2;L3C=3;GAU=4;L3A=5;
Ierr =0;
RSTATUS = zeros(20,1);
ISTATUS = zeros(20,1);
T1 = TIN;
T2 = TOUT;
Y(:,1)=Y0;
[T,Y,Ierr,RSTATUS,ISTATUS] = RungeKutta(N, T1,T2,Y0,RTOL,ATOL,RCNTRL,ICNTRL,RSTATUS,ISTATUS,Ierr,OdeFunction1,OdeJacobian1);





if(Ierr < 0)
    disp('Runge-kutta Unsuccessful exit at T=',num2str(TIN));
end
return

function [T,Y,Ierr,RSTATUS,ISTATUS]=RungeKutta(N,T1,T2,Y,RelTol,AbsTol,RCNTRL,...
    ICNTRL,ISTATUS,RSTATUS,Ierr,OdeFunction1,OdeJacobian1)
global ZERO  Roundoff SdirkError
T=T1;

if (ICNTRL(1)==0)
	    ITOL=1;
else
	    ITOL=0;
end

if(ICNTRL(9)==0)
    SdirkError = 0;
else
    SdirkError = 1;
end

switch(ICNTRL(2))
	case{0,1}
		Radau2A_Coefficients();
	case(2)
		Lobatto3C_Coefficients(); 
	case(3)
		Gauss_Coefficients();
	case(4)
		Radau1A_Coefficients();
	case(5)
		Lobatto3A_Coefficients();
	otherwise
		disp(['ICNTRL(2)=',num2str(ICNTRL(2))]);
		RK_ErrorMsg(-13,Tin,0,IERR);
		return;
end

if (ICNTRL(3)==0)
		Max_no_steps=200000;
else 
		Max_no_steps=ICNTRL(3);
        if(Max_no_steps<0)
            disp(['User-selected max no. of steps is -ve and is ',num2str(ICNTRL(3))]);
        end
end

if (ICNTRL(4)==0)
	NewtonMaxit=8;
else
	NewtonMaxit=ICNTRL(4);
	if NewtonMaxit<=0
		disp(['User-selected max no. of newton iterations is -ve and is ',num2str(ICNTRL(5))]);
		RK_ErrorMsg(-2, T, ZERO, IERR);
	end
end

if (ICNTRL(5)==0)
    StartNewton=1;
else
    StartNewton =0;
end

if(ICNTRL(10)==0)
    Gustafsson = 1;
else
    Gustafsson = 0;
end
Roundoff = eps;
if(RCNTRL(1)==ZERO)
    Hmin = ZERO;
else
    Hmin = min(abs(RCNTRL(1)),abs(T2 -T1));
end
if(RCNTRL(2)==ZERO)
    Hmax=abs(T2-T1);
else
	Hmax= min(abs(RCNTRL(2)),abs(T2-T1));
end
if RCNTRL(3)==0
	Hstart=0;
else
	Hstart = min(abs(RCNTRL(3)),abs(T2-T1));
end

%FacMin-- Lower bound on step decrease factor
if (RCNTRL(4) == 0)
	FacMin = 0.2;
else 
   	FacMin = RCNTRL(4);
end
%FacMax--Upper bound on step increase factor
if RCNTRL(5)==0
	FacMax=8.0;
else
	FacMax=RCNTRL(5);
end
%FacRej--step decrease factor after 2 consecutive rejections
if RCNTRL(6)==0
	FacRej=0.1;
else
	FacRej=RCNTRL(6);
end

%Facsafe:by which the new step is slightly smaller than the 
	%predicted value
if RCNTRL(7)==0
	FacSafe=0.9;
else
	FacSafe=RCNTRL(7);
end

if ( (FacMax < 1) ||( FacMin > 1) || (FacSafe <= 1.0e-03) || (FacSafe >= 1) )
	disp(['\n RCNTRL(5)=',num2str(RCNTRL(4)) ' RCNTRL[5]=',num2str(RCNTRL(5)) ' RCNTRL[6]=', num2str(RCNTRL(6)) ' RCNTRL[7]=',num2str(RCNTRL(7))]);
	RK_ErrorMsg(-4, T, ZERO, Ierr);
end

%ThetaMin: Decides whether the jacobian should be recomputed
if RCNTRL(8)==0
	ThetaMin=1.0e-03;
else
	ThetaMin=RCNTRL(8);
end
if (ThetaMin <= 0.0 || ThetaMin >= 1.0)
	disp(['RCNTRL[8]=', num2str(RCNTRL(8))]);
	RK_ErrorMsg(-5, Tin, ZERO, Ierr);
end

if RCNTRL(9)==0
	NewtonTol = 3.0e-02;
else 
	NewtonTol = RCNTRL(9);
	if NewtonTol<=Roundoff
		disp(['RCNTRL(9)=',num2str(RCNTRL(9))]);
		RK_ErrorMsg(-6, Tin, ZERO,Ierr);
	end
end

%Qmin and Qmax: If Qmin < Hnew/Hold < Qmax then step size is constant
if RCNTRL(10)==0
	Qmin=1;
else
	Qmin = RCNTRL(10);
end
if RCNTRL(11)==0
	Qmax=1.2;
else
	Qmax=RCNTRL(11);
end
if (Qmin > 1 || Qmax <1)
	disp(['Qmin=', num2str(RCNTRL(10)) ' Qmax=',num2str(RCNTRl(11))]);
	RK_ErrorMsg(-7, T, ZERO, Ierr);
end

%check if tolerances are reasonable
if ITOL==0
	if AbsTol(1)<=0 ||  RelTol(1) <= 10.0*Roundoff
		disp(['AbsTol=', num2str(AbsTol(1)) ' RelTol=',num2str(RelTol(1))]);
		RK_ErrorMsg(-8, T, ZERO, Ierr);
	end
else
	for i=1:N 
        if (( AbsTol(i) <= 0) || (RelTol(i) <= 10.0*Roundoff) ) 
            disp(['AbsTol(',num2str(i) ')= ',num2str(AbsTol(i))]);
			disp(['RelTol(',num2str(i) ')=  ',RelTol(i)]);
			RK_ErrorMsg(-8, T, ZERO, Ierr);
        end
	end 
end

if(Ierr < 0)
    return;
end
%Call the core Integrator
[T,Y,Ierr,ISTATUS,RSTATUS] = RK_Integrator(N,T1,T2,Y,Hstart,Hmin,...
                            Hmax,Roundoff,AbsTol,RelTol,...
                            ITOL,Max_no_steps,StartNewton,NewtonTol,ThetaMin,...
                            FacSafe,FacMin,FacMax,FacRej,Qmin,Qmax,NewtonMaxit,ISTATUS,RSTATUS,Gustafsson,Ierr,OdeFunction1,OdeJacobian1);
return;



function [T,Y,Ierr,ISTATUS,RSTATUS] = RK_Integrator(N,T1,T2,Y,Hstart,Hmin,Hmax,...
                                      Roundoff,AbsTol,RelTol,ITOL,Max_no_steps,...
                                      StartNewton,NewtonTol,ThetaMin,FacSafe,FacMin,...
                                      FacMax,FacRej,Qmin,Qmax,NewtonMaxit,ISTATUS,...
                                      RSTATUS,Gustafsson,Ierr,OdeFunction1,OdeJacobian1)
global SdirkError ZERO Nfun Njac Nsng  Nstp ONE rkMethod L3A rkBgam
global rkTheta rkGamma rkD rkELO Nacc Ntexit Nhacc Nhnew Hacc
T = T1;
CONT=zeros(N,3);
Tdirection = sign(T2-T);
H = min (max(abs(Hmin),abs(Hstart)),Hmax);
if(abs(H)<=10*Roundoff)
    H=1.0e-6;
end
H =sign(Tdirection)*abs(H);
Hold = H;
Reject = 0;
FirstStep = 1;
SkipJac = 0;
SkipLU = 0;
if((T+H*1.0001-T2)*Tdirection>=ZERO)
    H = T2 -T;
end
Nconsecutive = 0;
SCAL = RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y);
while((T2-T)*Tdirection-Roundoff>ZERO)
    F0 = OdeFunction1(T,Y);
    ISTATUS(Nfun)=ISTATUS(Nfun)+1;
    if(SkipLU == 0)
        if(SkipJac==0)
            FJAC = OdeJacobian1(T,Y);
            ISTATUS(Njac)=ISTATUS(Njac)+1;
        end
        [E1L,E1U,E2L,E2U,ISING,ISTATUS]=RK_Decomp(N,H,FJAC,ISTATUS);
        if(ISING ~= 0)
            ISTATUS(Nsng)=ISTATUS(Nsng)+1;
            Nconsecutive=Nconsecutive+1;
            if(Nconsecutive>=5)
                RK_ErrorMsg(-12,T,H,Ierr);
            end
            H=H*0.5;
            Reject = 1;
            SkipJac = 1;
            SkipLU = 0;
            continue;
        else
            Nconsecutive = 0;
        end
    end
    ISTATUS(Nstp)=ISTATUS(Nstp)+1;
    if (ISTATUS(Nstp)>Max_no_steps)
		disp(['Max number of time steps is =', num2str(Max_no_steps)]);
		RK_ErrorMsg(-9,T,H,IERR);
    end
    if (0.1*abs(H)<=abs(T)*Roundoff)
		RK_ErrorMsg(-10,T,H,Ierr);
    end
    %Loop for simplified newton iterations

	%Starting values for newton iterations
    if((FirstStep==1)||(StartNewton==0))
			Z1=zeros(N,1);
			Z2=zeros(N,1);
			Z3=zeros(N,1);
	else
		evaluate=2;
		[CONT, Z1, Z2, Z3]=RK_Interpolate(evaluate, N, H, Hold,Z1,Z2,Z3,CONT);
    end
    
    %/*~~~> Initializations for Newton iteration */
    
    NewtonDone = 0;
    Fac = 0.5;
    for NewtonIter=1:NewtonMaxit
        %Preparing RHS
        [R1,R2,R3]=RK_PrepareRHS(N,T,H,Y,F0,Z1,Z2,Z3,OdeFunction1);
        %Solve the linear systems
        [R1,R2,R3,ISTATUS] = RK_Solve(N,H,E1L,E1U,E2L,E2U,R1,R2,R3,ISTATUS);
        %Find NewtonIncrement
        NewtonIncrement=sqrt((RK_ErrorNorm(N,SCAL,R1)*RK_ErrorNorm(N,SCAL,R1)...
									+RK_ErrorNorm(N,SCAL,R2)*RK_ErrorNorm(N,SCAL,R2)...
									+RK_ErrorNorm(N,SCAL,R3)*RK_ErrorNorm(N,SCAL,R3))...
							/3.0);
        if(NewtonIter == 1)
            Theta = abs(ThetaMin);
            NewtonRate = 2.0;
        else
            Theta = NewtonIncrement/NewtonIncrementOld;
            if(Theta < 0.99)
                NewtonRate = Theta/(ONE-Theta);
            else
                break
            end
            if(NewtonIter<NewtonMaxit)
                NewtonPredictedErr = NewtonIncrement*Theta^(NewtonMaxit-NewtonIter)/(ONE - Theta);
                if(NewtonPredictedErr>=NewtonTol)
                    Qnewton = min(10.0,NewtonPredictedErr/NewtonTol);
                    Fac = 0.8*Qnewton^(-ONE/(1+NewtonMaxit-NewtonIter));
                    break
                end
            end
        end
        NewtonIncrementOld=max(NewtonIncrement,Roundoff);
        Z1=Z1-R1; 
        Z2=Z2-R2;
        Z3=Z3-R3;
        NewtonDone = (NewtonRate*NewtonIncrement<=NewtonTol);
        if(NewtonDone==1)
            break;
        end
        if(NewtonIter==NewtonMaxit)
            disp('Slow or no convergence in newton iterations');
			disp('Max no of newton iterations reached');
        end
    end
    if(NewtonDone == 0)
        H=Fac*H;
        Reject=1;
        SkipJac=1;
		SkipLU=0;
        continue;
    end
    %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    %/*~~~> SDIRK Stage */
    %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if(SdirkError==1)
        Z4=Z3;
        G=zeros(N,1);
        if(rkMethod ~= L3A)
            G=G+rkBgam(1)*H*F0;
        end
        G=G+Z1*rkTheta(1);
        G=G+Z2*rkTheta(2);
        G=G+Z3*rkTheta(3);
        NewtonDone = 0;
        Fac = 0.5;
        for NewtonIter=1:NewtonMaxit
            TMP=Y+Z4;
            R4=OdeFunction1(T+H,TMP);
            ISTATUS(Nfun)=ISTATUS(Nfun)+1;
            R4=R4-(rkGamma/H)*Z4;
            R4=R4+(rkGamma/H)*G;
            [R4,ISTATUS] = RK_KppSolve(E1L,E1U,R4,ISTATUS);
            NewtonIncrement = RK_ErrorNorm(N,SCAL,R4);
            if(NewtonIter==1)
                ThetaSD=abs(ThetaMin);
                NewtonRate = 2.0;
            else
                ThetaSD = NewtonIncrement/NewtonIncrementOld;
                if(ThetaSD<0.99)
                    NewtonRate = ThetaSD/(ONE - ThetaSD);
                    NewtonPredictedErr=NewtonIncrement*ThetaSD^((NewtonMaxit-NewtonIter)/(1-ThetaSD));
                    if(NewtonPredictedErr>=NewtonTol)
						Qnewton=min(10.0,NewtonPredictedErr/NewtonTol);
						Fac=0.8*Qnewton^(-1/(1+NewtonMaxit-NewtonIter));
						break;
                    end
                else
                    break;
                end
            end
            NewtonIncrementOld = NewtonIncrement;
            Z4 = Z4+R4;
            %Check error in NewtonIterations
            NewtonDone=(NewtonRate*NewtonIncrement<=NewtonTol);
            if(NewtonDone==1)
                break;
            end
        end
        if(NewtonDone == 0)
            H=Fac*H;
            Reject =1;
            SkipJac=1;
            SkipLU = 0;
            continue;
        end
    end
    %/*~~~> End of Simplified SDIRK Newton iterations */

    %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    %/*~~~> Error estimation */
    %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if(SdirkError==1)
        R4=zeros(N,1);
        if(rkMethod==L3A)
            R4=H*rkF(1)*F0;
            if(rkF(2)~=ZERO)
                R4=R4+Z1*rkF(1);
            end
            
            if(rkF(3)~=ZERO)
                R4=R4+Z2*rkF(2);
            end
            
            if(rkF(4)~=ZERO)
                R4=R4+Z3*rkF(3);
            end
            
            TMP = Y+Z4;
            R1 = OdeFunction1(T+H,TMP);
            R4=R4+H*rkBgam(5)*R1;
        else
            if(rkD(1)~=0)
                R4=R4+rkD(1)*Z1;
            end
            if(rkD(2)~=ZERO)
                R4=R4+rkD(2)*Z2;
            end
            if(rkD(3)~=ZERO)
                R4=R4+rkD(3)*Z3;
            end
            R4=R4-Z4;
        end
		Err=RK_ErrorNorm(N,SCAL,R4);
    else
        [Err,ISTATUS]=RK_ErrorEstimate(N,H,T,Y,F0,E1L,E1U,Z1,Z2,Z3,SCAL,...
							FirstStep,Reject,ISTATUS);
    end
    Fac=Err^(-1/rkELO)*min(FacSafe,(1+2*NewtonMaxit)/(NewtonIter+2*NewtonMaxit));
	Fac=min(FacMax,max(FacMin,Fac));
	Hnew=Fac*H;
    %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    %/*~~~> Accept/reject step */
    %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    %Accept
    if(Err < ONE)
        FirstStep =0;
        ISTATUS(Nacc)=ISTATUS(Nacc)+1;
        if(Gustafsson==ONE)
            if(ISTATUS(Nacc)>1)
                FacGus=FacSafe*(H/Hacc)*(Err*Err/ErrOld)^(-0.25);
				FacGus=min(FacMax,max(FacMin,FacGus));
				Fac=min(Fac,FacGus);
				Hnew=Fac*H;
            end
            Hacc=H;
			ErrOld=max(1.0e-02,Err);
        end
        Hold=H;
		T=T+H;
		if(rkD(1)~=0)
			Y=Y+Z1*rkD(1);
		end
		if(rkD(2)~=0)
			Y=Y+Z2*rkD(2);
		end;
		if(rkD(3)~=0)
			Y=Y+Z3*rkD(3);
		end;
        %Construct the solution quadratic interpolant Q(c_i)=Z_i,i=1:3
		if(StartNewton==1)
			evaluate=1;
			[CONT, Z1, Z2, Z3]=RK_Interpolate(evaluate, N, H, Hold,Z1,Z2,Z3,CONT);
		end
		SCAL=RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y);
		RSTATUS(Ntexit)=T;
		RSTATUS(Nhnew)=Hnew;
		RSTATUS(Nhacc)=H;
		Hnew=Tdirection*min(max(abs(Hnew),Hmin),Hmax);
        if(Reject==1)
			Hnew=Tdirection*min(abs(Hnew),abs(H));
        end;
		Reject=0;
        if((T + Hnew/Qmin - T2)*Tdirection>=0)
			H=T2-T;
		else
			Hratio=Hnew/H;
			SkipLU=((Theta<=ThetaMin) && (Hratio>=Qmin) && (Hratio<=Qmax));
			if(SkipLU==0)
				H=Hnew;
			end;
            
        end;
         SkipJac=0;
		
        else
			if((FirstStep==1) || (Reject==1))
				H=FacRej*H;
			else
				H=Hnew;
			end;
			Reject=1;
			SkipJac=1;
			SkipLU=0;
			if(ISTATUS(Nacc)>=1)
				ISTATUS(Nrej)=ISTATUS(Nrej)+1;
			end;
		end;
	end;
	Ierr=1;

return;

function [Err,ISTATUS]=RK_ErrorEstimate(N,H,T,Y,F0,E1L,E1U,Z1,Z2,Z3,SCAL,...
							FirstStep,Reject,ISTATUS)
global rkE rkMethod R1A GAU L3A
HrkE1=rkE(2)/H;
HrkE2=rkE(3)/H;
HrkE3=rkE(4)/H;
F2=HrkE1*Z1+HrkE2*Z2+HrkE3*Z3;
TMP = rkE(1)*F0+F2;
[TMP,ISTATUS]=RK_KppSolve(E1L,E1U,TMP,ISTATUS);
if((rkMethod == R1A) || (rkMethod == GAU) || (rkMethod == L3A))
	[TMP,ISTATUS]=RK_KppSolve(E1L,E1U,TMP,ISTATUS);
end

if (rkMethod == GAU)
	[TMP,ISTATUS]=RK_KppSolve(E1L,E1U,TMP,ISTATUS);
end
Err=RK_ErrorNorm(N,SCAL,TMP);
if(Err < 1)
    return;
end
if((FirstStep==1 )||(Reject==1))
    TMP=TMP+Y;
    F1=OdeFunction1(T,TMP);
    ISTATUS(Nfun)=ISTATUS(Nfun)+1;
    TMP=F1+F2;
    [TMP,ISTATUS]=RK_KppSolve(E1L,E1U,TMP,ISTATUS);
    Err=RK_ErrorNorm(N,SCAL,TMP);
end
return

function RK_ErrorMsg(Code, T, H, IERR)
	
	Code=IERR;
	disp('Forced to exit from RungeKutta due to the following error:');
	switch(Code)
		case(-1)
			disp('--> Improper value for maximal no of steps');
		case (-2)
      		disp('--> Improper value for maximal no of Newton iterations');
   		case (-3)
      		disp('--> Hmin/Hmax/Hstart must be positive');
   		case (-4)
      		disp('--> Improper values for FacMin/FacMax/FacSafe/FacRej');
   		case (-5)
      		disp('--> Improper value for ThetaMin');
   		case (-6)
      		disp('--> Newton stopping tolerance too small');
   		case (-7)
      		disp('--> Improper values for Qmin, Qmax');
   		case (-8)
      		disp('--> Tolerances are too small');
   		case (-9)
      		disp('--> No of steps exceeds maximum bound');
   		case (-10)
      		disp('--> Step size too small: (T + 10*H = T) or H < Roundoff');
   		case (-11)
      		disp('--> Matrix is repeatedly singular');
   		case (-12)
      		disp('--> Non-convergence of Newton iterations');
   		case (-13)
      		disp('--> Requested RK method not implemented');
   		otherwise
      		disp(['Unknown Error code:', num2str(Code)]);
	end;

   	disp(['T=',num2str(T) 'H=',num2str(H)]);
return;

function SCAL=RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y)
	if(ITOL==0)
        SCAL=1./(AbsTol(1)+RelTol(1)*abs(Y));
    else
        SCAL=1./(AbsTol+RelTol.*abs(Y));
	end;
return;

function ErrorNorm = RK_ErrorNorm(N,SCAL,DY)

    temp=DY.*DY.*SCAL.*SCAL;
    ErrorNorm = sum(temp);
ErrorNorm=max(sqrt(ErrorNorm/N),1.0e-10);

function [CONT, Z1, Z2, Z3]=RK_Interpolate(action, N, H, Hold,Z1,Z2,Z3,CONT)
	global rkC
    

	if(action==1) %Make
		den=(rkC(3)-rkC(2))*(rkC(2)-rkC(1))*(rkC(1)-rkC(3));
		for i=1:N
			CONT(i,1)=(-rkC(3)*rkC(3)*rkC(2)*Z1(i)...
					  +Z3(i)*rkC(2)*rkC(1)*rkC(1)...
	   				  +rkC(2)*rkC(2)*rkC(3)*Z1(i)...
	   				  -rkC(2)*rkC(2)*rkC(1)*Z3(i)...
	   				  +rkC(3)*rkC(3)*rkC(1)*Z2(i)...
	   		          -Z2(i)*rkC(3)*rkC(1)*rkC(1))/den-Z3(i);
			CONT(i,2) = -( rkC(1)*rkC(1)*(Z3(i)-Z2(i))...
	   		+ rkC(2)*rkC(2)*(Z1(i)-Z3(i))...
	   		+ rkC(3)*rkC(3)*(Z2(i)-Z1(i)) )/den;

			CONT(i,3) = ( rkC(1)*(Z3(i)-Z2(i))...
	   		+ rkC(2)*(Z1(i)-Z3(i))...
	   		+ rkC(3)*(Z2(i)-Z1(i)) )/den;
		end;
    end
    
    if (action==2) %Eval
		r=H/Hold;
		x1=1+rkC(1)*r;
		x2 = 1 + rkC(2)*r;
		x3 = 1 + rkC(3)*r;
		for i=1:N
		
			Z1(i) = CONT(i,1)+x1*(CONT(i,2)+x1*CONT(i,3));	
			Z2(i) = CONT(i,1)+x2*(CONT(i,2)+x2*CONT(i,3));	
			Z3(i) = CONT(i,1)+x3*(CONT(i,2)+x3*CONT(i,3));
		end;
    end;
return;



function [R1,R2,R3]=RK_PrepareRHS(N,T,H,Y,FO,Z1,Z2,Z3,OdeFunction1)
	global rkMethod 
	global rkA rkC 
    global L3A
	R1=Z1;
	R2=Z2;
	R3=Z3;
	if(rkMethod==L3A)
		R1=R1-H*rkA(1,1)*FO;
		R2=R2-H*rkA(2,1)*FO;
		R3=R3-H*rkA(3,1)*FO;
	end
	TMP=Y+Z1;
	F=OdeFunction1(T+rkC(1)*H,TMP);	
	R1=R1-H*rkA(1,1)*F;
	R2=R2-H*rkA(2,1)*F;
	R3=R3-H*rkA(3,1)*F;

	TMP=Y+Z2;
	F=OdeFunction1(T+rkC(2)*H,TMP);
	R1=R1-H*rkA(1,2)*F;
	R2=R2-H*rkA(2,2)*F;
	R3=R3-H*rkA(3,2)*F;

	TMP=Y+Z3;
	F=OdeFunction1(T+rkC(3)*H,TMP);
	R1=R1-H*rkA(1,3)*F;
	R2=R2-H*rkA(2,3)*F;
	R3=R3-H*rkA(3,3)*F;
return;

function [E1L,E1U,E2L,E2U,ISING,ISTATUS]=RK_Decomp(N, H, FJAC,ISTATUS)
	global Ndec 
	global rkGamma rkAlpha rkBeta 
    ISING =0;
	Gamma = rkGamma / H;
   	Alpha = rkAlpha / H;
   	Beta  = rkBeta / H;
	E1=Gamma*eye(N);
    E1=E1-FJAC;


	
	[E1L,E1U]=lu(E1);
    ISTATUS(Ndec)=ISTATUS(Ndec)+1;
	if(det(E1L)==0 || det(E1U)==0)
		ISING=1;
	end;
    
	if(ISING ~= 0)
		
		return;
	end;
    E2R=complex(Alpha,Beta)*eye(N)-FJAC;
	
	[E2L,E2U]=lu(E2R);
    ISTATUS(Ndec)=ISTATUS(Ndec)+1;
	if(abs(det(E2L))==0 || abs(det(E2U))==0)
		ISING=1;
	end;
	if(ISING ~= 0)
		return;
	end;
    
return

function [R1,R2,R3,ISTATUS]=RK_Solve(N,H,E1L,E1U,E2L,E2U,R1,R2,R3,ISTATUS)
	global Nsol  
	global rkT rkTinvAinv 
	for i=1:N
		x1=R1(i)/H;
		x2=R2(i)/H;
		x3=R3(i)/H;
		R1(i) = rkTinvAinv(1,1)*x1 + rkTinvAinv(1,2)*x2 + rkTinvAinv(1,3)*x3;
		R2(i) = rkTinvAinv(2,1)*x1 + rkTinvAinv(2,2)*x2 + rkTinvAinv(2,3)*x3;
		R3(i) = rkTinvAinv(3,1)*x1 + rkTinvAinv(3,2)*x2 + rkTinvAinv(3,3)*x3;
	end;
	tmp = E1L\R1;
	R1 =  E1U\tmp;
	BCR=R2;
	BCI=R3;
	a1=complex(BCR,BCI);
	tmp = E2L\a1;
	BC = E2U\tmp;
	R2=real(BC);
	R3=imag(BC);
	for i=1:N
		x1 = R1(i);
		x2 = R2(i);
		x3 = R3(i);
        R1(i) = rkT(1,1)*x1 + rkT(1,2)*x2 + rkT(1,3)*x3;
		R2(i) = rkT(2,1)*x1 + rkT(2,2)*x2 + rkT(2,3)*x3;
		R3(i) = rkT(3,1)*x1 + rkT(3,2)*x2 + rkT(3,3)*x3;
	end
	ISTATUS(Nsol)=ISTATUS(Nsol)+1;
return

function [R4,ISTATUS]=RK_KppSolve(E1L,E1U,R4,ISTATUS)
global Nsol
temp = E1L\R4;
R4=E1U\temp;
ISTATUS(Nsol)=ISTATUS(Nsol)+1;
return




function Radau2A_Coefficients ()

	global rkMethod SdirkError
	global rkT rkTinv rkTinvAinv rkAinvT rkA rkB rkC rkD rkE
	global rkBgam  rkTheta  rkGamma rkAlpha rkBeta rkELO
	global R2A 

	
	if(SdirkError==1)
		b0=0.2e-01;
	else
		b0=0.5e-01;
	end;
	rkMethod=R2A;
   	rkA(1,1) = 1.968154772236604258683861429918299e-01;
	rkA(1,2) = -6.55354258501983881085227825696087e-02;
	rkA(1,3) = 2.377097434822015242040823210718965e-02;
	rkA(2,1) = 3.944243147390872769974116714584975e-01;
	rkA(2,2) = 2.920734116652284630205027458970589e-01;
	rkA(2,3) = -4.154875212599793019818600988496743e-02;
	rkA(3,1) = 3.764030627004672750500754423692808e-01;
	rkA(3,2) = 5.124858261884216138388134465196080e-01;
	rkA(3,3) = 1.111111111111111111111111111111111e-01;

	rkB(1) = 3.764030627004672750500754423692808e-01;
	rkB(2) = 5.124858261884216138388134465196080e-01;
	rkB(3) = 1.111111111111111111111111111111111e-01;

	rkC(1) = 1.550510257216821901802715925294109e-01;
	rkC(2) = 6.449489742783178098197284074705891e-01;
	rkC(3) = 1;

	% New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j 
	rkD(1) =0;
	rkD(2) =0;
	rkD(3) =1;

	% Classical error estimator: */
	% H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */
	rkE(1) = 1*b0;
	rkE(2) = -10.04880939982741556246032950764708*b0;
	rkE(3) = 1.382142733160748895793662840980412*b0;
	rkE(4) = -.3333333333333333333333333333333333*b0;

    %	/* Sdirk error estimator */
	rkBgam(1) = b0;
	rkBgam(2) = .3764030627004672750500754423692807...
			-1.558078204724922382431975370686279*b0;
	rkBgam(3) = 0.8914115380582557157653087040196118*b0...
			+0.5124858261884216138388134465196077;
	rkBgam(4) = (-0.1637777184845662566367174924883037)...
			-0.3333333333333333333333333333333333*b0;
	rkBgam(5) = 0.2748888295956773677478286035994148;

	% H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
	rkTheta(1) = (-1.520677486405081647234271944611547)...
			-10.04880939982741556246032950764708*b0;
	rkTheta(2) = (2.070455145596436382729929151810376)...
			+1.382142733160748895793662840980413*b0;
	rkTheta(3) = -0.3333333333333333333333333333333333*b0...
			-.3744441479783868387391430179970741;

	%/* Local order of error estimator */
	if ( b0 == 0)
		rkELO  = 6.0;
	else
		rkELO  = 4.0;
	end;

	%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 % ~~~> Diagonalize the RK matrix:
	  %rkTinv * inv(rkA) * rkT =
	  %|  rkGamma      0           0     |
	  %|      0      rkAlpha   -rkBeta   |
	  %|      0      rkBeta     rkAlpha  |
	  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	rkGamma = 3.637834252744495732208418513577775;
	rkAlpha = 2.681082873627752133895790743211112;
	rkBeta  = 3.050430199247410569426377624787569;

	rkT(1,1) = 9.443876248897524148749007950641664e-02;
	rkT(1,2) = -1.412552950209542084279903838077973e-01;
	rkT(1,3) = -3.00291941051474244918611170890539e-02;
	rkT(2,1) = 2.502131229653333113765090675125018e-01;
	rkT(2,2) = 2.041293522937999319959908102983381e-01;
	rkT(2,3) = 3.829421127572619377954382335998733e-01;
	rkT(3,1) =  1;
	rkT(3,2) =  1;
	rkT(3,3) =  0;

	rkTinv(1,1) = 4.178718591551904727346462658512057;
	rkTinv(1,2) = 3.27682820761062387082533272429617e-01;
	rkTinv(1,3) = 5.233764454994495480399309159089876e-01;
	rkTinv(2,1) = -4.178718591551904727346462658512057;
	rkTinv(2,2) = -3.27682820761062387082533272429617e-01;
	rkTinv(2,3) = 4.766235545005504519600690840910124e-01;
	rkTinv(3,1) = -5.02872634945786875951247343139544e-01;
	rkTinv(3,2) = 2.571926949855605429186785353601676;
	rkTinv(3,3) = -5.960392048282249249688219110993024e-01;

	rkTinvAinv(1,1) = 1.520148562492775501049204957366528e+01;
	rkTinvAinv(1,2) = 1.192055789400527921212348994770778;
	rkTinvAinv(1,3) = 1.903956760517560343018332287285119;
	rkTinvAinv(2,1) = -9.669512977505946748632625374449567;
	rkTinvAinv(2,2) = -8.724028436822336183071773193986487;
	rkTinvAinv(2,3) = 3.096043239482439656981667712714881;
	rkTinvAinv(3,1) = -1.409513259499574544876303981551774e+01;
	rkTinvAinv(3,2) = 5.895975725255405108079130152868952;
	rkTinvAinv(3,3) = -1.441236197545344702389881889085515e-01;

	rkAinvT(1,1) = .3435525649691961614912493915818282;
	rkAinvT(1,2) = -.4703191128473198422370558694426832;
	rkAinvT(1,3) = .3503786597113668965366406634269080;
	rkAinvT(2,1) = .9102338692094599309122768354288852;
	rkAinvT(2,2) = 1.715425895757991796035292755937326;
	rkAinvT(2,3) = .4040171993145015239277111187301784;
	rkAinvT(3,1) = 3.637834252744495732208418513577775;
	rkAinvT(3,2) = 2.681082873627752133895790743211112;
	rkAinvT(3,3) = -3.050430199247410569426377624787569;
return;

function Lobatto3C_Coefficients()
	
	global rkMethod SdirkError
	global rkT rkTinv rkTinvAinv rkAinvT rkA rkB rkC rkD rkE
	global rkBgam rkBhat rkTheta rkGamma rkAlpha rkBeta rkELO
	global L3C 

	

	rkMethod=L3C;
	if(SdirkError==1)
		b0=0.2;
	else
		b0=0.5;
	end;
	rkA(1,1) = .1666666666666666666666666666666667;
	rkA(1,2) = -.3333333333333333333333333333333333;
	rkA(1,3) = .1666666666666666666666666666666667;
	rkA(2,1) = .1666666666666666666666666666666667;
	rkA(2,2) = .4166666666666666666666666666666667;
	rkA(2,3) = -.8333333333333333333333333333333333e-01;
	rkA(3,1) = .1666666666666666666666666666666667;
	rkA(3,2) = .6666666666666666666666666666666667;
	rkA(3,3) = .1666666666666666666666666666666667;

	rkB(1) = .1666666666666666666666666666666667;
	rkB(2) = .6666666666666666666666666666666667;
	rkB(3) = .1666666666666666666666666666666667;

	rkC(1) = 0;
	rkC(2) = 0.5;
	rkC(3) = 1;

	%/* Classical error estimator, embedded solution: */
	rkBhat(1) = b0;
	rkBhat(2) = .16666666666666666666666666666666667-b0;
	rkBhat(3) = .66666666666666666666666666666666667;
	rkBhat(4) = .16666666666666666666666666666666667;

%	/* New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j */
	rkD(1) = 0;
	rkD(2) = 0;
	rkD(3) = 1;

%	/* Classical error estimator: */
%	/* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */
	rkE(1) = .3808338772072650364017425226487022*b0;
	rkE(2) = (-1.142501631621795109205227567946107)*b0;
	rkE(3) = (-1.523335508829060145606970090594809)*b0;
	rkE(4) = .3808338772072650364017425226487022*b0;

%	/* Sdirk error estimator */
	rkBgam(1) = b0;
	rkBgam(2) = .1666666666666666666666666666666667-1.0*b0;
	rkBgam(3) = .6666666666666666666666666666666667;
	rkBgam(4) = (-.2141672105405983697350758559820354);
	rkBgam(5) = .3808338772072650364017425226487021;

%	/* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
	rkTheta(1) = -3.0*b0-.3808338772072650364017425226487021;
	rkTheta(2) = -4.0*b0+1.523335508829060145606970090594808;
	rkTheta(3) = (-.142501631621795109205227567946106)+b0;

%	/* Local order of error estimator */
	if (b0 == 0)
			rkELO  = 5.0;
	else
			rkELO  = 4.0;

	%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	%  ~~~> Diagonalize the RK matrix:
	%  rkTinv * inv(rkA) * rkT =
	%  |  rkGamma      0           0     |
	%  |      0      rkAlpha   -rkBeta   |
	%  |      0      rkBeta     rkAlpha  |
	%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	rkGamma = 2.625816818958466716011888933765284;
	rkAlpha = 1.687091590520766641994055533117359;
	rkBeta  = 2.508731754924880510838743672432351;

	rkT(1,2) = 1;
	rkT(1,2) = 1;
	rkT(1,3) = 0;
	rkT(2,1) = .4554100411010284672111720348287483;
	rkT(2,2) = -.6027050205505142336055860174143743;
	rkT(2,3) = -.4309321229203225731070721341350346;
	rkT(3,1) = 2.195823345445647152832799205549709;
	rkT(3,2) = -1.097911672722823576416399602774855;
	rkT(3,3) = .7850032632435902184104551358922130;

	rkTinv(1,1) = .4205559181381766909344950150991349;
	rkTinv(1,2) = .3488903392193734304046467270632057;
	rkTinv(1,3) = .1915253879645878102698098373933487;
	rkTinv(2,1) = .5794440818618233090655049849008650;
	rkTinv(2,2) = -.3488903392193734304046467270632057;
	rkTinv(2,3) = -.1915253879645878102698098373933487;
	rkTinv(3,1) = -.3659705575742745254721332009249516;
	rkTinv(3,2) = -1.463882230297098101888532803699806;
	rkTinv(3,3) = .4702733607340189781407813565524989;

	rkTinvAinv(1,1) = 1.104302803159744452668648155627548;
	rkTinvAinv(1,2) = .916122120694355522658740710823143;
	rkTinvAinv(1,3) = .5029105849749601702795812241441172;
	rkTinvAinv(2,1) = 1.895697196840255547331351844372453;
	rkTinvAinv(2,2) = 3.083877879305644477341259289176857;
	rkTinvAinv(2,3) = -1.502910584974960170279581224144117;
	rkTinvAinv(3,1) = .8362439183082935036129145574774502;
	rkTinvAinv(3,2) = -3.344975673233174014451658229909802;
	rkTinvAinv(3,3) = .312908409479233358005944466882642;

	rkAinvT(1,1) = 2.625816818958466716011888933765282;
	rkAinvT(1,2) = 1.687091590520766641994055533117358;
	rkAinvT(1,3) = -2.508731754924880510838743672432351;
	rkAinvT(2,1) = 1.195823345445647152832799205549710;
	rkAinvT(2,2) = -2.097911672722823576416399602774855;
	rkAinvT(2,3) = .7850032632435902184104551358922130;
	rkAinvT(3,1) = 5.765829871932827589653709477334136;
	rkAinvT(3,2) = .1170850640335862051731452613329320;
	rkAinvT(3,3) = 4.078738281412060947659653944216779;
end;
return;

	% /* Lobatto3C_Coefficients */
function Gauss_Coefficients()
	
	global rkMethod 
	global rkT rkTinv rkTinvAinv rkAinvT rkA rkB rkC rkD rkE
	global rkBgam rkBhat rkTheta rkGamma rkAlpha rkBeta rkELO
	global GAU 


   rkMethod = GAU;
   
   b0 = 0.1;

   %/* The coefficients of the Gauss method */
   rkA(1,1) = .1388888888888888888888888888888889;
   rkA(1,2) = -.359766675249389034563954710966045e-01;
   rkA(1,3) = .97894440153083260495800422294756e-02;
   rkA(2,1) = .3002631949808645924380249472131556;
   rkA(2,2) = .2222222222222222222222222222222222;
   rkA(2,3) = -.224854172030868146602471694353778e-01;
   rkA(3,1) = .2679883337624694517281977355483022;
   rkA(3,2) = .4804211119693833479008399155410489;
   rkA(3,3) = .1388888888888888888888888888888889;

   rkB(1) = .2777777777777777777777777777777778;
   rkB(2) = .4444444444444444444444444444444444;
   rkB(3) = .2777777777777777777777777777777778;

   rkC(1) = .1127016653792583114820734600217600;
   rkC(2) = .5000000000000000000000000000000000;
   rkC(3) = .8872983346207416885179265399782400;

%   /* Classical error estimator, embedded solution: */
   rkBhat(1) = b0;
   rkBhat(2) = -1.4788305577012361475298775666303999*b0...
   		 +.27777777777777777777777777777777778;
   rkBhat(3) = .44444444444444444444444444444444444...
    	         +.66666666666666666666666666666666667*b0;
   rkBhat(4) = -.18783610896543051913678910003626672*b0...
  	         +.27777777777777777777777777777777778;

 %  /* New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j */
   rkD(1) = (.1666666666666666666666666666666667e+01);
   rkD(2) = (-.1333333333333333333333333333333333e+01);
   rkD(3) = (.1666666666666666666666666666666667e+01);

   %/* Classical error estimator: */
   %/* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */
   rkE(1) = .2153144231161121782447335303806954*b0;
   rkE(2) = (-2.825278112319014084275808340593191)*b0;
   rkE(3) = .2870858974881495709929780405075939*b0;
   rkE(4) = (-.4558086256248162565397206448274867e-01)*b0;

   %/* Sdirk error estimator */
   rkBgam(1) = 0;
   rkBgam(2) = .2373339543355109188382583162660537;
   rkBgam(3) = .5879873931885192299409334646982414;
   rkBgam(4) = -.4063577064014232702392531134499046e-01;
   rkBgam(5) = .2153144231161121782447335303806955;

   %/* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta(1) = (-2.594040933093095272574031876464493);
   rkTheta(2) = 1.824611539036311947589425112250199;
   rkTheta(3) = .1856563166634371860478043996459493;

   %/* ELO = local order of classical error estimator */
   rkELO = 4.0;

   %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   %~~~> Diagonalize the RK matrix:
   % rkTinv * inv(rkA) * rkT =
   %           |  rkGamma      0           0     |
   %           |      0      rkAlpha   -rkBeta   |
   %           |      0      rkBeta     rkAlpha  |
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = 4.644370709252171185822941421408064;
   rkAlpha = 3.677814645373914407088529289295970;
   rkBeta  = 3.508761919567443321903661209182446;

   rkT(1,1) = .7215185205520017032081769924397664e-01;
   rkT(1,2) = -.8224123057363067064866206597516454e-01;
   rkT(1,3) = -.6012073861930850173085948921439054e-01;
   rkT(2,1) = .1188325787412778070708888193730294;
   rkT(2,2) = .5306509074206139504614411373957448e-01;
   rkT(2,3) = .3162050511322915732224862926182701;
   rkT(3,1) = 1;
   rkT(3,2) = 1;
   rkT(3,3) = 0;

   rkTinv(1,1) = 5.991698084937800775649580743981285;
   rkTinv(1,2) = 1.139214295155735444567002236934009;
   rkTinv(1,3) = .4323121137838583855696375901180497;
   rkTinv(2,1) = -5.991698084937800775649580743981285;
   rkTinv(2,2) = -1.139214295155735444567002236934009;
   rkTinv(2,3) = .5676878862161416144303624098819503;
   rkTinv(3,1) = -1.246213273586231410815571640493082;
   rkTinv(3,2) = 2.925559646192313662599230367054972;
   rkTinv(3,3) = -.2577352012734324923468722836888244;

   rkTinvAinv(1,1) = 27.82766708436744962047620566703329;
   rkTinvAinv(1,2) = 5.290933503982655311815946575100597;
   rkTinvAinv(1,3) = 2.007817718512643701322151051660114;
   rkTinvAinv(2,1) = (-17.66368928942422710690385180065675);
   rkTinvAinv(2,2) = (-14.45491129892587782538830044147713);
   rkTinvAinv(2,3) = 2.992182281487356298677848948339886;
   rkTinvAinv(3,1) = (-25.60678350282974256072419392007303);
   rkTinvAinv(3,2) = 6.762434375611708328910623303779923;
   rkTinvAinv(3,3) = 1.043979339483109825041215970036771;

   rkAinvT(1,1) = .3350999483034677402618981153470483;
   rkAinvT(1,2) = -.5134173605009692329246186488441294;
   rkAinvT(1,3) = .6745196507033116204327635673208923e-01;
   rkAinvT(2,1) = .5519025480108928886873752035738885;
   rkAinvT(2,2) = 1.304651810077110066076640761092008;
   rkAinvT(2,3) = .9767507983414134987545585703726984;
   rkAinvT(3,1) = 4.644370709252171185822941421408064;
   rkAinvT(3,2) = 3.677814645373914407088529289295970;
   rkAinvT(3,3) = -3.508761919567443321903661209182446;
return;
 %/* Gauss_Coefficients */

%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%	The coefficients of the 3-stage Gauss method
%	(given to ~30 accurate digits)
%	The parameter b3 can be chosen by the user
%	to tune the error estimator
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
function Radau1A_Coefficients()
	global rkMethod 
	global rkT rkTinv rkTinvAinv rkAinvT rkA rkB rkC rkD rkE
	global rkBgam rkBhat rkTheta rkGamma rkAlpha rkBeta rkELO
	global  R1A 

	global ZERO 


%   /* The coefficients of the Radau1A method */
   b0 = 0.1;
   rkMethod = R1A;
   rkA(1,1) =  .1111111111111111111111111111111111;
   rkA(1,2) = (-.1916383190435098943442935597058829);
   rkA(1,3) = (.8052720793239878323318244859477174e-01);
   rkA(2,1) = .1111111111111111111111111111111111;
   rkA(2,2) = .2920734116652284630205027458970589;
   rkA(2,3) = (-.481334970546573839513422644787591e-01);
   rkA(3,1) = .1111111111111111111111111111111111;
   rkA(3,2) = .5370223859435462728402311533676479;
   rkA(3,3) = .1968154772236604258683861429918299;

   rkB(1) = .1111111111111111111111111111111111;
   rkB(2) = .5124858261884216138388134465196080;
   rkB(3) = .3764030627004672750500754423692808;

   rkC(1) = 0;
   rkC(2) = .3550510257216821901802715925294109;
   rkC(3) = .8449489742783178098197284074705891;

  % /* Classical error estimator, embedded solution: */
   rkBhat(1) = b0;
   rkBhat(2) = .11111111111111111111111111111111111-b0;
   rkBhat(3) = .51248582618842161383881344651960810;
   rkBhat(4) = .37640306270046727505007544236928079;

   %/* New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j */
   rkD(1) = .3333333333333333333333333333333333;
   rkD(2) = -.8914115380582557157653087040196127;
   rkD(3) = 1.558078204724922382431975370686279;

   %/* Classical error estimator: */
   %/* H* Sum (b_j-bhat_j) f(Z_j) = H*E(0)*F(0) + Sum E_j Z_j */
   rkE(1) = .2748888295956773677478286035994148*b0;
   rkE(2) = -1.374444147978386838739143017997074*b0;
   rkE(3) = -1.335337922441686804550326197041126*b0;
   rkE(4) = .235782604058977333559011782643466*b0;

   %/* Sdirk error estimator */
   rkBgam(1) = ZERO;
   rkBgam(2) = (.1948150124588532186183490991130616e-01);
   rkBgam(3) = .7575249005733381398986810981093584;
   rkBgam(4) = (-.518952314149008295083446116200793e-01);
   rkBgam(5) = .2748888295956773677478286035994148;

   %/* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta(1) = (-1.224370034375505083904362087063351);
   rkTheta(2) = .9340045331532641409047527962010133;
   rkTheta(3) = .4656990124352088397561234800640929;

   %/* ELO = local order of classical error estimator */
   rkELO = 4.0;

   %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   %~~~> Diagonalize the RK matrix:
   % rkTinv * inv(rkA) * rkT =
   %           |  rkGamma      0           0     |
   %           |      0      rkAlpha   -rkBeta   |
   %           |      0      rkBeta     rkAlpha  |
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = 3.637834252744495732208418513577775;
   rkAlpha = 2.681082873627752133895790743211112;
   rkBeta  = 3.050430199247410569426377624787569;

   rkT(1,1) = .424293819848497965354371036408369;
   rkT(1,2) = -.3235571519651980681202894497035503;
   rkT(1,3) = -.522137786846287839586599927945048;
   rkT(2,1) = .57594609499806128896291585429339e-01;
   rkT(2,2) = .3148663231849760131614374283783e-02;
   rkT(2,3) = .452429247674359778577728510381731;
   rkT(3,1) = 1;
   rkT(3,2) = 1;
   rkT(3,3) = 0;

   rkTinv(1,1) = 1.233523612685027760114769983066164;
   rkTinv(1,2) = 1.423580134265707095505388133369554;
   rkTinv(1,3) = .3946330125758354736049045150429624;
   rkTinv(2,1) = -1.233523612685027760114769983066164;
   rkTinv(2,2) = -1.423580134265707095505388133369554;
   rkTinv(2,3) = .6053669874241645263950954849570376;
   rkTinv(3,1) = -.1484438963257383124456490049673414;
   rkTinv(3,2) = 2.038974794939896109682070471785315;
   rkTinv(3,3) = -.544501292892686735299355831692542e-01;

   rkTinvAinv(1,1) = 4.487354449794728738538663081025420;
   rkTinvAinv(1,2) = 5.178748573958397475446442544234494;
   rkTinvAinv(1,3) = 1.435609490412123627047824222335563;
   rkTinvAinv(2,1) = -2.854361287939276673073807031221493;
   rkTinvAinv(2,2) = -1.003648660720543859000994063139137e+01;
   rkTinvAinv(2,3) = 1.789135380979465422050817815017383;
   rkTinvAinv(3,1) = -4.160768067752685525282947313530352;
   rkTinvAinv(3,2) = 1.124128569859216916690209918405860;
   rkTinvAinv(3,3) = 1.700644430961823796581896350418417;

   rkAinvT(1,1) = 1.543510591072668287198054583233180;
   rkAinvT(1,2) = -2.460228411937788329157493833295004;
   rkAinvT(1,3) = -.412906170450356277003910443520499;
   rkAinvT(2,1) = .209519643211838264029272585946993;
   rkAinvT(2,2) = 1.388545667194387164417459732995766;
   rkAinvT(2,3) = 1.20339553005832004974976023130002;
   rkAinvT(3,1) = 3.637834252744495732208418513577775;
   rkAinvT(3,2) = 2.681082873627752133895790743211112;
   rkAinvT(3,3) = -3.050430199247410569426377624787569;
return;
% /* Radau1A_Coefficients */

%/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%	The coefficients of the 4-stage Lobatto-3A method
%	(given to ~30 accurate digits)
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
function Lobatto3A_Coefficients()
	
	global rkMethod 
	global rkT rkTinv rkTinvAinv rkAinvT rkA rkB rkC rkD rkE
	global rkBgam rkBhat rkTheta rkF rkGamma rkAlpha rkBeta rkELO
	global L3A


   %/* The coefficients of the Lobatto-3A method */
   rkMethod = L3A;

   rkA(1,1) = 0;
   rkA(1,2) = 0;
   rkA(1,3) = 0;
   rkA(1,4) = 0;
   rkA(2,1) = .11030056647916491413674311390609397;
   rkA(2,2) = .1896994335208350858632568860939060;
   rkA(2,3) = -.339073642291438837776604807792215e-01;
   rkA(2,4) = .1030056647916491413674311390609397e-01;
   rkA(3,1) = .73032766854168419196590219427239365e-01;
   rkA(3,2) = .4505740308958105504443271474458881;
   rkA(3,3) = .2269672331458315808034097805727606;
   rkA(3,4) = -.2696723314583158080340978057276063e-01;
   rkA(4,1) = .83333333333333333333333333333333333e-01;
   rkA(4,2) = .4166666666666666666666666666666667;
   rkA(4,3) = .4166666666666666666666666666666667;
   rkA(4,4) = .8333333333333333333333333333333333e-01;

   rkB(1) = .83333333333333333333333333333333333e-01;
   rkB(2) = .4166666666666666666666666666666667;
   rkB(3) = .4166666666666666666666666666666667;
   rkB(4) = .8333333333333333333333333333333333e-01;

   rkC(1) = 0;
   rkC(2) = .2763932022500210303590826331268724;
   rkC(3) = .7236067977499789696409173668731276;
   rkC(4) = 1;

   %/* New solution: H*Sum B_j*f(Z_j) = Sum D_j*Z_j */
   rkD(1) = 0;
   rkD(2) = 0;
   rkD(3) = 0;
   rkD(4) = 1;

   %/* Classical error estimator, embedded solution: */
   rkBhat(1) = .90909090909090909090909090909090909e-01;
   rkBhat(2) = .39972675774621371442114262372173276;
   rkBhat(3) = .43360657558711961891219070961160058;
   rkBhat(4) = .15151515151515151515151515151515152e-01;

   %/* Classical error estimator: */
   %/* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */

   rkE(1) = .1957403846510110711315759367097231e-01;
   rkE(2) = -.1986820345632580910316020806676438;
   rkE(3) = .1660586371214229125096727578826900;
   rkE(4) = -.9787019232550553556578796835486154e-01;

   %/* Sdirk error estimator: */
   rkF(1) = 0;
   rkF(2) = -.66535815876916686607437314126436349;
   rkF(3) = 1.7419302743497277572980407931678409;
   rkF(4) = -1.2918865386966730694684011822841728;

   %/* ELO = local order of classical error estimator */
   rkELO = 4.0;

   %/* Sdirk error estimator: */
   rkBgam(1) = .2950472755430528877214995073815946e-01;
   rkBgam(2) = .5370310883226113978352873633882769;
   rkBgam(3) = .2963022450107219354980459699450564;
   rkBgam(4) = -.7815248400375080035021681445218837e-01;
   rkBgam(5) = .2153144231161121782447335303806956;

   %/* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta(1) = 0.0;
   rkTheta(2) = -.6653581587691668660743731412643631;
   rkTheta(3) = 1.741930274349727757298040793167842;
   rkTheta(4) = -.291886538696673069468401182284174;


   %/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   %~~~> Diagonalize the RK matrix:
   % rkTinv * inv(rkA) * rkT =
   %          |  rkGamma      0           0     |
   %           |      0      rkAlpha   -rkBeta   |
   %           |      0      rkBeta     rkAlpha  |
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = 4.644370709252171185822941421408063;
   rkAlpha = 3.677814645373914407088529289295968;
   rkBeta  = 3.508761919567443321903661209182446;

   rkT(1,1) = .5303036326129938105898786144870856e-01;
   rkT(1,2) = -.7776129960563076320631956091016914e-01;
   rkT(1,3) = .6043307469475508514468017399717112e-02;
   rkT(2,1) = .2637242522173698467283726114649606;
   rkT(2,2) = .2193839918662961493126393244533346;
   rkT(2,3) = .3198765142300936188514264752235344;
   rkT(3,1) = 1;
   rkT(3,2) = 0;
   rkT(3,3) = 0;

   rkTinv(1,1) = 7.695032983257654470769069079238553;
   rkTinv(1,2) = -.1453793830957233720334601186354032;
   rkTinv(1,3) = .6302696746849084900422461036874826;
   rkTinv(2,1) = -7.695032983257654470769069079238553;
   rkTinv(2,2) = .1453793830957233720334601186354032;
   rkTinv(2,3) = .3697303253150915099577538963125174;
   rkTinv(3,1) = -1.066660885401270392058552736086173;
   rkTinv(3,2) = 3.146358406832537460764521760668932;
   rkTinv(3,3) = -.7732056038202974770406168510664738;

   rkTinvAinv(1,1) = 35.73858579417120341641749040405149;
   rkTinvAinv(1,2) = -.675195748578927863668368190236025;
   rkTinvAinv(1,3) = 2.927206016036483646751158874041632;
   rkTinvAinv(2,1) = -24.55824590667225493437162206039511;
   rkTinvAinv(2,2) = -10.50514413892002061837750015342036;
   rkTinvAinv(2,3) = 4.072793983963516353248841125958369;
   rkTinvAinv(3,1) = -30.92301972744621647251975054630589;
   rkTinvAinv(3,2) = 12.08182467154052413351908559269928;
   rkTinvAinv(3,3) = -1.546411207640594954081233702132946;

   rkAinvT(1,1) = .2462926658317812882584158369803835;
   rkAinvT(1,2) = -.2647871194157644619747121197289574;
   rkAinvT(1,3) = .2950720515900466654896406799284586;
   rkAinvT(2,1) = 1.224833192317784474576995878738004;
   rkAinvT(2,2) = 1.929224190340981580557006261869763;
   rkAinvT(2,3) = .4066803323234419988910915619080306;
   rkAinvT(3,1) = 4.644370709252171185822941421408064;
   rkAinvT(3,2) = 3.677814645373914407088529289295968;
   rkAinvT(3,3) = -3.508761919567443321903661209182446;
return;
		%/* Lobatto3A_Coefficients */

