# Safety Controller with CBF

Using Control Barrier Function for Shoulder Exoskeleton Robot

Now This CODE is not final version




# Version Check
03.25.2022 - Modified the closed form QP solver matrix. After that check it have same result with using QP solver package(CVXOPT)


# Theory Part

Control Barrier Fucntion

Barrier Function : Barrier function is a continuous function which value on a point increase to infinity as the point accproaches the boundary of the feasible region of an optimization problem
If C is the safety set
reciprocal barrier function : B(x) -> infinity as x -> boundary of C
zeroing barrier function : h(x) -> 0 as x -> boundary of C, I use this for this code.

Control Barrier Function 
Control brreier fucntion make control input which can make all states are in invariant safety sets with Optimization Porblem and Safety Constraints.

Solving Flow
First, Construct our control affine systems x_dot  = f(x) + g(x)u
Second, Construct barrier function h(x) which match with our safey conditions (Ex : Our safety set is made by joint limit and joint acceleration limit on human shoulder.)
Third, Construct safety constraint by h_dot = dh/dx * dx/dt = dh/dx * (f(x) + g(x)u)
dh/dx*f(x) = L_fh , dh/dx*g(x) = L_gh, alpha is the class kappa function.
Using this we can make constraints as, L_fh + L_gh >= -alpha(h(x))   
Fourth, make QP problem with safety constraints and control input constraints.



References : 
1) Ames, Aaron D., et al. "Control barrier function based quadratic programs for safety critical systems." IEEE Transactions on Automatic Control 62.8 (2016): 3861-3876.
2) https://dev10110.github.io/tech-notes/
3) https://en.wikipedia.org/wiki/Barrier_function
