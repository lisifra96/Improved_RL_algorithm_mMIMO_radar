# Reinforcement Learning based algorithm for massive MIMO radar with adaptive parameters

The present repository contains the Matlab code used to simulate a RL-based algorithm for a massive MIMO radar. The RL-based algorithm shapes the transmitted beampattern by selecting the weighting matrix of the transmitted waveforms based on the position of the detected targets. The ε and α parameters can be adaptively selected (default option). The quasi ε-greedy policy with target recovery is the default one. The user can compare the performance of the algorithm over one parameter, such as ε, α or the policy, while keeping all the others constant. 

## Offline mode

To speed up the algorithm the user can compute the W matrices that are solutions of the optimization problem associated with all the possible Omega sets and store them in a cube (Wcube). Then the algorithm doesn't have to solve the optimization problem at each iteration but can select the correct W matrix from the cube. This mode of operation is called "offline" mode. Before running the code in "offline" mode the user has to run the W_opt_generator script to generate the cube.

## UPDATE

The optimization problem in eq. (19) of the paper [1] and in eq. (6) of 
the paper [2] can be solved in closed form by using the procedure 
discussed in Sec. III.B of [3]. An extensive analysis of this point will be proposed in a paper in preparation.

The related Matlab function is `Closed_Form_W.m`. The files `W_opt_generator.m` and `getWfromTargetIndexes_online` have been updated accordingly. 
The updated `MonteCarlo_online.m` script exploits the parfor loop and the closed-form solution of the optimization problem to reduce the computational time. As a consequence, we recommend using the "online" mode. **The updated version of the code does not require the user to install cvx and Mosek anymore**

The function `Alg2v1.m` that solves the optimization problem using the cvx package is not used in the current version of the algorithm. If a user wants to change the optimization function and use cvx, he can update the objective function in `Alg2v1.m` and update the functions `W_opt_generator.m` and `getWfromTargetIndexes_online` accordingly.

[1] A. M. Ahmed, A. A. Ahmad, S. Fortunati, A. Sezgin, M. S. Greco and F. Gini, "A Reinforcement Learning Based Approach for Multitarget Detection in Massive MIMO Radar," in *IEEE Transactions on Aerospace and Electronic Systems*, vol. 57, no. 5, pp. 2622-2636, Oct. 2021.

[2] F. Lisi, S. Fortunati, M. S. Greco and F. Gini, "Enhancement of a State-of-the-Art RL-Based Detection Algorithm for Massive MIMO Radars," in *IEEE Transactions on Aerospace and Electronic Systems*, vol. 58, no. 6, pp. 5925-5931, Dec. 2022.

[3] P. Stoica, J. Li and Y. Xie, "On Probing Signal Design For MIMO Radar," in *IEEE Transactions on Signal Processing*, vol. 55, no. 8, pp. 4151-4161, Aug. 2007.
