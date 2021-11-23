# Reinforcement Learning based algorithm for massive MIMO radar with adaptive parameters

The present repository contains the Matlab code used to simulate a RL-based algorithm for a massive MIMO radar. The RL-based algorithm shapes the transmitted beampattern by selecting the weighting matrix of the transmitted waveforms based on the position of the detected targets. The ε and α parameters can be adaptively selected (default option). The quasi ε-greedy policy with target recovery is the default one. The user can compare the performance of the algorithm over one parameter, such as ε, α or the policy, while keeping all the others constant. 

## Mosek with CVX

In order to run the code the user has to download Mosek. <br />
All the instructions to install the software can be found at http://cvxr.com/cvx/doc/mosek.html

## Offline mode

To speed up the algorithm the user can compute the W matrices that are solution of the optimization problem associated to all the possible Omega sets and store them in a cube (Wcube). Then the algorithm doesn't have to solve the optimization problem at each iteration, but can select the correct W matrix from the cube. This mode of operation is called "offline" mode. Before running the code in "offline" mode the user has to run the W_opt_generator script to generate the cube.