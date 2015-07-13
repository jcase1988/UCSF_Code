function [params,Lik]=RL_behaviorfit_search(numruns,rewards,offered,choice,correct_choice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calls RL_behaviorfit.  It runs a search of
% parameters eta and b with random initialization.  
%
% params = optimized parameters
% Lik = minimum liklihood
%
%**NOTE** right now the number of eta and b parameters is hard-coded.  If you
% change the RL model to have a different number of eta and b's you need to
% change the search function to pass extra parameters to be optimized.
%
% eta = learning rate parameter -
%                       this controls how fast the model updates values based on feedback
% b = inverse temperature -
%                       this controls how consistenly the model choses the
%                       option with the higher value estimate vs. making a
%                       stochastic choice
%
% 
%  inputs:       numruns = number of random seeds - start with 100 - 500
%                   rewards = n x 1 array of feedback on each trial (in points)
%                   offered = n x 2 array of stimulus IDs (1 to 3), 2 for each trial (first option, second option).
%                   choice = n x 1 array of choices (1 = chose first stimulus, 2 = chose second stimulus).
%                   correct_choice = n x 1 array of stimulus IDs (1 to 3) of higher value stimulus in each pair of 'offered'.
%                               Note: if data includes a reversal, the 'correct_choice' for each pair will change at the block change.
% 
% elr 5/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lik=1000;
for k=1:numruns
    eta=rand;  %randomly initiate eta between 0 and 1
    b = rand*10; %randomly initiate b between 0 and 10
    [outparams,lik,exitflag] = fminsearch(@(x) RL_behaviorfit(x(1),x(2),rewards,offered,choice,correct_choice), [eta b]);
    if exitflag==1 && lik < Lik%if a solution was found that is better than the previous
        Lik = lik;
        params = outparams;
    end
end
