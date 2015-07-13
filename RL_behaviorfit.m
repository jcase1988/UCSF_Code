function [lik,Pr, statevals_cum]=RL_behaviorfit(eta,b,rewards,offered,choice,correct_choice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the probability of choosing the best option (Pr)
% and the value estimate for each stimulus (statevals_cum) at every trial, based on 
% free parameters eta and b.
%
% Called by RL_behaviorfit_search
%
% eta = learning rate parameter -
%                       this controls how fast the model updates values based on feedback
% b = inverse temperature -
%                       this controls how consistenly the model choses the
%                       option with the higher value estimate vs. making a
%                       stochastic choice
%
% other outputs:  lik = negative log liklihood of the softmax equation
% against actual choice behavior.  This is the parameter to be minimized in
% the optimization step
%
% other inputs: rewards = n x 1 array of feedback on each trial (in points)
%                   offered = n x 2 array of stimulus IDs (1 to 3), 2 for each trial (first option, second option).
%                   choice = n x 1 array of choices (1 = chose first stimulus, 2 = chose second stimulus).
%                   correct_choice = n x 1 array of stimulus IDs (1 to 3) of higher value stimulus in each pair of 'offered'.
%                               Note: if data includes a reversal, the 'correct_choice' for each pair will change at the block change.
% 
% elr 5/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


trialn=length(choice);
statevals=[0 0 0];  %initialize values at zero
statevals_cum=zeros(length(choice),3);
Pr=zeros(trialn,1);
llik = 0;

%make the variable offeredpair
offeredpair=NaN(size(choice));
for j=1:length(offeredpair)
    if ismember(1,offered(j,:)) && ismember(2,offered(j,:))
        offeredpair(j)=1;
    elseif ismember(1,offered(j,:)) && ismember(3,offered(j,:))
        offeredpair(j)=2;
    elseif ismember(2,offered(j,:)) && ismember(3,offered(j,:))
        offeredpair(j)=3;
    end
end

for j=1:trialn
    chosen=offered(j,choice(j));  %stim ID chosen 
    statevals_offered = statevals(offered(j,:));
    
    %calculate the probability of choosing the correct option
    Pr(j)=(exp(b*statevals(correct_choice(j))))/((exp(b*statevals(offered(j,1))))+(exp(b*statevals(offered(j,2)))));
    
    %compute the liklihood
    ebv = exp(b.*statevals_offered);
    %llik = llik + b*statevals(chosen) - log(sum(ebv)) - (eta^2 + b^2);
    llik = llik + b*statevals(chosen) - log(sum(ebv)) - 1*(eta + b/10);
    
    
    %%%
    %An alternative model with separate betas for each offer pair:
    %Pr(j)=(exp(b(offeredpair(j))*statevals(correct_choice(j))))/((exp(b(offeredpair(j))*statevals(offered(j,1))))+(exp(b(offeredpair(j))*statevals(offered(j,2)))));
    %ebv = exp(b(offeredpair(j)).*statevals_offered);
    %llik = llik + b(offeredpair(j))*statevals(chosen) - log(sum(ebv));
    %%%
    
    %calculte reward prediction error; output if desired
    RPE(j)=rewards(j) - statevals(chosen);
    
  
    %update the value estimates, store in statevals_cum for output if desired
    statevals(chosen) = statevals(chosen) + eta*(rewards(j) - statevals(chosen));
    statevals_cum(j,:)=statevals;
    
    %%%
    %
    % Below is an alternative model for value updating.  It has separate learning rates for 
    % positive and negative feedback.  To use this model, eta must have two
    % elements.
    %
    % Other alternative models could use separate etas for positive and
    % negative RPEs; separate betas for different offer pairs, etc.
    %
    % if rewards(j)>=0
    %       statevals(chosen) = statevals(chosen) + eta(1)*(rewards(j) - statevals(chosen));
    % elseif rewards(j)<0
    %       statevals(chosen) = statevals(chosen) + eta(2)*(rewards(j) - statevals(chosen));
    % end
    %
    %%%


end
lik=-llik;

