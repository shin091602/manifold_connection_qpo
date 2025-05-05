function [value,isterminal,direction] = odestop_hetero_1(t,x,mu)

value = x(1)-1+mu; % dot{r}, the value that we want to be zero
isterminal = 1; % if isterminal == 1, the calcuration will stop when the event occurs.
direction = 1; % 2021-04-27  commment by Kyo Tetsu:  if it's 0 here, Apoaposis will also be outputed.
end