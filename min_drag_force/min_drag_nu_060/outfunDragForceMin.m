function stop = outfunDragForceMin(x,optimValues,state)
stop = false;

switch state
    case 'iter'
        % Make updates to plot or guis as needed
        if optimValues.iteration ~=0
            shape_current = shape3Dadjoint(x);
            shape_current.printresults;
        end
    case 'interrupt'
        % Probably no action here. Check conditions to see whether optimization should quit.
    case 'init'
        % Setup for plots or guis
    case 'done'
        % Cleanup of plots, guis, or final plot
    otherwise
end

end