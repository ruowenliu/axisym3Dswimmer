function stop = outfunDragForceMin(x,optimValues,state)
stop = false;

global nu

switch state
    case 'iter'
        % Make updates to plot or guis as needed
        if optimValues.iteration ~=0

            shape_current = shape3Dadjoint(x);
            shape_current.printresults;

            fID = fopen(['./all_designvec_mindragforce_nu_' num2str(nu, '%.2e') '.txt'], 'a');
            fprintf(fID, '%.15f \n', x);
            fclose(fID);

        end
end

end