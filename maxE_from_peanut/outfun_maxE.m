function stop = outfun_maxE(x,optimValues,state)
stop = false;

global fixed_RN fixed_RS fixed_ZN fixed_ZS
global figK nu design_vec_iteration

switch state
    case 'iter'
        if optimValues.iteration ~=0
            design_vec = fixpolesreverse(x,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS);            
            design_vec_iteration = [design_vec_iteration, design_vec];
            shape_current = shape3Dmaxefficiency2(design_vec);
            shape_current.plotblack;
            shape_current.printresults;

            figK = figK+1;
            saveas(gcf, ['./iteration_' num2str(figK) '.png']);
            close all

            cc = 1;
            if (shape_current.rvol-nu) > cc
                stop = true;
                fprintf('\nStop due to nu too large from target, criterion %g. \n',cc);
            end
            if (nu-shape_current.rvol) > cc
                stop = true;
                fprintf('\nStop due to nu too small from target, criterion %g. \n',cc);
            end
        end
        if ~isempty(optimValues.directionalderivative)
            cc = 0.0001;
            if norm(optimValues.directionalderivative) < cc
                    stop = true;
                    fprintf('\nStop due to directionalderivative small, criterion %g. \n',cc);
            end
        end
        if ~isempty(optimValues.lssteplength)
            cc = 10;
            if optimValues.lssteplength > cc
                stop = true;
                fprintf('\nStop due to Step-size too large, criterion %g. \n',cc);
            end
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
