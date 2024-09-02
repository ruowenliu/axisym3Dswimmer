function stop = outfunMaxEffi(x,optimValues,state)
stop = false;

switch state
    case 'iter'
        if optimValues.iteration == 0
% %             previous_lssteplength = 0;
        else
            shape_current = shape3Dmaxefficiency(x);
            shape_current.printresults;
%             shape_current.plotshapesimple2;
%             drawnow
            %             cc = 0.05;
            %             if abs(shape_current.rvol-nu) > cc
            %                 stop = true;
            %                 fprintf('\n -> Stop because reduced volume is too far from target, criterion %g, current rvol %g. \n', cc, shape_current.rvol);
            %             end
        end
        %         if ~isempty(optimValues.directionalderivative)
        %             cc = 0.0001;
        %             if norm(optimValues.directionalderivative) < cc
        %                     stop = true;
        %                     fprintf('\n -> Stop due to directional derivative is less than criterion %g. \n',cc);
        %             end
        %         end
%         if ~isempty(optimValues.lssteplength)
%             if optimValues.lssteplength < previous_lssteplength
%                 stop = true;
%                 fprintf('\n --> Stop because Step-size decreases. -- \n');
%             else
%                 stop = false;
%                 previous_lssteplength = optimValues.lssteplength;
%             end
%         end
    case 'interrupt'
        % Probably no action here. Check conditions to see whether optimization should quit.
    case 'init'
        % Setup for plots or guis
    case 'done'
        % Cleanup of plots, guis, or final plot
    otherwise
end

end
