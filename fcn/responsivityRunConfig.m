function cfg = responsivityRunConfig(action, value)
%RESPONSIVITYRUNCONFIG Share output settings with script-based workflows.
% The exploratory scripts retain their own local switches when no override
% is active. reproduce_all uses this persistent configuration so that one
% entry point can control figure and table generation despite scripts
% clearing their workspaces.

    persistent override

    if nargin < 1 || isempty(action)
        action = "get";
    end

    switch lower(string(action))
        case "get"
            % No action required.
        case "set"
            if nargin < 2 || ~isstruct(value)
                error('responsivityRunConfig:InvalidConfiguration', ...
                    'The configuration override must be supplied as a struct.');
            end
            override = value;
        case "reset"
            override = [];
        otherwise
            error('responsivityRunConfig:UnknownAction', ...
                'Unknown configuration action: %s', action);
    end

    if isempty(override)
        cfg = struct( ...
            'OverrideOutputSwitches', false, ...
            'SaveFigures', false, ...
            'SaveTables', false, ...
            'CloseFiguresBeforeRun', true, ...
            'ValidationScenario', "");
    else
        cfg = override;
    end
end
