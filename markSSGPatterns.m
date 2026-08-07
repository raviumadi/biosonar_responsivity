function [T, eventTable, summary] = markSSGPatterns(T, ssgParams)
% markSSGPatterns
% Add SSG pattern labels to a call-level timing table.
%
% This is an analysis function, not part of the simulator. It detects stable
% doublets/triplets/quadruplets using within-group IPI CV and flanking IPI
% boundary criteria.
%
% Inputs:
%   T         - call-level table from simulateResponsivityCore
%   ssgParams - optional struct with fields:
%       groupSizes         default [2 3 4]
%       maxWithinCV        default 0.05
%       boundaryMultiplier default 1.20
%
% Outputs:
%   T          - input table with SSG label columns added
%   eventTable - one row per detected SSG event
%   summary    - diagnostic counts

    if nargin < 2 || isempty(ssgParams)
        ssgParams = struct();
    end
    ssgParams = setDefault(ssgParams, 'groupSizes', [2 3 4]);
    ssgParams = setDefault(ssgParams, 'maxWithinCV', 0.05);
    ssgParams = setDefault(ssgParams, 'boundaryMultiplier', 1.20);

    if ~ismember('IPI_ms', T.Properties.VariableNames)
        T.IPI_ms = 1000 * T.IPI_s;
    end

    T.IsSSG = zeros(height(T), 1);
    T.SSG_ID = zeros(height(T), 1);
    T.SSG_GroupSize_calls = zeros(height(T), 1);
    T.SSG_WithinCV = nan(height(T), 1);
    T.SSG_WithinMedianIPI_ms = nan(height(T), 1);
    T.SSG_LeftBoundaryRatio = nan(height(T), 1);
    T.SSG_RightBoundaryRatio = nan(height(T), 1);

    eventRows = struct([]);

    summary = struct();
    summary.StableCandidates = containers.Map('KeyType', 'double', 'ValueType', 'double');
    summary.LeftBoundaryPass = containers.Map('KeyType', 'double', 'ValueType', 'double');
    summary.RightBoundaryPass = containers.Map('KeyType', 'double', 'ValueType', 'double');
    summary.TrueEvents = containers.Map('KeyType', 'double', 'ValueType', 'double');

    for gs = ssgParams.groupSizes
        summary.StableCandidates(gs) = 0;
        summary.LeftBoundaryPass(gs) = 0;
        summary.RightBoundaryPass(gs) = 0;
        summary.TrueEvents(gs) = 0;
    end

    seqIDs = unique(T.SeqID);
    globalSSGID = 0;

    for s = 1:numel(seqIDs)
        sid = seqIDs(s);
        seqIdx = find(T.SeqID == sid);
        S = T(seqIdx, :);
        nCalls = height(S);
        assignedLocal = false(nCalls, 1);

        for groupSize = ssgParams.groupSizes
            if nCalls < groupSize + 2
                continue
            end

            for startCall = 2:(nCalls - groupSize)
                groupLocal = startCall:(startCall + groupSize - 1);

                if any(assignedLocal(groupLocal))
                    continue
                end

                groupIPI = S.IPI_ms(groupLocal);
                withinMedianIPI = median(groupIPI, 'omitnan');
                withinCV = std(groupIPI, 'omitnan') / withinMedianIPI;

                if withinCV > ssgParams.maxWithinCV
                    continue
                end

                summary.StableCandidates(groupSize) = summary.StableCandidates(groupSize) + 1;

                leftBoundaryIPI = S.IPI_ms(startCall - 1);
                rightBoundaryIPI = S.IPI_ms(startCall + groupSize);

                leftRatio = leftBoundaryIPI / withinMedianIPI;
                rightRatio = rightBoundaryIPI / withinMedianIPI;

                leftOK = leftRatio >= ssgParams.boundaryMultiplier;
                rightOK = rightRatio >= ssgParams.boundaryMultiplier;

                if leftOK
                    summary.LeftBoundaryPass(groupSize) = summary.LeftBoundaryPass(groupSize) + 1;
                end
                if rightOK
                    summary.RightBoundaryPass(groupSize) = summary.RightBoundaryPass(groupSize) + 1;
                end

                if ~(leftOK && rightOK)
                    continue
                end

                globalSSGID = globalSSGID + 1;
                summary.TrueEvents(groupSize) = summary.TrueEvents(groupSize) + 1;

                globalRows = seqIdx(groupLocal);
                T.IsSSG(globalRows) = 1;
                T.SSG_ID(globalRows) = globalSSGID;
                T.SSG_GroupSize_calls(globalRows) = groupSize;
                T.SSG_WithinCV(globalRows) = withinCV;
                T.SSG_WithinMedianIPI_ms(globalRows) = withinMedianIPI;
                T.SSG_LeftBoundaryRatio(globalRows) = leftRatio;
                T.SSG_RightBoundaryRatio(globalRows) = rightRatio;
                assignedLocal(groupLocal) = true;

                eventRow = makeEventRow(T(globalRows, :), S, sid, globalSSGID, ...
                    groupSize, startCall, groupLocal, withinMedianIPI, withinCV, ...
                    leftRatio, rightRatio);

                if isempty(eventRows)
                    eventRows = eventRow;
                else
                    eventRows(end + 1, 1) = eventRow; %#ok<AGROW>
                end
            end
        end
    end

    if isempty(eventRows)
        eventTable = table();
    else
        eventTable = struct2table(eventRows);
    end

    summary.NumCalls = height(T);
    summary.NumSSGCalls = sum(T.IsSSG == 1);
    summary.NumSSGEvents = height(eventTable);
    summary.PercentSSGCalls = 100 * mean(T.IsSSG == 1);
end

function eventRow = makeEventRow(E, S, seqID, ssgID, groupSize, startCall, ...
    groupLocal, stableIPI, withinCV, leftRatio, rightRatio)

    eventRow = struct();
    eventRow.SeqID = seqID;
    eventRow.SSG_ID = ssgID;
    eventRow.GroupSize = groupSize;
    eventRow.FirstCall = min(E.CallNumber);
    eventRow.LastCall = max(E.CallNumber);
    eventRow.SequenceLength_calls = height(S);
    eventRow.NormalisedPosition_percent = 100 * mean(groupLocal) / height(S);
    eventRow.StableIPI_ms = stableIPI;
    eventRow.WithinCV = withinCV;
    eventRow.LeftFlankRatio = leftRatio;
    eventRow.RightFlankRatio = rightRatio;
    eventRow.AnchorDistance_m = mean(E.AnchorDistance_m, 'omitnan');
    eventRow.NearestDistance_m = mean(E.NearestTargetDistance_m, 'omitnan');
    eventRow.AnchorSameWithinSSG = numel(unique(E.AnchorTargetID)) == 1;
    eventRow.AnchorSwitchWithinSSG = numel(unique(E.AnchorTargetID)) > 1;

    if ismember('ConditionName', E.Properties.VariableNames)
        eventRow.ConditionName = E.ConditionName(1);
    else
        eventRow.ConditionName = "unlabelled";
    end
end

function s = setDefault(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end
