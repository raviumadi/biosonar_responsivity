function result = simulateResponsivityCore(params, opts)
% simulateResponsivityCore
% Unified call-by-call simulator for the responsivity framework.
%
% The simulator generates raw timing, trajectory, anchor, call-duration, and
% wingbeat-state variables. It intentionally does not classify SSGs, buzzes,
% wingbeat synchrony, or feasibility regimes; those labels should be added by
% downstream analysis functions.
%
% Usage:
%   result = simulateResponsivityCore()
%   result = simulateResponsivityCore(params)
%   result = simulateResponsivityCore(params, opts)
%
% Main output:
%   result.calls   - call-level table, one row per emitted call interval
%   result.params  - resolved parameter struct
%   result.opts    - resolved option struct
%   result.summary - compact run summary

    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    params = applyDefaultParams(params);
    opts = applyDefaultOpts(opts);

    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    rows = struct([]);

    for seqID = 1:params.numSequences
        seqRows = simulateOneSequence(seqID, params, opts);
        if isempty(seqRows)
            continue
        end
        if isempty(rows)
            rows = seqRows;
        else
            rows = [rows; seqRows(:)]; %#ok<AGROW>
        end
    end

    if isempty(rows)
        calls = table();
    else
        calls = struct2table(rows);
    end

    result = struct();
    result.calls = calls;
    result.params = params;
    result.opts = opts;
    result.summary = makeSummary(calls);
end

function rows = simulateOneSequence(seqID, params, opts)
    c = params.c;
    cEps = 1e-9;

    batPos = [0 0 0];
    targetPos = initialiseTargets(params, opts);
    targetVel = zeros(opts.numTargets, 3);
    previousAnchorID = NaN;
    elapsedTime = 0;
    callOnsetTime = 0;
    previousTa = NaN;
    previousTcall = NaN;
    contractionHasStarted = false;
    lagState = 0;
    durationJitterState = 0;

    rows = struct([]);

    for callNumber = 1:params.maxCalls
        targetDistances = vecnorm(targetPos - batPos, 2, 2);
        nearestDistance = min(targetDistances);

        if nearestDistance <= params.interceptDistance_m
            break
        end
        if elapsedTime >= params.maxElapsedTime_s
            break
        end

        anchorTargetID = selectAnchor(targetDistances, previousAnchorID, opts);
        anchorPos = targetPos(anchorTargetID, :);
        anchorDistance = targetDistances(anchorTargetID);

        if anchorDistance > params.maxAnchorDistance_m
            break
        end

        if isnan(previousAnchorID)
            anchorSwitch = 0;
        else
            anchorSwitch = double(anchorTargetID ~= previousAnchorID);
        end

        batSpeed = sampleBatSpeed(params, opts);
        losToAnchor = anchorPos - batPos;
        losNorm = norm(losToAnchor);
        if losNorm <= 0
            break
        end
        losUnit = losToAnchor ./ losNorm;
        batVel = makeBatVelocity(batSpeed, losUnit, opts);

        targetVel = sampleTargetVelocities(targetVel, params, opts);
        anchorVel = targetVel(anchorTargetID, :);

        relativeVelocityVec = batVel - anchorVel;
        relativeClosingVelocity = dot(relativeVelocityVec, losUnit);
        relativeClosingVelocityForTiming = min(max(relativeClosingVelocity, 0), c - cEps);

        if strcmpi(opts.timingMode, "staticDelay")
            Tdelay = 2 * anchorDistance / c;
        elseif strcmpi(opts.timingMode, "motionAware")
            Tdelay = 2 * anchorDistance / (c + relativeClosingVelocityForTiming);
        else
            error('Unknown opts.timingMode: %s', string(opts.timingMode));
        end

        [Tcall, TcallDeterministic, TcallJitter, contractionHasStarted, durationJitterState] = ...
            computeCallDuration(callNumber, Tdelay, previousTa, contractionHasStarted, ...
            durationJitterState, params, opts);

        phi = computePhi(callNumber, Tcall, previousTcall, previousTa, params, opts);
        Tphi = phi * Tcall;
        Ta = Tdelay + Tphi;

        TbCore = params.kr * Ta;
        if opts.enforceResponseFloor
            TbCore = max(TbCore, params.responseFloor_s);
        end

        [lagValue, lagState] = computeLag(TbCore, lagState, opts);
        TbEffective = TbCore + lagValue;
        IPI = Ta + TbEffective;

        if opts.enforceMaxCallRate
            IPI = max(IPI, 1 / params.maxCallRate_Hz);
        end

        callRate = 1 / IPI;
        callOffsetTime = callOnsetTime + Tcall;
        echoStartTime = callOnsetTime + Tdelay;
        echoAnchorTime = callOnsetTime + Ta;
        nextCallOnsetTime = callOnsetTime + IPI;
        silentInterval = IPI - Tcall;

        wing = sampleWingbeatState(callOnsetTime, callRate, params, opts);

        batPosNext = batPos + batVel * IPI;
        targetPosNext = targetPos + targetVel * IPI;
        targetDistancesNext = vecnorm(targetPosNext - batPosNext, 2, 2);
        nearestDistanceNext = min(targetDistancesNext);
        anchorDistanceNext = targetDistancesNext(anchorTargetID);

        distanceChangePredicted = relativeClosingVelocityForTiming * IPI;
        distanceChangeObserved = anchorDistance - anchorDistanceNext;

        stopReason = "running";
        if nearestDistanceNext <= params.interceptDistance_m
            stopReason = "intercept_next";
        elseif nextCallOnsetTime >= params.maxElapsedTime_s
            stopReason = "max_time";
        elseif callNumber == params.maxCalls
            stopReason = "max_calls";
        elseif anchorDistanceNext > params.maxAnchorDistance_m
            stopReason = "anchor_too_far";
        end

        row = makeCallRow(seqID, callNumber, callOnsetTime, callOffsetTime, ...
            nextCallOnsetTime, Tcall, TcallDeterministic, TcallJitter, ...
            silentInterval, anchorTargetID, anchorSwitch, anchorDistance, ...
            nearestDistance, targetDistances, batPos, batVel, batSpeed, ...
            targetPos, targetVel, anchorVel, norm(anchorVel), ...
            relativeClosingVelocity, relativeClosingVelocityForTiming, ...
            Tdelay, Tphi, Ta, TbCore, lagValue, TbEffective, IPI, callRate, ...
            echoStartTime, echoAnchorTime, wing, distanceChangePredicted, ...
            distanceChangeObserved, stopReason, params, opts);

        if isempty(rows)
            rows = row;
        else
            rows(end + 1, 1) = row; %#ok<AGROW>
        end

        if stopReason ~= "running"
            break
        end

        previousAnchorID = anchorTargetID;
        previousTa = Ta;
        previousTcall = Tcall;
        elapsedTime = nextCallOnsetTime;
        callOnsetTime = nextCallOnsetTime;
        batPos = batPosNext;
        targetPos = targetPosNext;
    end
end

function params = applyDefaultParams(params)
    params = setDefault(params, 'c', 343);
    params = setDefault(params, 'kr', 5);
    params = setDefault(params, 'initialDistance_m', 5);
    params = setDefault(params, 'initialBatSpeed_m_s', 5);
    params = setDefault(params, 'initialCallDuration_s', 0.005);
    params = setDefault(params, 'minCallDuration_s', 0.0005);
    params = setDefault(params, 'callDurationGain', 1 / params.kr);
    params = setDefault(params, 'callDurationJitterSD_s', 0.0004);
    params = setDefault(params, 'callDurationJitterFrac', 0.05);
    params = setDefault(params, 'callDurationMargin_s', 0);
    params = setDefault(params, 'maxCalls', 250);
    params = setDefault(params, 'maxElapsedTime_s', 5);
    params = setDefault(params, 'maxAnchorDistance_m', 10);
    params = setDefault(params, 'interceptDistance_m', 0.15);
    params = setDefault(params, 'numSequences', 1);
    defaultMaxCallRate_Hz = (params.c + params.initialBatSpeed_m_s) / ...
        (2 * (1 + params.kr) * params.interceptDistance_m);
    params = setDefault(params, 'maxCallRate_Hz', defaultMaxCallRate_Hz);
    params = setDefault(params, 'responseFloor_s', 0.004);
    params = setDefault(params, 'initialWingbeatFreq_Hz', 8);
    params = setDefault(params, 'maxWingbeatFreq_Hz', 12);
    params = setDefault(params, 'dynamicWingbeatFreqMultiplier', 2);
    params = setDefault(params, 'initialWingbeatExcursion', 1);
    params = setDefault(params, 'minWingbeatExcursion', 0.1);
    params = setDefault(params, 'bandwidth_Hz', [45e3 90e3]);
    params = setDefault(params, 'cfFrequency_Hz', 70e3);
end

function opts = applyDefaultOpts(opts)
    opts = setDefault(opts, 'rngSeed', []);
    opts = setDefault(opts, 'geometryMode', "1D");
    opts = setDefault(opts, 'numTargets', 1);
    opts = setDefault(opts, 'initialTargetSeparation_m', 0.5);

    opts = setDefault(opts, 'targetMotion', false);
    opts = setDefault(opts, 'targetVelocityMode', "stationary");
    opts = setDefault(opts, 'targetVelocityScale', 0.5);
    opts = setDefault(opts, 'targetVelocityJitterFrac', 0.15);

    opts = setDefault(opts, 'batVelocityMode', "constant");
    opts = setDefault(opts, 'batVelocityJitterFrac', 0);
    opts = setDefault(opts, 'batSteeringMode', "anchor");

    opts = setDefault(opts, 'anchorMode', "single");
    opts = setDefault(opts, 'timingMode', "motionAware");
    opts = setDefault(opts, 'enforceMaxCallRate', false);
    opts = setDefault(opts, 'enforceResponseFloor', false);
    opts = setDefault(opts, 'phiMode', "constant");
    opts = setDefault(opts, 'phi', 0);
    opts = setDefault(opts, 'phiBounds', [0 1]);
    if ~isfield(opts, 'phiDuration') || ~isstruct(opts.phiDuration)
        opts.phiDuration = struct();
    end
    if ~isfield(opts.phiDuration, 'phiMin')
        opts.phiDuration.phiMin = 0;
    end
    if ~isfield(opts.phiDuration, 'phiMax')
        opts.phiDuration.phiMax = 1;
    end
    if ~isfield(opts.phiDuration, 'alpha')
        opts.phiDuration.alpha = 1;
    end
    if ~isfield(opts.phiDuration, 'usePreviousCall')
        opts.phiDuration.usePreviousCall = true;
    end

    opts = setDefault(opts, 'lag', struct());
    opts.lag = setDefault(opts.lag, 'enabled', false);
    opts.lag = setDefault(opts.lag, 'frac', 0.05);
    opts.lag = setDefault(opts.lag, 'sigma', 0.15);
    opts.lag = setDefault(opts.lag, 'rho', 0.70);
    opts.lag = setDefault(opts.lag, 'maxFrac', 0.25);

    opts = setDefault(opts, 'callDurationMode', "constant");
    opts = setDefault(opts, 'callDurationJitter', struct());
    opts.callDurationJitter = setDefault(opts.callDurationJitter, 'enabled', false);
    opts.callDurationJitter = setDefault(opts.callDurationJitter, 'mode', "additive");
    opts.callDurationJitter = setDefault(opts.callDurationJitter, 'rho', 0);

    opts = setDefault(opts, 'wing', struct());
    opts.wing = setDefault(opts.wing, 'enabled', false);
    opts.wing = setDefault(opts.wing, 'frequencyMode', "fixed");
    opts.wing = setDefault(opts.wing, 'excursionMode', "fixed");
    opts.wing = setDefault(opts.wing, 'phaseMode', "continuous");

    opts = setDefault(opts, 'call', struct());
    opts.call = setDefault(opts.call, 'type', "FM");
end

function s = setDefault(s, fieldName, value)
    if ~isfield(s, fieldName) || isempty(s.(fieldName))
        s.(fieldName) = value;
    end
end

function targetPos = initialiseTargets(params, opts)
    targetPos = nan(opts.numTargets, 3);
    targetPos(1, :) = [params.initialDistance_m 0 0];

    for k = 2:opts.numTargets
        if strcmpi(opts.geometryMode, "3D")
            offsetMag = opts.initialTargetSeparation_m * rand();
            offsetDir = randomUnitVector3D();
            targetPos(k, :) = targetPos(1, :) + offsetMag * offsetDir;
        else
            offset = opts.initialTargetSeparation_m * (k - 1);
            targetPos(k, :) = [params.initialDistance_m + offset 0 0];
        end
    end
end

function anchorTargetID = selectAnchor(targetDistances, previousAnchorID, opts)
    switch lower(char(opts.anchorMode))
        case {'single', 'target1'}
            anchorTargetID = 1;
        case 'target2'
            if numel(targetDistances) < 2
                error('opts.anchorMode is target2, but opts.numTargets < 2.');
            end
            anchorTargetID = 2;
        case 'nearest'
            [~, anchorTargetID] = min(targetDistances);
        case 'random'
            anchorTargetID = randi(numel(targetDistances));
        case 'hysteretic'
            if isnan(previousAnchorID) || rand() < 0.15
                [~, anchorTargetID] = min(targetDistances);
            else
                anchorTargetID = previousAnchorID;
            end
        otherwise
            error('Unknown opts.anchorMode: %s', string(opts.anchorMode));
    end
end

function batSpeed = sampleBatSpeed(params, opts)
    batSpeed = params.initialBatSpeed_m_s;
    if strcmpi(opts.batVelocityMode, "jittered")
        batSpeed = batSpeed * (1 + opts.batVelocityJitterFrac * randn());
    end
    batSpeed = max(batSpeed, 0.01);
end

function batVel = makeBatVelocity(batSpeed, losUnit, opts)
    switch lower(char(opts.batSteeringMode))
        case {'anchor', 'target1'}
            batDir = losUnit;
        case 'fixeddirection'
            batDir = [1 0 0];
        otherwise
            error('Unknown opts.batSteeringMode: %s', string(opts.batSteeringMode));
    end
    batVel = batSpeed * batDir;
end

function targetVel = sampleTargetVelocities(previousTargetVel, params, opts)
    targetVel = previousTargetVel;

    if ~opts.targetMotion || strcmpi(opts.targetVelocityMode, "stationary")
        targetVel(:) = 0;
        return
    end

    nTargets = size(previousTargetVel, 1);
    baseSpeed = opts.targetVelocityScale * params.initialBatSpeed_m_s;

    switch lower(char(opts.targetVelocityMode))
        case 'constant'
            if all(previousTargetVel(:) == 0)
                for k = 1:nTargets
                    targetVel(k, :) = baseSpeed * randomMotionDirection(opts);
                end
            end
        case 'stochastic'
            for k = 1:nTargets
                speed = baseSpeed * (1 + opts.targetVelocityJitterFrac * randn());
                speed = max(speed, 0.01);
                targetVel(k, :) = speed * randomMotionDirection(opts);
            end
        otherwise
            error('Unknown opts.targetVelocityMode: %s', string(opts.targetVelocityMode));
    end
end

function direction = randomMotionDirection(opts)
    if strcmpi(opts.geometryMode, "3D")
        direction = randomUnitVector3D();
    else
        if rand() < 0.5
            direction = [1 0 0];
        else
            direction = [-1 0 0];
        end
    end
end

function [Tcall, TcallDet, TcallJitter, contractionHasStarted, jitterState] = ...
    computeCallDuration(callNumber, Tdelay, previousTa, contractionHasStarted, ...
    jitterState, params, opts)

    mode = lower(char(opts.callDurationMode));
    switch mode
        case 'constant'
            TcallDet = params.initialCallDuration_s;
        case {'previousta', 'decrement'}
            if callNumber == 1 || isnan(previousTa)
                TcallDet = params.initialCallDuration_s;
            else
                TcallDet = params.callDurationGain * previousTa;
            end
        case 'threshold'
            thresholdActive = params.initialCallDuration_s >= ...
                (Tdelay + params.callDurationMargin_s);
            contractionHasStarted = contractionHasStarted || thresholdActive;
            if contractionHasStarted
                TcallDet = params.callDurationGain * Tdelay;
            else
                TcallDet = params.initialCallDuration_s;
            end
        case 'function'
            if ~isfield(opts, 'callDurationFcn') || isempty(opts.callDurationFcn)
                error('opts.callDurationMode is function, but opts.callDurationFcn is missing.');
            end
            state = struct('callNumber', callNumber, 'Tdelay_s', Tdelay, ...
                'previousTa_s', previousTa, 'contractionHasStarted', contractionHasStarted);
            TcallDet = opts.callDurationFcn(state, params, opts);
        otherwise
            error('Unknown opts.callDurationMode: %s', string(opts.callDurationMode));
    end

    TcallDet = boundValue(TcallDet, params.minCallDuration_s, params.initialCallDuration_s);

    TcallJitter = 0;
    if opts.callDurationJitter.enabled
        rho = opts.callDurationJitter.rho;
        if rho > 0
            jitterState = rho * jitterState + sqrt(max(0, 1 - rho^2)) * randn();
            eta1 = jitterState;
        else
            eta1 = randn();
        end

        switch lower(char(opts.callDurationJitter.mode))
            case 'additive'
                TcallJitter = params.callDurationJitterSD_s * eta1;
            case 'proportional'
                TcallJitter = params.callDurationJitterFrac * TcallDet * eta1;
            case 'both'
                TcallJitter = params.callDurationJitterSD_s * eta1 + ...
                    params.callDurationJitterFrac * TcallDet * randn();
            otherwise
                error('Unknown opts.callDurationJitter.mode: %s', ...
                    string(opts.callDurationJitter.mode));
        end
    end

    Tcall = boundValue(TcallDet + TcallJitter, ...
        params.minCallDuration_s, params.initialCallDuration_s);
    TcallJitter = Tcall - TcallDet;
end

function phi = computePhi(callNumber, Tcall, previousTcall, previousTa, params, opts) %#ok<INUSD>
    switch lower(char(opts.phiMode))
        case 'constant'
            phi = opts.phi;
        case 'schedule'
            if ~isfield(opts, 'phiSchedule') || isempty(opts.phiSchedule)
                error('opts.phiMode is schedule, but opts.phiSchedule is missing.');
            end
            idx = min(callNumber, numel(opts.phiSchedule));
            phi = opts.phiSchedule(idx);
        case 'function'
            if ~isfield(opts, 'phiFcn') || isempty(opts.phiFcn)
                error('opts.phiMode is function, but opts.phiFcn is missing.');
            end
            phi = opts.phiFcn(callNumber, params, opts);
        case 'durationcoupled'
            phiMin = opts.phiDuration.phiMin;
            phiMax = opts.phiDuration.phiMax;
            alpha = opts.phiDuration.alpha;
            if opts.phiDuration.usePreviousCall
                refCall = previousTcall;
                if callNumber == 1 || isnan(refCall)
                    refCall = params.initialCallDuration_s;
                end
            else
                refCall = Tcall;
            end
            denom = max(params.initialCallDuration_s - params.minCallDuration_s, eps);
            u = (refCall - params.minCallDuration_s) / denom;
            u = boundValue(u, 0, 1);
            phi = phiMin + (phiMax - phiMin) * (u .^ alpha);
        otherwise
            error('Unknown opts.phiMode: %s', string(opts.phiMode));
    end
    if numel(opts.phiBounds) ~= 2 || opts.phiBounds(2) < opts.phiBounds(1)
        error('opts.phiBounds must be an increasing two-element vector.');
    end
    phi = boundValue(phi, opts.phiBounds(1), opts.phiBounds(2));
end

function [lagValue, lagState] = computeLag(TbCore, lagState, opts)
    if ~opts.lag.enabled
        lagValue = 0;
        return
    end

    meanLag = opts.lag.frac * TbCore;
    lagSD = opts.lag.sigma * meanLag;
    lagState = opts.lag.rho * lagState + lagSD * randn();
    lagCandidate = max(0, meanLag + lagState);
    lagCap = opts.lag.maxFrac * TbCore;
    lagValue = min(lagCandidate, lagCap);
end

function wing = sampleWingbeatState(callOnsetTime, callRate, params, opts)
    wing = struct();

    if ~opts.wing.enabled
        wing.WingbeatFreq_Hz = NaN;
        wing.WingbeatPeriod_s = NaN;
        wing.WingbeatPhase_rad = NaN;
        wing.WingbeatPhase_frac = NaN;
        wing.WingbeatExcursion = NaN;
        wing.WingbeatAngularPosition = NaN;
        wing.WingbeatCycleNumber = NaN;
        wing.WingCondition = "none";
        return
    end

    switch lower(char(opts.wing.frequencyMode))
        case 'fixed'
            wingFreq = params.initialWingbeatFreq_Hz;
        case 'dynamic'
            if strcmpi(opts.wing.excursionMode, "dynamic")
                wingFreqCeiling = params.dynamicWingbeatFreqMultiplier * params.maxWingbeatFreq_Hz;
            else
                wingFreqCeiling = params.maxWingbeatFreq_Hz;
            end
            wingFreq = min(callRate, wingFreqCeiling);
        otherwise
            error('Unknown opts.wing.frequencyMode: %s', string(opts.wing.frequencyMode));
    end

    switch lower(char(opts.wing.excursionMode))
        case 'fixed'
            wingExcursion = params.initialWingbeatExcursion;
        case 'dynamic'
            wingExcursion = params.initialWingbeatExcursion * ...
                (params.maxWingbeatFreq_Hz / max(wingFreq, eps));
            wingExcursion = min(wingExcursion, params.initialWingbeatExcursion);
            wingExcursion = max(wingExcursion, params.minWingbeatExcursion);
        otherwise
            error('Unknown opts.wing.excursionMode: %s', string(opts.wing.excursionMode));
    end

    wingPeriod = 1 / wingFreq;
    wingPhaseFrac = mod(callOnsetTime, wingPeriod) / wingPeriod;
    wingPhaseRad = 2 * pi * wingPhaseFrac;
    wingAngularPosition = wingExcursion * sin(wingPhaseRad);
    wingCycleNumber = floor(callOnsetTime / wingPeriod);

    wing.WingbeatFreq_Hz = wingFreq;
    wing.WingbeatPeriod_s = wingPeriod;
    wing.WingbeatPhase_rad = wingPhaseRad;
    wing.WingbeatPhase_frac = wingPhaseFrac;
    wing.WingbeatExcursion = wingExcursion;
    wing.WingbeatAngularPosition = wingAngularPosition;
    wing.WingbeatCycleNumber = wingCycleNumber;
    wing.WingCondition = string(opts.wing.frequencyMode) + "_" + string(opts.wing.excursionMode);
end

function row = makeCallRow(seqID, callNumber, callOnsetTime, callOffsetTime, ...
    nextCallOnsetTime, Tcall, TcallDet, TcallJitter, silentInterval, ...
    anchorTargetID, anchorSwitch, anchorDistance, nearestDistance, ...
    targetDistances, batPos, batVel, batSpeed, targetPos, targetVel, ...
    anchorVel, anchorSpeed, relativeClosingVelocity, relativeClosingVelocityForTiming, ...
    Tdelay, Tphi, Ta, TbCore, lagValue, TbEffective, IPI, callRate, ...
    echoStartTime, echoAnchorTime, wing, distanceChangePredicted, ...
    distanceChangeObserved, stopReason, params, opts)

    row = struct();
    row.SeqID = seqID;
    row.CallNumber = callNumber;
    row.CallOnsetTime_s = callOnsetTime;
    row.CallOffsetTime_s = callOffsetTime;
    row.NextCallOnsetTime_s = nextCallOnsetTime;
    row.Tcall_s = Tcall;
    row.TcallDeterministic_s = TcallDet;
    row.TcallJitter_s = TcallJitter;
    row.CallDurationMode = string(opts.callDurationMode);
    row.CallDurationJitterMode = string(opts.callDurationJitter.mode);
    row.SilentInterval_s = silentInterval;

    row.AnchorTargetID = anchorTargetID;
    row.AnchorMode = string(opts.anchorMode);
    row.AnchorSwitch = anchorSwitch;
    row.AnchorDistance_m = anchorDistance;
    row.NearestTargetDistance_m = nearestDistance;

    row.BatX_m = batPos(1);
    row.BatY_m = batPos(2);
    row.BatZ_m = batPos(3);
    row.BatVx_m_s = batVel(1);
    row.BatVy_m_s = batVel(2);
    row.BatVz_m_s = batVel(3);
    row.BatSpeed_m_s = batSpeed;

    row.AnchorTargetVx_m_s = anchorVel(1);
    row.AnchorTargetVy_m_s = anchorVel(2);
    row.AnchorTargetVz_m_s = anchorVel(3);
    row.AnchorTargetSpeed_m_s = anchorSpeed;
    row.RelativeClosingVelocity_m_s = relativeClosingVelocity;
    row.RelativeClosingVelocityForTiming_m_s = relativeClosingVelocityForTiming;

    row.Tdelay_s = Tdelay;
    row.Tphi_s = Tphi;
    row.Ta_s = Ta;
    row.TbCore_s = TbCore;
    row.Lag_s = lagValue;
    row.TbEffective_s = TbEffective;
    row.IPI_s = IPI;
    row.CallRate_Hz = callRate;
    row.EchoStartTime_s = echoStartTime;
    row.EchoAnchorTime_s = echoAnchorTime;

    row.WingbeatFreq_Hz = wing.WingbeatFreq_Hz;
    row.WingbeatPeriod_s = wing.WingbeatPeriod_s;
    row.WingbeatPhase_rad = wing.WingbeatPhase_rad;
    row.WingbeatPhase_frac = wing.WingbeatPhase_frac;
    row.WingbeatExcursion = wing.WingbeatExcursion;
    row.WingbeatAngularPosition = wing.WingbeatAngularPosition;
    row.WingbeatCycleNumber = wing.WingbeatCycleNumber;
    row.WingCondition = wing.WingCondition;

    row.DistanceChangePredicted_m = distanceChangePredicted;
    row.DistanceChangeObserved_m = distanceChangeObserved;
    row.StopReason = string(stopReason);

    row.CallType = string(opts.call.type);
    row.StartFreq_Hz = min(params.bandwidth_Hz);
    row.EndFreq_Hz = max(params.bandwidth_Hz);
    row.CFreq_Hz = params.cfFrequency_Hz;
    row.BandwidthLow_Hz = min(params.bandwidth_Hz);
    row.BandwidthHigh_Hz = max(params.bandwidth_Hz);

    for k = 1:opts.numTargets
        row.(sprintf('Target%dDistance_m', k)) = targetDistances(k);
        row.(sprintf('Target%dX_m', k)) = targetPos(k, 1);
        row.(sprintf('Target%dY_m', k)) = targetPos(k, 2);
        row.(sprintf('Target%dZ_m', k)) = targetPos(k, 3);
        row.(sprintf('Target%dVx_m_s', k)) = targetVel(k, 1);
        row.(sprintf('Target%dVy_m_s', k)) = targetVel(k, 2);
        row.(sprintf('Target%dVz_m_s', k)) = targetVel(k, 3);
    end
end

function summary = makeSummary(calls)
    summary = struct();

    if isempty(calls)
        summary.NumSequences = 0;
        summary.NumCalls = 0;
        summary.MeanIPI_s = NaN;
        summary.MedianIPI_s = NaN;
        summary.MeanCallRate_Hz = NaN;
        summary.NumAnchorSwitches = 0;
        summary.PercentAnchorSwitches = NaN;
        summary.StopReasonCounts = table();
        return
    end

    summary.NumSequences = numel(unique(calls.SeqID));
    summary.NumCalls = height(calls);
    summary.MeanIPI_s = mean(calls.IPI_s, 'omitnan');
    summary.MedianIPI_s = median(calls.IPI_s, 'omitnan');
    summary.MeanCallRate_Hz = mean(calls.CallRate_Hz, 'omitnan');
    summary.NumAnchorSwitches = sum(calls.AnchorSwitch, 'omitnan');
    possibleSwitches = max(0, height(calls) - summary.NumSequences);
    if possibleSwitches > 0
        summary.PercentAnchorSwitches = 100 * summary.NumAnchorSwitches / possibleSwitches;
    else
        summary.PercentAnchorSwitches = NaN;
    end

    seqIDs = unique(calls.SeqID);
    terminalStopReason = strings(numel(seqIDs), 1);
    for ii = 1:numel(seqIDs)
        idx = find(calls.SeqID == seqIDs(ii));
        terminalStopReason(ii) = calls.StopReason(idx(end));
    end

    [stopGroups, stopNames] = findgroups(terminalStopReason);
    stopCounts = splitapply(@numel, terminalStopReason, stopGroups);
    summary.StopReasonCounts = table(stopNames, stopCounts, ...
        'VariableNames', {'StopReason', 'NumSequences'});
end

function value = boundValue(value, lowerBound, upperBound)
    value = min(max(value, lowerBound), upperBound);
end

function u = randomUnitVector3D()
    v = randn(1, 3);
    nv = norm(v);
    if nv == 0
        u = [1 0 0];
    else
        u = v ./ nv;
    end
end
