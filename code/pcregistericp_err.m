function [tform, movingReg, rmse, err] = pcregistericp_err(moving, fixed, varargin)
%PCREGISTERICP Register two point clouds using ICP algorithm.
%   tform = PCREGISTERICP(moving, fixed) returns the rigid transformation
%   that registers the moving point cloud with the fixed point cloud.
%   moving and fixed are pointCloud objects. tform is a rigid3d object
%   that describes the rigid 3-D transform. The rigid transformation
%   between moving and fixed are estimated using the Iterative Closest
%   Point (ICP) algorithm. Consider downsampling point clouds using
%   pcdownsample before using PCREGISTERICP to improve accuracy and
%   efficiency of registration.
%
%   [tform, movingReg] = PCREGISTERICP(moving, fixed) additionally
%   returns the transformed  point cloud, movingReg, that is aligned with
%   the fixed point cloud.
%
%   [..., rmse] = PCREGISTERICP(moving, fixed) additionally returns the
%   root mean squared error of the Euclidean distance between the inlier
%   points of the  aligned point clouds.
%
%   [...] = PCREGISTERICP(...,Name, Value) specifies additional
%   name-value pairs described below:
%
%   'Metric'            A string used to specify the metric of the
%                       minimization function. The ICP algorithm minimizes
%                       the distance between the two point clouds according
%                       to the given metric. Valid strings are
%                       'pointToPoint' and 'pointToPlane'. Setting Metric
%                       to 'pointToPlane' can reduce the number of
%                       iterations to process. However, this metric
%                       requires extra algorithmic steps within each
%                       iteration. The 'pointToPlane' metric helps
%                       registration of planar surfaces.
%
%                       Default: 'pointToPoint'
%
%   'Extrapolate'       A boolean to turn on/off the extrapolation step
%                       that traces out a path in the registration state
%                       space, described in the original paper of ICP by
%                       Besl and McKay (1992). This may reduce the number
%                       of iterations to converge.
%
%                       Default: false
%
%   'InlierRatio'       A scalar to specify the percentage of inliers.
%                       During an ICP iteration, every point in the moving
%                       point cloud is matched to its nearest neighbor in
%                       the fixed point cloud. The pair of matched points
%                       is considered as an inlier if its Euclidean
%                       distance falls into the given percentage of the
%                       distribution of matching distance. By default, all
%                       matching pairs are used.
%
%                       Default: 1
%
%   'InitialTransform'  A rigid3d object to specify the initial rigid
%                       transformation. This is useful when a coarse
%                       estimation can be provided externally.
%
%                       Default: rigid3d()
%
%   'MaxIterations'     A positive integer to specify the maximum number
%                       of iterations before ICP stops.
%
%                       Default: 20
%
%   'Tolerance'         A 2-element vector, [Tdiff, Rdiff], to specify
%                       the tolerance of absolute difference in translation
%                       and rotation estimated in consecutive ICP
%                       iterations. Tdiff measures the Euclidean distance
%                       between two translation vectors, while Rdiff
%                       measures the angular difference in degrees. The
%                       algorithm stops when the average difference between
%                       estimated rigid transformations in the three most
%                       recent consecutive iterations falls below the
%                       specified tolerance value.
%
%                       Default: [0.01, 0.5]
%
%   'Verbose'           Set true to display progress information.
%
%                       Default: false
%
%   Notes
%   -----
%   - The registration algorithm is based on Iterative Closest Point (ICP)
%     algorithm, which is an iterative process. Best performance might
%     require adjusting the options for different data.
%
%   - Point cloud normals are required by the registration algorithm when
%    'pointToPlane' metric is chosen. If the Normal property of the second
%    input is empty, the function fills it.
%
%   - When the Normal property is filled automatically, the number of
%     points, K, to fit local plane is set to 6. This default may not work
%     under all circumstances. If the registration with 'pointToPlane'
%     metric fails, consider calling pcnormals function with a custom value
%     of K.
%
%   Class Support
%   -------------
%   moving and fixed must be pointCloud objects.
%
%   Example: Align two point clouds
%   -------------------------------
%   ptCloud = pcread('teapot.ply');
%   figure
%   pcshow(ptCloud)
%   title('Teapot')
%
%   % Create a transform object with 30 degree rotation along z-axis and
%   % translation [5, 5, 10]
%   theta = pi/6;
%   rot   = [cos(theta) sin(theta) 0; ...
%           -sin(theta) cos(theta) 0; ...
%                    0          0  1];
%   trans = [5 5 10];
%
%   tform1 = rigid3d(rot, trans);
%
%   % Transform the point cloud
%   ptCloudTformed = pctransform(ptCloud, tform1);
%
%   figure
%   pcshow(ptCloudTformed)
%   title('Transformed Teapot')
%
%   % Apply the rigid registration
%   [tform, ptCloudReg] = PCREGISTERICP(ptCloudTformed, ptCloud, 'Extrapolate', true);
%
%   % Visualize the alignment
%   pcshowpair(ptCloud, ptCloudReg)
%
%   % Compare the result with the true transformation
%   disp(tform1.T);
%   tform2 = invert(tform);
%   disp(tform2.T);
%
% See also pointCloud, pcregisterndt, pcregistercpd, pctransform, rigid3d,
%          pcdownsample, pcshowpair, pcdenoise, pcmerge, pcshow

% Copyright 2014-2021 The MathWorks, Inc.
%
% References
% ----------
% Besl, Paul J.; N.D. McKay (1992). "A Method for Registration of 3-D
% Shapes". IEEE Trans. on Pattern Analysis and Machine Intelligence (Los
% Alamitos, CA, USA: IEEE Computer Society) 14 (2): 239-256.
%
% Chen, Yang; Gerard Medioni (1991). "Object modelling by registration of
% multiple range images". Image Vision Comput. (Newton, MA, USA:
% Butterworth-Heinemann): 145-155

%#codegen

% Validate inputs
[metric, doExtrapolate, inlierRatio, maxIterations, tolerance, ...
    initialTransform, verbose, useDegree] = parseInputs(moving, fixed, varargin{:});

if isSimMode
    printer = vision.internal.MessagePrinter.configure(verbose);
end

isPointToPoint = strcmpi(metric, 'PointToPoint');
useAllMatches  = inlierRatio == 1;

[ptCloudA, ptCloudB] = preparePointClouds(moving, fixed, isPointToPoint);

Rs = zeros(3, 3, maxIterations+1);
Ts = zeros(3, maxIterations+1);

% Quaternion and translation vector
qs = [ones(1, maxIterations+1); zeros(6, maxIterations+1)];

% The difference of quaternion and translation vector in consecutive
% iterations
dq = zeros(7, maxIterations+1);

% The angle between quaternion and translation vectors in consecutive
% iterations
dTheta = zeros(maxIterations+1, 1);

% Inlier RMSE
err = zeros(maxIterations+1, 1);

% Apply the initial condition.
% We use pre-multiplication format in this algorithm.
Rs(:,:,1) = initialTransform.Rotation';
Ts(:,1)   = initialTransform.Translation';
qs(:,1)   = [vision.internal.quaternion.rotationToQuaternion(Rs(:,:,1)); Ts(:,1)];

locA = ptCloudA.Location;
if qs(1) ~= 0 || any(qs(2:end,1))
    locA = rigidTransform(ptCloudA.Location, Rs(:,:,1), Ts(:,1));
end

stopIteration   = maxIterations;
numPointsA      = ptCloudA.Count;
upperBound      = max(1, round(inlierRatio * numPointsA));

% Start ICP iterations
for i = 1 : maxIterations
        if verbose
            printer.linebreak;
            printer.print('--------------------------------------------\n');
            printer.printMessage('vision:pointcloud:icpIteration',i);
            printer.printMessageNoReturn('vision:pointcloud:findCorrespondenceStart');
        end

    if useAllMatches
        % Find correspondences
        inlierIndicesA = (1 : upperBound).';
        [inlierIndicesB, inlierDist] = multiQueryKNNSearchImpl(ptCloudB, locA, 1);
    else
        % Find correspondences
        [indices, dists] = multiQueryKNNSearchImpl(ptCloudB, locA, 1);

        % Remove outliers
        keepInlierA = false(numPointsA, 1);
        [~, idx] = mink(dists, upperBound, 2);
        keepInlierA(idx) = true;
        inlierIndicesA = find(keepInlierA);
        if isSimMode()
            inlierIndicesB = indices(keepInlierA);
            inlierDist     = dists(keepInlierA);
        else
            inlierIndicesB = indices(keepInlierA');
            inlierDist = dists(keepInlierA');
        end
    end

    coder.internal.errorIf(numel(inlierIndicesA) < 3, 'vision:pointcloud:notEnoughPoints');

    if i == 1
        err(i) = sqrt(sum(inlierDist(:))/numel(inlierDist));
    end

    if verbose
        printer.printMessage('vision:pointcloud:stepCompleted');
        printer.printMessageNoReturn('vision:pointcloud:estimateTransformStart');
    end

    % Estimate transformation given correspondences
    if isPointToPoint
        [R, T] = vision.internal.calibration.computeRigidTransform(...
            locA(inlierIndicesA, :), ...
            ptCloudB.Location(inlierIndicesB, :));

    else % PointToPlane
        [R, T] = pointToPlaneMetric(locA(inlierIndicesA, :), ...
            ptCloudB.Location(inlierIndicesB, :), ptCloudB.Normal(inlierIndicesB, :));
    end

    % Bad correspondence may lead to singular matrix
    coder.internal.errorIf(any(isnan(T)) || any(isnan(R), 'all') , 'vision:pointcloud:singularMatrix');

    % Update the total transformation
    Rnew = R * Rs(:,:,i);
    Tnew = R * Ts(:,i) + T;
    Rs(:,:,i+1) = Rnew;
    Ts(:,i+1)   = Tnew;
    
    if verbose
        printer.printMessage('vision:pointcloud:stepCompleted');
    end

    % RMSE
    locA = rigidTransform(ptCloudA.Location, Rnew, Tnew);
    squaredError = sum((locA(inlierIndicesA, :) - ptCloudB.Location(inlierIndicesB, :)).^2, 2);
    err(i+1) = sqrt(sum(squaredError)/numel(squaredError));

    % Convert to vector representation
    qs(:,i+1) = [vision.internal.quaternion.rotationToQuaternion(Rnew); Tnew];

    % With extrapolation, we might be able to converge faster
    if doExtrapolate
        if verbose
            printer.printMessageNoReturn('vision:pointcloud:updateTransformStart');
        end
        extrapolateInTransformSpace;
        if verbose
            printer.printMessage('vision:pointcloud:stepCompleted');
        end
    end

    % Check convergence
    % Compute the mean difference in R/T from the recent three iterations.
    [dR, dT, rdiff, tdiff] = getChangesInTransformation;

        if verbose
            if useDegree
                rdiff = rdiff*180/pi;
            end
            printer.printMessage('vision:pointcloud:checkConverge', ...
                num2str(tdiff), num2str(rdiff), num2str(err(i+1)));
        end

    % Stop ICP if it already converges
    if dT <= tolerance(1) && dR <= tolerance(2)
        stopIteration = i;
        break;
    end
end

% Make the R to be orthogonal as much as possible
R = Rs(:,:,stopIteration+1)';
[U, ~, V] = svd(R);
R = U * V';

tform = rigid3d(R, Ts(:, stopIteration+1)');

rmse = err(stopIteration+1);


    if verbose
        printer.linebreak;
        printer.print('--------------------------------------------\n');
        printer.printMessage('vision:pointcloud:icpSummary',stopIteration, num2str(rmse));
    end

if nargout >= 2
    movingReg = pctransform(moving, tform);
end

%----------------------------------------------------------------------
% Nested function to perform extrapolation
% Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes.
% IEEE Transactions on pattern analysis and machine intelligence, p245.
%----------------------------------------------------------------------
    function extrapolateInTransformSpace
        dq(:,i+1) = qs(:,i+1) - qs(:,i);
        n1 = norm(dq(:,i));
        n2 = norm(dq(:,i+1));
        dTheta(i+1) = (180/pi)*acos(dot(dq(:,i),dq(:,i+1))/(n1*n2));

        angleThreshold = 10;
        scaleFactor    = 25;
        if i > 2 && dTheta(i+1) < angleThreshold && dTheta(i) < angleThreshold
            d = [err(i+1), err(i), err(i-1)];
            v = [0, -n2, -n1-n2];
            vmax = scaleFactor * n2;
            dv = extrapolate(v,d,vmax);
            if dv ~= 0
                q = qs(:,i+1) + dv * dq(:,i+1)/n2;
                q(1:4) = q(1:4)/norm(q(1:4));
                % Update transformation and data
                qs(:,i+1) = q;
                Rs(:,:,i+1) = vision.internal.quaternion.quaternionToRotation(q(1:4));
                Ts(:,i+1) = q(5:7);
                locA = rigidTransform(ptCloudA.Location, Rs(:,:,i+1), Ts(:,i+1));
            end
        end
    end

%----------------------------------------------------------------------
% Nested function to compute the changes in rotation and translation
%----------------------------------------------------------------------
    function [dR, dT, rdiff, tdiff] = getChangesInTransformation
        dR = 0;
        dT = 0;
        rdiff = 0;
        tdiff = 0;
        count = 0;
        for k = max(i-2,1):i
            % Rotation difference in radians
            rdiff = acos(min(1, max(-1, dot(qs(1:4,k),qs(1:4,k+1))/(norm(qs(1:4,k))*norm(qs(1:4,k+1))))));
            % Euclidean difference
            tdiff = sqrt(sum((Ts(:,k)-Ts(:,k+1)).^2));
            dR = dR + rdiff;
            dT = dT + tdiff;
            count = count + 1;
        end
        dT = dT/count;
        dR = dR/count;
    end
end

%--------------------------------------------------------------------------
function [ptCloudA, ptCloudB] = preparePointClouds(moving, fixed, isPointToPoint)

% Unorganized M-by-3 data
ptCloudA                        = removeInvalidPoints(moving);
[ptCloudB, validPtCloudIndices] = removeInvalidPoints(fixed);

% At least three points are needed to determine a 3-D transformation
coder.internal.errorIf( ptCloudA.Count < 3 || ptCloudB.Count < 3, 'vision:pointcloud:notEnoughPoints');

% Normal vector is needed for PointToPlane metric
if ~isPointToPoint
    % Compute the unit normal vector if it is not provided.
    if isempty(fixed.Normal)
        fixedCount = fixed.Count;
        % Use 6 neighboring points to estimate a normal vector. You may use
        % pcnormals with customized parameter to compute normals upfront.
        fixed.Normal = surfaceNormalImpl(fixed, 6);
        ptCloudB.Normal = [fixed.Normal(validPtCloudIndices), ...
            fixed.Normal(validPtCloudIndices + fixedCount), ...
            fixed.Normal(validPtCloudIndices + fixedCount * 2)];
    end

    % Remove points if their normals are invalid
    validIndices = all(isfinite(ptCloudB.Normal), 2);
    if nnz(validIndices) < ptCloudB.Count
        [loc, ~, nv] = subsetImpl(ptCloudB, validIndices);
        ptCloudB = pointCloud(loc, 'Normal', nv);
        coder.internal.errorIf(ptCloudB.Count < 3, 'vision:pointcloud:notEnoughPoints');
    end
end

end

%--------------------------------------------------------------------------
% Parameter validation
%--------------------------------------------------------------------------
function [metric, doExtrapolate, inlierRatio, maxIterations, tolerance, ...
    initTformParsed, verbose, useDegree] = parseInputs(moving, fixed, varargin)

funcName = mfilename;
validateattributes(moving, {'pointCloud'}, {'scalar'}, funcName, 'moving');
validateattributes(fixed,  {'pointCloud'}, {'scalar'}, funcName, 'fixed');

if isSimMode()
    [metric, doExtrapolate, inlierRatio, maxIterations, tolerance, ...
        initTform, verbose, useDegree] = vision.internal.pc.parseICPOptionsSim(varargin{:});
else
    [metric, doExtrapolate, inlierRatio, maxIterations, tolerance, ...
        initTform, verbose, useDegree] = parseInputsCG(varargin{:});
end



% Convert from degree to radian internally
if useDegree
    tolerance(2) = tolerance(2)*pi/180;
end

if isa(initTform, 'affine3d')
    initTformParsed = rigid3d(initTform.T);
else
    initTformParsed = initTform;
end
end

%--------------------------------------------------------------------------
function [metric, doExtrapolate, inlierRatio, maxIterations, tolerance, ...
    initTform, verbose, useDegree] = parseInputsCG(varargin)
defaults = struct(...
    'Metric', 'PointToPoint', ...
    'Extrapolate',  false, ...
    'InlierRatio', 1.0,...
    'MaxIterations', 20,...
    'Tolerance', [0.01, 0.5],...
    'InitialTransform', rigid3d(),...
    'Verbose', false, ...
    'UseDegree', true);

pvPairs = struct(...
    'Metric', uint32(0),...
    'Extrapolate', uint32(0),...
    'InlierRatio', uint32(0),...
    'MaxIterations', uint32(0),...
    'Tolerance', uint32(0),...
    'InitialTransform', uint32(0),...
    'Verbose', uint32(0),...
    'UseDegree', uint32(0));

properties =  struct( ...
    'CaseSensitivity', false, ...
    'StructExpand',    true, ...
    'PartialMatching', false);

optarg = coder.internal.parseParameterInputs(pvPairs, properties, varargin{:});

metric = coder.internal.getParameterValue(optarg.Metric, defaults.Metric, varargin{:});
validateattributes(metric, {'char','string'},{'scalartext'});
validatestring(metric, {'PointToPoint', 'PointToPlane'}, mfilename, 'Metric');

doExtrapolate   = coder.internal.getParameterValue(optarg.Extrapolate, defaults.Extrapolate, varargin{:});
validateattributes(doExtrapolate, {'logical'}, {'scalar','nonempty'});

inlierRatio     = coder.internal.getParameterValue(optarg.InlierRatio, defaults.InlierRatio, varargin{:});
validateattributes(inlierRatio, {'single', 'double'}, {'real','nonempty','scalar','>',0,'<=',1});

maxIterations   = coder.internal.getParameterValue(optarg.MaxIterations, defaults.MaxIterations, varargin{:});
validateattributes(maxIterations, {'single', 'double'}, {'real','scalar','integer','positive'})

tolerance       = coder.internal.getParameterValue(optarg.Tolerance, defaults.Tolerance, varargin{:});
validateattributes(tolerance, {'single', 'double'}, {'real','nonnegative','numel', 2});

useDegree       = coder.internal.getParameterValue(optarg.UseDegree, defaults.UseDegree, varargin{:});
validateattributes(useDegree, {'logical'}, {'scalar','nonempty'});

% verbose is not applicable for generated code, but however validated to
% maintain consistency
verboseIn = coder.internal.getParameterValue(optarg.Verbose, defaults.Verbose, varargin{:});
validateattributes(verboseIn, {'logical'}, {'scalar','nonempty'});

verbose = false;

initTform = coder.internal.getParameterValue(optarg.InitialTransform, defaults.InitialTransform, varargin{:});
validateattributes(initTform, {'rigid3d', 'affine3d'}, {'scalar'})
coder.internal.errorIf(~(vision.internal.isRigidTransform(initTform.T)),'vision:pointcloud:rigidTransformOnly');

end

function B = rigidTransform(A, R, T)
B = A * R' + T';
end

%--------------------------------------------------------------------------
% Solve the following minimization problem:
%       min_{R, T} sum(|dot(R*p+T-q,nv)|^2)
%
% p, q, nv are all N-by-3 matrix, and nv is the unit normal at q
%
% Here the problem is solved by linear approximation to the rotation matrix
% when the angle is small.
%--------------------------------------------------------------------------
function [R, T] = pointToPlaneMetric(p, q, nv)

% Set up the linear system
cn = [cross(p,nv,2),nv];
C = cn'*cn;
qp = q-p;
b =  [...
    sum(qp .* cn(:,1) .* nv, 'all');
    sum(qp .* cn(:,2) .* nv, 'all');
    sum(qp .* cn(:,3) .* nv, 'all');
    sum(qp .* cn(:,4) .* nv, 'all');
    sum(qp .* cn(:,5) .* nv, 'all');
    sum(qp .* cn(:,6) .* nv, 'all')];

% X is [alpha, beta, gamma, Tx, Ty, Tz]
X = C\b;

cx = cos(X(1));
cy = cos(X(2));
cz = cos(X(3));
sx = sin(X(1));
sy = sin(X(2));
sz = sin(X(3));

R = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
    cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
    -sy,          sx*cy,          cx*cy];

T = X(4:6);
end

%--------------------------------------------------------------------------
% Extrapolation in quaternion space. Details are found in:
% Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes.
% IEEE Transactions on pattern analysis and machine intelligence, 239-256.
%--------------------------------------------------------------------------
function dv = extrapolate(v,d,vmax)
p1 = polyfit(v,d,1);    % linear fit
p2 = polyfit(v,d,2);    % parabolic fit
v1 = -p1(2)/p1(1);      % linear zero crossing point
v2 = -p2(2)/(2*p2(1));  % polynomial top point

if (issorted([0 v2 v1 vmax]) || issorted([0 v2 vmax v1]))
    % Parabolic update
    dv = v2;
elseif (issorted([0 v1 v2 vmax]) || issorted([0 v1 vmax v2])...
        || (v2 < 0 && issorted([0 v1 vmax])))
    % Line update
    dv = v1;
elseif (v1 > vmax && v2 > vmax)
    % Maximum update
    dv = vmax;
else
    % No extrapolation
    dv = 0;
end
end

%--------------------------------------------------------------------------
% Codegen support flag
%--------------------------------------------------------------------------
function flag = isSimMode()
flag = isempty(coder.target);
end