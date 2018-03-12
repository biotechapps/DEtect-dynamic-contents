function [sys,X0] = PHYCOB(data,Order,varargin)
%N4SID  Estimate state-space model using a subspace method.
%     
%   SYS = N4SID(DATA, NX)
%     estimates a state-space model of order NX using time- or
%     frequency-domain data DATA. SYS is an IDSS object which 
%     encapsulates an identified state-space model:
%             dx(t) = A x(t) + B u(t) + K e(t)
%              y(t) = C x(t) + D u(t) + e(t)
%     where A, B, C, D and K are state-space matrices, u(t) is the input,
%     y(t) the output, e(t) the disturbance and x(t) the vector of NX
%     states. See help on IDSS for more information. All the entries of A,
%     B, C and K are considered free parameters. D is fixed to zero values
%     (no feedthrough) by default except for static systems. Use the
%     "Form", "Feedthrough" and "DisturbanceModel" name/value pairs to
%     modify the default behavior.
%
%     DATA: Time domain estimation data must be specified as an IDDATA
%     object. Frequency domain data can mean either the measured frequency
%     response (FRF) or the frequency domain input-output signals which are
%     related to the corresponding time domain signals via the Fourier
%     Transform. Use an IDFRD or FRD object to specify the frequency
%     response data and an IDDATA object with Domain = 'Frequency' for the
%     frequency domain signals. By default, the sample time of SYS matches
%     that of the estimation data. For estimating a continuous-time model,
%     specify 'Ts'/0 as name-value pair, or use continuous-time frequency
%     domain data.
%
%     NX: Specify order NX as a positive integer. You can also specify a
%     vector as value of NX, as in NX = 1:15, which creates a plot from
%     where you can choose a suitable model order. The plot shows the Hankel
%     singular values for models of different orders. States with
%     relatively small Hankel singular values can be safely discarded. A
%     default choice is suggested in the plot. 
%     You can also specify NX = 'best', as in N4SID(DATA,'best'), in which
%     case the optimal order is chosen automatically in the 1:10 range.
%
%  SYS = N4SID(DATA, NX, 'Name1', Value1, 'Name2', Value2,...)
%     specifies additional attributes of the model structure as name-value
%     pairs. The following name-value pairs are available:
%     'Ts':   Specify the sample time of the model. The value must be 0
%             (default) or equal to the sample time of data (DATA.Ts). Use
%             'Ts' to create continuous-time models, as in:
%             SYS = N4SID(DATA,2,'Ts',0) 
%     'Form': Specify the structure of the A, B and C matrices
%             in some common forms:
%                  'free': All entries of A, B, C are treated as free and
%                          updated during estimation. Default. 
%             'companion': Companion form of the model where the characteristic
%                          polynomial appears in the rightmost column.
%                 'modal': Modal decomposition form where the state matrix
%                          A is block diagonal, each block corresponding to
%                          a cluster of nearby modes.
%             'canonical': A, B, and C are parameterized in an observability
%                          canonical form, as described in [1]. 
%    'Feedthrough': Denote presence of direct feedthrough from the input
%             u(t) to the output y(t). Specify as a logical vector of
%             length equal to the number of inputs (Nu). By default
%             feedthrough is assumed to be absent (false(1,Nu)) unless the
%             model has no states in which case it is considered to be
%             present for all inputs (true(1,Nu)).
%    'DisturbanceModel': Specify if the noise component of the model should
%             be estimated or not. Set the value to one of:
%                 'none': Noise component is not estimated. The value of
%                         the K matrix is fixed to zero value.
%             'estimate': The K matrix is treated as a free parameter.
%    'InputDelay': Specify input delay as a double vector of length equal
%             to number of inputs. In discrete-time case, delay must be
%             specified in number of lags (multiples of DATA.Ts).
%
%  SYS = N4SID(DATA, NX, ... OPTIONS)
%     specifies estimation options that configure the initial states,
%     estimation objective, and subspace algorithm related choices to be
%     used for estimation. Use the "n4sidOptions" command to configure
%     OPTIONS. 
%
%  [SYS, X0] = N4SID(DATA, ...) 
%    also returns the value of the initial states computed during the
%    estimation of SYS. This value is also stored as a part of the model's
%    estimation report. If DATA contains multiple experiments, X0 is a
%    matrix with one column of values for each data experiment.
%
% Reference:
%  [1] Lennart Ljung: System Identification - Theory For the User, Appendix
%  4A, pp 132-134, 2nd ed, PTR Prentice Hall, Upper Saddle River, N.J., 1999.
%
% See also N4SIDOPTIONS, IDSS, SSEST, TFEST, PROCEST, POLYEST, IDDATA,
% IDFRD, DYNAMICSYSTEM/CANON, IDGREY, PEM, DYNAMICSYSTEM/HSVD.

%   M. Viberg, 8-13-1992, T. McKelvey, L. Ljung 9-26-1993.
%   Rewritten; L. Ljung 8-3-2000, updated R. Singh 04/2010.

%  Copyright 1986-2016 The MathWorks, Inc.

% Note: n4sid(data, model) is undocumented; it provides meta data but no
% parameter specs.

narginchk(1,Inf)

%% Parse inputs and convert orders, if specified, into an idss template.
if nargin<2,Order = 'best'; end
% Set estimation data name.
if (isa(data,'iddata') || isa(data,'frd')) && isempty(data.Name)
   data.Name = inputname(1);
end

% Validate input arguments.
try
   [sys, EstimData, Order] = localProcessInputs(data, Order, varargin{:});
catch E
   throw(E)
end

Options = getDefaultOptions(sys);
Disp = ~strcmpi(Options.Display,'off');
if Disp
   W = Options.ProgressWindow;
   Str = ctrlMsgUtils.message('Ident:estimation:msgDispSsest1');
   idDisplayEstimationInfo('Intro',{Str,' '},W);
   W.STOP = true;
end

%% Perform estimation. Options are contained in sys.
try
   [sys, X0] = n4sid_(sys, EstimData, Order);
   if isequal(sys,[]), return, end % GUI order selection mode
catch E
   if Disp
      S{1} = sprintf('<font color="red">%s</font>',E.message);
      S{2} = ctrlMsgUtils.message('Ident:estimation:msgAbortEstimation');
      idDisplayEstimationInfo('Error',S,W);
   end
   throw(E)
end

%% Reconcile metadata between model and data.
sys = copyEstimationDataMetaData(sys, EstimData);

%---------------------------------------------------------------------------
function [sys, data, Order] = localProcessInputs(data, Order, varargin)
% Process input arguments of the N4SID command.
%
% Syntaxes supported:
%  model = n4sid(data, orders, options) -> orders can be 'best', integer (array).
%  model = n4sid(data, model, options)
%  model = n4sid(data, <orders/model>, PV pairs) -> no options
%  model = n4sid(data, 'nx', orders, PV pairs)
%  model = n4sid(data, orders, ...GUIdata) -> undocumented
%  model = n4sid(data, orders, ... 'trace') -> deprecated; 'trace' can appear anywhere
%  model = n4sid(data, orders, 'trace', aux, dkx, maxsize, Ts);

% Note:
% Old incompatibility: If 'trace' is third input argument and there are
% more input arguments, processing may throw error.
% InitialState values will be mapped to 'zero', 'estimate' or numeric
% vector.

Nymod = NaN; Numod = NaN; Tsmod = NaN; nk = []; Tsumod = '';
% Replace (data, 'nx', N, ..) with (data, N, ..)
if strcmpi(Order,'nx')
   if nargin==2
      ctrlMsgUtils.error('Ident:estimation:n4sidCheck1')
   else
      Order = varargin{1};
      varargin = varargin(2:end);
   end
end

% Validate order
InitModel = false; CharOrder = false;
if isempty(Order) || strcmpi(Order,'best')
   Order = 1:10;
   CharOrder = true;
elseif isnumeric(Order)
   Order = Order(:).';
   if ~idpack.isNonnegIntVector(Order)
      ctrlMsgUtils.error('Ident:estimation:ssInvalidModelOrder')
   end
elseif isa(Order,'lti') && ~isa(Order, 'frd')
   if ndims(Order)>2 %#ok<ISMAT>
      ctrlMsgUtils.error('Ident:general:unsupportedCommandForArrays','n4sid')
   end
   sys = idss(Order);
   InitModel = true;
   [Nymod, Numod] = size(sys);
   Tsmod = sys.Ts;
   ISB = getDefaultISB(sys);
   TU = sys.TimeUnit;
   Order = order(sys);
else
   ctrlMsgUtils.error('Ident:general:InvalidSyntax','n4sid','n4sid')
end

% Handle deprecated syntax n4sid(data, orders, 'trace', ny, aux, dkx, maxsize, Ts)
I = idpack.findOptionInList('trace',varargin,2);

DispV5 = false; V5Syntax = false;
if numel(I)>1
   ctrlMsgUtils.error('Ident:general:InvalidSyntax','n4sid','n4sid')
elseif ~isempty(I) || (~isempty(varargin) && ~ischar(varargin{1}))
   V5Syntax = true;
   if ~isempty(I) && length(varargin)>I && ...
         any(strncmpi(varargin{I+1},{'on','off'},min(2,length(varargin{I+1}))))
      % trace specified as PV pair
      varargin{I} = 'display';
   elseif ~isempty(I) && (I==1 || (I>1 && ~ischar(varargin{1})))
      % Old syntax; 'trace' appears by itself.
      varargin(I) = [];
      DispV5 = ~isempty(I);      
      % Regenerate varargin
      [Nymod, Tsmod, nk, varargin] = localIntepretV5Args(varargin);
   elseif isempty(I) && all(cellfun(@(x_)isnumeric(x_),varargin))
      % n4sid(z,Nx,Ny,auxord,dkx,...);
      [Nymod, Tsmod, nk, varargin] = localIntepretV5Args(varargin);
   else
      V5Syntax = false;
   end
end

% Determine Ts
Tsi = find(strcmpi(varargin,'Ts'));
if ~isempty(Tsi)
   if length(varargin)<=Tsi(end)
      ctrlMsgUtils.error('Ident:general:CompleteOptionsValuePairs','n4sid','n4sid')
   else
      Tsmod = ltipack.utValidateTs(varargin{Tsi(end)+1});
      varargin([Tsi, Tsi+1]) = [];
   end
   if InitModel, sys.Ts = Tsmod; end
end
Tsui = idpack.findOptionInList('TimeUnit',varargin);
if ~isempty(Tsui)
   if length(varargin)<=Tsui(end)
      ctrlMsgUtils.error('Ident:general:CompleteOptionsValuePairs','n4sid','n4sid')
   else
      Tsumod = varargin{Tsui(end)+1};
      varargin([Tsui, Tsui+1]) = [];
   end
   if InitModel, sys.TimeUnit = Tsumod; end
end

DoubleData = isnumeric(data) && ismatrix(data) && size(data,1)>size(data,2);
FRdata = false;
if isa(data,'iddata') || isa(data,'frd')
   if isa(data, 'frd')
      if size(data,2)==0
         ctrlMsgUtils.error('Ident:estimation:FRDataTimeSeriesModel');
      elseif ndims(data)>2 %#ok<ISMAT>
         error(message('Ident:estimation:frdArrayForEstimation'))
      end
      data = idfrd(data);
      FRdata = true;
      [Nydat, Nudat] = size(data);
      Tsdat = data.Ts;
      Ncaps = length(data.Frequency);
   else
      [~,Nydat, Nudat] = size(data);
      Tsdat = pvget(data, 'Ts'); Tsdat = Tsdat{1};
      Ncaps = size(data,1);
   end   
elseif DoubleData
   if ~isnan(Nymod)
      % Either initial model or obsolete option was used to define Ny
      Nydat = Nymod;
   else
      Nydat = 1;
      if size(data,2)>1
         ctrlMsgUtils.warning('Ident:general:doubleDataNyAmbiguity')
      elseif size(data,2)==0
         ctrlMsgUtils.error('Ident:estimation:noOutputChannel')
      end
   end
   Nudat = size(data,2)-Nydat;
   Ncaps = size(data,1);
   if Nudat<0
      ctrlMsgUtils.error('Ident:general:modelDataDimMismatch')
   end
   Tsdat = 1;
   data = iddata(data(:,1:Nydat), data(:,Nydat+1:end), Tsdat); % Sample time could change
   if InitModel
      data.InterSample = ISB;
      data.TimeUnit = TU;
   end
else
   ctrlMsgUtils.error('Ident:general:InvalidSyntax','n4sid','n4sid')
end

if isnan(Nymod), Nymod = Nydat; end
if isnan(Numod), Numod = Nudat; end
if Nymod==0
   error(message('Ident:estimation:noOutputChannel'))
elseif Nymod~=Nydat || Numod~=Nudat
   ctrlMsgUtils.error('Ident:general:modelDataDimMismatch')
end

if V5Syntax
   % Fix nk.
   if isscalar(nk)
      nk = repmat(nk,[1,Numod]);
   elseif Tsmod==0 && any(nk>1)
      ctrlMsgUtils.warning('Ident:estimation:CTModelNkVal')
      nk = min(nk,1);
   end
   varargin = [varargin, {'nk', nk}];
   
   % Update display for old syntax
    if DispV5, varargin = [varargin, {'Display','on'}]; end 
end

if ~InitModel
   % Create template IDSS model. D and X0 are fixed at zero values while K
   % is zero but free.
   if isnan(Tsmod)
      % Ts was not specified; make it same as first data Ts.
      Tsmod = Tsdat;
   end
   
   nxx = min(Order);
   c = [eye(Nymod),zeros(Nymod,nxx)]; c = c(:,1:nxx);
   RS = rng;
   S = pmodel.ss(randn(nxx), ones(nxx,Numod), c, ...
      zeros(Nymod, Numod),zeros(nxx, Nymod));
   S.k.Free = true;  % estimate disturbance model by default
   S.d.Free = nxx==0; % feedthrough = 'off' by default unless static
   SSData = idpack.ssdata(S,Tsmod);
   SSData.X0.Value = zeros(nxx,1); % estimate initial states by default
   sys = idss.make(SSData,[Nymod,Numod]);
   sys = cacheDefaultRandStream(sys,RS);
   if ~isempty(Tsumod)
      sys.TimeUnit = Tsumod;
   elseif ~DoubleData
      sys.TimeUnit = data.TimeUnit;
   end
   
   % Check data size.   
   if CharOrder && (min(Ncaps) <= ...
         (ceil(10/Nymod)+1)*(1+Nymod+Numod) + ...
         (1+Nymod+Numod)*ceil((10-Nymod+1)/(Nymod+Numod)))
      ctrlMsgUtils.error('Ident:estimation:tooFewSamplesOrderTest')
   end
else
   % Refresh rand stream
   sys = cacheDefaultRandStream(sys);
end

if DoubleData
   Tsdat = abs(Tsmod);
   data.Ts = Tsdat;
   data.TimeUnit = sys.TimeUnit;
end

% Handle 'Feedthrough' which is not a current or obsolete IDSS property.
% Since user could pass a mix of 'nk' and Feedthrough' specifications in
% any order, the feedthrough specification must be replaced by 'nk' in
% place.
FtInd = idpack.findOptionInList('Feedthrough',varargin,2);
if ~isempty(FtInd)
   if length(varargin)<=FtInd(end)
      ctrlMsgUtils.error('Ident:general:CompleteOptionsValuePairs','n4sid','n4sid')
   else
      FtValue = varargin{FtInd(end)+1}(:).';
      if (~islogical(FtValue) && ~isequal(logical(FtValue),FtValue)) || ...
            (~isscalar(FtValue) && length(FtValue)~=Numod)
         ctrlMsgUtils.error('Ident:estimation:IncorrectFeedthruValue')
      end
      varargin([FtInd(1:end-1),FtInd(1:end-1)+1]) = [];
      if isscalar(FtValue), FtValue = FtValue(ones(1,Numod)); end
      FtValue = logical(FtValue);
      varargin(FtInd(end):FtInd(end)+1) = {'nk',double(~FtValue)};
   end
end

% Handle 'SSParameterization'/'Form' specification separately. Setting it
% as a regular model property will cause an unnecessary transformation to
% canonical form.
SSInd1 = idpack.findOptionInList('SSParameterization',varargin,2);
SSInd2 = idpack.findOptionInList('Form',varargin,2);
SSInd = [SSInd1, SSInd2]; 
if ~isempty(SSInd)
   if length(varargin)<=SSInd(end)
      ctrlMsgUtils.error('Ident:general:CompleteOptionsValuePairs','n4sid','n4sid')
   else
      if strcmpi(varargin{SSInd(end)+1},'c')
         SSValue = 'canonical'; % b.c.
      else
         SSValue = ltipack.matchKey(varargin{SSInd(end)+1}, ...
            {'free', 'canonical', 'structured', 'modal', 'companion'});
      end
      if isempty(SSValue)
         ctrlMsgUtils.error('Ident:estimation:ssestFormCheck','n4sid')
      else
         varargin([SSInd,SSInd+1]) = [];
      end
   end
end

% Construct estimation options and SET PV pairs.
[sys, data] = idpack.utPrepareEstimationOptions(sys, data, varargin, ...
   n4sidOptions, 'n4sid', ~InitModel);

Options = getDefaultOptions(sys);
ut = Options.Utility;

if ~isempty(SSInd)
   ut.SSParameterization = SSValue;
   % ut.CanonicalIndices = sys.CanonicalIndices; % required for conversion to canonical form.
end
% Transform initial state specification.
optinit = Options.IC_;

% Use zero i.c. if idfrd data.
if FRdata && ~strcmpi(optinit,'zero')
%    if ~strcmpi(optinit,'auto')
%       ctrlMsgUtils.warning('Ident:estimation:X0EstForIDFRD1')
%    end
   optinit = 'zero';
end

% N4SID can only handle 'zero' or 'estimate' i.c.
if strcmpi(optinit,'auto')
   optinit = 'estimate';
elseif strcmpi(optinit,'backcast')
   optinit = 'zero'; % (b.c.)
end

Options.IC_ = optinit;

% Update interactive vs automatic order selection
ut.Interactive = ~CharOrder;
Options.Utility = ut;

sys = setDefaultOptions(sys,Options);

%--------------------------------------------------------------------------
function [ny,Ts,nk,arg] = localIntepretV5Args(arg)
% Translate old (v5) n4sid syntax.

ctrlMsgUtils.warning('Ident:estimation:n4sidOldSyntax')

try
   ni = length(arg);
   auxord = []; Kest = 'Estimate'; estX = 'Auto'; MaxSize = 'Auto'; Ts = 1;
   nk = 0; ny = 1;
   
   if ni>0
      ny = arg{1};
      if ~idpack.isPosIntScalar(ny)
         ctrlMsgUtils.error('Ident:general:PosIntScalarRequired','ny')
      end
      if ni>1
         arg2 = arg{2};
         if isempty(arg2)
            auxord = 'Auto';
         else
            if isrow(arg2)
               arg2 = arg2';
            end
            auxord = repmat(arg2,[1 3]);
         end
         
         if ni>2
            dkx = arg{3};
            if dkx(1), nk = 0; else nk = 1; end
            if ~dkx(2), Kest = 'None'; end
            if dkx(3), estX = 'Estimate'; else estX = 'Zero'; end
            
            if length(dkx)>3, nk = dkx(4:end); end
            if ni>3
               MaxSize = arg{4};
               if ni>4
                  Ts = arg{5};
               end
            end
         end
      end
   end
   
   arg = {'N4Horizon',auxord,'InitialState',estX,'DisturbanceModel',Kest,...
      'MaxSize',MaxSize};
catch %#ok<CTCH>
   ctrlMsgUtils.error('Ident:general:InvalidSyntax','n4sid','n4sid')
end
