classdef WecSystemModel < handle
    %WECSYSTEMMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    
    % ---------------------------------------------------------------------
    % Private properties
    properties (SetAccess = private)
        mass
        kHyd
        bGen
        Ainf
        addMass
        hydFile
        radIrf
        radFreq
        feIrf
        feFreq 
        simulinkModel = 'wec_model';
        ssRad
    end
    
    properties (SetAccess = private, GetAccess = private)
        verify_results
    end
    
    
    % ---------------------------------------------------------------------
    % Public Methods 
    %
    
    methods (Access = public)
        % Constructor
        function obj = WecSystemModel(hydFilename, varargin)
            % WecModel = WecSystemModel(hydFilename) uses default values
            %
            % WecModel = WecSystemModel(hydFilename, options) where options
            % is a structure will set parameters to values in options. See
            % below for available options fields
            %
            % Permissable options, default value listed in ()
            %
            % .mass  (140075.5) - dry mass in [kg] of wec
            % .kHyd  (789740) - hydrostatic force coefficient
            % .bGen  (1e5) - constant generator damping ratio in N-s/m
            % .feIrf.t (-25:0.05:25) - time span used for fe IRF
            % .radIrf.t (-0:0.05:20) - time span for radiation IRF
            
            obj.hydFile = hydFilename;  
            
            
            def.mass = 140076.5;
            def.kHyd = 1025 * 9.81 * pi * 25;
            def.feIrf.t = -25:0.05:25;
            def.radIrf.t = 0:0.05:20;
            def.bGen = 1e5;
            
            % if option was given, set the corresponding property to the
            % value. Otherwise use the default value
            fields = fieldnames(def);
            if nargin == 2
                for ii = 1:numel(fields)
                    if isfield(varargin{1}, fields{ii})
                        obj.(fields{ii}) = varargin{1}.(fields{ii});
                    else
                        obj.(fields{ii}) = def.(fields{ii});
                    end
                end
            else % no input = set all defaults
                for ii = 1:numel(fields)
                    obj.(fields{ii}) = def.(fields{ii});
                end
            end
            
            
            % Load hydrdynamic data file
                     
            [hydParms] = obj.read_hyd_file(obj.hydFile);
            obj.feFreq = [hydParms.freq hydParms.dfkMag hydParms.dfkPhase];
            obj.radFreq = [hydParms.freq hydParms.radDamp];
            obj.addMass = [hydParms.freq hydParms.addMass];
            obj.Ainf    = obj.addMass(end,2);
            clear hydParms
            


            % Construct my impulse response functions
            obj.feIrf  = calc_fe_irf(obj);
            obj.radIrf = calc_radiation_irf(obj);
            
            % TODO
            % Construct my state space model including radiation force
            % stuff
            
            
        end % WecSystemModel
       
        
        % Run simulation
        function logsout = run_simulation(obj, eta, dt)
            % calculate excitation
            [fe, feTime] = calc_excitation(obj, eta, dt);
            
            % ramp in fe
            fe(1:200) = linspace(0,1,200) .* fe(1:200);
            
            % Convert ratiation irf to a trapezoidal coefficients to be
            % used in simulation
            n  = length(obj.radIrf.t) - 1;
            tr = obj.radIrf.t(end) - obj.radIrf.t(1);
            radIrfDt = tr / n;
            radIrfCoeff = (tr/(2*n)) .* ...
                [ 1 repmat(2,1,length(obj.radIrf.irf)-2) 1 ] .* ...
                obj.radIrf.irf;
         
            
            % load and run simulink model
            
            if bdIsLoaded(obj.simulinkModel) == true
                close_system(obj.simulinkModel, true);
            end
            load_system(obj.simulinkModel)
            % assign values to simulink workspace
            hws = get_param(obj.simulinkModel, 'modelworkspace');
            hws.assignin('m', obj.mass);
            hws.assignin('b_gen', obj.bGen);
            hws.assignin('kHyd', obj.kHyd);
            hws.assignin('Ainf', obj.Ainf);
            hws.assignin('feSeries', [feTime' fe']);
            hws.assignin('feDt', dt);
            hws.assignin('radIrfDt', radIrfDt);
            hws.assignin('radIrfCoeff', radIrfCoeff);            
            % run model
            set_param(obj.simulinkModel, 'StopTime', num2str(feTime(end)));
            sim(obj.simulinkModel);             
            close_system(obj.simulinkModel, false);
        end % run_simulation
        

        
        function varargout = disp_irf(obj, identifier)
            switch identifier
                case 'rad'
                    t = obj.radIrf.t;
                    x = obj.radIrf.irf;
                case 'fe'
                    t = obj.feIrf.t;
                    x = obj.feIrf.irf;
                otherwise
                    error('bad identifier')
            end
            h = figure;
            set(h,'color','w')
            
            plot(t, x)
            xlabel('Time (sec)')
            ylabel('Impulse response')
            grid on
            
            if nargout == 1
               varargout{1} = h;
            end
        end % disp_irf()

        function SS = construct_state_space_model(obj)
            % State vector x = [Xrad z zdot]
            % Output y = [z zDot]
            
            % check to see if the approximation has been done already
            if isempty(obj.ssRad)
                warning(['State Space Radiation Approximation not calculated.' ...
                    'calculating using default values'])
                state_space_radiation_approx(obj, 3); 
            end
            
            % Build the state space model
            n = size(obj.ssRad.A, 1);
            a = (obj.mass + obj.Ainf);
            SS.A = [ obj.ssRad.A zeros(n,1) obj.ssRad.B ;
                     zeros(1,n+1) 1 ;
                     -obj.ssRad.C ./ a -obj.kHyd / a -obj.bGen/a ];
            SS.B = [zeros(n+1, 1) ; 1/a];
            SS.C = [zeros(2,n)  eye(2) ];
            SS.Bd = [];
            % clear the old verification results
            obj.verify_results = [];
        end % construct_state_space_model()
        
        function set_bgen(obj, bGen)
            obj.bGen = bGen;
        end
        
    end % methods (public)
    
    % ---------------------------------------------------------------------
    % Private Methods stored in seperate files
    %
    % these methods are stored in the class folder seperately
    methods (Access = private)
        radIrf = calc_radiation_irf(obj, varargin);
        feIrf  = calc_fe_irf(obj, varargin);
        
        verify_freq_response(obj, freqs);  
    end % methods (Access = private)
    
    
    % ---------------------------------------------------------------------
    % Public Methods stored in seperate files
    % 
    % These methods are stored in other files but accessable to the user
    methods (Access = public)
        [varargout] = state_space_radiation_approx(obj, modelOrder, varargin);
        simResults = run_state_space_simulation(obj, eta, dt, varargin)
        verify_state_space_model(obj, varargin);
        [fe, feTime, varargout] = calc_excitation(obj, eta, dt);
        [varargout] = radiation_force_balanced_reduction(obj, tol, varargin)
    end
    
    
    
    % ---------------------------------------------------------------------
    % Static Methods
    
    methods(Static)
        [hydParms] = read_hyd_file(filename);
        
        function obj = set_options(options)
            % takes the input optiosn structure and sets all of these
            % options. If it was empty it returns the default values
            
            % Specify default
            def.mass = 140076.5;
            def.kHyd = 1025 * 9.81 * pi * 25;
            def.feIrf.t = -25:0.05:25;
            def.radIrf.t = 0:0.05:20;
            def.bGen = 1e5;
            
            % if option was given, set the corresponding property to the
            % value. Otherwise use the default value
            fields = fieldnames(def);
            if nargin == 1
                for ii = 1:numel(fields)
                    if isfield(options, fields{ii})
                        obj.(fields{ii}) = options.(fields{ii});
                    else
                        obj.(fields{ii}) = def.(fields{ii});
                    end
                end
            else % no input = set all defaults
                for ii = 1:numel(fields)
                    obj.(fields{ii}) = def.(fields{ii});
                end
            end
        end % set_options()
         
    end    
    
end

