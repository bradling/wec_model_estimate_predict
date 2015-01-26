classdef WecSystemModel < handle
    %WECSYSTEMMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    
    % TODO:
    %  Figure out a good way to specify mass, kHyd, bGen,     
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
    
    methods (Access = public)
        % Constructor
        function obj = WecSystemModel(hydFilename, varargin)
            % WecModel = WecSystemModel(hydFilename)
            %
            % WecModel = WecSystemModel(hydFilename, bGen) also sets the
            % constant generator damping value.
            
            obj.bGen = 0;
            if nargin == 2
                obj.bGen = 1e5;
            end
            
            % Load hydrdynamic data file
            obj.hydFile = hydFilename;           
            [hydParms] = read_hyd_file(obj);
            obj.feFreq = [hydParms.freq hydParms.dfkMag hydParms.dfkPhase];
            obj.radFreq = [hydParms.freq hydParms.radDamp];
            obj.addMass = [hydParms.freq hydParms.addMass];
            obj.Ainf    = obj.addMass(end,2);
            clear hydParms
            
            % these are hardcoded for now. Need to make this dynamic
            obj.mass = 140076.5;
            obj.kHyd = 1025 * 9.81 * pi * 25;

            % Construct my impulse response functions
            obj.feIrf  = calc_fe_irf(obj);
            obj.radIrf = calc_radiation_irf(obj);
            
            % TODO
            % Construct my state space model including radiation force
            % stuff
            
            
        end
       
        
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
        end
        

        
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
        end
    end
    
    methods (Access = public)
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
            
            % clear the old verification results
            obj.verify_results = [];
        end
    end
    
    % these methods are stored in the class folder seperately
    methods (Access = private)
        [hydParms] = read_hyd_file(obj);
        radIrf = calc_radiation_irf(obj, varargin);
        feIrf  = calc_fe_irf(obj, varargin);
        [fe, feTime] = calc_excitation(obj, eta, dt);
        verify_freq_response(obj, freqs);
    end
    
    % These methods are stored in other files but accessable to the user
    methods (Access = public)
        [varargout] = state_space_radiation_approx(obj, modelOrder, varargin);
        simResults = run_state_space_simulation(obj, eta, dt, varargin)
        verify_state_space_model(obj, varargin);
    end
    
    
end

