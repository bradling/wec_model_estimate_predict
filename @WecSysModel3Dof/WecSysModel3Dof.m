classdef WecSysModel3Dof < handle
    %WECSYSTEMMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    
    % ---------------------------------------------------------------------
    % Private properties
    properties (SetAccess = private)
        inertia
        kHyd
        bGen
        Ainf
        addMass
        radIrf
        radFreq
        freq
        feIrf
        feMag
        fePhase
        ssRad
        ssRadApproxOrder
    end
    
    properties (SetAccess = public)
        fcnPto = [];
        fcnMooring = [];
        userData
    end
    
    properties (SetAccess = private, GetAccess = private)
        verify_results
    end
    
    
    % ---------------------------------------------------------------------
    % Public Methods 
    %
    
    methods 
        % Constructor
        function obj = WecSysModel3Dof(hydFilename, varargin)
            % WecModel = WecSystemModel3Dof(ansys_ls_file) uses default values
            %
            % WecModel = WecSystemModel3Dof(ansys_ls_file, options) uses
            %             options in options structure
            %
            % Options structure fields
            % .bGen  (1e5) - constant generator damping ratio in N-s/m
            % .feIrf.dt = 0.05;
            % .feIrf.t (-20:0.05:20) - time span used for fe IRF
            % 
            % .radIrf.dt = 0.05;
            % .radIrf.t (-0:0.05:20) - time span for radiation IRF
                        
            if nargin > 0 
            
            def.feIrf.dt  = 0.05;
            def.feIrf.t   = -20:def.feIrf.dt:20;
            def.radIrf.dt = 0.05;
            def.radIrf.t  = 0:def.radIrf.dt:20;
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
            
            ext = hydFilename(end-2:end);
            switch ext
                case 'LIS'
                    [hydParms] = WecSysModel3Dof.read_ansys_output(hydFilename);
                case 'mat'
                    hydParms = load(hydFilename);
            end
            
            
            obj.freq    = hydParms.freq;
            obj.feMag   = hydParms.dfkMag(:, [1 3 5]);
            obj.fePhase = hydParms.dfkPhase(:, [1 3 5]);
            obj.radFreq = hydParms.RadDamp(:, [1 3 5 7 8 11]);
            obj.addMass = hydParms.AddMass(:, [1 3 5 7 8 11]);
            obj.Ainf    = hydParms.AddMass(end, [1 3 5 7 8 11]);
            obj.kHyd    = hydParms.K([1 3 5], [1 3 5]);
            obj.inertia = hydParms.inertia([1 3 5], [1 3 5]);
            clear hydParms
            


            % Construct my impulse response functions
            obj.feIrf.irf  = nan(size(obj.feMag,2), size(obj.feIrf.t,2));
            
            for ii = 1:size(obj.feIrf.irf,1)
                [~, obj.feIrf.irf(ii,:)]  = calc_fe_irf(obj, obj.feMag(:,ii), obj.fePhase(:,ii));
                
            end
            
            obj.radIrf.irf = nan(size(obj.radFreq,2), size(obj.radIrf.t,2));
            for ii = 1:size(obj.radFreq,2)
                [~, obj.radIrf.irf(ii,:)] = calc_radiation_irf(obj, obj.radFreq(:,ii));
            end
            % TODO
            % Construct my state space model including radiation force
            % stuff
            
            end
            
            
        end % WecSystemModel
        
        function varargout = construct_rad_approx(obj, order, varargin)
            [obj.ssRad, obj.ssRadApproxOrder] = rad_ss_approx(obj, order);
            
            if nargin > 2
                [h] = obj.plot_rad_approx();
                if nargout > 0, varargout{1} = h; end
            else
                if nargout > 0, varargout{1} = []; end
            end
        end
       
        function varargout = disp_freq(obj)
            % generates frequency response data plots
            
            h1 = figure;
            set(h1, 'color', 'w')
            set(h1, 'units', 'inches')
            set(h1, 'position', [1 1 5 5])
            set(h1, 'name', 'Excitation Frequency Response')
            
            subplot(2,1,1)
            plot(obj.freq, obj.feMag'./1e3)
            grid on
            xlabel('Frequency (rad/s)')
            ylabel('Magnitude (kN/m)')
            legend('Surge', 'Heave', 'Pitch', 'location', 'northeast')
            
            subplot(2,1,2)
            plot(obj.freq, obj.fePhase')
            grid on
            xlabel('Frequency (rad/s)')
            ylabel('Phase (rad)')
            
            
            h2 = figure;
            set(h2, 'color', 'w')
            set(h2, 'units', 'inches')
            set(h2, 'position', [1.3 0.7 5 5])
            set(h2, 'name', 'Radiation Frequency Response')
            
            subplot(2,1,1)
            plot(obj.freq, obj.radFreq'./1e3)
            grid on
            xlabel('Frequency (rad/s)')
            ylabel('Radiation Mag (kN-s/m)')
            title('Radiation Damping Coefficient')
            legend('Surge', 'Heave', 'Pitch', 'Surge-Heave', ...
                'Surge-Pitch', 'Heave-Pitch', 'location', 'northwest')
            
%             h3 = figure;
%             set(h3, 'color', 'w')
%             set(h3, 'units', 'inches')
%             set(h3, 'position', [1.6 0.4 5 3.2])
%             set(h3, 'name', 'Added Mass Frequency Response')
            
            subplot(2,1,2)
            plot(obj.freq, obj.addMass'./1e3)
            grid on
            title('Added Mass Coefficient')
            xlabel('Frequency (rad/s)')
            ylabel('Added Mass (kN-s^2/m)')
            legend('Surge', 'Heave', 'Pitch', 'Surge-Heave', ...
                'Surge-Pitch', 'Heave-Pitch', 'location', 'northwest')
            
            if nargout > 0, varargout{1} = [h1 h2]; end
            
        end %disp_freq()
        
        function varargout = disp_fe_irf(obj)
            h = figure;
            set(h, 'color', 'w')
            set(h, 'units', 'inches')
            set(h, 'position', [1.3 0.6 5 5.5])
            set(h, 'name', 'Fe Impulse Responses')
            
            subplot(3,1,1)
            plot(obj.feIrf.t, obj.feIrf.irf(1,:)./1e3, 'r')
            ylabel('Fe Irf (kN/m)')
            grid on
            title('Surge')
            
            subplot(3,1,2)
            plot(obj.feIrf.t, obj.feIrf.irf(2,:)./1e3, 'g')
            ylabel('Fe Irf (kN/m)')
            grid on
            title('Heave')
            
            subplot(3,1,3)
            plot(obj.feIrf.t, obj.feIrf.irf(3,:)./1e3, 'b')
            xlabel('Time (sec)')
            ylabel('Fe Irf (kN/m)')
            grid on
            title('Pitch')
            
            if nargout > 0, varargout{1} = h; end
        end % disp_fe_irf()

        function ss = construct_state_space_model(obj)
            % ss = obj.construct_state_space_model()
            %
            % ss is an output strcutre containing a, b, c, d state space
            % matrics in the following form:
            %
            % zdot = ss.a*z + ss.b*u + ss.d*w,
            %    y = ss.c*z
            %
            % NOTES: - model is in continuous time
            %        - mooring forces and PTO forces are not currently
            %          included, assumed to be a part of the input u
            
            % First lets build the radiation submodel
            m = obj.ssRadApproxOrder;
            Ar = zeros(6*m);
            jj = 1;
            for ii = 1:m:(6*m)-1
                Ar(ii:ii+m-1,ii:ii+m-1) = obj.ssRad{jj}.a;
                jj=jj+1;
            end
            Br = [obj.ssRad{1}.b zeros(m,1)     zeros(m,1) ; 
                  zeros(m,1)     obj.ssRad{2}.b zeros(m,1) ;
                  zeros(m,1)     zeros(m,1)     obj.ssRad{3}.b ;
                  obj.ssRad{4}.b zeros(m,1)     zeros(m,1) ; 
                  obj.ssRad{5}.b zeros(m,1)     zeros(m,1) ;
                  zeros(m,1)     obj.ssRad{6}.b zeros(m,1) ];
            Cr = [obj.ssRad{1}.c zeros(1,m)     zeros(1,m)     obj.ssRad{4}.c obj.ssRad{5}.c zeros(1,m) ;
                  zeros(1,m)     obj.ssRad{2}.c zeros(1,m)     obj.ssRad{4}.c zeros(1,m)     obj.ssRad{6}.c ;
                  zeros(1,m)     zeros(1,m)     obj.ssRad{3}.c zeros(1,m)     obj.ssRad{5}.c obj.ssRad{6}.c ];
            Dr = [(obj.ssRad{1}.d+obj.ssRad{4}.d+obj.ssRad{5}.d) 0               0 ;
                  (obj.ssRad{4}.d + obj.ssRad{5}.d)              obj.ssRad{2}.d  0 ;
                   obj.ssRad{5}.d                                obj.ssRad{6}.d  obj.ssRad{3}.d];
               
            % Now combine the other things with the radiation model for
            % full state space model
            
            %TODO - Add Generator force and mooring force components.
            
            Ainf = [obj.Ainf(1) obj.Ainf(4) obj.Ainf(5) ;
                    obj.Ainf(4) obj.Ainf(2) obj.Ainf(6) ;
                    obj.Ainf(5) obj.Ainf(6) obj.Ainf(3) ];
            Jinv = inv(obj.inertia + Ainf);
            ss.a = [-(Jinv * Dr) -(Jinv * obj.kHyd) -(Jinv * Cr)   ; 
                     eye(3)        zeros(3)           zeros(3,6*m) ; 
                     Br            zeros(6*m, 3)      Ar           ];
            ss.b = [Jinv ; zeros(3+6*m, 3) ];
            ss.d = [Jinv ; zeros(3+6*m, 3) ];
            ss.c = [eye(6) zeros(6, 6*m)];
            
        end % construct_state_space_model()

        function set.fcnMooring(obj, fcn)
            % obj.fcnMooring = @fcn, wher @fcn is a function handle with
            % the following form:
            %
            % Fm = @fcn(s,sdot), where Fm is the column vector of mooring
            % forces, s and sdot are column vectors for position and
            % velocity of the device.
            
            if ~isa(fcn, 'function_handle')
                error('Must be a function handle');
            end
            obj.fcnMooring = fcn;
        end
        
        function set.fcnPto(obj, fcn)
            % obj.fcnPto = @fcn, wher @fcn is a function handle with
            % the following form:
            %
            % F = @fcn(s,sdot), where Fpto is the column vector of PTO
            % forces, s and sdot are column vectors for position and
            % velocity of the device.
            
            if ~isa(fcn, 'function_handle')
                error('Must be a function handle');
            end
            obj.fcnPto = fcn;
        end
        
        function newObj = copy(obj)
            newObj = WecSysModel3Dof();
            props = properties(obj);
            for ii = 1:numel(props)
                newObj.(props{ii}) = obj.(props{ii});
            end
        end
        
        function verify_model(obj)
            % verify model results
            omega = linspace(obj.freq(1), obj.freq(end), floor(length(obj.freq)/2));
            dt = 0.1;
            t = 0:dt:1000;
            feMag = nan(3,length(omega));
            rao = nan(3, length(omega));
            for ii = 1:length(omega)
                fprintf('%.0f of %.0f\n',ii,length(omega));
                eta = sin(omega(ii).*t);
                
                results = obj.run_sim(eta,dt);
                idxEnd = floor(0.8.*length(results.t)):(length(results.t));
                feMag(1,ii) = max(results.fe(1,idxEnd));
                feMag(2,ii) = max(results.fe(2,idxEnd));
                feMag(3,ii) = max(results.fe(3,idxEnd));
                rao(1,ii) = max(abs(results.s(1,idxEnd)));
                rao(2,ii) = max(abs(results.s(2,idxEnd)));
                rao(3,ii) = max(abs(results.s(3,idxEnd)));
                
                
                
            end
            
            figure
            plot(obj.freq, obj.feMag(:,1), 'r-', ...
                 obj.freq, obj.feMag(:,2), 'g-', ...
                 obj.freq, obj.feMag(:,3), 'b-')
            hold on
            plot(omega, feMag(1,:), 'rx', ...
                 omega, feMag(2,:), 'go', ...
                 omega, feMag(3,:), 'bd')
            hold off
            legend('Surge','Heave','Pitch','location','best')
            xlabel('Freq (rad/s)')
            ylabel('Excitation Magnitude Response')
            grid on
            
            
        end
        
    end % methods (public)
    
    % ---------------------------------------------------------------------
    % Private Methods stored in seperate files
    %
    % these methods are stored in the class folder seperately
    methods (Access = private)
        [t, irf] = calc_radiation_irf(obj, mag);
        [t, irf]  = calc_fe_irf(obj, mag, phase);
        [ssRad, order] = rad_ss_approx(obj, order)
        
%        verify_freq_response(obj, freqs);  
    end % methods (Access = private)
    
    
    % ---------------------------------------------------------------------
    % Public Methods stored in seperate files
    % 
    % These methods are stored in other files but accessable to the user
    methods (Access = public)
        [varargout] = plot_rad_approx(obj)
        
        results = run_sim(obj, eta, dt);
        
%        [varargout] = state_space_radiation_approx(obj, modelOrder, varargin);
%        simResults = run_state_space_simulation(obj, eta, dt, varargin)
%        verify_state_space_model(obj, varargin);
        [fe, feTime, varargout] = calc_excitation(obj, eta, dt);
%        [varargout] = radiation_force_balanced_reduction(obj, tol, varargin)
    end
    
    
    
    % ---------------------------------------------------------------------
    % Static Methods
    
    methods(Static)
        hydParms = read_ansys_output(filename)
        
        function obj = set_options(options)
            % takes the input optiosn structure and sets all of these
            % options. If it was empty it returns the default values
            
            % Specify default
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

