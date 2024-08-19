classdef RefTaylorState < TaylorState
    %REFTAYLORSTATE A subclass representing a reference Taylor state
    %   A reference Taylor state in a domain Omega bounded by m toroidal 
    %   surfaces is a solution B to the BVP
    % 
    %          curl B = k B           in Omega
    %           n . B = n . B_0       on the boundary of Omega
    %     \oint_{S_i} = Phi_i(B_0)    i = 1, ..., m
    % 
    %   where n is the outward unit normal to Omega, S_i is a topologically 
    %   distinct cross-sectional surface through Omega, B_0 is a magnetic 
    %   field that satisfies the PDE, and Phi_i is the flux by B_0 through
    %   S_i. A uniqueness theorem says that the solution B equals B_0. This
    %   class is used for testing. 
    % 
    %   NOTE: Currently assumes that the cross-sectional surface is 
    %   parallel to the xz-plane
    %   
    %   TODO: 
    %     1. Allow for other cross-sections
    %     2. Compute flux here (rather than req. precompute)?

    properties
        B0 % solution to curl B = k B
        xs_nodes % cross-section (XS) quadrature nodes 
        xs_weights % cross-section quadrature weights 
        quad_opts_xs_taylor % near quad. for grad/curl Sk for xs_nodes 
        quad_opts_xs_laphelm % near quad. for Sk for xs_nodes
        eps_xs_taylor % tolerance for grad and curl Sk 
        eps_xs_laphelm % tolerance for Sk
    end

    methods
        function obj = RefTaylorState(domain,domparams,zk,flux,B0, ...
                xs_nodes,xs_weights,tols)
            %REFTAYLORSTATE Construct an instance of a RefTaylorState
            %   Arguments:
            %     domain [surfacemesh] toroidal boundary surface
            %     domparams [int(3)] parameters describing surface
            %     zk [double complex] Beltrami parameter
            %     flux [double] XS flux due to B0 (precomputed)
            %     B0 [surfacefunv] solution to curl B = k B
            %     xs_nodes [double(3,x)] quad. nodes for XS
            %     xs_weights [double(x)] quad. weights for XS
            %     tols [double or double(3)] quad. and GMRES tolerances
            [tolsr, tolsc] = size(tols);
            if isscalar(tols)
                supertols = tols;
            elseif tolsr == 5 && tolsc == 1 || tolsr == 1 && tolsc == 5
                supertols = tols(1:3);
            else
                error(['Invalid call to RefTaylorState constructor. ' ...
                    'Seventh argument should be 1 or 5 tolerance(s). '])
            end
            obj@TaylorState(domain,domparams,zk,flux,supertols);
            if isscalar(tols)
                obj.eps_xs_taylor = tols;
                obj.eps_xs_laphelm = tols;
            else
                obj.eps_xs_taylor = tols(4);
                obj.eps_xs_laphelm = tols(5);
            end
            
            if isa(B0, 'surfacefunv')
                obj.B0 = B0;
            else
                error(['Invalid call to RefTaylorState constructor. ' ...
                    'Fourth argument should be a surfacefunv ' ...
                    'corresponding with a solution to curl B = k B. '])
            end

            [xs_nodes_r, xs_nodes_c] = size(xs_nodes);
            [xs_weights_r, xs_weights_c] = size(xs_weights);
            if xs_nodes_r == 3 
                obj.xs_nodes = xs_nodes;
                n_nodes = xs_nodes_c;
            elseif xs_nodes_c == 3
                obj.xs_nodes = xs_nodes.';
                n_nodes = xs_nodes_r;
            end
            if xs_weights_r == 1 && xs_weights_c == n_nodes
                obj.xs_weights = xs_weights;
            elseif xs_weights_c == 1 && xs_weights_r == n_nodes
                obj.xs_weights = xs_weights;
            else
                error(['Invalid call to RefTaylorState constructor. ' ...
                    'Fifth argument, xs_nodes, and sixth argument, ' ...
                    'xs_weights, are arrays with incompatible sizes. '])
            end
        end

        function obj = get_quad_corr_xs_taylor(obj,varargin)
            %GET_QUAD_CORR_XS_TAYLOR Get near quad. for grad/curl Sk in XS
            %   Argument (optional)
            %     opts [struct]
            %       opts.format - Storage format for sparse matrices
            %         'rsc' - row sparse compressed format
            %         'csc' - column sparse compressed format
            %         'sparse' - sparse matrix format
            %       opts.quadtype - currently only 'ggq' supported
            format = 'rsc';
            if nargin > 1 && isa(varargin{1},'struct')
                opts = varargin{1};
                if ~isfield(opts,'format')
                    opts.format = format;
                else
                    format = opts.format;
                end
            else
                opts = [];
                opts.format = format;
            end
            
            targinfo = [];
            targinfo.r = obj.xs_nodes;
            if abs(obj.zk) < eps
                Q = taylor.static.get_quadrature_correction(obj.surf, ...
                    obj.eps_taylor,targinfo,opts);
            else
                Q = taylor.dynamic.get_quadrature_correction(obj.surf, ...
                    obj.zk,obj.eps_taylor,targinfo,opts);
            end
            
            obj.quad_opts_xs_taylor = [];
            obj.quad_opts_xs_taylor.format = format;
            obj.quad_opts_xs_taylor.precomp_quadrature = Q;
        end

        function obj = get_quad_corr_xs_laphelm(obj,varargin)
            %GET_QUAD_CORR_XS_LAPHELM Get near quad. for Sk in XS
            %   Argument (optional)
            %     opts [struct]
            %       opts.format - Storage format for sparse matrices
            %         'rsc' - row sparse compressed format
            %         'csc' - column sparse compressed format
            %         'sparse' - sparse matrix format
            %       opts.quadtype - currently only 'ggq' supported
            format = 'rsc';
            if nargin > 1 && isa(varargin{1},'struct')
                opts = varargin{1};
                if ~isfield(opts,'format')
                    opts.format = format;
                else
                    format = opts.format;
                end
            else
                opts = [];
                opts.format = format;
            end
            
            targinfo = [];
            targinfo.r = obj.xs_nodes;
            if abs(obj.zk) < eps
                Q = lap3d.dirichlet.get_quadrature_correction(obj.surf, ...
                    obj.eps_laphelm,[1.0,0],targinfo,opts);
            else
                Q = helm3d.dirichlet.get_quadrature_correction(obj.surf, ...
                    obj.eps_laphelm,obj.zk,[1.0,0],targinfo,opts);
            end
            
            obj.quad_opts_xs_laphelm = [];
            obj.quad_opts_xs_laphelm.format = format;
            obj.quad_opts_xs_laphelm.precomp_quadrature = Q;
        end

        function wfunc = solveforW(obj,varargin)
            %SOLVEFORW Helper function for GMRES when solving for sigma
            %          and alpha
            %   Optional argument:
            %     time [boolean, default false]: report GMRES timings
            if nargin > 1
                time = varargin{1};
            else 
                time = false;
            end
            nB0 = dot(obj.vn,obj.B0);
            b = surfacefun_to_array(nB0,obj.dom,obj.surf);
            if time
                t1 = tic;
            end
            [W,~,~,iter] = gmres(@(s) TaylorState.gmresA(s,obj.dom,obj.surf, ...
                obj.zk,obj.eps_taylor,obj.eps_laphelm, ...
                obj.quad_opts_taylor,obj.quad_opts_laphelm),b,[],obj.eps_gmres,50);
            if time
                t2 = toc(t1);
                fprintf('\tGMRES for A11*D = A12: %f s / %d iter. = %f s\n', ...
                    t2, iter(2), t2/iter(2))
            end
            wfunc = array_to_surfacefun(W,obj.dom,obj.surf);
        end

        function obj = compute_sigma_alpha(obj,varargin)
            %COMPUTE_SIGMA_ALPHA Solve integral equation representation 
            %                    for sigma and alpha, overriding the
            %                    TaylorState version
            %   Optional argument:
            %     time [boolean, default false]: report GMRES timings
            if nargin > 1
                time = varargin{1};
            else
                time = false;
            end
            dfunc = solveforD(obj,time);
            wfunc = solveforW(obj,time);
            fluxsigmaD = RefTaylorState.mtxfluxsigmanontaylor(obj.surf, ...
                obj.dom,obj.xs_nodes,obj.xs_weights,dfunc,obj.zk,obj.eps_taylor, ...
                obj.eps_laphelm,obj.surf,obj.quad_opts_xs_taylor, ...
                obj.quad_opts_xs_laphelm);
            fluxsigmaW = RefTaylorState.mtxfluxsigmanontaylor(obj.surf, ...
                obj.dom,obj.xs_nodes,obj.xs_weights,wfunc,obj.zk,obj.eps_taylor, ...
                obj.eps_laphelm,obj.surf,obj.quad_opts_xs_taylor, ...
                obj.quad_opts_xs_laphelm);
            fluxalpha = RefTaylorState.mtxfluxalphanontaylor(obj.surf, ...
                obj.dom,obj.xs_nodes,obj.xs_weights,obj.mH,obj.zk,obj.eps_taylor, ...
                obj.eps_laphelm,obj.surf,obj.quad_opts_xs_taylor, ...
                obj.quad_opts_xs_laphelm);
            obj.alpha = -1i*(obj.flux - fluxsigmaW)...
                /(-fluxsigmaD + fluxalpha);
            obj.sigma = 1i*obj.alpha.*dfunc - wfunc;
        end

        function obj = solve(obj,varargin)
            %SOLVE Solve the integral equation formulation of the 
            %      reference Taylor state BVP, overriding the TaylorState
            %      version
            % 
            %   Optional argument:
            %     time [boolean, default false]: print timings 
            if nargin > 1
                time = varargin{1};
            else
                time = false;
            end
            
            if time; t1 = tic; end
            obj = obj.compute_mH();
            if time
                t2 = toc(t1);
                fprintf('compute mH: %f s\n',t2)
                t1 = tic;
            end
            obj = obj.get_quad_corr_laphelm();
            if time
                t2 = toc(t1);
                fprintf('get Laplace/Helmholtz quad. corr.: %f s\n',t2)
                t1 = tic;
            end
            obj = obj.get_quad_corr_taylor();
            if time
                t2 = toc(t1);
                fprintf('get +taylor routine quad. corr.: %f s\n',t2)
                t1 = tic;
            end
            obj = obj.get_quad_corr_xs_laphelm();
            if time
                t2 = toc(t1);
                fprintf('get XS Laplace/Helmholtz quad. corr.: %f s\n',t2)
                t1 = tic;
            end
            obj = obj.get_quad_corr_taylor();
            if time
                t2 = toc(t1);
                fprintf('get XS +taylor routine quad. corr.: %f s\n',t2)
                t1 = tic;
            end
            obj = obj.compute_sigma_alpha(time);
            if time
                t2 = toc(t1);
                fprintf('compute sigma and alpha: %f s\n',t2)
                t1 = tic;
            end
            obj = obj.compute_m();
            if time; t2 = toc(t1); fprintf('compute m: %f s\n',t2); end
        end
    end

    methods (Static)
        fluxsigma = mtxfluxsigmanontaylor(S,dom,nodes,weights,sigma, ...
            zk,epstaylor,epslh,varargin)
        fluxalpha = mtxfluxalphanontaylor(S,dom,nodes,weights,mH,zk, ...
            epstaylor,epslh,varargin)
    end
end