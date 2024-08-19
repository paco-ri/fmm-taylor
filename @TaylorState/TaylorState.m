classdef TaylorState
    %TAYLORSTATE A class describing a Taylor state with parameter k
    %   A Taylor state in a domain Omega bounded by m toroidal surfaces is 
    %   a solution to the BVP
    % 
    %          curl B = k B     in Omega
    %           n . B = 0       on the boundary of Omega
    %     \oint_{S_i} = Phi_i   i = 1, ..., m
    % 
    %   where n is the outward unit normal to Omega, S_i is a topologically 
    %   distinct cross-sectional surface through Omega, and Phi_i is a 
    %   prescribed flux.

    properties
        dom % surface as a surfacefun.surfacemesh
        surf % surface as an fmm3dbie.surfer
        domparams % parameters describing surface
        zk % Beltrami parameter 
        flux % cross-sectional flux
        vn % outward unit normal vector on surface 
        mH % surface harmonic vector field on surface 
        quad_opts_taylor % near quadrature correction for grad and curl Sk
        quad_opts_laphelm % near quad. corr. for Sk
        sigma % scalar density in integral equation solution
        alpha % mH coefficient in integral equation solution 
        m0 % sigma-dependent part of vector density m 
        m % vector density in integral equation solution 
        eps_gmres % GMRES tolerance 
        eps_taylor % tolerance for grad and curl Sk 
        eps_laphelm % tolerance for Sk
    end

    methods
        function obj = TaylorState(domain,domparams,zk,flux,tols)
            %TAYLORSTATE Construct an instance of a TaylorState
            %   Arguments:
            %     domain [surfacemesh] toroidal boundary surface
            %     domparams [int(3)] parameters describing surface
            %     zk [double complex] Beltrami parameter
            %     flux [double] cross-sectional flux
            %     tols [double or double(3)] quad. and GMRES tolerances

            if isa(domain, 'surfacemesh')
                obj.dom = domain;
                surf = surfer.surfacemesh_to_surfer(domain);
                obj.surf = surf;
                obj.vn = normal(domain);
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'First argument should be a surfacemesh.'])
            end

            if isnumeric(domparams)
                obj.domparams = domparams;
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Second argument should be an array of three ' ...
                    'integers.'])
            end

            if isnumeric(zk)
                obj.zk = complex(zk);
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Third argument should be the Beltrami parameter.'])
            end

            if isnumeric(flux)
                obj.flux = flux;
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Fourth argument should be the cross-sectional flux.'])
            end
            
            if isnumeric(tols)
                if isscalar(tols)
                    obj.eps_gmres = tols;
                    obj.eps_taylor = tols;
                    obj.eps_laphelm = tols;
                else
                    [tolsr, tolsc] = size(tols);
                    if tolsr == 3 && tolsc == 1 || tolsr == 1 && tolsc == 3
                        obj.eps_gmres = tols(1);
                        obj.eps_taylor = tols(2);
                        obj.eps_laphelm = tols(3);
                    else
                    error(['Invalid call to TaylorState constructor. ' ...
                        'Fifth argument should be 1 or 3 tolerance(s).'])
                    end
                end
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Fifth argument should be 1 or 3 tolerance(s).'])
            end
        end

        function obj = compute_mH(obj,varargin)
            %COMPUTE_MH Compute surface harmonic vector field 
            %   Argument (optional)
            %     dummy [surfacefunv] arbitrary surface vector field 
            if nargin > 1 && isa(varargin{2},'surfacefunv')
                dummy = varargin{1};
            else
                sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
                cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
                phihat = surfacefunv(@(x,y,z) -sinphi(x,y,z), ...
                     @(x,y,z) cosphi(x,y,z), ...
                     @(x,y,z) 0.*z, obj.dom);
                dummy = cross(obj.vn, phihat); 
            end
                
            [~, ~, vH] = hodge(dummy);
            obj.mH = vH + 1i.*cross(obj.vn,vH);
        end

        function obj = get_quad_corr_taylor(obj,varargin)
            %GET_QUAD_CORR_TAYLOR Get near quad. corr. for grad / curl Sk 
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
            
            if abs(obj.zk) < eps
                Q = taylor.static.get_quadrature_correction(obj.surf, ...
                    obj.eps_taylor,obj.surf,opts);
            else
                Q = taylor.dynamic.get_quadrature_correction(obj.surf, ...
                    obj.zk,obj.eps_taylor,obj.surf,opts);
            end
            
            obj.quad_opts_taylor = [];
            obj.quad_opts_taylor.format = format;
            obj.quad_opts_taylor.precomp_quadrature = Q;
        end

        function obj = get_quad_corr_laphelm(obj,varargin)
            %GET_QUAD_CORR_LAPHELM Get near quad. corr. for Sk 
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
            
            if abs(obj.zk) < eps
                Q = lap3d.dirichlet.get_quadrature_correction(obj.surf, ...
                    obj.eps_laphelm,[1.0,0],obj.surf,opts);
            else
                Q = helm3d.dirichlet.get_quadrature_correction(obj.surf, ...
                    obj.eps_laphelm,obj.zk,[1.0,0],obj.surf,opts);
            end
            
            obj.quad_opts_laphelm = [];
            obj.quad_opts_laphelm.format = format;
            obj.quad_opts_laphelm.precomp_quadrature = Q;
        end

        function dfunc = solveforD(obj,varargin)
            %SOLVEFORD Helper function for GMRES when solving for sigma
            %          and alpha
            %   Optional argument:
            %     time [boolean, default false]: report GMRES timings
            if nargin > 1
                time = varargin{1};
            else 
                time = false;
            end
            Balpha = TaylorState.mtxBalpha(obj.surf,obj.dom,obj.mH,obj.zk, ...
                obj.eps_taylor,obj.eps_laphelm,obj.surf, ...
                obj.quad_opts_taylor,obj.quad_opts_laphelm);
            b = surfacefun_to_array(Balpha,obj.dom,obj.surf);
            if time
                t1 = tic;
            end
            [D,~,~,iter] = gmres(@(s) TaylorState.gmresA(s,obj.dom,obj.surf, ...
                obj.zk,obj.eps_taylor,obj.eps_laphelm, ...
                obj.quad_opts_taylor,obj.quad_opts_laphelm),b,[],obj.eps_gmres,50);
            if time
                t2 = toc(t1);
                fprintf('\tGMRES for A11*D = A12: %f s / %d iter. = %f s\n', ...
                    t2, iter(2), t2/iter(2))
            end
            dfunc = array_to_surfacefun(D,obj.dom,obj.surf);
        end

        function obj = compute_sigma_alpha(obj,varargin)
            %COMPUTE_SIGMA_ALPHA Solve integral equation representation 
            %                    for sigma and alpha
            %   Optional argument:
            %     time [boolean, default false]: report GMRES timings
            if nargin > 1
                time = varargin{1};
            else
                time = false;
            end
            dfunc = solveforD(obj,time);
            fluxsigmaD = TaylorState.mtxfluxsigma(obj.surf,obj.dom, ...
                obj.domparams,dfunc,obj.zk,obj.eps_taylor, ...
                obj.eps_laphelm,obj.surf,obj.quad_opts_taylor, ...
                obj.quad_opts_laphelm);
            fluxalpha = TaylorState.mtxfluxalpha(obj.surf,obj.dom, ...
                obj.domparams,obj.mH,obj.zk,obj.eps_taylor, ...
                obj.eps_laphelm,obj.surf,obj.quad_opts_taylor, ...
                obj.quad_opts_laphelm);
            obj.alpha = obj.flux/(1i*fluxsigmaD + fluxalpha);
            obj.sigma = 1i*obj.alpha.*dfunc;
        end

        function obj = compute_m0(obj,varargin)
            %COMPUTE_M0 Compute the sigma-dependent term of the vector
            %           density m
            %   Optional argument:
            %     np [int] additional polynomial order for oversampling
            oversample = false;
            if nargin > 1
                oversample = true;
                np = varargin{2};
            end

            pdo = [];
            pdo.lap = 1;

            % resample
            n = size(obj.sigma.vals{1,1},1);
            if oversample
                sigma1 = resample(obj.sigma,n+np);
            else
                sigma1 = obj.sigma;
            end
            dom1 = sigma1.domain;

            % solve lap(u) = sigma
            L = surfaceop(dom1, pdo, sigma1);
            L.rankdef = true;
            u = L.solve();

            if oversample
                vn1 = normal(dom1);
            else
                vn1 = obj.vn;
            end
            obj.m0 = 1i.*obj.zk.*(grad(u) + 1i.*cross(vn1, grad(u)));

            if oversample
                obj.m0 = resample(obj.m0,n);
            end
        end

        function obj = compute_m(obj)
            %COMPUTE_M Compute the vector density m
            obj = obj.compute_m0();
            obj.m = obj.m0 + obj.alpha.*obj.mH;
        end

        function obj = solve(obj,varargin)
            %SOLVE Solve the integral equation formulation of the Taylor
            %      state BVP
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
                fprintf('get +taylor route quad. corr.: %f s\n',t2)
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

        function B = surface_B(obj)
            %SURFACE_B Return the surface magnetic field as a surfacefunv 
            sigmavals = surfacefun_to_array(obj.sigma,obj.dom,obj.surf);
            mvals = surfacefun_to_array(obj.m,obj.dom,obj.surf);
            if abs(obj.zk) < eps
                gradSksigma = taylor.static.eval_gradS0(obj.surf, ...
                    sigmavals.',obj.eps_taylor,obj.surf, ...
                    obj.quad_opts_taylor);
                gradSksigma = array_to_surfacefun(gradSksigma.', ...
                    obj.dom,obj.surf);
                curlSkm = taylor.static.eval_curlS0(obj.surf,mvals.', ...
                    obj.eps_taylor,obj.surf,obj.quad_opts_taylor);
                curlSkm = array_to_surfacefun(curlSkm.',obj.dom,obj.surf);

                B = -obj.sigma.*obj.vn./2 + obj.m./2 - gradSksigma + ...
                    1i.*curlSkm;
            else
                gradSksigma = taylor.dynamic.eval_gradSk(obj.surf, ...
                    obj.zk,sigmavals.',obj.eps_taylor,obj.surf, ...
                    obj.quad_opts_taylor);
                gradSksigma = array_to_surfacefun(gradSksigma.', ...
                    obj.dom,obj.surf);
                curlSkm = taylor.dynamic.eval_curlSk(obj.surf,obj.zk, ...
                    mvals.',obj.eps_taylor,obj.surf,obj.quad_opts_taylor);
                curlSkm = array_to_surfacefun(curlSkm.',obj.dom,obj.surf);
            
                dpars = [1.0,0];
                Skm1 = helm3d.dirichlet.eval(obj.surf,mvals(:,1), ...
                    obj.surf,obj.eps_laphelm,obj.zk,dpars, ...
                    obj.quad_opts_laphelm);
                Skm2 = helm3d.dirichlet.eval(obj.surf,mvals(:,2), ...
                    obj.surf,obj.eps_laphelm,obj.zk,dpars, ...
                    obj.quad_opts_laphelm);
                Skm3 = helm3d.dirichlet.eval(obj.surf,mvals(:,3), ...
                    obj.surf,obj.eps_laphelm,obj.zk,dpars, ...
                    obj.quad_opts_laphelm);
                Skm = [Skm1 Skm2 Skm3];
                Skm = array_to_surfacefun(Skm,obj.dom,obj.surf);
            
                B = -obj.sigma.*obj.vn./2 + obj.m./2 + 1i.*obj.zk.*Skm ...
                    - gradSksigma + 1i.*curlSkm;
            end
        end

        function B = interior_B(obj,intpts)
            %INTERIOR_B Compute the magnetic field at points inside the
            %           domain
            %   Argument:
            %     intpt [double(3,*)] interior points
            if ~isnumeric(intpts)
                error(['Invalid call to interior_B(). Argument should' ...
                    'be an array of interior points. '])
            end
            [intptsm, intptsn] = size(intpts);
            if intptsm ~= 3 
                if intptsn == 3
                    intpts = intpts.';
                else
                    error(['Invalid call to interior_B(). intpts is' ...
                        'an invalid size. '])
                end
            end

            targinfo = [];
            targinfo.r = intpts;

            opts_int = [];
            opts_int.format = 'rsc';
            if abs(obj.zk) < eps
                Q = taylor.static.get_quadrature_correction(obj.surf, ...
                    obj.eps_taylor,targinfo,opts_int);
            else
                Q = taylor.dynamic.get_quadrature_correction(obj.surf, ...
                    obj.zk,obj.eps_taylor,targinfo,opts_int);
            end
            opts_int.precomp_quadrature = Q;

            sigmavals = surfacefun_to_array(obj.sigma,obj.dom,obj.surf);
            mvals = surfacefun_to_array(obj.m,obj.dom,obj.surf);
            if abs(obj.zk) < eps
                gradSksigma = taylor.static.eval_gradS0(obj.surf, ...
                    sigmavals.',obj.eps_taylor,targinfo,opts_int);
                curlSkm = taylor.static.eval_curlS0(obj.surf,mvals.', ...
                    obj.eps_taylor,targinfo,opts_int);
            
                B = -gradSksigma + 1i.*curlSkm;
            else
                gradSksigma = taylor.dynamic.eval_gradSk(obj.surf, ...
                    obj.zk,sigmavals.',obj.eps_taylor,targinfo,opts_int);
                curlSkm = taylor.dynamic.eval_curlSk(obj.surf, ...
                    obj.zk,mvals.',obj.eps_taylor,targinfo,opts_int);
            
                dpars = [1.0, 0.0];
                opts_int_lh = [];
                opts_int_lh.format = 'rsc';
                Qlh = helm3d.dirichlet.get_quadrature_correction( ...
                    obj.surf,obj.eps_laphelm,obj.zk,dpars,targinfo, ...
                    opts_int_lh);
                opts_int_lh.precomp_quadrature = Qlh;

                Skm1 = helm3d.dirichlet.eval(obj.surf,mvals(:,1), ...
                    targinfo,obj.eps_laphelm,obj.zk,dpars,opts_int_lh);
                Skm2 = helm3d.dirichlet.eval(obj.surf,mvals(:,2), ...
                    targinfo,obj.eps_laphelm,obj.zk,dpars,opts_int_lh);
                Skm3 = helm3d.dirichlet.eval(obj.surf,mvals(:,3), ...
                    targinfo,obj.eps_laphelm,obj.zk,dpars,opts_int_lh);
                Skm = [Skm1 Skm2 Skm3];
                Skm = Skm.';
            
                B = 1i.*obj.zk.*Skm - gradSksigma + 1i.*curlSkm;
            end
        end

        function [errB, curlB, kB] = fd_test(obj,intpt,h)
            %FD_TEST Run a quick test that curl B = k B is satisfied using
            %        finite differences
            %   Arguments:
            %     intpt [double(3)] point in interior of volume where 
            %                       condition is checked
            %     h [double] finite difference increment
            if ~isnumeric(intpt)
                error(['Invalid call to interior_B(). First argument ' ...
                    'should be an interior point. '])
            end
            [intptm, intptn] = size(intpt);
            if intptn == 1 
                if intptm == 3
                    intpt = intpt.';
                else
                    error(['Invalid call to interior_B(). intpt is ' ...
                        'an invalid size. '])
                end
            elseif intptn == 3
                if intptm ~= 1
                    error(['Invalid call to interior_B(). intpt is ' ...
                        'an invalid size. '])
                end
            else
                error(['Invalid call to interior_B(). intpt is ' ...
                    'an invalid size. '])
            end

            intpts = [intpt;
                intpt + [h 0 0];
                intpt - [h 0 0];
                intpt + [0 h 0];
                intpt - [0 h 0];
                intpt + [0 0 h];
                intpt - [0 0 h]];
            intpts = intpts.';
            
            Bint = obj.interior_B(intpts);
            Bint = Bint.';
            curlB = [Bint(4,3)-Bint(5,3)-Bint(6,2)+Bint(7,2); % DyBz-DzBy
                Bint(6,1)-Bint(7,1)-Bint(2,3)+Bint(3,3); % DzBx-DxBz
                Bint(2,2)-Bint(3,2)-Bint(4,1)+Bint(5,1)]; % DxBy-DyBx
            curlB = curlB./(2*h);
            kB = obj.zk.*Bint(1,:).';
            errB = norm(curlB-kB);
        end
    end

    methods (Static)
        Bsigma = mtxBsigma(S,dom,sigma,zk,epstaylor,epslh,varargin)
        A = gmresA(s,dom,S,zk,epstaylor,epslh,opts,optslh)
        Balpha = mtxBalpha(S,dom,mH,zk,epstaylor,epslh,varargin)
        fluxsigma = mtxfluxsigma(S,dom,domparams,sigma,zk,epstaylor,epslh,varargin)
        fluxalpha = mtxfluxalpha(S,dom,domparams,mH,zk,epstaylor,epslh,varargin)
    end
end