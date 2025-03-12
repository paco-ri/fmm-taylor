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
        domain % Domain object
        zk % Beltrami parameter 
        flux % cross-sectional flux
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
            
            obj.domain = Domain(domain,domparams);

            if isnumeric(zk)
                obj.zk = complex(zk);
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Third argument should be the Beltrami parameter.'])
            end

            if isnumeric(flux)
                if obj.domain.nsurfaces == length(flux)
                    obj.flux = flux;
                end
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Fourth argument should be the cross-sectional ' ...
                    'flux(es). There should be as many fluxes as ' ...
                    'there are surfaces. '])
            end
            
            if isnumeric(tols)
                if isscalar(tols)
                    obj.eps_gmres = tols;
                    obj.eps_taylor = tols;
                    obj.eps_laphelm = tols;
                else
                    if length(tols)==3
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
            
            nsurf = obj.domain.nsurfaces;
            obj.quad_opts_taylor = cell(nsurf);
            for i = 1:nsurf
                for j = 1:nsurf
                    if abs(obj.zk) < eps
                        Q = taylor.static.get_quadrature_correction(obj.domain.surf{i}, ...
                            obj.eps_taylor,obj.domain.surf{j},opts);
                    else
                        Q = taylor.dynamic.get_quadrature_correction(obj.domain.surf{i}, ...
                            obj.zk,obj.eps_taylor,obj.domain.surf{j},opts);
                    end
                    obj.quad_opts_taylor{i,j} = [];
                    obj.quad_opts_taylor{i,j}.format = format;
                    obj.quad_opts_taylor{i,j}.precomp_quadrature = Q;                    
                end
            end
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
            
            nsurf = obj.domain.nsurfaces;
            obj.quad_opts_laphelm = cell(nsurf);
            for i = 1:nsurf
                for j = 1:nsurf
                    if abs(obj.zk) < eps
                        Q = lap3d.dirichlet.get_quadrature_correction(obj.domain.surf{i}, ...
                            obj.eps_laphelm,[1.0,0],obj.domain.surf{j},opts);
                    else
                        Q = helm3d.dirichlet.get_quadrature_correction(obj.domain.surf{i}, ...
                            obj.eps_laphelm,obj.zk,[1.0,0],obj.domain.surf{j},opts);
                    end
                    obj.quad_opts_laphelm{i,j} = [];
                    obj.quad_opts_laphelm{i,j}.format = format;
                    obj.quad_opts_laphelm{i,j}.precomp_quadrature = Q;
                end
            end
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
            
            nsurf = obj.domain.nsurfaces;
            Balpha = TaylorState.mtxBalpha(obj.domain,obj.zk, ...
                obj.eps_taylor,obj.eps_laphelm,obj.domain.surf, ...
                obj.quad_opts_taylor,obj.quad_opts_laphelm);
            b = zeros(nsurf*obj.domain.nptspersurf,nsurf);
            for i = 1:nsurf
                for j = 1:nsurf
                    inds = obj.domain.nptspersurf*(i-1)+1:obj.domain.nptspersurf*i;
                    b(inds,j) = surfacefun_to_array(Balpha{i,j},...
                        obj.domain.dom{i},obj.domain.surf{i});
                end
            end
            if time
                t1 = tic;
            end
            D = zeros(nsurf*obj.domain.nptspersurf,nsurf);
            totiter = 0;
            for i = 1:nsurf
                [D(:,i),~,~,iter,resvec] = gmres(@(s) TaylorState.gmresA(s, ...
                    obj.domain,obj.zk,obj.eps_taylor,obj.eps_laphelm, ...
                    obj.quad_opts_taylor,obj.quad_opts_laphelm), b(:,i), [], ...
                    obj.eps_gmres, 150);
                totiter = totiter+iter(2);
            end
            if time
                t2 = toc(t1);
                fprintf(['\t(GMRES for A11*D = A12 (part of "compute ' ...
                    'sigma and alpha"): %f s / %d iter. = %f s)\n'], ...
                    t2, totiter, t2/totiter)
            end

            dfunc = cell(nsurf);
            for i = 1:nsurf
                inds = obj.domain.nptspersurf*(i-1)+1:obj.domain.nptspersurf*i;
                for j = 1:nsurf
                    dfunc{i,j} = array_to_surfacefun(D(inds,j),obj.domain.dom{i},obj.domain.surf{i});
                end
            end
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

            if obj.domain.nsurfaces == 1
                fluxsigmaD = TaylorState.mtxfluxsigma(obj.domain,dfunc{1}, ...
                    obj.zk,obj.eps_taylor,obj.eps_laphelm, ...
                    obj.domain.surf{1},obj.quad_opts_taylor{1}, ...
                    obj.quad_opts_laphelm{1});
                fluxalpha = TaylorState.mtxfluxalpha( ...
                    obj.domain,obj.zk,obj.eps_taylor,obj.eps_laphelm, ...
                    obj.domain.surf{1},obj.quad_opts_taylor{1}, ...
                    obj.quad_opts_laphelm{1});
            else
                fluxsigmaD = TaylorState.mtxfluxsigma(obj.domain,dfunc,obj.zk, ...
                    obj.eps_taylor,obj.eps_laphelm,obj.domain.surf, ...
                    obj.quad_opts_taylor,obj.quad_opts_laphelm);
                fluxalpha = TaylorState.mtxfluxalpha(obj.domain,obj.zk, ...
                    obj.eps_taylor,obj.eps_laphelm,obj.domain.surf, ...
                    obj.quad_opts_taylor,obj.quad_opts_laphelm);
            end

            if obj.domain.nsurfaces == 1
                obj.alpha = obj.flux/(1i*fluxsigmaD + fluxalpha);
                obj.sigma{1} = 1i*obj.alpha.*dfunc{1};
            else
                A = 1i*fluxsigmaD + fluxalpha;
                obj.alpha = A\obj.flux.';
                obj.sigma = cell(1,2);
                obj.sigma{1} = 1i*(obj.alpha(1).*dfunc{1,1} + obj.alpha(2).*dfunc{1,2}); 
                obj.sigma{2} = 1i*(obj.alpha(1).*dfunc{2,1} + obj.alpha(2).*dfunc{2,2});
            end
        end

        function obj = compute_m0(obj,varargin)
            %COMPUTE_M0 Compute the sigma-dependent term of the vector
            %           density m
            %   Optional argument:
            %     np [int] additional polynomial order for oversampling

            % oversample = false;
            % if nargin > 1
            %     oversample = true;
            %     np = varargin{2};
            % end
            
            % deal with oversampling later 
            
            nsurf = obj.domain.nsurfaces;
            obj.m0 = cell(1,nsurf);
            % for j = 1:obj.nsurfaces
            %     if oversample
            %         obj.m0{j} = TaylorState.debyem0(obj.sigma{j},obj.zk,np);
            %     else
            %         obj.m0{j} = TaylorState.debyem0(obj.sigma{j},obj.zk);
            %     end
            % end
            for j = 1:nsurf
                obj.m0{j} = TaylorState.debyem0(obj.sigma{j},obj.zk, ...
                    obj.domain.L{j},obj.domain.vn{j});
            end
        end

        function obj = compute_m(obj)
            %COMPUTE_M Compute the vector density m
            obj = obj.compute_m0();
            obj.m = cell(1,2);
            for i = 1:obj.domain.nsurfaces
                obj.m{i} = obj.m0{i} + obj.alpha(i).*obj.domain.mH{i};
            end
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

            nsurf = obj.domain.nsurfaces;
            B = cell(1,nsurf);
            for j = 1:nsurf
                sigmavals = surfacefun_to_array(obj.sigma{j}, ...
                    obj.domain.dom{j},obj.domain.surf{j});
                mvals = surfacefun_to_array(obj.m{j},obj.domain.dom{j}, ...
                    obj.domain.surf{j});
                if abs(obj.zk) < eps
                    gradSksigma = taylor.static.eval_gradS0( ...
                        obj.domain.surf{j},sigmavals.',obj.eps_taylor, ...
                        obj.domain.surf{j},obj.quad_opts_taylor);
                    gradSksigma = array_to_surfacefun(gradSksigma.', ...
                        obj.domain.dom{j},obj.domain.surf{j});
                    curlSkm = taylor.static.eval_curlS0(obj.domain.surf{j}, ...
                        mvals.',obj.eps_taylor,obj.domain.surf{j}, ...
                        obj.quad_opts_taylor);
                    curlSkm = array_to_surfacefun(curlSkm.',obj.domain.dom{j}, ...
                        obj.domain.surf{j});
    
                    B{j} = -obj.sigma{j}.*obj.domain.vn{j}./2 + obj.m{j}./2 ...
                        - gradSksigma + 1i.*curlSkm;

                    if nsurf == 2
                        % need to add contribution of other surface
                        sigmavalsother = surfacefun_to_array(obj.sigma{3-j}, ...
                            obj.domain.dom{3-j},obj.domain.surf{3-j});
                        mvalsother = surfacefun_to_array(obj.m{3-j},obj.domain.dom{3-j}, ...
                            obj.domain.surf{3-j});
                        gradSksigmaother = taylor.static.eval_gradS0( ...
                            obj.domain.surf{3-j},sigmavalsother.',obj.eps_taylor, ...
                            obj.domain.surf{j});
                        gradSksigmaother = array_to_surfacefun( ...
                            gradSksigmaother.',obj.domain.dom{j},obj.domain.surf{j});
                        curlSkmother = taylor.static.eval_curlS0( ...
                            obj.domain.surf{3-j},mvalsother.',obj.eps_taylor, ...
                            obj.domain.surf{j});
                        curlSkmother = array_to_surfacefun( ...
                            curlSkmother.',obj.domain.dom{j},obj.domain.surf{j});
    
                        B{j} = B{j} - gradSksigmaother + 1i.*curlSkmother;
                    end
                else
                    gradSksigma = taylor.dynamic.eval_gradSk( ...
                        obj.domain.surf{j},obj.zk,sigmavals.',obj.eps_taylor, ...
                        obj.domain.surf{j},obj.quad_opts_taylor);
                    gradSksigma = array_to_surfacefun(gradSksigma.', ...
                        obj.domain.dom{j},obj.domain.surf{j});
                    curlSkm = taylor.dynamic.eval_curlSk(obj.domain.surf{j}, ...
                        obj.zk,mvals.',obj.eps_taylor,obj.domain.surf{j}, ...
                        obj.quad_opts_taylor);
                    curlSkm = array_to_surfacefun(curlSkm.',obj.domain.dom{j}, ...
                        obj.domain.surf{j});
                
                    dpars = [1.0,0];
                    Skm1 = helm3d.dirichlet.eval(obj.domain.surf{j}, ...
                        mvals(:,1),obj.domain.surf{j},obj.eps_laphelm,obj.zk, ...
                        dpars,obj.quad_opts_laphelm);
                    Skm2 = helm3d.dirichlet.eval(obj.domain.surf{j}, ...
                        mvals(:,2),obj.domain.surf{j},obj.eps_laphelm,obj.zk, ...
                        dpars,obj.quad_opts_laphelm);
                    Skm3 = helm3d.dirichlet.eval(obj.domain.surf{j}, ...
                        mvals(:,3),obj.domain.surf{j},obj.eps_laphelm,obj.zk, ...
                        dpars,obj.quad_opts_laphelm);
                    Skm = [Skm1 Skm2 Skm3];
                    Skm = array_to_surfacefun(Skm,obj.domain.dom{j},obj.domain.surf{j});
                
                    B{j} = -obj.sigma{j}.*obj.domain.vn{j}./2 + obj.m{j}./2 ... 
                        + 1i.*obj.zk.*Skm - gradSksigma + 1i.*curlSkm;

                    if nsurf == 2
                        % need to add contribution of other surface
                        sigmavalsother = surfacefun_to_array(obj.sigma{3-j}, ...
                            obj.domain.dom{3-j},obj.domain.surf{3-j});
                        mvalsother = surfacefun_to_array(obj.m{3-j},obj.domain.dom{3-j}, ...
                            obj.domain.surf{3-j});
                        Skmother = taylor.helper.helm_dir_vec_eval( ...
                            obj.domain.surf{3-j},mvalsother.',obj.domain.surf{j}, ...
                            obj.eps_laphelm,obj.zk,dpars);
                        Skmother = array_to_surfacefun(Skmother.', ...
                            obj.domain.dom{j},obj.domain.surf{j});
                        gradSksigmaother = taylor.dynamic.eval_gradSk( ...
                            obj.domain.surf{3-j},obj.zk,sigmavalsother.', ...
                            obj.eps_taylor,obj.domain.surf{j});
                        gradSksigmaother = array_to_surfacefun( ...
                            gradSksigmaother.',obj.domain.dom{j},obj.domain.surf{j});
                        curlSkmother = taylor.dynamic.eval_curlSk( ...
                            obj.domain.surf{3-j},obj.zk,mvalsother.', ...
                            obj.eps_taylor,obj.domain.surf{j});
                        curlSkmother = array_to_surfacefun( ...
                            curlSkmother.',obj.domain.dom{j},obj.domain.surf{j});
    
                        B{j} = B{j} + 1i.*obj.zk.*Skmother ...
                            - gradSksigmaother + 1i.*curlSkmother;
                    end
                end
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
            
            B = 0;
            for j = 1:obj.domain.nsurfaces
                opts_int = [];
                opts_int.format = 'rsc';
                if abs(obj.zk) < eps
                    Q = taylor.static.get_quadrature_correction( ...
                        obj.domain.surf{j},obj.eps_taylor,targinfo,opts_int);
                else
                    Q = taylor.dynamic.get_quadrature_correction( ...
                        obj.domain.surf{j},obj.zk,obj.eps_taylor,targinfo, ...
                        opts_int);
                end
                opts_int.precomp_quadrature = Q;

                sigmavals = surfacefun_to_array(obj.sigma{j}, ...
                    obj.domain.dom{j},obj.domain.surf{j});
                mvals = surfacefun_to_array(obj.m{j},obj.domain.dom{j}, ...
                    obj.domain.surf{j});
                if abs(obj.zk) < eps
                    gradSksigma = taylor.static.eval_gradS0( ...
                        obj.domain.surf{j},sigmavals.',obj.eps_taylor, ...
                        targinfo,opts_int);
                    curlSkm = taylor.static.eval_curlS0(obj.domain.surf{j}, ...
                        mvals.',obj.eps_taylor,targinfo,opts_int);
                
                    B = B - gradSksigma + 1i.*curlSkm;
                else
                    gradSksigma = taylor.dynamic.eval_gradSk( ...
                        obj.domain.surf{j},obj.zk,sigmavals.',obj.eps_taylor, ...
                        targinfo,opts_int);
                    curlSkm = taylor.dynamic.eval_curlSk(obj.domain.surf{j}, ...
                        obj.zk,mvals.',obj.eps_taylor,targinfo,opts_int);
                
                    dpars = [1.0, 0.0];
                    opts_int_lh = [];
                    opts_int_lh.format = 'rsc';
                    Qlh = helm3d.dirichlet.get_quadrature_correction( ...
                        obj.domain.surf{j},obj.eps_laphelm,obj.zk,dpars, ...
                        targinfo,opts_int_lh);
                    opts_int_lh.precomp_quadrature = Qlh;
    
                    Skm1 = helm3d.dirichlet.eval(obj.domain.surf{j}, ...
                        mvals(:,1),targinfo,obj.eps_laphelm,obj.zk, ...
                        dpars,opts_int_lh);
                    Skm2 = helm3d.dirichlet.eval(obj.domain.surf{j}, ...
                        mvals(:,2),targinfo,obj.eps_laphelm,obj.zk, ...
                        dpars,opts_int_lh);
                    Skm3 = helm3d.dirichlet.eval(obj.domain.surf{j}, ...
                        mvals(:,3),targinfo,obj.eps_laphelm,obj.zk, ...
                        dpars,opts_int_lh);
                    Skm = [Skm1 Skm2 Skm3];
                    Skm = Skm.';
                
                    B = B + 1i.*obj.zk.*Skm - gradSksigma + 1i.*curlSkm;
                end
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
        Bsigma = mtxBsigma(domain,zk,epstaylor,epslh,varargin)
        A = gmresA(s,domain,zk,epstaylor,epslh,opts,optslh)
        Balpha = mtxBalpha(domain,zk,epstaylor,epslh,varargin)
        fluxsigma = mtxfluxsigma(domain,zk,epstaylor,epslh,varargin)
        fluxalpha = mtxfluxalpha(domain,epstaylor,epslh,varargin)
        fluxsigmanontaylor = mtxfluxsigmanontaylor(domain,nodes,weights, ...
            sigma,zk,epstaylor,epslh,varargin)
        fluxalphanontaylor = mtxfluxalphanontaylor(domain,nodes,weights, ...
            mH,zk,epstaylor,epslh,varargin)
        m0 = debyem0(sigma,lambda,L,vn)
        integrala = intacyc(f,n,nv)
        integralb = intbcyc(f,n,nu)
        [u, v, w, curlfree, divfree] = hodge_inward(f)
        [qnodes, qweights] = toroidalfluxquadaxi(nr,nt,ro,ao,ri,ai)
        [qnodes, qweights] = poloidalfluxquadaxi(nr,np,ro,ao,ri,ai)
    end
end