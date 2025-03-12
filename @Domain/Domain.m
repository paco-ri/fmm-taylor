classdef Domain
    % A ``Domain`` object represents the surface :math:`\Gamma`, which has 
    % one or two connected components, each of which is a toroidal surface.
    % This class is essentially a wrapper for classes that represent 
    % surfaces in ``fmm3dbie`` and  ``surfacefun``, along with some 
    % additional information.
    % 
    % :param domain: the surface as a ``surfacefun.@surfacemesh`` or a cell 
    %   array of these
    % :type domain: :class:`surfacefun.surfacemesh`
    % :param domparams: parameters describing the surface: number of 
    %   Chebyshev nodes in each dimension on each quadrilateral patch, 
    %   number of patches in toroidal direction, number of patches in 
    %   poloidal direction
    % :type domparams: list of ints

    properties
        nsurfaces % number of nested toroidal surfaces (1 or 2)
        dom % surface as a :class``surfacefun.surfacemesh``
        surf % surface as an :class``fmm3dbie.surfer``

        % parameters describing the surface: number of 
        %   Chebyshev nodes in each dimension on each quadrilateral patch, 
        %   number of patches in toroidal direction, number of patches in 
        %   poloidal direction
        domparams

        nptspersurf % number of points on each surface
        vn % outward unit normal vector on surface 
        mH % surface harmonic vector field on surface 
        L % surfaceop for Laplace-Beltrami operator on surf
        aquad % A-cycle quadrature
        bquad % B-cycle quadrature
    end

    methods
        function obj = Domain(domain,domparams)
            % Construct an instance of a Domain (see parameters above)
            if isa(domain, 'surfacemesh')
                obj.nsurfaces = 1;
                obj.dom = {domain};
                surf = surfer.surfacemesh_to_surfer(domain);
                obj.surf = {surf};
                obj.vn = {normal(domain)};
            elseif isscalar(domain)
                if isa(domain{1}, 'surfacemesh')
                    obj.nsurfaces = 1;
                    obj.dom = domain;
                    surf = surfer.surfacemesh_to_surfer(domain{1});
                    obj.surf = {surf};
                    obj.vn = {normal(domain{1})};
                end
            elseif length(domain) == 2
                if isa(domain{1}, 'surfacemesh') && isa(domain{2}, 'surfacemesh')
                    obj.nsurfaces = 2;
                    obj.dom = domain;
                    surf = cell(1,2);
                    surf{1} = surfer.surfacemesh_to_surfer(domain{1});
                    surf{2} = surfer.surfacemesh_to_surfer(domain{2});
                    obj.surf = surf;
                    vn = cell(1,2);
                    vn{1} = normal(obj.dom{1});
                    vn{2} = -normal(obj.dom{2});
                    obj.vn = vn;
                end
            else
                error(['Invalid call to Domain constructor. ' ...
                    'First argument should be a surfacemesh or a cell ' ...
                    'array with two surfacemeshes.'])
            end

            if isnumeric(domparams)
                obj.domparams = domparams;
                obj.nptspersurf = obj.domparams(1)^2*obj.domparams(2)...
                *obj.domparams(3);
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Second argument should be an array of three ' ...
                    'integers.'])
            end

            obj.L = cell(1,obj.nsurfaces);
            pdo = [];
            pdo.lap = 1;
            for i = 1:obj.nsurfaces
                obj.L{i} = surfaceop(obj.dom{i}, pdo);
                obj.L{i}.rankdef = true;
                obj.L{i}.build();
            end

            obj = obj.compute_mH();
            
            obj.aquad = cell(1,obj.nsurfaces);
            if obj.nsurfaces == 1
                [x,xv,w] = Domain.acycquad(obj.dom{1},domparams);
                obj.aquad{1} = [];
                obj.aquad{1}.x = x;
                obj.aquad{1}.xv = xv;
                obj.aquad{1}.w = w;
            else
                obj.bquad = cell(1,2);
                for i = 1:2
                    [x,xv,w] = Domain.acycquad(domain{i},domparams);
                    obj.aquad{i} = [];
                    obj.aquad{i}.x = x;
                    obj.aquad{i}.xv = xv;
                    obj.aquad{i}.w = w;
                    [x,xu,w] = Domain.bcycquad(domain{i},domparams);
                    obj.bquad{i} = [];
                    obj.bquad{i}.x = x;
                    obj.bquad{i}.xu = xu;
                    obj.bquad{i}.w = w;
                end
            end
           
        end

        function obj = compute_mH(obj)
            %COMPUTE_MH Compute surface harmonic vector field
            sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
            cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
            obj.mH = cell(1,obj.nsurfaces);
            for i = 1:obj.nsurfaces
                phihat = surfacefunv(@(x,y,z) -sinphi(x,y,z), ...
                     @(x,y,z) cosphi(x,y,z), ...
                     @(x,y,z) 0.*z, obj.dom{i});
                dummy = cross(obj.vn{i}, phihat);
                    
                if i == 1
                    [~, ~, vH] = hodge(dummy);
                else
                    [~, ~, vH] = Domain.hodge_inward(dummy);
                end
                obj.mH{i} = vH + 1i.*cross(obj.vn{i},vH);
            end
        end
    end

    methods (Static)
        [x,xv,w] = acycquad(dom,domparams);
        [x,xu,w] = bcycquad(dom,domparams);
        [u, v, w, curlfree, divfree] = hodge_inward(f)
    end
end