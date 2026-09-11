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

        nptspersurf % number of points on each surface (vector, one per surface)
        vn % outward unit normal vector on surface
        mH % surface harmonic vector field on surface
        L % surfaceop for Laplace-Beltrami operator on surf
        aquad % A-cycle quadrature
        bquad % B-cycle quadrature

        % Quadforest and patch-to-quadforest map for adaptively refined
        % meshes, one cell entry per surface. Empty (or all-zero Morton
        % codes) means the mesh is a uniform nu-by-nv grid.
        qf
        p2q
    end

    methods
        function obj = Domain(domain,domparams,varargin)
            % Construct an instance of a Domain (see parameters above)
            if isa(domain, 'Domain')
                obj = domain;
                return
            end
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
                obj.nptspersurf = zeros(1,obj.nsurfaces);
                for i = 1:obj.nsurfaces
                    obj.nptspersurf(i) = obj.domparams(1)^2 ...
                        *length(obj.dom{i}.x);
                end
            else
                error(['Invalid call to TaylorState constructor. ' ...
                    'Second argument should be an array of three ' ...
                    'integers.'])
            end

            % Optional quadforest / patch-to-quadforest map for adaptive
            % meshes: Domain(domain,domparams,qf,p2q).
            obj.qf = cell(1,obj.nsurfaces);
            obj.p2q = cell(1,obj.nsurfaces);
            if numel(varargin) >= 2 && ~isempty(varargin{1})
                qf_in = varargin{1};
                p2q_in = varargin{2};
                if ~iscell(qf_in), qf_in = {qf_in}; end
                if ~iscell(p2q_in), p2q_in = {p2q_in}; end
                if length(qf_in) ~= obj.nsurfaces ...
                        || length(p2q_in) ~= obj.nsurfaces
                    error(['Invalid call to Domain constructor. There ' ...
                        'should be as many quadforests and patch maps ' ...
                        'as there are surfaces.'])
                end
                for i = 1:obj.nsurfaces
                    if ~isempty(p2q_in{i}) ...
                            && size(p2q_in{i},1) ~= length(obj.dom{i}.x)
                        error(['Invalid call to Domain constructor. ' ...
                            'p2q for surface %d has %d rows but the ' ...
                            'mesh has %d patches.'], i, ...
                            size(p2q_in{i},1), length(obj.dom{i}.x))
                    end
                end
                obj.qf = qf_in;
                obj.p2q = p2q_in;
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
            
            % Compute A- and B-cycle quadrature for uniform surface meshes
            obj.aquad = cell(1,obj.nsurfaces);
            obj.bquad = cell(1,obj.nsurfaces);
            for i = 1:obj.nsurfaces
                if obj.is_adaptive(i)
                    continue
                end
                [x,xv,w] = Domain.acycquad(obj.dom{i},domparams);
                obj.aquad{i} = [];
                obj.aquad{i}.x = x;
                obj.aquad{i}.xv = xv;
                obj.aquad{i}.w = w;
                if obj.nsurfaces > 1
                    [x,xu,w] = Domain.bcycquad(obj.dom{i},domparams);
                    obj.bquad{i} = [];
                    obj.bquad{i}.x = x;
                    obj.bquad{i}.xu = xu;
                    obj.bquad{i}.w = w;
                end
            end
           
        end

        function tf = is_adaptive(obj,varargin)
            %IS_ADAPTIVE True if a surface carries a nontrivial quadforest
            %   tf = obj.is_adaptive()  -> true if ANY surface is adaptive
            %   tf = obj.is_adaptive(i) -> true if surface i is adaptive
            %
            %   A p2q whose Morton codes are all zero describes an
            %   unrefined mesh.
            if nargin > 1
                inds = varargin{1};
            else
                inds = 1:obj.nsurfaces;
            end
            tf = false;
            for i = inds
                if ~isempty(obj.p2q) && numel(obj.p2q) >= i ...
                        && ~isempty(obj.p2q{i}) ...
                        && any(obj.p2q{i}(:,3) ~= 0)
                    tf = true;
                    return
                end
            end
        end

        function off = blockoffsets(obj)
            %BLOCKOFFSETS Offsets into a vector stacking all surfaces
            %   Returns a vector of length nsurfaces+1 such that the block
            %   representing a function on surface i occupies indices 
            %   off(i)+1:off(i+1) (necessary because the surfaces may have
            %   different numbers of patches).
            off = [0 cumsum(obj.nptspersurf(:).')];
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
                % obj.mH{i} = vH + 1i.*cross(obj.vn{i},vH);
                obj.mH{i} = vH + times(1i, cross(obj.vn{i}, vH));
            end
        end
    end

    methods (Static)
        [x,xv,w] = acycquad(dom,domparams);
        [x,xu,w] = bcycquad(dom,domparams);
        [u, v, w, curlfree, divfree] = hodge_inward(f)
    end
end