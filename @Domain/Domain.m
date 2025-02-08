classdef Domain
    %DOMAIN A class describing a domain in which Taylor states are solved
    %   This class is a wrapper class for the surfacemesh and surfer 
    %   classes 

    properties
        nsurfaces % number of nested toroidal surfaces (1 or 2)
        dom % surface as a surfacefun.surfacemesh
        nfine % high polynomial order for upsampled domain
        domfine % dom with polynomial order nfine
        surf % surface as an fmm3dbie.surfer
        domparams % parameters describing surface
        nptspersurf % number of points on each surface
        vn % outward unit normal vector on surface 
        vnfine % vn on domfine
        mH % surface harmonic vector field on surface 
        L % surfaceop for Laplace-Beltrami operator on surf
        aquad % A-cycle quadrature
        bquad % B-cycle quadrature
    end

    methods
        function obj = Domain(domain,domparams,nfine)
            %UNTITLED3 Construct an instance of a Domain
            %   Arguments: see above
            if isa(domain, 'surfacemesh')
                obj.nsurfaces = 1;
                obj.dom = {domain};
                obj.domfine = {resample(domain,nfine)};
                surf = surfer.surfacemesh_to_surfer(domain);
                obj.surf = {surf};
                obj.vn = {normal(domain)};
                obj.vnfine = {normal(obj.domfine)};
            elseif isscalar(domain)
                if isa(domain{1}, 'surfacemesh')
                    obj.nsurfaces = 1;
                    obj.dom = domain;
                    obj.domfine = {resample(domain{1},nfine)};
                    surf = surfer.surfacemesh_to_surfer(domain{1});
                    obj.surf = {surf};
                    obj.vn = {normal(domain{1})};
                    obj.vn = {normal(obj.domfine{1})};
                end
            elseif length(domain) == 2
                if isa(domain{1}, 'surfacemesh') && isa(domain{2}, 'surfacemesh')
                    obj.nsurfaces = 2;
                    obj.dom = domain;
                    obj.domfine = {resample(domain{1},nfine), ...
                        resample(domain{2},nfine)};
                    surf = cell(1,2);
                    surf{1} = surfer.surfacemesh_to_surfer(domain{1});
                    surf{2} = surfer.surfacemesh_to_surfer(domain{2});
                    obj.surf = surf;
                    vn = cell(1,2);
                    vn{1} = normal(domain{1});
                    vn{2} = -1.*normal(domain{2}); % flip inner normal
                    vnfine = cell(1,2);
                    vnfine{1} = normal(obj.domfine{1});
                    vnfine{2} = -normal(obj.domfine{2});
                    obj.vn = vn;
                    obj.vnfine = vnfine;
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
                     @(x,y,z) 0.*z, obj.domfine{i}); % obj.dom{i});
                % dummy = cross(obj.vn{i}, phihat); 
                dummy = cross(obj.vnfine{i}, phihat);
                    
                if i == 1
                    [~, ~, vH] = hodge(dummy);
                else
                    [~, ~, vH] = TaylorState.hodge_inward(dummy);
                end
                % obj.mH{i} = vH + 1i.*cross(obj.vn{i},vH);
                mHfine = vH + 1i.*cross(obj.vnfine{i},vH);
                obj.mH{i} = surfacefunv(obj.dom{i});
                for j = 1:3
                    obj.mH{i}.components{j} = ...
                        resample(mHfine.components{j},obj.domparams(1));
                end
            end
        end
    end

    methods (Static)
        [x,xv,w] = acycquad(dom,domparams);
        [x,xu,w] = bcycquad(dom,domparams);
    end
end