%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef class_expRadon_matlab < handle
    properties (SetAccess = public, Hidden = false)
        mu;
        N;
    end
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods (Access = private)
        %% Constructor - Create a new C++ class instance
        function this = class_expRadon_matlab(varargin)
            this.objectHandle = class_interface_mex('new', varargin{:});
        end
    end
    methods (Static)
        function singleObj = getInstance(varargin)
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = class_expRadon_matlab(varargin{:});
            end
            singleObj = localObj;
        end
    end
    methods
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            class_interface_mex('delete', this.objectHandle);
        end
        function varargout = eq2us(this, varargin)
            [varargout{1:nargout}] = class_interface_mex('eq2us', this.objectHandle, varargin{:});
        end
        function varargout = us2eq(this, varargin)
            [varargout{1:nargout}] = class_interface_mex('us2eq', this.objectHandle, varargin{:});
        end        
        function varargout = expfft1d(this, varargin)
            [varargout{1:nargout}] = class_interface_mex('expfft1d', this.objectHandle, varargin{:});
        end
        function varargout = expifft1d(this, varargin)
            [varargout{1:nargout}] = class_interface_mex('expifft1d', this.objectHandle, varargin{:});
        end
        function varargout = set_grids(this, varargin)
            [varargout{1:nargout}] = class_interface_mex('set_grids', this.objectHandle, varargin{:});
        end        
    end
end
