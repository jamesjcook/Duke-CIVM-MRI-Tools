% SPARSE_GRIDDER
% This class calculates and stores two sparse matrices (forward and
% inverse) for the fast calculation of gridding and inverse gridding. 
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef SparseGridder
    properties
        gridSize;
        overgridSize;
        dataSize;
        overgridding;
        nThreads;
        Atrans;
        A;
    end
    methods
        % Constructor
        % Pass traj in -0.5:0.5 units
        function obj = SparseGridder(traj,output_matrix_size,kern_sig,kern_ext,overgridding,nThreads)
            obj.gridSize = ceil(output_matrix_size);
            obj.overgridSize = ceil(output_matrix_size.*overgridding);
            obj.overgridding = obj.overgridSize./obj.gridSize;
            obj.dataSize = [1 size(traj,2)]; % needs to be 1xtotSampleCount
            
            % scale trajectories into overgridded pixel units
            nDims = length(output_matrix_size);
            for iDim = 1:nDims
                traj(iDim,:) = (traj(iDim,:)+0.5)*obj.overgridSize(iDim);%mem rise 30->31
                kern_sig(iDim) = kern_sig(iDim)*obj.overgridding(iDim);
                kern_ext(iDim) = kern_ext(iDim)*obj.overgridding(iDim);
                
                % Check that no trajectory exceeds bounds
                invalidTraj = (traj(iDim,:) > obj.overgridSize(iDim)) | (traj(iDim,:)<1);
                if(any(invalidTraj))
                    error('ERROR: trajectory must be within matrix bounds!');
                end
            end
            
            % Convert to porper classes
            obj.overgridSize = uint32(obj.overgridSize);
            obj.gridSize = uint32(obj.gridSize);
            obj.nThreads = uint32(nThreads);
            
            % Requires that traj be in pixel units not -0.5:0.5!
            obj.Atrans = mex_thread_calcSparseGridMatrix(double(traj), obj.overgridSize,...
                double(kern_sig), double(kern_ext), obj.nThreads); % mem rise from 31->79->50  OR 7->51->24
            obj.A = obj.Atrans';% mem rise 50->73 OR 24->47
        end
        
        
        function varargout = grid(varargin) %cart =  nonCart
            % If user gave output, assume its preallocated and use in-place
            % calculations
            obj = varargin{1};
            nonCart = varargin{2};
            if((nargin == 3) && (nargout == 0))
                % Use in-place
                cart = double(varargin{3});
            else
                % Allocate a new cartesian matrix
                if(isreal(nonCart))
                    cart = zeros(obj.overgridSize);
                else
                    cart = complex(zeros(obj.overgridSize));
                end
            end
            
            mex_thread_sparseMultiply(obj.A,nonCart,cart,obj.nThreads);
            
            % If not using in-place calculations, give back output
            if(nargout == 1)
                varargout{1} = cart;
            end
        end
        
        function varargout = ungrid(varargin) %nonCart = cart
            % If user gave output, assume its preallocated and use in-place
            % calculations
            obj = varargin{1};
            cart = varargin{2};
            if((nargin == 3) && (nargout == 0))
                % Use in-place
                nonCart = double(varargin{3});
            else
                % Allocate a new non-cartesian matrix
                if(isreal(cart))
                    nonCart = zeros(obj.dataSize);
                else
                    nonCart = complex(zeros(obj.dataSize));
                end
            end
            
            mex_thread_sparseMultiply(obj.Atrans,cart,nonCart,obj.nThreads);
            
            % If not using in-place calculations, give back output
            if(nargout == 1)
                varargout{1} = nonCart;
            end
        end
        
        %         function cart = grid(obj, nonCart)
        %             if(isreal(nonCart))
        %                 cart = zeros(obj.overgridSize);
        %             else
        %                 cart = complex(zeros(obj.overgridSize));
        %             end
        %
        %             mex_thread_sparseMultiply(obj.A,nonCart,cart,obj.nThreads);
        %         end
        %
        %         function nonCart = ungrid(obj, cart) %nonCartData, cartData
        %             if(isreal(cart))
        %                 nonCart = zeros(obj.dataSize);
        %             else
        %                 nonCart = complex(zeros(obj.dataSize));
        %             end
        %             mex_thread_sparseMultiply(obj.Atrans,cart,nonCart,obj.nThreads);
        %         end
    end
end