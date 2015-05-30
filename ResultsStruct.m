classdef ResultsStruct
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        eta
        tTrim
        etaTrim
        bulkStats = struct('Hs', nan, 'Ts', nan, 'dt', nan)
        modelHandle
        estimator
        predictor
        s
        sDot
        sNoisy
        sDotNoisy
        fe
        feHat
        feStar
        userData
    end
    
    methods
        function obj = ResultsStruct(varargin)
            % obj = ResultsStruct(t, eta)
            if nargin > 1
                obj.t = varargin{1};
                obj.eta = varargin{2};
            end
        end
        
        function obj = add_noise(obj, sVar, sDotVar)
            if size(sVar,1) < size(sVar,2); sVar = sVar'; end
            if size(sDotVar,1) < size(sDotVar,2); sDotVar = sDotVar'; end
            obj.sNoisy = obj.s + repmat(sVar, 1, size(obj.s,2)) .* ...
                randn(size(obj.s));
            obj.sDotNoisy = obj.sDot + repmat(sDotVar, 1, size(obj.sDot,2)) .* ...
                randn(size(obj.sDot));
        end
       
        function [rsq, mse] = calc_estimation_results(obj)
            rsq = nan(3,1);
            mse = nan(3,1);
            % if we have an array of results, return results of all element
            % compbined
            if numel(obj) > 1
                % figure out how much data we have
                lenData = 0;
                for ii = 1:numel(obj)
                    lenData = lenData + size(obj(ii).fe,2);
                end
                % fill fe and feHat
                fe = nan(3,lenData);
                feHat = nan(3,lenData);
                idxStart = 1;
                for ii = 1:numel(obj)
                    idxEnd = idxStart-1 + size(obj(ii).fe,2);
                    fe(:,idxStart:idxEnd) = obj(ii).fe;
                    feHat(:,idxStart:idxEnd) = obj(ii).feHat;
                    idxStart = idxEnd + 1;
                end
            else
                % operating on a single element, fe and feHat are passed
                % through
                fe = obj.fe;
                feHat = obj.feHat;
            end
            % calculate mse and rsq
            for ii = 1:size(obj(1).s,1)
                mse(ii) = mean( (fe(ii, :) - feHat(ii, :)).^2);
                rsq(ii) = regression_analysis(fe(ii,:), feHat(ii,:));
            end
            
        end
        
        
        
        function [rsq, mse] = calc_prediction_results(obj)
            rsq = nan(size(obj(1).feStar,1), size(obj(1).feStar,3));
            mse = nan(size(obj(1).feStar,1), size(obj(1).feStar,3));
            if numel(obj) > 1
                % first lets initialize my arrays
                lenData = 0;
                for ii = 1:numel(obj)
                    lenData = lenData + size(obj(1).fe, 2);
                end
                % now fill fe and feStar
                fe = nan(3,lenData);
                feStar = nan(3,lenData,size(obj(1).feStar, 3));
                idxStart = 1;
                idxThereAll = false(lenData,1);
                for ii = 1:numel(obj)
                    idxEnd = idxStart-1 + size(obj(ii).fe,2);
                    idx = idxStart:idxEnd;
                    fe(:,idx) = obj(ii).fe;
                    feStar(:,idx,:) = obj(ii).feStar(:, idx, :);
                    idxNan = find(~isnan(feStar(1,idx,1)));
                    idxThereAll(idxNan(1000:end)) = true;
                end
                % convert this index from logicals to actual indices
                idxThereAll = find(idxThereAll);
                for ii = 1:size(feStar, 1)
                    for jj = 1:size(feStar, 3)
                        rsq(ii, jj) = regression_analysis(...
                                        squeeze(feStar(ii,idxThereAll,jj)), ...
                                        squeeze(fe(ii,idxThereAll+jj-1)));
                        mse(ii, jj) = mean(...
                                       (squeeze(feStar(ii,idxThereAll,jj)) - ...
                                        squeeze(fe(ii,idxThereAll+jj-1))).^2);
                    end
                end
            else
                for ii = 1:size(obj.feStar, 1);
                    for jj = 1:size(obj.feStar, 3)
                        idx = find(~isnan(obj.feStar(ii,:,jj) ));
                        idx = idx(1000:end);
                        rsq(ii, jj) = regression_analysis(...
                                        squeeze(obj.feStar(ii,idx,jj)), ...
                                        squeeze(obj.fe(ii,idx+jj-1)));
                        mse(ii, jj) = mean(...
                                       (squeeze(obj.feStar(ii,idx,jj)) - ...
                                        squeeze(obj.fe(ii,idx+jj-1))).^2);
                    end
                end
            end
        end
        

    end
    
end

