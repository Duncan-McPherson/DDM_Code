%  The DDM Kernel - Digital and Dynamic Manufacturing Kernel
%
%  This implementation a kernel combining both exp and sine elements to
%  better model cutting forces for feed drive identification.  
%
%  [KMN] = DK(XN, XM, omega, sig1, sig2)
% 
%    XN, XM:       Inputs, XM can never be larger than XN
%    omega:        Frequency of cutting force
%    sig1, sgi2:   Length scale for sine and decay respectivly
%
%   For further details, send email to jdmcpher@uvic.ca
% 2023 - Duncan McPherson
function KMN = DK(mode,XN,XM,hyp)
%     %Attempt to make Rasmussen style code
%     omega = 1/55;
%     % no default
%     if nargin < 1, error('Mode cannot be empty.'); end  
%     % default
%     if nargin < 2, par = []; end                                           
%     varargout = cell(max(1, nargout), 1); 
%     % allocate mem for output
%     if nargin < 4, varargout{1} = covMaha(mode,par); return, end
%     % covariance and derivative
%     k = @(d2) exp(-(sin(pi.*d2./omega).^2)/(2*(sig2^2)))*exp(-d2/(2*(sig1^2)));
%     dk = @(d2,k) (-1/2)*k;         
%     [varargout{:}] = covMaha(mode, par, k, dk, varargin{:});

    %Set Hyperparameters
    sig1 = cell2mat(hyp(1));

    if mode == "LP"
        %Set other hyperparameters
        sig2 = cell2mat(hyp(2));
        omega = cell2mat(hyp(3));

        %Un/comment this for LP kernel code
        KMN1 = exp(-(sin(pi.*pdist2(XN(:,end),XM(:,end))./omega).^2)/(2*(sig2^2)));
        KMN2 = exp(-(pdist2(XN(:,end),XM(:,end)).^2)/(2*(sig1^2)));
        KMN = KMN1.*KMN2;

    elseif mode == "SE"
        %Un/comment this for SE kernel code
        KMN = exp(-sq_dist(XN'/sig1,XM'/sig1)/2);

    elseif mode == "SP"
        %Set other hyperparameters
        sig2 = cell2mat(hyp(2));
        omega = cell2mat(hyp(3));
        P = cell2mat(hyp(4));
        %Un/comment this for making semiparametric kernel
        KMN1 = exp(-(sin(pi.*pdist2(XN(:,end),XM(:,end))./omega).^2)/(2*(sig2^2)));
        KMN2 = exp(-(pdist2(XN(:,end),XM(:,end)).^2)/(2*(sig1^2)));
        KMN = XM*P*XM' + KMN1.*KMN2;
    end
end