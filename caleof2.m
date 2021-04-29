function [pc,expvar,maps,rmaps,corrs] = caleof2(M,M_lat,M_lon,method,N,plots)

% [PC,EV,MAPS,RMAPS,CORRS] = CALEOF2(M,LAT,LON,METHOD,N)
%
% This function conducts an EOF analysis with several options.
%
% Input arguments...
%
% M:      Data matrix <time x lat x lon>
% LAT:    Vector array containing latitudes (must be <# of lats x 1>)
% LON:    Vector array containing longitudes (order of dimensions does not matter)
% METHOD: (1) EIG ("classic") ...somewhat slow
%         (2) EIGS (rapid "classic") ...very fast
%         (3) SVD ...extremely slow
%         (4) SVDS (rapid SVD) ...quite fast
%         Based on a simple test in the equatorial Pacific,
%         all four methods produce very similar results.
% N:      Number of desired EOFs
% plots:  0 if none, 1 if most, 2 if all including spectra (slower)
%
% Outputs...
%
% PC:    Principal component time series (or expansion coefficients) <N x time>
% EV:    Variance explained by each EOF <N>
% MAPS:  EOFs in map form <lat x lon x N>
% RMAPS: Rotated EOFs in map for <lat x lon x N>
% CORRS: Homogeneous correlation maps (maps of the correlation of the
%        preprocessed input data with each PC) <lat x lon x N>
%
% Guillaume MAZE - LPO/LMD - March 2004
% Revised July 2006, gmaze@univ-brest.fr
% Modified February 2011, Kristopher B. Karnauskas
% Modified April 2019, Marilena Oltmanns
%
% ****************************************


% Get the dimensions of the original input data set
% for reshaping the EOFs back into maps at the end
[nt,ny,nx]=size(M);

% Some preprocessing of the data
disp('Beginning preprocessing steps...')

% Reshape the input matrix for EOF calculations
% from <time x lat x lon> to <time x "station">
M=permute(M,[2 3 1]);
M=reshape(M,size(M,1)*size(M,2),size(M,3));
M=permute(M,[2 1]);

% EOF cannot be done on data w/ NaN values. Change NaNs to a constant.
land=isnan(M);
M(land)=0;

% Remove time mean (i.e., make anomalies)
M=detrend(M,'constant');

% Remove linear trend
%M=detrend(M);

% Filter data with "ideal" bandpass filter
% interval=[(1/500) (1/50)];
% Mf=zeros(size(M));
% for station=1:size(M,2)
%     M_ts=timeseries(M(:,station));
%     M_ts_filtered=idealfilter(M_ts,interval,'pass');
%     Mf(:,station)=M_ts_filtered.Data;
% end
% M=Mf;

% Filter data with recursive low-pass Butterworth filter
% M=filtrage2(M,'low',9,50);

% Weight data by the square root of the cosine of latitude
% this step adds noticeable time to the process
phi=repmat(M_lat',nt,nx);
M=M.*sqrt(cosd(phi));

disp('Preprocessing complete.')

% Need to save a copy of the preprocessed data
% for rotation and correlation/variance maps
M_orig=M;

% Get dimensions of the reshaped/preprocessed data
[n p]=size(M);

% Temporal covariance is p*p matrix, that why max EOF computable is p, 
% so we perform a test on parameter N:
% (Cannot have more EOFs than "stations")
if(N>p)
 disp('Requested N too large. Reducing.')
 N = p; 
end

switch method

    case 1 % CLASSIC METHOD
%================================================================
disp('Beginning EOF calculations with method 1: EIG ("classic")') 

% Transform the data matrix in the correct form (map*time) for eig
F = M';

% Covariance Matrix (inner product over space = covariance in time)
R = F * F';

% Eigenanalysis of the covariance matrix R
[E,L] = eig(R);

% Get PC by projecting eigenvectors on original data
Z = E'*F;

% Make them clear for output
for iN=1:N
   e(iN,:) = squeeze( E(:,p-(iN-1)) )';
   pc(iN,:) = squeeze( Z(p-(iN-1),:) );
end

% Amount of explained variance (at 0.1%)
dsum = diag(L)./trace(L);
for iN=1:N
   expvar(iN)=fix((dsum(p-(iN-1))*100/sum(dsum))*10)/10;
end

% Make EOFs in map form
maps=reshape(permute(e,[2 1]),ny,nx,N);

    case 2 % RAPID CLASSIC METHOD 
%================================================================
disp('Beginning EOF calculations with method 2: EIGS ("rapid classic")')

F = M;

% Covariance Matrix
if n >= p
   R = F' * F;
else 
   R = F * F';
end

% Eigen analysis of the square covariance matrix
[E,L] = eigs(R,N);
if n < p
  E = F' * E;
  sq = [sqrt(diag(L))+eps]';
  sq = sq(ones(1,p),:);
  E = E ./ sq;
end

% Get PC by projecting eigenvectors on original data
if n >= p
   Z = (F*E)';
else
   Z =  E'*F';
end

% Make them clear for output
for iN=1:N
    e(iN,:) = squeeze( E(:,iN) )';
   pc(iN,:) = squeeze( Z(iN,:) );
end

% Amount of variance explained a 0.1 pres et en %
% dsum=diag(L)./trace(L);
dsum=diag(L)/trace(R);
for iN=1:N
   % expvar(iN)=fix((dsum(iN)*100/sum(dsum))*10)/10;
   expvar(iN)=dsum(iN)*100;
end

% Make EOFs in map form
maps=reshape(permute(e,[2 1]),ny,nx,N);

    case 3 % SVD METHOD
%================================================================
disp('Beginning EOF calculations with method 3: SVD')

% Assume that M is (time*map) matrix
[n p]=size(M);

F = M;

% Form the covariance matrix:
R = F'*F;

% Find eigenvectors and singular values
[C,L,CC] = svd(R);
% Eigenvectors are in CC and the squared diagonal values of L
% are the eigenvalues of the temporal covariance matrix R=F'*F

% find the PC corresponding to eigenvalue
PC = F*CC;

% Make them clear for output
for iN=1:N
    e(iN,:) = squeeze( CC(:,iN) )';
   pc(iN,:) = squeeze( PC(:,iN) )';
end

% Amount of variance explained at 0.1%
dsum=diag(L)./trace(L);
if length(dsum)<N % L was not squared
  dsum = [dsum ;zeros(N-length(dsum),1)];
end
for iN = 1 : N
   expvar(iN)=fix( ( dsum(iN)*100/sum(dsum) )*10 ) /10;
end

% Make EOFs in map form
maps=reshape(permute(e,[2 1]),ny,nx,N);

    case 4 % FAST SVD METHOD
%================================================================
disp('Beginning EOF calculations with method 4: SVDS')

% Assume that M is (time*map) matrix
[n p]=size(M);

F = M;

% Form the covariance matrix:
R = F' * F;

% Find eigenvectors and singular values
[C,L,CC,flag] = svds(R,N);
% Eigenvectors are in CC and the squared diagonal values of L
% are the eigenvalues of the temporal covariance matrix R=F'*F
% (Sometimes, CC stops for nul eigenvector, then we need to fill to reach N)
if size(CC,2)<N
  CC = [CC  zeros(size(CC,1),N-size(CC,2)+1)];
end

% find the PC corresponding to eigenvalue
PC = F*CC;
% Which is similar to: C*L

% Make them clear for output
for iN=1:N
    e(iN,:) = squeeze( CC(:,iN) )';
   pc(iN,:) = squeeze( PC(:,iN) )';
end

% Amount of variance explained a 0.1 pres et en %
dsum=diag(L)./trace(L);
if length(dsum)<N % L was not squared
  dsum = [dsum ;zeros(N-length(dsum),1)];
end
for iN=1:N
   expvar(iN)=fix( ( dsum(iN)*100/sum(dsum) )*10 ) /10;
end

% Make EOFs in map form
maps=reshape(permute(e,[2 1]),ny,nx,N);

end % switch method

disp('EOF calculations complete.')

% Rotation (make an option, if nargin=4 (a number), use that many
% top EOFs to form the truncated basis and proceed to calc and plot rotated EOFs)
disp('Performing rotation...')

% Form truncated basis (first few modes) and reconstruct data
% on that basis rather than using full Z
eof1=permute(e(1,:),[2 1]);
eof2=permute(e(2,:),[2 1]);
eof3=permute(e(3,:),[2 1]);
eof4=permute(e(4,:),[2 1]);
E=[ eof1 eof2 eof3 eof4 ];
A=[ M_orig*eof1 M_orig*eof2 M_orig*eof3 M_orig*eof4 ];
ZT=A*E';

U=E;
B=ZT*U;
Unew=zeros(size(U));
%tol=1e-10; % tolerance parameter, set to desired convergence level
tol=1e-10; % tolerance parameter, set to desired convergence level
limit=1;
while limit>tol
    D=diag(mean(B.^2));
    C=B.^3;
    V=ZT'*(C-B*D);
    [w,k]=eig(V'*V);

    % The following is according to 7.26 in Preisendorfer (and no the w are not wrong)
    % This formulation of k^(-1/2) is numerically sensitive, and can be rewritten:
    %     kin=diag(diag(k.^(-1/2)));
    %     kin(~finite(kin))=zeros(size(kin(~finite(kin))));
    %     Ginv=w*kin*w';

    Ginv=w*k^(-1/2)*w';
    Unew=V*Ginv;
    limit=max(max(Unew-U));
    Bold=B;
    U=Unew;
    B=ZT*U;
end

rmaps=reshape(U,ny,nx,4);

disp('Rotation complete.')

% calculate correlations for homogeneous correlation maps
maps_orig=reshape(permute(M_orig,[2 1]),ny,nx,nt);
corrs=zeros(ny,nx,N);

for fn=1:N
    for y=1:ny
        for x=1:nx
            r=corrcoef(maps_orig(y,x,:),pc(fn,:));
            corrs(y,x,fn)=(r(1,2));
        end
    end
end

if plots>0

% Plot some essential results
disp('Plotting some results...')

% EOF maps/eigenvectors
figure(1)
for fn=1:N
    subplot(N,1,fn)
    contourf(M_lon,M_lat,(squeeze(maps(:,:,fn))))
    colorbar
%     hold on
%     set(gca,'FontSize',14)
%     [coast_x,coast_y]=getcoast;
%     hold on
%     plot(coast_x,coast_y,'k')
%     xlabel('Longitude (\circE)')
%     ylabel('Latitude (\circN)')
%     hold off
end

% Rotated EOF maps
figure(2)
for fn=1:4
    subplot(N,1,fn)
    contourf(M_lon,M_lat,rmaps(:,:,fn))
    colorbar
end

% PC time series/expansion coefficients
figure(3)
for fn=1:N
    subplot(N,1,fn)
    plot(pc(fn,:))
end

% Explained variance
figure(4)
bar(expvar(1:N))
ylim([0 100])

if plots==2
% Power spectra
figure(5)
for fn=1:N

    [pow,freq,per,sigpow]=dospec(pc(fn,:),1,95);

    subplot(N,1,fn)
    loglog(freq,pow,'k')
    hold on
    loglog(freq,sigpow,'r')
    ylabel('Power');
    xlabel('Frequency');
end
end % this ends the if plots=2

% Homogeneous correlation maps
% (correlation of original (preprocessed) input data w/ each PC)
figure(6)
for fn=1:N
    subplot(N,1,fn)
    contourf(M_lon,M_lat,corrs(:,:,fn))
    colorbar
end

end % this ends the if plots>0

disp('All done!')