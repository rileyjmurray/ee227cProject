function [x,esq,j] = kmeanlbg_solo(d,k)
%KMEANLBG Vector quantisation using the Linde-Buzo-Gray algorithm [X,ESQ,J]=(D,K)
%
%Inputs:
% D contains data vectors (one per row)
% K is number of centres required
%
%Outputs:
% X is output row vectors (K rows)
% ESQ is mean square error
% J indicates which centre each data vector belongs to
%
%  Implements LBG K-means algorithm:
% Linde, Y., A. Buzo, and R. M. Gray,
% "An Algorithm for vector quantiser design,"
% IEEE Trans Communications, vol. 28, pp.84-95, Jan 1980.


%      Copyright (C) Mike Brookes 1998
%      Version: $Id: kmeanlbg.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc=size(d,2);
[x,esq,j]=v_kmeans_solo(d,1);
m=1;
while m<k
   n=min(m,k-m);
   m=m+n;
   e=1e-4*sqrt(esq)*rand(1,nc);
   [x,esq,j]=v_kmeans_solo(d,m,[x(1:n,:)+e(ones(n,1),:); x(1:n,:)-e(ones(n,1),:); x(n+1:m-n,:)]);
end

function [x,g,j,gg] = v_kmeans_solo(d,k,x0,l)
%KMEANS Vector quantisation using K-means algorithm [X,ESQ,J]=(D,K,X0,L)
memsize=voicebox_solo('memsize'); 
[n_vks,p] = size(d);
nb=min(n_vks,max(1,floor(memsize/(8*p*k))));    % block size for testing data points
nl=ceil(n_vks/nb);                  % number of blocks
if nargin<4
    l=300;                  % very large max iteration count
    if nargin<3
        x0='f';             % use 'f' initialization mode
    end
end
if ischar(x0)
    if k<n_vks
        if any(x0)=='p'                  % Initialize using a random partition
            ix=ceil(rand(1,n_vks)*k);       % allocate to random clusters
            ix(rnsubset_solo(k,n_vks))=1:k;      % but force at least one point per cluster
            x=zeros(k,p);
            for i=1:k
                x(i,:)=mean(d(ix==i,:),1);
            end
        else                                % Forgy initialization: choose k random points [default] 
            x=d(rnsubset_solo(k,n_vks),:);         % sample k centres without replacement
        end
    else
        x=d(mod((1:k)-1,n_vks)+1,:);    % just include all points several times
    end
else
    x=x0;
end
m_vks=zeros(n_vks,1);           % minimum distance to a centre
j=zeros(n_vks,1);           % index of closest centre
gg=zeros(l,1);
wp=ones(1,p);
kk=1:p;
kk=kk(ones(n_vks,1),:);
kk=kk(:);

if l>0
    for ll=1:l                 % loop until x==y causes a break
        
        % find closest centre to each data point [m(:),j(:)] = distance, index
        
        ix=1;
        jx=n_vks-nl*nb;
        for il=1:nl
            jx=jx+nb;        % increment upper limit
            ii=ix:jx;
            z = disteusq_solo(d(ii,:),x,'x');
            [m_vks(ii),j(ii)] = min(z,[],2);
            ix=jx+1;
        end
        y = x;              % save old centre list
        
        % calculate new centres as the mean of their assigned data values (or zero for unused centres)
        
        nd=full(sparse(j,1,1,k,1));         % number of points allocated to each centre
        md=max(nd,1);                       % remove zeros
        jj=j(:,wp);
        x=full(sparse(jj(:),kk,d(:),k,p))./md(:,wp);    % calculate the new means 
        fx=find(nd==0);
        
        % if any centres are unused, assign them to data values that are not exactly on centres
        % choose randomly if there are more such points than needed
        
        if ~isempty(fx)
            q=find(m_vks~=0);
            if length(q)<=length(fx)
                x(fx(1:length(q)),:)=d(q,:);
            else
                if length(fx)>1
                    [~,ri]=sort(rand(length(q),1));
                    x(fx,:)=d(q(ri(1:length(fx))),:);
                else
                    x(fx,:) = d(q(ceil(rand(1)*length(q))),:);
                end
            end
        end
        
        % quit if the centres are unchanged
        
        gg(ll)=sum(m_vks,1);
        if x==y
            break
        end
    end
    gg=gg(1:ll)/n_vks;
%     ll % *** DEBUG ***
%     gg' % *** DEBUG ***
    g=gg(end);
else            % if l==0 then just calculate G and J (but rename as X and G)
    ix=1;
    jx=n_vks-nl*nb;
    for il=1:nl
        jx=jx+nb;        % increment upper limit
        ii=ix:jx;
        z = disteusq_solo(d(ii,:),x,'x');
        [m_vks(ii),j(ii)] = min(z,[],2);
        ix=jx+1;
    end
    x=sum(m_vks,1)/n_vks;
    g=j;
end
end

function y=voicebox_solo(f,v)
%VOICEBOX  set global parameters for Voicebox functions Y=(FIELD,VAL)
persistent PP
if isempty(PP)

    % System-dependent directory paths and constants

    PP.dir_temp='F:\TEMP';                      % directory for storing temporary files
    PP.dir_data='E:\dmb\data\speech';           % default directory to preappend to speech data file names
    PP.shorten='C:\bin\shorten.exe';            % location of shorten executable
    PP.flac='C:\bin\flac.exe';                  % location of flac executable
    PP.sfsbin='F:\Program Files\SFS\Program';   % location of Speech Filing Sysytem binaries
    PP.sfssuffix='.exe';                        % suffix for Speech Filing Sysytem binaries
    PP.memsize=50e6;                            % Maximum amount of temporary memory to use (Bytes)

    % DYPSA glottal closure identifier

    PP.dy_cpfrac=0.3;           % presumed closed phase fraction of larynx cycle
    PP.dy_cproj=0.2;            % cost of projected candidate
    PP.dy_cspurt=-0.45;         % cost of a talkspurt
    PP.dy_dopsp=1;              % Use phase slope projection (1) or not (0)?
    PP.dy_ewdly=0.0008;         % window delay for energy cost function term [~ energy peak delay from closure] (sec)
    PP.dy_ewlen=0.003;          % window length for energy cost function term (sec)
    PP.dy_ewtaper=0.001;        % taper length for energy cost function window (sec)
    PP.dy_fwlen=0.00045;        % window length used to smooth group delay (sec)
    PP.dy_fxmax=500;            % max larynx frequency (Hz)
    PP.dy_fxmin=50;             % min larynx frequency (Hz)
    PP.dy_fxminf=60;            % min larynx frequency (Hz) [used for Frobenius norm only]
    PP.dy_gwlen=0.0030;         % group delay evaluation window length (sec)
    PP.dy_lpcdur=0.020;         % lpc analysis frame length (sec)
    PP.dy_lpcn=2;               % lpc additional poles
    PP.dy_lpcnf=0.001;          % lpc poles per Hz (1/Hz)
    PP.dy_lpcstep=0.010;        % lpc analysis step (sec)
    PP.dy_nbest=5;              % Number of NBest paths to keep
    PP.dy_preemph=50;           % pre-emphasis filter frequency (Hz) (to avoid preemphasis, make this very large)
    PP.dy_spitch=0.2;           % scale factor for pitch deviation cost
    PP.dy_wener=0.3;            % DP energy weighting
    PP.dy_wpitch=0.5;           % DP pitch weighting
    PP.dy_wslope=0.1;           % DP group delay slope weighting
    PP.dy_wxcorr=0.8;           % DP cross correlation weighting
    PP.dy_xwlen=0.01;           % cross-correlation length for waveform similarity (sec)

    % RAPT pitch tracker

    PP.rapt_f0min=50;           % Min F0 (Hz)
    PP.rapt_f0max=500;          % Max F0 (Hz)
    PP.rapt_tframe=0.01;        % frame size (s)
    PP.rapt_tlpw=0.005;         % low pass filter window size (s)
    PP.rapt_tcorw=0.0075;       % correlation window size (s)
    PP.rapt_candtr=0.3;         % minimum peak in NCCF
    PP.rapt_lagwt=0.3;          % linear lag taper factor
    PP.rapt_freqwt=0.02;        % cost factor for F0 change
    PP.rapt_vtranc=0.005;       % fixed voice-state transition cost
    PP.rapt_vtrac=0.5;          % delta amplitude modulated transition cost
    PP.rapt_vtrsc=0.5;          % delta spectrum modulated transition cost
    PP.rapt_vobias=0.0;         % bias to encourage voiced hypotheses
    PP.rapt_doublec=0.35;       % cost of exact doubling or halving
    PP.rapt_absnoise=0;         % absolute rms noise level
    PP.rapt_relnoise=2;         % rms noise level relative to noise floor
    PP.rapt_signoise=0.001;     % ratio of peak signal rms to noise floor (0.001 = 60dB)
    PP.rapt_ncands=20;          % max hypotheses at each frame
    PP.rapt_trms=0.03;                      % window length for rms measurement
    PP.rapt_dtrms=0.02;                     % window spacing for rms measurement
    PP.rapt_preemph=-7000;                  % s-plane position of preemphasis zero
    PP.rapt_nfullag=7;                      % number of full lags to try (must be odd)
    
    % now see if an environment variable has been set
    
    vbenv=winenvar_solo('VOICEBOX');
    if exist(vbenv,'file');     % update with locally defined values if defined
        run(vbenv)
    end

    % now check some of the key values for validity

    if exist(PP.dir_temp,'dir')~=7        % check that temp directory exists
        PP.dir_temp = winenvar_solo('temp');     % else use windows temp directory
    end

    [fnp,~,~]=fileparts(mfilename('fullpath'));
    if exist(PP.shorten)~=2        % check that shorten executable exists
        PP.shorten=fullfile(fnp,'shorten.exe'); % next try local directory
        if exist(PP.shorten)~=2        % check if it exists in local directory
            PP.shorten='shorten.exe'; % finally assume it is on the search path
        end
    end

    if exist(PP.flac)~=2        % check that flac executable exists
        PP.flac=fullfile(fnp,'flac.exe'); % next try local directory
        if exist(PP.flac)~=2        % check if it exists in local directory
            PP.flac='flac.exe'; % finally assume it is on the search path
        end
    end

end
if nargin==0
    if nargout==0
        % list all fields
        nn=sort(fieldnames(PP));
        cnn=char(nn);
        fprintf('%d Voicebox parameters:\n',length(nn));

        for i=1:length(nn);
            if ischar(PP.(nn{i}))
                fmt='  %s = %s\n';
            else
                fmt='  %s = %g\n';
            end
            fprintf(fmt,cnn(i,:),PP.(nn{i}));
        end
    else
        y=PP;
    end
elseif nargin==1
    if isfield(PP,f)
        y=PP.(f);
    else
        y=[];
    end
else
    if isfield(PP,f)
        PP.(f)=v;
        y=PP;
    else
        error('''%s'' is not a valid voicebox field name',f);
    end
end
end

function m = rnsubset_solo(k,n)
%RNSUBSET choose k distinct random integers from 1:n M=(K,N)
if k>n
    error('rnsubset: k must be <= n');
end
% We use two algorithms according to the values of k and n
[~,e_rs]=log2(n);
if k>0.03*n*(e_rs-1)
[~,m]=sort(rand(1,n)); % for large k, just do a random permutation
else
    v=ceil(rand(1,k).*(n:-1:n-k+1));
    m=1:n;
    for i=1:k
        j_rs=v(i)+i-1;
        x_rs=m(i);
        m(i)=m(j_rs);
        m(j_rs)=x_rs;
    end
end
m=m(1:k);
end

function d=winenvar_solo(n)
%WINENVAR get windows environment variable [D]=(N)
p=['%',n,'%'];
[~,d]=system(['echo ',p]);
while d(end)<=' ';
    d(end)=[];
end
if strcmp(d,p)
    d=[];
end
end

function d=disteusq_solo(x,y,mode,w)
%DISTEUSQ calculate euclidean, squared euclidean or mahanalobis distance D=(X,Y,MODE,W)
[nx,p]=size(x); ny=size(y,1);
if nargin<3 || isempty(mode), mode='0'; end
if any(mode=='d') || (mode~='x' && nx==ny)

    % Do pairwise distance calculation

    nx=min(nx,ny);
    z=double(x(1:nx,:))-double(y(1:nx,:));
    if nargin<4
        d=sum(z.*conj(z),2);
    elseif min(size(w))==1
        wv=w(:).';
        d=sum(z.*wv(ones(size(z,1),1),:).*conj(z),2);
    else
        d=sum(z*w.*conj(z),2);
    end
else
    
    % Calculate full distance matrix
    
    if p>1
        
        % x and y are matrices
        
        if nargin<4
            z=permute(double(x(:,:,ones(1,ny))),[1 3 2])-permute(double(y(:,:,ones(1,nx))),[3 1 2]);
            d=sum(z.*conj(z),3);
        else
            nxy=nx*ny;
            z=reshape(permute(double(x(:,:,ones(1,ny))),[1 3 2])-permute(double(y(:,:,ones(1,nx))),[3 1 2]),nxy,p);
            if min(size(w))==1
                wv=w(:).';
                d=reshape(sum(z.*wv(ones(nxy,1),:).*conj(z),2),nx,ny);
            else
                d=reshape(sum(z*w.*conj(z),2),nx,ny);
            end
        end
    else
        
        % x and y are vectors
        
        z=double(x(:,ones(1,ny)))-double(y(:,ones(1,nx))).';
        if nargin<4
            d=z.*conj(z);
        else
            d=w*z.*conj(z);
        end
    end
end
if any(mode=='s')
    d=sqrt(d);
end
end

end