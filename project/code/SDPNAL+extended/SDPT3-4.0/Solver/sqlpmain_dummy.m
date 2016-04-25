%%*****************************************************************
%% sqlp: main solver 
%%*****************************************************************
%% SDPT3: version 4.0
%% Copyright (c) 1997 by
%% Kim-Chuan Toh, Michael J. Todd, Reha H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************

  function [blk,At,C,b] = sqlpmain_dummy(blk,At,C,b,par,parbarrier,X0,y0,Z0);

   global spdensity smallblkdim  printlevel msg
   global solve_ok  use_LU  exist_analytic_term numpertdiagschur  
   global schurfun  schurfun_par 
%%
   matlabversion = par.matlabversion;
   vers          = par.vers;
   predcorr      = par.predcorr;
   gam           = par.gam; 
   expon         = par.expon;
   gaptol        = par.gaptol;
   inftol        = par.inftol;
   steptol       = par.steptol;
   maxit         = par.maxit;
   printlevel    = par.printlevel;
   stoplevel     = par.stoplevel;
   scale_data    = par.scale_data;
   spdensity     = par.spdensity;
   rmdepconstr   = par.rmdepconstr;
   smallblkdim   = par.smallblkdim;
   schurfun      = par.schurfun;
   schurfun_par  = par.schurfun_par;
   ublksize      = par.ublksize; 
%%
   tstart = clock; 
   X = X0; y = y0; Z = Z0; 
   for p = 1:size(blk,1)
      if strcmp(blk{p,1},'u'); Z{p} = zeros(blk{p,2},1); end
   end
%%
%%-----------------------------------------
%% convert unrestricted blk to linear blk. 
%%-----------------------------------------
%%
   convertlen = 0; 
   [blk,At,C,X,Z,u2lblk,ublkidx] = sqlpu2lblk(blk,At,C,X,Z,par,convertlen);
   for p = 1:size(blk,1) 
      pblk = blk(p,:); 
      if (u2lblk(p) == 1) 
         n = 2*blk{p,2}; 
         blk{p,1} = 'l';  blk{p,2} = n;
         parbarrier{p} = zeros(1,n);
         At{p} = [At{p}; -At{p}];  
         tau = max(1,norm(C{p})); 
         C{p} = [C{p}; -C{p}]; 
         msg = 'convert ublk to lblk'; 
         if (printlevel); fprintf(' *** %s',msg); end
         b2 = 1 + abs(b');  
         normCtmp = 1+norm(C{p});
         normAtmp = 1+sqrt(sum(At{p}.*At{p}));
         if (n > 1000)
            const = sqrt(n); 
         else
	    const = n; 
         end
         if (par.startpoint == 1)
            X{p} = const* max([1,b2./normAtmp]) *ones(n,1); 
            Z{p} = const* max([1,normAtmp/sqrt(n),normCtmp/sqrt(n)]) *ones(n,1);
            X{p} = X{p}.*(1+1e-10*randmat(n,1,0,'u')); 
            Z{p} = Z{p}.*(1+1e-10*randmat(n,1,0,'u')); 
	 else
            const = max(abs(X{p})) + 100; 
            X{p} = [X{p}+const; const*ones(n/2,1)]; 
            %%old: const = 100; Z{p} = [const*ones(n/2,1); const*ones(n/2,1)];
            Z{p} = [abs(Z0{p}); abs(Z0{p})] + 1e-4; 
         end
      end
   end
%%-----------------------------------------
%% check whether {A1,...,Am} is 
%% linearly independent. 
%%-----------------------------------------
%%
   m0 = length(b); 
   [At,b,y,indeprows,par.depconstr,feasible,par.AAt] = ...
    checkdepconstr(blk,At,b,y,rmdepconstr);
   if (~feasible)
      obj = []; X = cell(size(blk,1),1); y = []; Z = cell(size(blk,1),1); 
      runhist = [];      
      msg = 'SQLP is not feasible'; 
      if (printlevel); fprintf('\n %s \n',msg); end
      return;
   end
   par.normAAt = norm(par.AAt,'fro'); 
%%
%%-----------------------------------------
%% scale SQLP data. Note: must be done only 
%% after checkdepconstr
%%-----------------------------------------
%%
   normA2 = 1+ops(At,'norm'); 
   normb2 = 1+norm(b); 
   normC2 = 1+ops(C,'norm'); 
   normX0 = 1+ops(X0,'norm'); 
   normZ0 = 1+ops(Z0,'norm'); 
   if (scale_data)
      [At,C,b,normA,normC,normb,X,y,Z] = scaling(blk,At,C,b,X,y,Z);
   else
      normA = 1; normC = 1; normb = 1; 
   end 
end
