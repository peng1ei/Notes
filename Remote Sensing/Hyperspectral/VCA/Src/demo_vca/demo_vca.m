demo_number=2
verbose='on';

switch demo_number 
  case 1
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % abundances（丰度） on strips
         % 3 endmembers (p=3) mixed in the image 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         a = [1.0 0.5 0.2 0.0 0.2 0.0;
              0.0 0.3 0.7 1.0 0.5 0.0;
              0.0 0.2 0.1 0.0 0.3 1.0];
         [p nt]=size(a);
         Lines = 36;
         Columns = 36;
         N = Lines * Columns; % number of pixels
         dim_t = Columns/nt;
         s_o = reshape(repmat( reshape(a',[1 p*nt]),[dim_t*Lines 1]),[N p]);
         
  case 2
         p=3;		%  number of endmembers

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % abundances with Dirichlet distributions（狄利克雷分布） 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Lines = 36;
         Columns = 36;
         N = Lines * Columns; % number of pixels
         s_o = dirichlet_rnd([ones(1,p)/p],N);

  case 3
   
         load cuprite_ref
         N = Lines * Columns; % number of pixels
         p = 6
         A = zeros(L,p);
         s = zeros(p,N);
         x_n = x;clear x;
otherwise         
         error('unknown demo number');
end        
     

if demo_number == 1 | demo_number == 2
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % signatures from USGS Library 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         load signatures
         if strcmp (verbose, 'on'),
            %fprintf('Available signatures:\n%s',names');
            fprintf('Selected signatures:\n%s',names(1:p,:)');
         end
         [L p_max]=size(A);
         A = A(:,1:p);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % illumination fluctuation
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         q=10;
         s = s_o .* random('beta',q,1,N,p);
         x = A * s' ;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % adding noise 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         SNR = 30; %dB
         varianc = sum(x(:).^2)/10^(SNR/10) /L/N ;
         n = sqrt(varianc)*randn([L N]);
         x_n = x + n;
         
end % if demo_number == 1 | demo_number == 2


         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Unmixing Procedure
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         if strcmp (verbose, 'on'), fprintf(1,'Unmixing with VCA algorithm\n'); end
         
         [A_est, location, y] = vca(x_n,...
                                    'Endmembers',p,...
                                    'SNR',SNR,...
                                    'verbose','on');
         s_est = (pinv(A_est) * y)';

         
switch demo_number 
  case { 1 , 2 }
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Permute Results
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
            CRD = corrcoef([A A_est]);
            D = abs(CRD(p+1:2*p,1:p));  
            % permute results
            perm_mtx = zeros(p,p);
            aux=zeros(p,1);
            for i=1:p
                [ld cd]=find(max(D(:))==D); 
                ld=ld(1);cd=cd(1); % in the case os more than one maximum
                perm_mtx(ld,cd)=1; 
                D(:,cd)=aux; D(ld,:)=aux';
            end
            A_est = A_est * perm_mtx;
            s_est = s_est * perm_mtx;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Measures
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ns = diag(s'*s);
            ns_est = diag(s_est'*s_est);
            ang_beta = 180/pi*acos( diag(s'*s_est) ./ sqrt(ns.*ns_est) );
            E_beta = mean(ang_beta.^2)^.5;%  
            
            nA = diag(A'*A);
            nA_est = diag(A_est'*A_est);
            ang_teta = 180/pi*acos( diag(A'*A_est) ./ sqrt(nA .* nA_est) );
            E_teta = mean(ang_teta.^2)^.5;
            
            pA = A./(repmat(sum(A),[L 1]));
            qA = A_est./(repmat(sum(A_est),[L 1])); 
            qA = abs(qA); %qA.*(qA>=0); % warning qA clouf d be < 0 !!!
            SID = sum(pA.*log(pA./qA)) + sum(qA.*log(qA./pA));
            E_SID = mean(SID.^2)^.5;
            
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Show results in Figure
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            coef_A = sqrt(nA_est./nA);	            
            A_est  = A_est ./ repmat(coef_A',[L 1]);
            coef_s = sqrt(ns_est./ns);
            s_est  = s_est ./ repmat(coef_s',[N 1]);

            if strcmp (verbose, 'on'),
               SNR=10*log10(sum(x(:).^2)/sum(n(:).^2));
               fprintf(1,'DEMO %d >\t SNR= %gdB\t E_sid= %.3g\tE_teta= %.3g篭t E_beta= %g篭n', demo_number,SNR,E_SID,E_teta,E_beta);
            end
            
            
            S_est = reshape(s_est,[Lines Columns p]);
            S = reshape(s,[Lines Columns p]);

            figure(1);colormap gray
            for i=1:p 
                   subplot(p,5,5*i-2),
                      imagesc(S(:,:,i));axis('off','equal');
                      if i==1, title('True abundance');end
                   subplot(p,5,5*i-3),
                      plot(1:L,A_est(:,i),'r');
                      if i==1, title('Abundance estimate');end
                   subplot(p,5,5*i-4),
                      plot(1:L,A(:,i),'b');
                      if i==1, title('True signature');end
                   subplot(p,5,5*i-1),
                      imagesc(S_est(:,:,i));axis('off','equal');
                      if i==1, title('Siganture estimate');end
                   subplot(p,5,5*i),
                      plot(1:Columns,S_est(1,:,i),'m',1:Columns,S(1,:,i),'b');
                      if i==1, title('Line of true and estimate abundance');end
            end
 
           figure(2)
           band_i=50;
           band_j=150;
           plot(x_n(band_i,:),x_n(band_j,:),'k.','Markersize',4),
           hold on;
           plot(A(band_i,[1:p 1]),A(band_j,[1:p 1]),'bx-',...
                A_est(band_i,[1:p 1]),A_est(band_j,[1:p 1]),'ro--','LineWidth',2);
           axis([0 1 0 1]);axis square;ylabel('channel 150');xlabel('channel 50');
           title('Reflectance');
           legend('data points','true','estimated')  
           hold off;
           
           figure(3)
           band_i=50;
           band_j=150;
           plot(y(band_i,:),y(band_j,:),'k.','Markersize',4),
           hold on;
           plot(A(band_i,[1:p 1]),A(band_j,[1:p 1]),'bx-',...
                A_est(band_i,[1:p 1]),A_est(band_j,[1:p 1]),'ro--','LineWidth',2);
           axis([0 1 0 1]);axis square;ylabel('channel 150');xlabel('channel 50');
           title('Reflectance');
           legend('projected data points','true','estimated')  
           hold off;

  case 3
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % show components
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
          S_est = reshape(s_est,[Lines Columns p]);
          a_est=nan*ones(224,p);
          a_est(BANDS,:)= A_est;     
          
          n=ceil(sqrt(p));
          for i= 1:p
              figure(1)
                subplot(n,n,i);
                imagesc(S_est(:,:,i));
                colormap gray
                axis('off','equal')
                eval(['title(''component ' int2str(i) ''')']);
              figure(2)  
                subplot(n,n,i);
                plot(wavlen,a_est(:,i),'b-','Linewidth',2);axis([400 2500 0 1]);
                eval(['title(''component ' int2str(i) ''')']);
                xlabel('wavelength (\mum)'); ylabel('reflectance (%)');
          end
            
end  % swhich