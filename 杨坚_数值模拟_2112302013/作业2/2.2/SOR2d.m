% SOR does not converge when double periodic b.c are implemented

SOR_ni_max=600;%max no. of iteration
SOR_c=1.950; %optimun value of the SOR method, between [1.5~2.0]
SOR_onemc=1-SOR_c;

slice=zeros(N-1);
SOR_ni_max=40;
er=[];

   psi_guess=rand(N); % psi^{n}
    psi_guess(1,:)=0;  psi_guess(N,:)=0; psi_guess(:,1)=0; psi_guess(:,N)=0; % Dirrichlet b.c.
  
    psi=zeros(N); % SOR psi^{n+1}

    nt=0; %iteration number
    while nt<SOR_ni_max
        nt=nt+1;
        
        for jj=2:N-1
            
%             if(jj==1)
%                 jjm=N; %jj-1
%             else
                jjm=jj-1;
%             end
%             
%             if(jj==N)
%                 jjp=1;
%             else
                jjp=jj+1;
%             end
            
            
            for kk=2:N-1
                
%                 if(kk==1)
%                     kkm=N;
%                 else
                    kkm=kk-1;
%                 end
%                 if(kk==N)
%                     kkp=1;
%                 else
                    kkp=kk+1;
%                 end
                
                % interior points
                psi(kk,jj)=0.25*SOR_c*(psi(kkm,jj)+psi_guess(kkp,jj)+psi(kk,jjm)+psi_guess(kk,jjp)-dx*dx*zeta(kk,jj))    +   .....
                                  SOR_onemc*psi_guess(kk,jj);
                
            end
        end


        

% %         diagnose the error with respect the zeta_true
 if mod( nt,50)==0
        slice = (psi(ix+1,jy) +psi(ix-1,jy) +psi(ix,jy+1) +psi(ix,jy-1)- 4*psi(ix,jy))*r_dx2; % zeta, relative vorticy
        slice = slice-zeta(2:end-1,2:end-1);
        
        slice=log10(abs( slice));
        er=[er max(slice(:))];
         if er(end)<-15
%              disp([num2str(nt),' SOR iters'])
            break
      
        end
 end

% psi(1,:)=psi(N-1,:);
%         psi(N,:)=psi(2,:);       
% psi(:,1)=psi(:,N-1);
%         psi(:,N)=psi(:,2);
%         psi(1,1)=psi(N-1,N-1);
%         psi(N,N)=psi(2,2);
%          psi(1,N)=psi(N-1,2);
%         psi(N,1)=psi(2,N-1);
        
        psi_guess=psi;
    end
    
% % % % %     for diagnostics    % % % %
%     figure;
%     plot(er,      'b.')
%     title('log10(error in zeta0)')
%     
%     figure;
%     imagesc(slice);colorbar
%     title('log10(diff error)')

