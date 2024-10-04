 function [pnew]=addrmp(pold,nn,nr,scale)
        % [pnew]=addrmp(pold,nn,nr,scale)
        % pold= old precision matrix
        % nn= number of creating new link
        % nr = number of removing link
        % scale = scale of new links scale*(rand(0->1)-0.5)
        %
        % Instead of adding random edges we alter existing
        % relationships, look in upper/lower trinagle for relationships
        % instead of just shifting these we now either create or destory
        % relationships setting precision matrix to be either zero or
        % non-zero...ie sparsity is varied at breakpoints
        %
        % Note: unlike shift, this just randomly adds and removes
        % relationships entirely, instead of adjusting their value
        %
        % output is checked to be positive definite
        %
        % Alex Gibberd March 2013
        
        %% Init
        p=size(pold,1);
    
        list=[];
        p=size(pold,1);
        
        for l=1:p-1
            for j=l+1:p
            if(pold(l,j)~=0)
                list=[list,[l;j]];  % get old edges
            end
            end
        end
        
        %% Error reporting
        len=length(list)
        if(nr>len)
           pnew=pold;
           disp('Not enough edges exist. No operation applied...decrease nr');
           return
        end
        
        if(nr>len)
           pnew=pold;
           disp('Too many edges already exist. No operation applied...decrease nn');
           return
        end
        
         pnew=[];    % Init new precision matrix, start with empty
        % Check for positive definitness
        % Perform final check to see if is actually positive semi-definite
        if(~isempty(pnew))
            [~,E]=chol(pnew);   % Check for cholesky decomposition
                if(E)
                     [pnew]=updateEdges2(pold,list,nr,nn,scale);
                else
                    return
                end
        else        %First run, generate first attempt 
                 [pnew]=updateEdges2(pold,list,nr,nn,scale);
        end
        
        %% Different updateEdges function to shiftp
     function [pnew]=updateEdges(pold,list,nr,nn,scale)
        pnew=pold;  % Used to do just this and add/remove according to edge permutation
        
         %% Select edges for removal
        sel=randperm(length(list));
        selected_edges=sel(1:nr);   % edges to select
        remove_list=list(:,selected_edges);
        
        %% Remove Edges
        for i=1:length(remove_list)
           delta=pold(remove_list(1,i),remove_list(2,i));
            pnew(remove_list(1,i),remove_list(2,i))=pold(remove_list(1,i),remove_list(2,i))-delta;
            pnew(remove_list(2,i),remove_list(1,i))=pold(remove_list(2,i),remove_list(1,i))-delta;    % make symmetric

            %% Try to make positive semi-definite...
            % Note no absolute value required...this should be ok as long
            % as pold is positive semi-definite....dunno about this
            % actually..
            
            pnew(remove_list(1,i),remove_list(1,i))=pold(remove_list(1,i),remove_list(1,i))-abs(delta);
            pnew(remove_list(2,i),remove_list(2,i))=pold(remove_list(2,i),remove_list(2,i))-abs(delta); 
        end
        
        %% Select new edges
        % Make sure they are not the same as those that previously existed
        new_list=[];
        master_list=list;
        while (length(new_list)<nn)
            si=1; sj=1; % to intialise loop below
            % While condition for pair is not met keep drawing pairs
        offij=[];
        offji=[];

            while (sj==si || ~isempty(offij) || ~isempty(offji))
                si=1+round(rand(1)*(p-1));
                sj=1+round(rand(1)*(p-1));
                offij=find(ismember([si,sj],master_list','rows'),1);
                offji=find(ismember([sj,si],master_list','rows'),1);
            end
            new_list=[new_list,[si;sj]];
            % Update master list
            master_list=[master_list,[si;sj]];

        end
        
        %% Create new Edges
        
        for i=1:size(new_list,2)
%            new=scale*(0.5-rand(1)); %Old version
           new=scale*2*(0.5-rand(1));   % New random edge generator 15/10/2014
           while(abs(new)<scale/2)
               new=scale*2*(rand(1)-0.5);   % Keep picking until accepted
           end
            pnew(new_list(1,i),new_list(2,i))=-new;
            pnew(new_list(2,i),new_list(1,i))=-new;    % make symmetric

            %% Try to make positive semi-definite...
            % Note no absolute value required...this should be ok as long
            % as pold is positive semi-definite.
            
            pnew(new_list(1,i),new_list(1,i))=pnew(new_list(1,i),new_list(1,i))+abs(new);
            pnew(new_list(2,i),new_list(2,i))=pnew(new_list(2,i),new_list(2,i))+abs(new); 
        end
        
     end
     %% Different updateEdges function to shiftp
     function [pnew]=updateEdges2(pold,list,nr,nn,scale)
         % Alternate edge updating procedure to ensure that variance doesnt
         % change dramatically either when updating
         % 16/10/2014
         
         pnew=0.5.*eye(p);
         
         %% Add all old edges - those to be removed
         
        %pnew=pold;  % Used to do just this and add/remove according to edge permutation
        
         %% Select edges for removal
        sel=randperm(length(list));
        selected_edges=sel(1:nr);   % edges to select
        remove_list=list(:,selected_edges);
        
        active_list=setdiff(list',remove_list','rows');
        active_list=active_list';
        
        
        %% Add active edges from old stuff
        for i=1:length(active_list)
           delta=pold(active_list(1,i),active_list(2,i));
            pnew(active_list(1,i),active_list(2,i))=pold(active_list(1,i),active_list(2,i))+delta;
            pnew(active_list(2,i),active_list(1,i))=pold(active_list(2,i),active_list(1,i))+delta;    % make symmetric

            %% Try to make positive semi-definite...
            % Note no absolute value required...this should be ok as long
            % as pold is positive semi-definite....dunno about this
            % actually..
            
            pnew(active_list(1,i),active_list(1,i))=pnew(active_list(1,i),active_list(1,i))+abs(delta);
            pnew(active_list(2,i),active_list(2,i))=pnew(active_list(2,i),active_list(2,i))+abs(delta); 
        end
        
        %% Select new edges
        % Make sure they are not the same as those that previously existed
        new_list=[];
        master_list=list;
        while (length(new_list)<nn)
            si=1; sj=1; % to intialise loop below
            % While condition for pair is not met keep drawing pairs
        offij=[];
        offji=[];

            while (sj==si || ~isempty(offij) || ~isempty(offji))
                si=1+round(rand(1)*(p-1));
                sj=1+round(rand(1)*(p-1));
                offij=find(ismember([si,sj],master_list','rows'),1);
                offji=find(ismember([sj,si],master_list','rows'),1);
            end
            new_list=[new_list,[si;sj]];
            % Update master list
            master_list=[master_list,[si;sj]];

        end
        
        %% Create new Edges
        
        for i=1:size(new_list,2)
%            new=scale*(0.5-rand(1)); %Old version
           new=scale*2*(0.5-rand(1));   % New random edge generator 15/10/2014
           while(abs(new)<scale/2)
               new=scale*2*(rand(1)-0.5);   % Keep picking until accepted
           end
            pnew(new_list(1,i),new_list(2,i))=-new;
            pnew(new_list(2,i),new_list(1,i))=-new;    % make symmetric

            %% Try to make positive semi-definite...
            % Note no absolute value required...this should be ok as long
            % as pold is positive semi-definite.
            
            pnew(new_list(1,i),new_list(1,i))=pnew(new_list(1,i),new_list(1,i))+abs(new);
            pnew(new_list(2,i),new_list(2,i))=pnew(new_list(2,i),new_list(2,i))+abs(new); 
        end
        
        
        %% If we want to normalise variance to 1 do below (16/10/2014)
        sigma=inv(pnew);  % calculate inverse of new sigma
        
        % Create vector of inverse variances 
        v=diag(sigma);
        v=1./sqrt(v);
        V=diag(v);
        % Divide by variance associated with edge updates
        sigma=V*sigma*V';
        
        pnew=inv(sigma);
        
        % Prune for really small values in the inverse (numerical??)
        pnew(find(abs(pnew)<1e-12))=0;
     end
 
 end
