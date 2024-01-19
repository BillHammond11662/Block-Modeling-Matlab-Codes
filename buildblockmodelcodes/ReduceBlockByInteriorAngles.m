function [blocks,faults]=ReduceBlockByInteriorAngles(blocks,nodes,faults,bnds,Amin,Amax,optverb)
%  [blocks,faults]=ReduceBlockByInteriorAngles(blocks,nodes,faults,bnds,Amin,Amax,optverb)
%
% Routine combines block with adjacent block along boundary  
% that is not a fault or model boundary. 

[~,~,~,~,~,angmax,angmin,~] = blockgeoms(blocks,nodes);

% Eliminate blocks with very small or large interior angles

bnames = fieldnames(blocks);
bnamesout = bnames(angmin<Amin | angmax>Amax);

for b=bnamesout'
    
    if optverb
       disp(['Trying to remove block #: ' char(b) ]);
    end
    
    % compute interior angles
    
    if isfield(blocks,char(b))
        nn =eval(['blocks.' char(b)]);
        if nn(1)==nn(end)
            nn(end)=[];
        end
        
        angs  = blockangles(nodes(nn,1),nodes(nn,2));  % one angle per node

        isfs=false(size(nn));
        L=nan(size(nn));
        seglist=[];
        for j=1:length(nn)
            lon1=nodes(nn(j),1);
            lat1=nodes(nn(j),2);
            if j==length(nn)
                lon2=nodes(nn(1),1);
                lat2=nodes(nn(1),2);
            else
                lon2=nodes(nn(j+1),1);
                lat2=nodes(nn(j+1),2);
            end
            [L(j),~] = distance(lat1,lon1,lat2,lon2);
            
            % is this segment a boundary segment?
            if ((lon1==bnds(1) && lon2==bnds(1)) || ...
                    (lon1==bnds(2) && lon2==bnds(2)) || ...
                    (lat1==bnds(3) && lat2==bnds(3)) || ...
                    (lat1==bnds(4) && lat2==bnds(4)) )
                isfs(j)=true;
            end
            
            % is this segment a fault segment?
            fnum=[];
            if j==length(nn)
                fnum=nodes2fault(faults,nn(j),nn(1));
            else
                fnum=nodes2fault(faults,nn(j),nn(j+1));
            end
            if ~isempty(fnum)
                isfs(j)=true;
            end
            
            if j<length(nn)
               jends = [j j+1];
            else
               jends = [length(nn) 1];
            end
            
            if isfs(j)==false
                if any(angs(jends)<Amin) || any(angs(jends)>Amax)
                    seglist=[seglist; j];
                end
            end
        end
        
        if isempty(seglist)
            if optverb
               disp(['   Not removing block #: ' char(b) '.  No candidate segments.']);
            end
        else
            if optverb
               disp(['   Working on removing block #: ' char(b)]);
            end
            
            % Go through the list of segment candidates and test their
            % exceedence.  
            exceed0=max([max([(angs-Amax) (Amin-angs)]) 0]);
            exceed=zeros(size(seglist));
            for k=1:length(seglist)
                sn=seglist(k);
                
                if sn==length(nn)
                    segnds=[nn(sn);nn(1)];
                else
                    segnds=[nn(sn);nn(sn+1)];
                end
                
                [blocks2,~,blocksout,faultsout]=RemoveBlock(blocks,nodes,faults,segnds);
                
                % don't allow this if faults are removed or if no block was
                % removed owing to error
                if (isempty(faultsout) && ~isempty(blocksout))
                    
                    bnames2 = fieldnames(blocks2);
                    blockschanged = setdiff(bnames2,bnames);
                    
                    for q=1:length(blockschanged)
                        
                        % check block angles of blocks that changed.
                        nn2=eval(['blocks2.' char(blockschanged(q))]);
                        
                        angs2  = blockangles(nodes(nn2,1),nodes(nn2,2));  % one angle per node
                        
                        exceed(k)= max([max([(angs2-Amax) (Amin-angs2)]) 0 exceed(k)]);
                        
                    end
                end
                
            end
            
            % Take the lowest exceedence if it is lower than the original.
            % i.e., if removing that segments results in improvement of the model (global) in terms
            % of exceedence
            
            if (min(exceed)<=exceed0 && min(exceed>0))
                
                [~,jmin]=min(exceed);
                sn=seglist(jmin);
                              
                if sn==length(nn)
                    segnds=[nn(sn);nn(1)];
                else
                    segnds=[nn(sn);nn(sn+1)];
                end
                
                [blocks,faults,~,~]=RemoveBlock(blocks,nodes,faults,segnds);                
                
            else
                if optverb
                   disp(['Warning: Not combining block ' char(b) ' Exceedence not reduced or removes fault.']);
                end
            end
                
        end
        
    else
        if optverb
           disp(['Cannot find block ' char(b) ', may have been combined already']);
        end
    end
end

%%  perform imporant checks on model

CheckFaults(nodes,blocks,faults);

blocks=CheckBlocks(blocks);


 
