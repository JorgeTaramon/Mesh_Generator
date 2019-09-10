function [MESH] = uzawa(MESH,STIFFNESS,F_le,SETTINGS)
error('not coded yet')
% #########################################################################
% NEED TO BE CODED!!!
% #########################################################################

% 5. Panelty force (Uzawa loop)
% Before Assembly, we add the penalty terms in, so they'll be added
% automatically to the right place when using "sparse"

% choose lambda
lambda = 1e8;
p_old = p - 1;
% add new items to KK1, KKi and KKj
iter =0;
penalty_force_temp_x = zeros(nbnd,1);
penalty_force_temp_y = zeros(nbnd,1);

while (max(abs(p(:)-p_old(:)))>1e-6 && iter < 10)
    
    p_old = p;
    x_0 = segBC(pseg,1);
    y_0 = segBC(pseg,2);
    a_temp = segBC(pseg,5);
    b_temp = segBC(pseg,6);
    penalty_force_temp_x = penalty_force_temp_x + lambda*a_temp.*(-a_temp.*x_0 - b_temp.*y_0);
    penalty_force_temp_y = penalty_force_temp_y + lambda*b_temp.*(-a_temp.*x_0 - b_temp.*y_0);
    %
    KK1_penalty = zeros(2,2,nbnd);
    KK1_penalty(1,1,:) = a_temp.^2;
    KK1_penalty(1,2,:) = a_temp.*b_temp;
    KK1_penalty(2,1,:) = a_temp.*b_temp;
    KK1_penalty(2,2,:) = b_temp.^2;
    KK1_penalty = -lambda*KK1_penalty;
    %
    KKi_penalty = zeros(2,2,nbnd);
    KKi_penalty(1,1,:) = 2*[nfix+1:nfix+nbnd]-1;
    KKi_penalty(1,2,:) = KKi_penalty(1,1,:);
    KKi_penalty(2,1,:) = 2*[nfix+1:nfix+nbnd];
    KKi_penalty(2,2,:) = KKi_penalty(2,1,:);
    %
    KKj_penalty = zeros(2,2,nbnd);
    KKj_penalty(1,1,:) = 2*[nfix+1:nfix+nbnd]-1;
    KKj_penalty(1,2,:) = 2*[nfix+1:nfix+nbnd];
    KKj_penalty(2,1,:) = KKj_penalty(1,1,:);
    KKj_penalty(2,2,:) = KKj_penalty(1,2,:);
    %
    penalty_force_vector_bnd = [penalty_force_temp_x, penalty_force_temp_y]';
    penalty_force_vector     = [zeros(2*nfix,1); penalty_force_vector_bnd(:); zeros(2*nint,1)];
    
    %     difference = max(abs(penalty_force_vector-penalty_force_vector_old))  %this line is just for debug
    
    % 6. Assemble into global stiffness (part of Uzawa loop)
    
    
    KK_0 = sparse2(KKi(:),KKj(:),KK1(:)); % MATLAB's sparse matrix maker accumulates all element KK's
    
    
    KKi_with_penalty = [KKi(:);KKi_CrossBars(:);KKi_penalty(:)];
    KKj_with_penalty = [KKj(:);KKj_CrossBars(:);KKj_penalty(:)];
    KK1_with_penalty = [KK1(:);KK1_CrossBars(:);KK1_penalty(:)];
    KK = sparse2(KKi_with_penalty(:),KKj_with_penalty(:),KK1_with_penalty(:)); % MATLAB's sparse matrix maker accumulates all element KK's
    % 5. 1st and last nodes are prescribed
    % free = 2*nfix+2*nbnd+1:2*npt;
    % fix  = 1:2*nfix+2*nbnd;
    free = 2*nfix+1:2*npt;
    fix  = 1:2*nfix;
    ptemp = p';
    ptemp = ptemp(:);
    
    % rhs = -KK(:,fix)*p(fix); % = 0 - KK*pfix
    %rhs = -KK(:,fix)*ptemp(fix); % = 0 - KK*pfix
    rhs = -KK(:,fix)*ptemp(fix) + penalty_force_vector + F_le ; % = 0 - KK*pfix
    KKfree   = KK(free,free);
    
    ptemp(free)  = KKfree\rhs(free);
    pfree = ptemp(free);
    pfree = reshape(pfree,2,[])';
    force = KK_0*ptemp - F_le;
    p=reshape(ptemp,2,[])';
    pbnd = p(nfix+1:nfix+nbnd,:);
    pint = p(nfix+nbnd+1:end,:);
    
    %vectorized version for 'Uzawa' penalty force update
    penalty_force_temp_x = penalty_force_temp_x + lambda*a_temp.*(a_temp.*pbnd(:,1) + b_temp.*pbnd(:,2));
    penalty_force_temp_y = penalty_force_temp_y + lambda*b_temp.*(a_temp.*pbnd(:,1) + b_temp.*pbnd(:,2));
    penalty_force_vector_bnd = [penalty_force_temp_x, penalty_force_temp_y]';
    penalty_force_vector     = [zeros(2*nfix,1); penalty_force_vector_bnd(:); zeros(2*nint,1)];
    %-------------------- old version -------------------------------------
    %     % 'Uzawa'
    %     segstart = 2*nfix+1;
    %     for i=1:nseg
    %         a_temp=segBC(i,5);
    %         b_temp=segBC(i,6);
    %         nodes_this_seg = find(pseg(:)==i);
    %         num_nodes_this_seg = size(nodes_this_seg,1);
    %         for j=1:num_nodes_this_seg    % needs change
    %             penalty_force_vector(segstart) = penalty_force_vector(segstart) + lambda*a_temp*(a_temp*ptemp(segstart) + b_temp*ptemp(segstart+1));
    %             penalty_force_vector(segstart+1) = penalty_force_vector(segstart+1) + lambda*b_temp*(a_temp*ptemp(segstart) + b_temp*ptemp(segstart+1));
    %             segstart = segstart+2;
    %         end
    %     end
    %     penalty_force_vector_old = penalty_force_vector; %for debug
    %     a_temp = segBC(pseg,5);  % for debug
    %     b_temp = segBC(pseg,6);  % for debug
    %----------------------------------------------------------------------
    iter = iter+1;
end
end