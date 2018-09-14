%----------------------------------------------------------------------------------------------------------------
% Authors: Hassan Talebi, Ute Mueller, Raimon Tolosana-Delgado
% Paper: "Joint simulation of compositional and categorical data via direct sampling technique – Application to
%         improve mineral resource confidence"
% Journal: Computer & Geosciences
%----------------------------------------------------------------------------------------------------------------


%%%--------------------------
% Simulation parameters
%%%--------------------------
n_cat = 1;                                         % number of categorical variables
n_cont = 4;                                        % number of continuous variables
n = 12;                                            % maximum number of points in the data event
f = 0.8;                                           % maximum fraction of the training image to be scanned
t = 0.025;                                         % overall distance threshold (between 0 and 1)
n_realiz = 1;                                      % number of realizations
postpone=1;                                        % times to postpone the simulation of a node when a close pattern not found
nx=250;                                            % number of nodes in x direction
ny=250;                                            % number of nodes in y direction
nz=1;                                              % number of nodes in z direction
dimx=1;                                            % cell dimension in x direction
dimy=1;                                            % cell dimension in y direction
dimz=1;                                            % cell dimension in z direction
expo_w=2;                                          % exponent of the inverse distance weights
cond_w=5;                                          % weight of hard data vs simulated data
alpha_w=0.5;                                       % mixing coefficient 
var_w=[0.25,0.25,0.25,0.25];                       % weight of each variable (first categorical variables 
rx=50;                                             % maximum search distances along x direction
ry=50;                                             % maximum search distances along y direction
rz=5;                                              % maximum search distances along z direction
denoise=1;                                         % noise removal (resimulating via fully informed nodes)
ilr_trans = true;                                  % transform compositional data to ilr space
plt=true;                                          % plot


% loading the simulation grid: the training data should be transformed to a regular grid. Emty nodes are filled by "NaN"
% first columns are the categorical variables followed by the continuous variables
% the last column is the identification code: 1: training data
%                                             2: to be simulated
%                                             NaN: masked nodes


Grid_TD=load('Grid.txt');
Grid_raw=Grid_TD;
final_simul=nan(nx*ny*nz,n_cat+n_cont,n_realiz);

if ilr_trans
    [ilr,clr,contrast_matrix] = ilr (Grid_TD(Grid_TD(:,end)==1,(n_cat+1):(n_cat+n_cont)));
    Grid_TD(:,end-1)=[];
    Grid_TD(Grid_TD(:,end)==1,(n_cat+1):(n_cat+n_cont-1))=ilr;
    n_cont=n_cont-1;
end

Grid_simul=repmat( Grid_TD,1,1,n_realiz);


% colsest distances between patterns
  bestmin_dist = nan(size(find(Grid_TD(:,end)==2),1)*(size(Grid_TD,2)-1),n_realiz,denoise+1);

% number of tries to find the best pattern
  nbtries = nan(size(find(Grid_TD(:,end)==2),1)*(size(Grid_TD,2)-1),n_realiz,denoise+1);

% list of training nodes in the grid
  list_TD = find(Grid_TD(:,end)==1);

% list of nodes to be simulated
  list_sim = repmat(find(Grid_TD(:,end)==2),size(Grid_TD,2)-1,1);

% square of range for each continuous variable
  cont_norm = (abs(max(Grid_TD(:,n_cat+1:end-1))-min(Grid_TD(:,n_cat+1:end-1)))).^2;

if plt
    coord_plot=IDtoCoord((1:size(Grid_TD,1))',nx,ny,dimx,dimy,dimz);
end


for ii=1:n_realiz
    
    fprintf('Realization number %0.1i \n' , ii)
    Grid=Grid_TD;
    
    
    for kk=1:(denoise+1)
        
        if kk>1
            fprintf('denoising level %0.1i \n' , kk-1)
        end
        progress=10;
        Grid_simul(:,:,ii) = Grid_TD;
    
        % defining a fully random path for simulation (number of nodes to be simulated * number of variables)
        path_sim = (list_sim(randperm(length(list_sim))))';
        postpond=ones(size(path_sim));
        
        for jj=1:(postpone+1)
            
            %looping simulation nodes
            for simnod = 1:size(path_sim,2)
                
                if postpond(simnod)==0
                    continue
                end
                
                %finding where we are in the simulation path
                SimCoord = IDtoCoord(path_sim(simnod),nx,ny,dimx,dimy,dimz);
                
                %finding distance with all training points and previously simulated points
                all_cond_pts = find(Grid(:,end)==1 | Grid(:,end)==3);
                all_cond_pts(all_cond_pts==path_sim(simnod))=[];
                
                CondCoord    = IDtoCoord(all_cond_pts,nx,ny,dimx,dimy,dimz);
                
                % anisotropic search
                CondCoord    = CondCoord(find(sum(abs(CondCoord - repmat(SimCoord,size(CondCoord,1),1)) > ...
                    repmat([rx ry rz],size(CondCoord,1),1),2)<0.1),:);
                
                d = sqrt(sum((CondCoord - repmat(SimCoord,size(CondCoord,1),1)).^2,2));
                
                %sorting
                [d,s] = sort(d);
                CondCoord = CondCoord (s',:);
                
                %taking the closest n or less
                if size(CondCoord,1) < n
                    informed_nodes_counter = size(CondCoord,1);
                else
                    informed_nodes_counter = n;
                end
                
                CondCoord = CondCoord (1:informed_nodes_counter,:);
                all_cond_pts = CoordtoID(CondCoord,nx,ny,nz,dimx,dimy,dimz);
                data_eve_vec = CondCoord - repmat (SimCoord,size(CondCoord,1),1);
                data_event_sim = Grid(all_cond_pts,:);
                
                weights_inv_dist = 1./(d(1:size(data_event_sim,1)).^expo_w);
                weights_cond = data_event_sim(:,end); weights_cond(weights_cond==1)=cond_w; weights_cond(weights_cond==3)=1;
                
                path_ti = (list_TD(randperm(length(list_TD))))';
                
                %scanning ti
                mindist = inf;  %initial best distance is set to inf. Updated with every best distance encountered
                nb_of_tries = ceil(size(path_ti,2)*f); %number of tries in the ti
                
                for i = 1:nb_of_tries
                    
                    % selecting a node in the training data and building the training pattern
                    tiCoord = IDtoCoord(path_ti(i),nx,ny,dimx,dimy,dimz);
                    [all_training_pts,out_win] = CoordtoID((repmat(tiCoord,size(data_eve_vec,1),1)+data_eve_vec),nx,ny,nz,dimx,dimy,dimz);
                    data_event_ti = Grid_TD(all_training_pts,:);
                    data_event_ti(out_win,:) = nan;
                    
                    % reducing the order of statistics
                    if sum(~isnan(data_event_ti(:,1)))<n/2
                        continue
                    end
                    
                    % mixing the weights
                    weights = alpha_w * weights_inv_dist(~isnan(data_event_ti(:,1)),:)/sum( weights_inv_dist(~isnan(data_event_ti(:,1)),:)) + ...
                        (1-alpha_w) * weights_cond(~isnan(data_event_ti(:,1)),:)/sum( weights_cond(~isnan(data_event_ti(:,1)),:));
                    
                    % measuring the distance between two patterns
                    distance = nan(1,n_cat+n_cont);
                    distance(1:n_cat) = sum(repmat(weights,1,n_cat).*(data_event_sim(~isnan(data_event_ti(:,1)),1:n_cat)~=...
                        data_event_ti(~isnan(data_event_ti(:,1)),1:n_cat)));
                    distance(n_cat+1:end) =sqrt( sum((repmat(weights,1,n_cont).*((data_event_sim(~isnan(data_event_ti(:,1)),n_cat+1:end-1)-...
                        data_event_ti(~isnan(data_event_ti(:,1)),n_cat+1:end-1)).^2))./repmat(cont_norm,size(weights,1),1)));
                    distance = sum(var_w .* distance);
                    
                    %checking if the distance is under the minimum distance found so far
                    if distance < mindist
                        mindist = distance;
                        bestpoint = path_ti(i);
                    end
                    
                    % if distance is under t, break the loop and accept the best point
                    if mindist <= t
                        break
                    end
                    
                end
                
                % after searching all the training data if a close pattern was not
                % found, postpone the simulation of this node
                
                if mindist > t && postpond(simnod)<(postpone+1)
                    postpond(simnod)=postpond(simnod)+1;
                    continue
                else
                    
                    % recording number of tries for finding the pattern
                    nbtries(simnod,ii,kk) = i;
                    
                    % recording the best min distance for this node
                    bestmin_dist(simnod,ii,kk) = mindist;
                    
                    % after finding the pattern, paste one variable randomly
                    n_simul=find(isnan(Grid_simul(path_sim(simnod),1:end-1,ii)));
                    n_simul=n_simul(randperm(length(n_simul)));
                    Grid_simul(path_sim(simnod),n_simul(1),ii)=Grid_TD(bestpoint,n_simul(1));
                    Grid(path_sim(simnod),n_simul(1))=Grid_TD(bestpoint,n_simul(1));
                    
                    
                    % if all the variables have been simulated, use this node for conditioning
                    if sum(isnan(Grid(path_sim(simnod),1:end-1)))==0
                        Grid(path_sim(simnod),end)=3;
                    end
                    
                    postpond(simnod)=0;
                      
                end
                
                if round(length(find(postpond==0))/length(postpond)*100)>=progress
                    disp(['  ',num2str(progress),'% completed']);
                    progress = progress+10;
                end          
            end       
        end
    end
    
    final_simul(:,1:n_cat,ii)=Grid_simul(:,1:n_cat,ii);
    final_simul(Grid_simul(:,end)==2,(n_cat+1):(n_cat+n_cont+1),ii) = b_ilr (Grid_simul(Grid_simul(:,end)==2,(n_cat+1):(end-1),ii),contrast_matrix,100);
    xlswrite(num2str(ii),final_simul(:,:,ii));
end

    % plot the training data
    if plt
        for i=1:(n_cat+n_cont+1)
            % training data
            subplot(2,(n_cat+n_cont+1),i),scatter(coord_plot(Grid_raw(:,end)==1,1),coord_plot(Grid_raw(:,end)==1,2),...
                2,Grid_raw(Grid_raw(:,end)==1,i),'fill');colormap jet;axis equal tight;
            % first simulation
            subplot(2,(n_cat+n_cont+1),(n_cat+n_cont+i+1)),scatter(coord_plot(Grid_raw(:,end)==2,1),coord_plot(Grid_raw(:,end)==2,2),...
                2,final_simul(Grid_raw(:,end)==2,i,1),'fill');colormap jet;axis equal tight;
            hold on
        end         
    end




