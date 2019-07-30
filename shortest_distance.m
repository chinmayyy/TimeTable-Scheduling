%%
n = 7  %number of vertices

%{
A = randi(2,n,n) -1 ;
A = A - tril(A, -1) + triu(A, 1)';
%adj = G.Adj
A = A .* (1 - eye(n,n))
%}

A = [0 1 1 0 0 0 0; 1 0 1 1 0 0 0; 1 1 0 1 0 0 0; 0 1 1 0 1 1 0; 0 0 0 1 0 0 1; 0 0 0 1 0 0 1; 0 0 0 0 1 1 0]

heatmap(A)
title('Adjacency Matrix')

G_try = digraph(A)

dist= ones(n,n) ./ 0.00000001
for i = 1:n
    for j = 1:n
        if i == j
            dist(i,j) = 0
        end
    end
end

paths = eye(n,n)

figure
plot(G_try)
title('Graph')

%{
for i = 1:n
    for j = 1:n
        neighbors = A(j,:)
        neighbor_index = find(neighbors)
        for k = neighbor_index
            x = [i,j]
            y = [i,k]
            if dist(y(1,1),y(1,2)) > dist(x(1,1),x(1,2)) + 1   %dist[y] > dist[x] +1 
                dist(y(1,1),y(1,2)) = dist(x(1,1),x(1,2)) +1 
                paths(y(1,1),y(1,2)) = paths(x(1,1),x(1,2)) 
            elseif dist(y(1,1),y(1,2)) == dist(x(1,1),x(1,2)) +1 
                paths(y(1,1),y(1,2)) = paths(y(1,1),y(1,2)) + paths(x(1,1),x(1,2)) 
            end
        end
    end
end
%}

for i = 1:n
    src = i
    visited = zeros(1,n)
    q = [src]
    visited(src) = 1
    while isempty(q) == 0
        curr = q(end)
        q = q(1:end-1)
        neighbors = A(curr,:)
        neighbor_index = find(neighbors)
        for j = neighbor_index
            if visited(j) == 0
                q = [j q]
                visited(j) = 1
            end
            if dist(i,j) > dist(i,curr) + 1
                dist(i,j) = dist(i,curr) + 1
                paths(i,j) = paths (i,curr)
            elseif dist(i,j) == dist(i,curr) + 1
                paths(i,j) = paths(i,j) + paths(i,curr)
            end
        end
    end
end

%{        
for i = 1:n
    src = i
    visited = zeros (1,n)
    q = [src]
    visited(src) = 1
    if isempty(q)==0
        curr = q(1,end)
        neighbors = A(curr,:)
        neighbor_index = find(neighbors)
        for k = neighbor_index
            if visited(k) == 0
                q = [k q]
                visited(k) = 1
            end
            if dist(i,k) > dist(i,curr) + 1
                dist(i,k) = dist(i,curr) + 1
                paths(i,k) = paths(i, curr)
            elseif dist(i,k) == dist(i,curr) + 1
                paths(i,k) = paths(i,k) + paths(i, curr)
            end
        end
    end
end
%}

dist
paths

%{
dist = dist - tril(dist,-1) + triu(dist,1)'
paths = paths - tril(paths,-1) + triu(paths,1)'
%}

%% FUNCTIONS %%
function [G]=erdosRenyi(nv,p,Kreg)
%Funciton [G]=edosRenyi(nv,p,Kreg) generates a random graph based on
%the Erdos and Renyi algoritm where all possible pairs of 'nv' nodes are
%connected with probability 'p'. 
%
% Inputs:
%   nv - number of nodes 
%   p  - rewiring probability
%   Kreg - initial node degree of for regular graph (use 1 or even numbers)
%
% Output:
%   G is a structure inplemented as data structure in this as well as other
%   graph theory algorithms.
%   G.Adj   - is the adjacency matrix (1 for connected nodes, 0 otherwise).
%   G.x and G.y -   are row vectors of size nv wiht the (x,y) coordinates of
%                   each node of G.
%   G.nv    - number of vertices in G
%   G.ne    - number of edges in G
%
%Created by Pablo Blinder. blinderp@bgu.ac.il
%
%Last update 25/01/2005
%build regular lattice 
A=sparse(nv,nv);
Kreg=fix(abs(Kreg)/2);Kreg=(Kreg<1)+Kreg;
for k=1:Kreg
    A=sparse(A+diag(ones(1,length(diag(A,k))),k)+diag(ones(1,length(diag(A,nv-k))),nv-k));
end
ne0=nnz(A);
%find connected pairs
[v1,v2]=find(A);
% P=permPairs(nv);%my version is faster
Dis=(rand(length(v1),1)<=p);%pairs to disconnect
A(v1(Dis),v2(Dis))=0;
vDis=unique([v1(Dis),v2(Dis)]);%disconnected vertices
nDis=ne0-nnz(A);sum(Dis);
%cycle trough disconnected pairs
disconPairs=[v1(Dis),v2(Dis)];
for n=1:nDis
    %choose one of the vertices from the disconnected pair
    i=ceil(rand*size(disconPairs,1));
    j=logical(1+rand>0.5);
    vDisToRec=disconPairs(i,j);
    %find non adjacent vertices and reconnect
    adj=[find(A(:,vDisToRec)) ; find(A(vDisToRec,:))'];
    nonAdj=setdiff(1:nv,adj);
    vToRec=nonAdj(ceil(rand*length(nonAdj)));
    S=sort([vDisToRec vToRec]);
    A(S(1),S(2))=1;
end
[x,y]=getNodeCoordinates(nv);
%make adjacency matrix symetric
A=A+fliplr((flipud(triu(A))));
G=struct('Adj',A,'x',x','y',y','nv',nv,'ne',nnz(A));
end
function [x,y]=getNodeCoordinates(nv)
%Adapted from circle.m by Zhenhai Wang <zhenhai@ieee.org>. For more details
%see under  MATLAB Central >  File Exchange > Graphics > Specialized
%Plot and Graph Types > Draw a circle.
center=[0,0];
theta=linspace(0,2*pi,nv+1);
rho=ones(1,nv+1);%fit radius and nv
[X,Y] = pol2cart(theta',rho');
X=X+center(1);
Y=Y+center(2);
x=X(1:end-1)*10;
y=Y(1:end-1)*10;
end
function P=permPairs(N)
%Produces all pairs of pairs from 1 to N.
%It is ~30% to 50% faster that nchoosek(1:N,2).
%Created by Pablo 02/12/2003
ini_i=1;
ini_j=ini_i+1;
r0=1;
P=[];
for i=ini_i:N-1
    lj=N-ini_i;
    P(:,r0:lj+r0-1)=[ones(1,lj)*i;ini_j:N];
    r0=r0+lj;
    ini_i=ini_i+1;
    ini_j=ini_i+1;
end
P=P';
end
function plotGraphBasic(G,markerSize,addText)
%function plotGraph(G) plots graph G.
% Inputs:
%   G is a structure inplemented as data structure in this as well as other
%   graph theory algorithms.
%   G.Adj   - is the adjacency matrix (1 for connected nodes, 0 otherwise).
%   G.x and G.y -   are row vectors of size nv wiht the (x,y) coordinates of
%                   each node of G.
%   G.nv    - number of vertices in G
%   G.ne    - number of edges in G
%
%   markerSize  -  controls the size of each node in the graph
%   addText - toggles text display (1 - on, 0 off).
%   Note: The color of each node is computed based on the
%         its degree. 
%
%Created by Pablo Blinder. blinderp@bgu.ac.il
%
%Last updated 25/01/2005
%generate plot. Decompose to single lines for more detailed formatting
figure;
[XX,YY]=gplot(G.Adj,[G.x' G.y'],'k-');
i=~isnan(XX);
XX=XX(i);YY=YY(i);
XX=reshape(XX,2,length(XX)/2);
YY=reshape(YY,2,length(YY)/2);
hLines=line(XX,YY);
set(hLines,'color','k');
hold on;
kv=full(diag(G.Adj*G.Adj));
kvGroups=unique(setdiff(kv,0));
nGroups=length(kvGroups);
map=jet(max(kvGroups)); 
kv(kv<1)=1;%scale lowest to first 
Pv=num2cell(map(kv,:),2);
if kvGroups==1; kvGroups=2; end %Safeguard aginst single values
set(gca,'Clim',[1 max(kvGroups)]);
Pn(1)={'MarkerFaceColor'};
% Now draw the plot, one line per point.
h = [];
for i=1:G.nv
    h = [h;plot(G.x(i),G.y(i),'ko')];
end
ht=[];
ti=1;
if addText
    for i=1:G.nv
        if ti
            ht = [ht;text(G.x(i)+0.1*G.x(i),G.y(i)+0.1*G.y(i),num2str(i))];
            ti=0;
        else
            ti=1;
        end
    end
end
set(h,'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerSize',markerSize,Pn,Pv);
set(gca,'Visible','Off','YDir','reverse');
colormap(map);
hc=colorbar;
set(hc,'FontSize',8,'FontW','Demi')
set(hc,'Visible','off')
set(gcf,'Color','w')
end
function [searchNodes,iterations] = bfs(source,target,names,startNode)
%    BFS performs breadth first search on graph with source and target
%    vectors.
%     
%    Syntax:     
%          
%    [searchNodes,iterations] = bfs(source,target,names,startNode)
%    [searchNodes,iterations] = bfs(source,target,startNode)
%    
%    Inputs:
%
%    source = Vector or cell array containing starting node of each of the edge.
%    target = Vector or cell array containing ending node of each of the edge.
%    names = Cell array containing string names of each of the node.
%    startNode = Initial node in the graph.
%    
%    Outputs:
%
%    path = Cell array containing search path.
%    iterations = Table containing bfs iteration summary.
% 
%    Example 01:
%     
%    s = {'A','A','A','B','B','C'};
%    t = {'B','C','D','E','F','G'};
%    [searchNodes,iterations] = bfs(s,t,'A')
%     
%    Example 02:
%     
%    s = [1 1 1 2 2 3];
%    t = [2 3 4 5 6 7];
%    names = {'A','B','C','D','E','F','G'};
%    [searchNodes,iterations] = bfs(s,t,names,'A')
%     
%    Example 03:
%     
%    s = [1 1 1 2 2 3];
%    t = [2 3 4 5 6 7];
%    [searchNodes,iterations] = bfs(s,t,1)
%     
%    Coded by Ali Asghar Manjotho
%    Lecturer, CSE-MUET
%    Email: ali.manjotho.ali@gmail.com
    iterations = table;
    
    % Refactor source, target, names & startNode vectors to numbers
    
    % If names argument missing 
    if (nargin<4)
        
        % Third argument (i.e. names) is starting Node
        startingNode = names;
        
        [s,t,n,sNode] = refactor(source,target,startingNode);
        
    else
        
        % Fourth argument (i.e. startNode) is starting Node
        startingNode = startNode;
        
        [s,t,n,sNode] = refactor(source,target,names,startingNode);        
        
    end
    
    % Get all unique nodes from source and target vectors
    uniqueNodes = getNodes(s,t);
     
    % Initialize visited list and queue
    visited = [];
    queue = [];
     
    % Set starting node as current node and add it in to visited list
    currentNode = sNode;
    visited = [visited sNode];
          
    % Local variables to track iteration number
    iteration = 1;   
    
    % Update Iterations table
    iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,'Starting Node')];
    
    
    % Repeat until queue is empty
    while(~isempty(currentNode))
        
        % Get all childs of current node
        childs = getChilds(s,t,currentNode);
        
            
        for i=1:length(childs)
           
            
            % If new unvisited child found add it in queue and visited list
            if(length(find(visited==childs(i)))==0)
                queue = [queue childs(i)];
                visited = [visited childs(i)];
                
                % Increase iteration number
                iteration = iteration + 1;
                iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,strcat('Unvisited node found�',n(childs(i))))];
            end          
        end
         
              
         
        % If no new child found for current node then remove first item
        % from queue and make it as current node
        if (length(queue)>0)            
            currentNode = queue(1);
            queue(1) = [];
            
            
            comments = strcat('Dequeue�',n(currentNode),' from queue');
            
            
            iteration = iteration + 1;
            iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,comments)];
        else
            currentNode = [];
        end
        
               
         
    end
    
    
    % Update iteration table for last state of queue, visited list, current
    % node and search path when after queue is empty
    if(iteration > 0)
        
        iteration = iteration +1;
        
        iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,'BFS Converged')];
    end
    
    searchNodes = n(visited);
    iterations.Properties.VariableNames = {'Iteration' 'CurrentNode' 'Queue' 'Visited' 'Comments'};
    
        
end
function childs = getChilds(source,target,node)
    
    childs = sort(target(find(source==node)));
    
end
function nodes = getNodes(s,t)
    nodes = unique(horzcat(s,t));
end
function [s,t,n,sn] = refactor(source,target,names,startNode)
    % If names argument missing 
    if (nargin<4)
        
        % Third argument (i.e. names) is starting Node
        sn = names;
    else
        
        % Fourth argument (i.e. startNode) is starting Node
        sn = startNode;  
        
    end
    
    % Get all unique nodes
    uNodes = unique(horzcat(source,target));
        
        
    
    % If source and target are cell arrays
    if(iscell(source) && iscell(target))
    
        % If names argument missing
        if(nargin<4)
            n = uNodes;
        else
            n = names;
        end
        
        
        % Get unique nodes cell array
        uNodes = unique(horzcat(source,target));
        s = [];
        t = [];
        % Populate source and target with equivalent numeric values
        for i=1:length(source)
            [sFound,sIndex] = ismember(source(i),uNodes);
            [tFound,tIndex] = ismember(target(i),uNodes);
            s = [s sIndex];
            t = [t tIndex];
        end
            
        
        
        
    else
        
        s = source;
        t = target;
        
        % If names argument missing
        if(nargin<4)    
            
            uNodes = unique(horzcat(source,target));
            n = cell(1,length(uNodes));
            
            
            for i=1:length(uNodes)
                n{i} = num2str(uNodes(i));
            end
            
        else
            n = names;
        end
    end
    
    
    
    % If starting node is not a number
    if(~isnumeric(sn))
        sn = find(ismember(n,sn));
        
    end
end
function tableIteration = updateTable(s,t,n,currentNode,queue,visited,iteration,comments)
   
    uniqueNodes = getNodes(s,t);
    % Display current queue
    queueStr = '[';
    
    for i=length(queue):-1:1
        if(i==1)
            queueStr = strcat(queueStr,sprintf('%s',char(n(find(uniqueNodes==queue(i))))));
        else
            queueStr = strcat(queueStr,sprintf('%s�',char(n(find(uniqueNodes==queue(i))))));
        end
    end
    queueStr = strcat(queueStr,sprintf(']'));
    
    
    % Display current visited list     
    visitedStr = '[';
    
    for i=1:length(visited)
        if(i==length(visited))
            visitedStr = strcat(visitedStr,sprintf('%s',char(n(find(uniqueNodes==visited(i))))));
        else
            visitedStr = strcat(visitedStr,sprintf('%s�',char(n(find(uniqueNodes==visited(i))))));
        end
    end
    
    visitedStr = strcat(visitedStr,sprintf(']'));
    % Display current node
    if(~isempty(currentNode))
        node = n(currentNode);        
        currentNodeStr = sprintf('[%s]',node{1,1}(1,1));
    else
        currentNodeStr = sprintf('[]');
    end
    array = {iteration currentNodeStr queueStr visitedStr comments};
    tableIteration = cell2table(array);
end
