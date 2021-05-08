%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is for calculating the stiffness of a Snapology-based cubic
%%% unit cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all
clf
% -----------------------------------------------------------------------------------------
% Input the defect case we want to calculate the stiffness
defectCase='a.1';
% -----------------------------------------------------------------------------------------
% Basic Info
L=1;
t=0.01*L;
v=1/3;
%%% Bar's stiffness
NormKst=t*(L^2-v*L^2)/(2*(L^2)*(1-v^2));
NormKsh=t*v*((2*L^2)^(3/2))/(2*(L^3)*(1-v^2));
NormKh=(0.441*L*t^3*(1/t)^(1/3))/(24*(1-v^2));
NormKb=(0.441*t^3*(L/t)^(1/3))/(12*(1-v^2));
%%%
Stiff1=1; Stiff2=1;
% F=-1.5*10e-2;  % maximum stiffness' loading
F=-8*10e-3;  % minimum stiffness' loading
scaleV=1;
nameFolder=[pwd,'\state#1',defectCase,'\'];   % Find the folder with data
PriNodes=L*scaleV*load(strcat(nameFolder,'mnode_state#1.txt'));  
extrudedEdge=load((strcat(nameFolder,'medge_state#1.txt')));
extrudedFace=load((strcat(nameFolder,'mface_state#1.txt')));
HingeEx=load((strcat(nameFolder,'mhinge_state#1.txt')));
LatVinfo=load((strcat(nameFolder,'mlatInfo_state#1.txt')));
nodeHingeEx=HingeEx;
initialNodes=PriNodes;
% Find lattice vector 
initialLatticeVector=[initialNodes(LatVinfo(1,1),:)-initialNodes(LatVinfo(1,2),:); initialNodes(LatVinfo(2,1),:)-initialNodes(LatVinfo(2,2),:);
    initialNodes(LatVinfo(3,1),:)-initialNodes(LatVinfo(3,2),:);];   
initialV=abs(dot(cross(initialLatticeVector(1,:),initialLatticeVector(2,:)),initialLatticeVector(3,:)));
extrudedFaceCell=mat2cell(extrudedFace, zeros(1,size(extrudedFace,1))+1 ,4);
%  For merging overplapping nodes
% [PriNodes,extrudedEdge,nodeHingeEx,extrudedFaceCell]=mergeNodeAdvance(PriNodes,extrudedEdge,nodeHingeEx,extrudedFaceCell);

% -----------------------------------------------------------------------------------------
% Plot the undeformed state
Polyhedron.edge=extrudedEdge;      Polyhedron.node=PriNodes;
viewPoint=[139  20];
Polyhedron.face=extrudedFaceCell;
figure(1)
Plot(Polyhedron,viewPoint,'text','Snapology',1)
%%%%
extrudedFace=cell2mat(extrudedFaceCell);
% For recording the min stiffness % For recording the max stiffness
minStiffness=1;
maxStiffness=0;
[ShearingSprings]=FindShearingSprings(extrudedFace);
[BendingSprings]=FindBendingSprings(extrudedFace);
[StretchingSprings]=FindStretchingSprings(extrudedFace);
% -----------------------------------------------------------------------------------------


% -----------------------------------------------------------------------------------------
% Calculating stiffness 
tick=0.1; % the density of the plot, the smaller the tick, the high resolution of directional stiffness plot
for theta=0:tick:pi   
for gama=0:tick:pi  
% theta=0*pi/180;
% gama=0*pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Basic Imformation
nodes=PriNodes;
Sizenode=size(nodes,1);
% PlotShapeEXternalFace(extrudedFace,nodes);
latticeVector=initialLatticeVector;

[nodes]=RotateUnitCell(nodes,gama,theta);
[latticeVector]=RotateUnitCell(latticeVector,gama,theta);
% PlotShapeEXternalFace(extrudedFace,nodes);   %% For examing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Construction of Jst %%%%%%%%%
[Jst]=buildJst(StretchingSprings,nodes);
Kst=NormKst*eye(size(StretchingSprings,1));  %% L=?
globaJst=Jst'*Kst*Jst;
augment1=zeros(size(globaJst,1),size(latticeVector,2));
augment2=zeros(size(globaJst,1)+size(latticeVector,2),size(latticeVector,2));
globaJst=[[globaJst,augment1]',augment2]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Construction of Jsh %%%%%%%%%
[Jsh]=buildJsh(ShearingSprings,nodes);
Ksh=NormKsh*eye(size(ShearingSprings,1));   %% L=?
globaJsh=Jsh'*Ksh*Jsh;
augment1=zeros(size(globaJsh,1),size(latticeVector,2));
augment2=zeros(size(globaJsh,1)+size(latticeVector,2),size(latticeVector,2));
globaJsh=[[globaJsh,augment1]',augment2]';

%%%%%%%%% Construction of Jh %%%%%%%%%
Jh=zeros(size(nodeHingeEx,1),size(nodes,1)*3);
thetaHinge=zeros(size(nodeHingeEx,1));
for i=1:size(nodeHingeEx,1)
    nodes(nodeHingeEx(i,:),:);
    index(1:3:12)=3*nodeHingeEx(i,:)-2;
    index(2:3:12)=3*nodeHingeEx(i,:)-1;
    index(3:3:12)=3*nodeHingeEx(i,:);
    [Jh(i,index),thetaHinge(i)]=JacobianHinge(nodes(nodeHingeEx(i,:),:));
end

Kh=NormKh*eye(size(nodeHingeEx,1)); %% L=?
globaJh=Jh'*Kh*Jh;
augment1=zeros(size(globaJh,1),size(latticeVector,2));
augment2=zeros(size(globaJh,1)+size(latticeVector,2),size(latticeVector,2));
globaJh=[[globaJh,augment1]',augment2]';

%%%%%%%%% Construction of Jb %%%%%%%%%
Jb=zeros(size(BendingSprings,1),size(nodes,1)*3);
thetaHinge=zeros(size(BendingSprings,1));
for i=1:size(BendingSprings,1)
    nodes(BendingSprings(i,:),:);
    index(1:3:12)=3*BendingSprings(i,:)-2;
    index(2:3:12)=3*BendingSprings(i,:)-1;
    index(3:3:12)=3*BendingSprings(i,:);
    [Jb(i,index),thetaHinge(i)]=JacobianHinge(nodes(BendingSprings(i,:),:));
end
Kb=NormKb*eye(size(BendingSprings,1)); %% L=?
globaJb=Jb'*Kb*Jb;
augment1=zeros(size(globaJb,1),size(latticeVector,2));
augment2=zeros(size(globaJb,1)+size(latticeVector,2),size(latticeVector,2));
globaJb=[[globaJb,augment1]',augment2]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PBC and Transformation Matrix
[PBCinfo]=findPBCinfo(nodes,latticeVector);
%%% This matrix includes all the information of all PBCs: ...
% 1. Row represents a individual PBC.
% 2. First column represents the slave node, second column represents the master node,
%    third column represents which lattice vector is used in this PBC.

[T]=buildTransformationMatrix2D(PBCinfo,latticeVector,Sizenode);

%%%%%%%%% Calculate stiffness %%%%%%%%%
%%% Delete the rows and columns corresponding to the nodes with zero displacement
Fixnode=1;
Fixednode=[3*Fixnode-2,3*Fixnode-1,3*Fixnode;];  
R=zeros(3*size(nodes,1)+3,1);    %%% applied forces
R(3*size(nodes,1)+1)=F;               %%% add the force on H11
globaK=globaJst+globaJsh+globaJh+globaJb;               %globaJst+globaJsh+globaJh+globaJb; %+globaJst+globaJsh;                 %globaKh+globaJst;    %%% Correct option : choose globaKst
%globaK=globaJst;
%%%%%%%%% Add Transformation Matrix (Add PBC)
utimateK=T'*globaK*T;
%utimateKrref=rref(utimateK);
utimateR=T'*R;
%%%%%%%%% Delete rows and columns
utimateKafterBC=utimateK;
utimateKafterBC(Fixednode,:)=[];
utimateKafterBC(:,Fixednode)=[];
utimateRafterBC=utimateR;
utimateRafterBC(Fixednode,:)=[];
UafterBC=(pinv(utimateKafterBC)*utimateRafterBC)';
Dispmet=[0,0,0,UafterBC(:,1:4),UafterBC(:,5:end)]';   
OriDispmet=T*Dispmet;
OriDispmet([size(R,1)-2,size(R,1)-1,size(R,1)])=[];
Denodes=zeros(size(nodes,1),3);
for i=1:size(nodes,1)
    Denodes(i,:)=nodes(i,:)+[OriDispmet(3*i-2),OriDispmet(3*i-1),OriDispmet(3*i);];
end
H11=UafterBC(size(UafterBC,2)-2);
H22=UafterBC(size(UafterBC,2)-1);
H33=UafterBC(size(UafterBC,2));
Stiffness(Stiff1,Stiff2)=(F/initialV)/UafterBC(size(UafterBC,2)-2); % The magnitude of force is divided by the prejection of displacement on the unit vector of force
if Stiffness(Stiff1,Stiff2)<minStiffness
minStiffness=Stiffness(Stiff1,Stiff2);
minGama=gama;
minTheta=theta;
end
if Stiffness(Stiff1,Stiff2)>maxStiffness
maxStiffness=Stiffness(Stiff1,Stiff2);
maxGama=gama;
maxTheta=theta;
end
Stiff1=Stiff1+1;
 end
Stiff1=1;
Stiff2=Stiff2+1;
end

% -----------------------------------------------------------------------------------------
% Plot the directional stiffness
[StiffX,StiffY]=meshgrid(0:tick:pi);
[thetaplot,gamaplot]=meshgrid(0:0.01:pi);
gridSitffness = griddata(StiffX,StiffY,Stiffness,thetaplot,gamaplot,'cubic');
figure(1)
mesh(thetaplot,gamaplot,gridSitffness);
hold on
axis([0,pi    0,pi   0,9e-2]);
set(gcf,'color','w')
set(gca,'XTick',(0:pi/2:pi)); set(gca,'xtickLabel',{'0','\pi/2','\pi'},'FontSize',20,'FontWeight','bold');
set(gca,'YTick',(0:pi/2:pi)); set(gca,'YtickLabel',{'0','\pi/2','\pi'},'FontSize',20,'FontWeight','bold');
xlabel('\theta','FontSize',20,'FontWeight','bold'),ylabel('\gamma','FontSize',20,'FontWeight','bold');
axis on
axis square
grid off
view(0,90)
colorbar
maxStiffValue=maxStiffness;
caxis([0,maxStiffValue])
colorbar('YTick', [0 maxStiffValue],'YTickLabel',{'0','13\times10^{-3}'})
% nameFig=strcat('exnode3','fig');
% saveas(gcf, nameFig)
% -----------------------------------------------------------------------------------------
% Plot the deformed configuration 
figure(2)
dPolyhedron.edge=extrudedEdge;      dPolyhedron.node=Denodes;
viewPoint=[-127  22];
dPolyhedron.face=extrudedFaceCell;
Plot(dPolyhedron,viewPoint,'noxt','Snapology',1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------------
%%%%%%% Sub Function %%%%%%%%%%%%%%%%
function [PBCinfo]=findPBCinfo(nodes,laV)
n=1;
for i=1:size(laV,1)
    for j=1:size(nodes,1)
        for k=1:size(nodes,1)
            if j~=k
                judge=nodes(j,:)-nodes(k,:)-laV(i,:);
                if norm(judge)<=1e-8
                    PBCinfo(n,:)=[j k i];
                    n=n+1;
                end
            end
        end
    end
end
end

function [MaxStiffArray, MinStiffArray]= FindMaxMinStiff (Stiffness)
[MaxStiffness,~] = max(Stiffness(:));
% [Imax_row, Imax_col] = ind2sub(size(Stiffness),Imax);  % Imax is the ~
m=1; n=1;
[MinStiffness,~] = min(Stiffness(:));
% [Imin_row, Imin_col] = ind2sub(size(Stiffness),Imin);  % Imin is the ~
for i=1:size(Stiffness,1) % Rows mean the gama angle
    for j=1:size(Stiffness,2)  % Columns mean theta angle
        if  abs(Stiffness(i,j)-MaxStiffness)<9e-8
            MaxStiffArray(m,:)=[j*180/size(Stiffness,2),i*180/size(Stiffness,1),Stiffness(i,j);];
            m=m+1;
        end
        if  abs(Stiffness(i,j)-MinStiffness)<9e-8
            MinStiffArray(n,:)=[j*180/size(Stiffness,2),i*180/size(Stiffness,1),Stiffness(i,j);];
            n=n+1;
        end
    end
end
plot3(MaxStiffArray(:,1)*pi/180,MaxStiffArray(:,2)*pi/180,10*MaxStiffArray(:,3),'o','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
plot3(MinStiffArray(:,1)*pi/180,MinStiffArray(:,2)*pi/180,MinStiffArray(:,3)+10e-3,'o','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','w')
hold on
end

function [ShearingSprings]=FindShearingSprings(Face)
ShearingSprings=zeros(2*size(Face,1),2);
for i=1:size(Face,1)
    ShearingSprings(2*i-1,:)=[Face(i,1),Face(i,3)];
    ShearingSprings(2*i,:)=[Face(i,2),Face(i,4)];
end
end

function [BendingSprings]=FindBendingSprings(Face)
BendingSprings=zeros(size(Face,1),4);
for i=1:size(Face,1)
    BendingSprings(i,:)=[Face(i,3),Face(i,1),Face(i,4),Face(i,2)];
end
end


function [mnode,medge,mhinge,mface]=mergeNodeAdvance(node,edge,hinge,face)
% This version will delete the overplapping nodes and update the
% information in the edge and face
nM=1; % the cell matrix index colloect the information of overlapping node
for a=1:size(node,1)
    for b=1:size(node,1)
        if a~=b
            distance=sqrt((node(b,1)-node(a,1))^2+(node(b,2)-node(a,2))^2+(node(b,3)-node(a,3))^2);
            if distance<=10e-4
                mNode(nM,:)=[a b];
                nM=nM+1;
            end
        end
    end
end
mNodeUp=unique(sort(mNode,2),'rows');  % Keep the key information of deleteing nodes, the second column is the nodes needed to be deleted
% Keep only one nodes when one position has more than two nodes
deRo=[];
for i=1:size(mNodeUp,1)
    if ismember(mNodeUp(i,1),mNodeUp(:,2))==1
        if ismember(mNodeUp(i,2),mNodeUp(:,2))==1
            deRo=[deRo;i;];
        end
    end
end
mNodeUp(deRo,:)=[];
slave=mNodeUp(:,2);
% Replace the index of slave overlapped nodes
edgeUp=edge;
hingeUp=hinge;
faceUp=face;
for i=1:size(slave,1)
    edgeUp(edge==slave(i))=mNodeUp(i,1);  % update the element in the edges
    hingeUp(hinge==slave(i))=mNodeUp(i,1);
    for j=1:size(face,1)
        faceUp{j}(face{j}==slave(i))=mNodeUp(i,1);
    end
end
% edgeUp=unique(sort(edgeUp,2),'rows');  % delete redundant edge constraints, but did not find a way to effiently deleted the
% Delete overlapping node
nindex = linspace(1,size(node,1),size(node,1));
nindex(mNodeUp(:,2))=[];
node(mNodeUp(:,2),:)=[];
% Update the new index to the edge and face
edgeUp2=edgeUp;
hingeUp2=hingeUp;
faceUp2=faceUp;
for i=1:size(nindex,2)
    edgeUp2(edgeUp==nindex(i))=i;  % update the element in the edges
    hingeUp2(hingeUp==nindex(i))=i;
    for j=1:size(faceUp,1)
        faceUp2{j}(faceUp{j}==nindex(i))=i;
    end
end
% edgeUp2=unique(sort(edgeUp2,2),'rows');
mnode=node;  medge=edgeUp2;  mface=faceUp2; mhinge= hingeUp2;% mface is not uniqued
end

function [J,t]=JacobianHinge(p0)   % This function is writtern by Dr. Bas Overvelde
%Node coordinates
p=reshape(p0',12,1);  %p is 4 nodes' coodinates
%Vectors
a=p(4:6)-p(1:3);    %Rotation axis
b=p(7:9)-p(1:3);    %Vector on face 1
c=p(10:12)-p(1:3);  %Vector on face 2

da=[-1 0 0
    0 -1 0
    0 0 -1
    1 0 0
    0 1 0
    0 0 1
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0]';
db=[-1 0  0
    0  -1 0
    0  0 -1
    0  0  0
    0  0  0
    0  0  0
    1  0  0
    0  1  0
    0  0  1
    0  0  0
    0  0  0
    0  0  0]';
dc=[-1 0  0
    0  -1 0
    0  0 -1
    0  0  0
    0  0  0
    0  0  0
    0  0  0
    0  0  0
    0  0  0
    1  0  0
    0  1  0
    0  0  1]';

ka=norm(a);
kab=sqrt((a'*a)*(b'*b)-(a'*b)^2);   %%cab
kca=sqrt((c'*c)*(a'*a)-(c'*a)^2);   %%cca
na=a/ka;
nab=crossvector(a,b)/kab;
nca=crossvector(c,a)/kca;

detf=na'*(crossvector(nab,nca));
dotf=nab'*nca;
t = atan2(detf,dotf);

dna=da/ka-a*(a'*da)/ka^3;

dna2=[-1/ka+(a(1)^2)/ka^3  a(2)*a(1)/ka^3  a(3)*a(1)/ka^3
    a(1)*a(2)/ka^3  -1/ka+(a(2)^2)/ka^3  a(3)*a(2)/ka^3
    a(1)*a(3)/ka^3  a(2)*a(3)/ka^3  -1/ka+(a(3)^2)/ka^3
    1/ka-(a(1)^2)/ka^3  -a(2)*a(1)/ka^3  -a(3)*a(1)/ka^3
    -a(1)*a(2)/ka^3  +1/ka-(a(2)^2)/ka^3  -a(3)*a(2)/ka^3
    -a(1)*a(3)/ka^3  -a(2)*a(3)/ka^3  +1/ka-(a(3)^2)/ka^3
    0  0  0
    0  0  0
    0  0  0
    0  0  0
    0  0  0
    0  0  0]';
max(max(dna-dna2));
dkab=1/(kab)*((da'*a)*(b'*b)+(a'*a)*(db'*b)-(a'*b)*((da'*b)+(db'*a)));
dkca=1/(kca)*((dc'*c)*(a'*a)+(c'*c)*(da'*a)-(c'*a)*((dc'*a)+(da'*c)));
dnab=1/kab^2*((crossvector(da,b)+crossvector(a,db))*kab-crossvector(a,b)*dkab');
dnca=1/kca^2*((crossvector(dc,a)+crossvector(c,da))*kca-crossvector(c,a)*dkca');

ddetf=dna'*(crossvector(nab,nca))+(na'*(crossvector(dnab,nca)+crossvector(nab,dnca)))';
ddotf=dnab'*nca+dnca'*nab;

J = (-dotf*ddetf+detf*ddotf)/(detf^2+dotf^2);
end

function [Jst]=buildJst(extrudedEdge,DeNode)
Jst=zeros(size(extrudedEdge,1),3*size(DeNode,1));

for i=1:size(extrudedEdge,1)
    x1=DeNode(extrudedEdge(i,1),1); x2=DeNode(extrudedEdge(i,2),1);
    y1=DeNode(extrudedEdge(i,1),2); y2=DeNode(extrudedEdge(i,2),2);
    z1=DeNode(extrudedEdge(i,1),3); z2=DeNode(extrudedEdge(i,2),3);
    patialDe1=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*x1-2*x2);   %%x1
    patialDe2=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*y1-2*y2);   %%y1
    patialDe3=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*z1-2*z2);   %%z1
    patialDe4=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*x2-2*x1);   %%x2
    patialDe5=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*y2-2*y1);   %%y2
    patialDe6=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*z2-2*z1);   %%z2
    
    Jst(i,3*extrudedEdge(i,1)-2)=patialDe1; Jst(i,3*extrudedEdge(i,1)-1)=patialDe2; Jst(i,3*extrudedEdge(i,1))=patialDe3;
    Jst(i,3*extrudedEdge(i,2)-2)=patialDe4; Jst(i,3*extrudedEdge(i,2)-1)=patialDe5; Jst(i,3*extrudedEdge(i,2))=patialDe6;
    
end
end
function [StretchingSprings]=FindStretchingSprings(extrudedFace)
StretchingSprings=zeros(size(extrudedFace,1)*size(extrudedFace,2),2);
j=1;
for i=1:size(extrudedFace,1)
    StretchingSprings(j,:)=[extrudedFace(i,1),extrudedFace(i,2);];
    StretchingSprings(j+1,:)=[extrudedFace(i,2),extrudedFace(i,3);];
    StretchingSprings(j+2,:)=[extrudedFace(i,3),extrudedFace(i,4);];
    StretchingSprings(j+3,:)=[extrudedFace(i,4),extrudedFace(i,1);];
    j=j+4;
end
end
function [Jsh]=buildJsh(ShearingSprings,DeNode)
Jsh=zeros(size(ShearingSprings,1),3*size(DeNode,1));
for i=1:size(ShearingSprings,1)
    x1=DeNode(ShearingSprings(i,1),1); x2=DeNode(ShearingSprings(i,2),1);
    y1=DeNode(ShearingSprings(i,1),2); y2=DeNode(ShearingSprings(i,2),2);
    z1=DeNode(ShearingSprings(i,1),3); z2=DeNode(ShearingSprings(i,2),3);
    patialDe1=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*x1-2*x2);   %%x1
    patialDe2=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*y1-2*y2);   %%y1
    patialDe3=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*z1-2*z2);   %%z1
    patialDe4=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*x2-2*x1);   %%x2
    patialDe5=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*y2-2*y1);   %%y2
    patialDe6=(0.5*((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(-0.5))*(2*z2-2*z1);   %%z2
    
    Jsh(i,3*ShearingSprings(i,1)-2)=patialDe1; Jsh(i,3*ShearingSprings(i,1)-1)=patialDe2; Jsh(i,3*ShearingSprings(i,1))=patialDe3;
    Jsh(i,3*ShearingSprings(i,2)-2)=patialDe4; Jsh(i,3*ShearingSprings(i,2)-1)=patialDe5; Jsh(i,3*ShearingSprings(i,2))=patialDe6;
    
end
end

function [T]=buildTransformationMatrix2D(PBCinfo,latticeVector,sizenodes)
DimeH=3;   %%% The dimension of displacement gradient
Trows=3*sizenodes+DimeH;   %%% the number of row of the Transformation matrix should equal the number of slave nodes*3 + DimH
Tcolumns=3*sizenodes+DimeH-3*size(PBCinfo,1);
T=zeros(Trows,Tcolumns);
%%%% For constructing T conveniently, use Umaster as a index
Uslave=zeros(Trows,1);
for j=1:Trows
    Uslave(j)=j;
end
Umaster=Uslave;
DeleteRows=zeros(3*size(PBCinfo,1),1);
for j=1:size(PBCinfo,1)
    DeleteRows(3*j-2)=3*PBCinfo(j,1)-2;
    DeleteRows(3*j-1)=3*PBCinfo(j,1)-1;
    DeleteRows(3*j)=3*PBCinfo(j,1);
end
Umaster(DeleteRows,:)=[];  %%% Obtain the information of rows have been deleted.
%%%% Start constructing T matrix
for t=1:(size(T,1)/3)            %%% t is for constructing every row of T  ;;  all the coefficient 2 is for 2D cases
    for PBCinfoIndex=1:size(PBCinfo,1)        %%% PBCinfoIndex is for identify the PBC when constructing every row of T
        switch t
            case PBCinfo(PBCinfoIndex,1)     %%% t is one of slave nodes  %PBCinfo=[6,1,1;5,2,1;3,8,2;4,7,2;];
                for k=1:size(Umaster,1)   %%% this loop is for constructing each column of the matrix
                    if 3*PBCinfo(PBCinfoIndex,2)-2==Umaster(k)   %%% identify the column when construct a PBC corresponding the t row
                        T(3*t-2,k)=1;
                        T(3*t-2,Tcolumns-2)=latticeVector(PBCinfo(PBCinfoIndex,3),1);
                        T(3*t-1,k+1)=1;
                        T(3*t-1,Tcolumns-1)=latticeVector(PBCinfo(PBCinfoIndex,3),2);
                        T(3*t,k+2)=1;
                        T(3*t,Tcolumns)=latticeVector(PBCinfo(PBCinfoIndex,3),3);
                    end
                end
            otherwise                    %%% t is one of master nodes
                for k=1:size(Umaster,1)
                    if 3*t-2==Umaster(k)
                        T(3*t-2,k)=1;
                        T(3*t-1,k+1)=1;
                        T(3*t,k+2)=1;
                    end
                end
        end
    end
end
end

function [FictitiousNode]=RotateUnitCell(FictitiousNode,gama,theta)
Rz=[cos(gama),-sin(gama),0;sin(gama),cos(gama),0;0,0,1;];
yprime1=Rz*[0,1,0;]';
xprime1=Rz*[1,0,0;]';
zprime1=[0,0,1]';
xprime2=cos(theta)*(xprime1-dot(xprime1,yprime1)*yprime1)+sin(theta)*cross(yprime1,xprime1)+dot(xprime1,yprime1)*yprime1;
zprime2=cos(theta)*(zprime1-dot(zprime1,yprime1)*yprime1)+sin(theta)*cross(yprime1,zprime1)+dot(zprime1,yprime1)*yprime1;
yprime2=yprime1;

FictitiousNode=FictitiousNode(:,1)*xprime2'+FictitiousNode(:,2)*yprime2'+FictitiousNode(:,3)*zprime2';
end

function w=crossvector(u,v)
w=[u(2,:).*v(3,:)-u(3,:).*v(2,:);-u(1,:).*v(3,:)+u(3,:).*v(1,:);u(1,:)*v(2,:)-u(2,:)*v(1,:)];
end


function Plot(Polyhedron,viewPoint,Text,object,pic)
str={};
if strcmp(Text,'text')==1
    for i=1:1000
        s=num2str(i);
        str{i}=s;
    end
end
% Loop from i-th polyhedron
for i=1:length(Polyhedron)
    if strcmp(object,'interPolyhedra')==1
        facecolor=[199,21,133]/256;
    elseif strcmp(object,'Snapology')==1
        excolor='m';
        Lightpoint=2*[cosd(viewPoint(2))*sind(viewPoint(1)), -cosd(viewPoint(2))*cosd(viewPoint(1)), sind(viewPoint(2))];
    end
    if strcmp(Text,'text')==1
        for kk=1:size(Polyhedron(i).node,1)  % text the nodes
            text(Polyhedron(i).node(kk,1),Polyhedron(i).node(kk,2),Polyhedron(i).node(kk,3),str{kk},'fontsize',10);
        end
    end
    for j=1:length(Polyhedron(i).face)
        if strcmp(Text,'text')==1
            indexX=sum(Polyhedron(i).node(:,1))/size(Polyhedron(i).node,1);
            indexY=sum(Polyhedron(i).node(:,2))/size(Polyhedron(i).node,1);
            indexZ=sum(Polyhedron(i).node(:,3))/size(Polyhedron(i).node,1);
            text(indexX, indexY, indexZ,str{i},'fontsize',14,'Color','b');  % indicate j-th polyhedron this is
            faceCenter=sum(Polyhedron(i).node([Polyhedron(i).face{j}(1:end)],:))/length(Polyhedron(i).face{j});
            text( faceCenter(1), faceCenter(2), faceCenter(3),str{j},'fontsize',12,'Color','m');
        end
        if strcmp(object,'Snapology')==1
            patch('Faces',Polyhedron(i).face{j},'Vertices',Polyhedron(i).node,'FaceColor',excolor,'facealpha',0.75,'LineWidth',0.88); %
        end
    end
end
% End loop
view(viewPoint(1),viewPoint(2))
axis equal
axis on
grid off
if strcmp(object,'Snapology')==1
    light('color','w','style','infinite','position',Lightpoint)
    lighting flat
end
material dull
set(gcf,'color','w')
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
set(gca,'ztick',[],'zticklabel',[])
% xlabel('x','FontSize',40,'FontWeight','bold'),ylabel('y','FontSize',40,'FontWeight','bold'),zlabel('z','FontSize',40,'FontWeight','bold')
set(gca,'xcolor','r')
set(gca,'LineWidth',1)
end